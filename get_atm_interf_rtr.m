function [dt, da, dg, Ht, Ha, Hg, N, de, der, rtr] = get_atm_interf_rtr (e, H, pos, date, rtr, opt)
% GET_ATM_INTERF_RTR:  Closed-form interferometric atmospheric delay, using ray-tracer results.
% 
% SYNTAX:
%    [dt, da, dg, Ht, Ha, Hg, N, de, der, rtr] = get_atm_interf_rtr (e, H, rtr, opt);
%
% INPUT: 
%    e: [vector] satellite elevation angle (in degrees)
%    H: [vector] reflector depth or antenna height above reflector (in meters)
%
% OUTPUT: 
%    dt: [vector] total delay (in meters)
%    da: [vector] along-path delay (in meters)
%    dg: [vector] geometric delay (in meters)
%    Ht: [vector] total atmospheric altimetry correction (in meters)
%    Ha: [vector] along-path atmospheric altimetry correction (in meters)
%    Hg: [vector] geometric atmospheric altimetry correction (in meters)
%    N: [vector] layer average refractivity (N=n-1, unitless)
%    de: [vector] satellite elevation angle bending (in degrees)
%    der: [vector] rate of change of elevation bending w.r.t. elevation angle (in degrees per degree)
%    rtr: [struct] raytracer data (see raytrace_interf_series.m)
% 
% OPTIONAL INPUT: 
%    pos: [vector] position, in geodetic coordinates (latitude, longitude, altitude; in degrees, degrees, meters)
%    date: [vector] date (year month day, in this order)
%    rtr: [struct] raytracer data (see raytrace_interf_series.m)
%    rtr.target: [structure] raytracer input target; as output by set_target_interf.m
%    rtr.atm: [structure] raytracer input atmospheric model; as output by setup_atm_which.m
%    rtr.opt: [structure] raytracer input options/settings; as output by set_opt.m
%    rtr.result: [struct] raytracer output propagation quantities (see raytrace_interf_series.m)
%    rtr.extra: [struct] raytracer output extra results (see raytrace_interf_series.m)
%    opt: [struct] closed formula options
%    opt.n_zenith_only: [scalar, boolean] index of refraction (n) -- use zenithal values only? defaults to false
%    opt.n_field_name: index of refraction (n) -- field name for average value [char] defaults to 'logavg' (see raytrace_interf_more.m)
%    opt.rtr: [structure] supplementary raytracer input options/settings, takes precedence over rtr.opt.
%    opt.* (see get_atm_interf_met.m for other options)
%
% EXAMPLE:
%    e = 45;  % degrees
%    H = 10;  % meters
%    [dt, ~, ~, Ht] = get_atm_interf_rtr (e, H)

    if (nargin < 3),  pos = [];  end
    if (nargin < 4),  date = [];  end
    if (nargin < 5),  rtr = [];  end
    if (nargin < 6),  opt = [];  end
    if isempty(pos),  pos  = [0 0 0];  end
    if isempty(date), date = [2000 0 0];  end
    if isempty(rtr),  rtr  = struct();  end
    opt = structmergenonempty(get_opt_default(), opt);
    
    [N, de, der] = get_atm_rtr (e, H, pos, date, rtr, opt);

    if (nargout < 4)
        [dt, da, dg] = get_atm_interf_gen (e, H, N, de, der);
    else
        [dt, da, dg, Ht, Ha, Hg] = get_atm_interf_gen (e, H, N, de, der);
    end
end

%%
function opt = get_opt_default ()
    opt = struct();
    opt.numerical_noend = true;
    opt.n_zenith_only = false;
    opt.n_field_name = 'logavg';
    opt.height_sat_factor = 1;
    opt.atm_planar = false;
    opt.rtr = struct();
    opt.rtr.tol = 1e-6;
    %opt.rtr.quad_routine = 'adaptive';
    opt.rtr.interf_convergence_check_simplified = false;
    %opt.rtr.interf_approach = 'rigorous';
    opt.rtr.interf_approach = 'rectilinear-mixed';
    opt.rtr.interf_elev_input = 'geometric';
    opt.rtr.interf_series_mesh_elev_rh = false;    
    opt.rtr.interf_series_use_input_elev_geom_for_height = false;
    opt.rtr.interf_refraction_index = true;
end

%%
function [N, de, der] = get_atm_rtr (e, H, pos, date, rtr, opt)
    rtr = get_atm_rtr_aux (e, H, pos, date, rtr, opt);
    get_field2 = @(f, f2) arrayfun(@(x) x.(f).(f2), rtr.extra);    

    % refractivity:
    n = get_field2('n', opt.n_field_name);
    N = n-1;  clear n
    if opt.n_zenith_only
        idx = (e == 90);
          assert(sum(idx)==1)
        N = N(idx);
    end

    % bending angle:
    eg = get_field2('elev','geom');
    ea = get_field2('elev','appar');
    de = ea - eg;
    %myassert(e, eg)  % WRONG! it's only max(abs(e-eg))<rtr.opt.tol/dd_de;

    % bending angle rate w.r.t. elevation angle:
    if opt.rtr.interf_series_use_input_elev_geom_for_height
        e2 = e;
    else
        % (see raytrace_interf_series.m for details)
        e2 = eg;
    end
    der = gradient_all(de, e2, [], opt.numerical_noend);
end

%%
function rtr = get_atm_rtr_aux (e, H, pos, date, rtr, opt2)
    if ~isfieldempty(rtr, 'result') ...
    && ~isfieldempty(rtr, 'extra')
        return;
    end

    if isfieldempty(rtr, 'opt')
        rtr.opt = set_opt();
    end
    rtr.opt = structmergenonempty(rtr.opt, opt2.rtr);
    if opt2.height_sat_factor~=1
      rtr.opt.height_sat = rtr.opt.height_sat*opt2.height_sat_factor;  % EXPERIMENTAL
      %rtr.opt.nominal_radius_earth = 1000*rtr.opt.nominal_radius_earth;  % EXPERIMENTAL
    end
    
    if isfieldempty(rtr, 'target')        
        rtr.target = set_target (pos, [], [], mydatenum(date));        
    end
    
    if isfieldempty(rtr, 'atm')
        rtr.atm = setup_atm_3d_cira (rtr.target, [],  rtr.opt);
        if opt2.atm_planar
            rtr.atm = setup_atm_tanplane (rtr.target, rtr.atm, rtr.opt);  % EXPERIMENTAL
        else
            rtr.atm = setup_atm_sphosc (rtr.target, rtr.atm, rtr.opt);
        end
    end
    
    [rtr.result, ~, ~, ~, rtr.extra] = raytrace_interf_series (...
        rtr.target, rtr.atm, rtr.opt, e, H);
end

