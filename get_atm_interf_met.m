function [dt, da, dg, Ht, Ha, Hg, N, de, der] = get_atm_interf_met (e, H, P, T, s, r, opt)
% GET_ATM_INTERF_MET:  Closed-form interferometric atmospheric delay, using in-situ meteorological data.
% 
% SYNTAX:
%    [dt, da, dg, Ht, Ha, Hg, N, de] = get_atm_interf_met (e, H, P, T, s);
%
% INPUT: 
%    e: [vector] satellite elevation angle (in degrees)
%    H: [vector] reflector depth or antenna height above reflector (in meters)
%    P: [scalar] pressure at the antenna (in pascals)
%    T: [scalar] temperature at the antenna (in kelvin)
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
%
% OPTIONAL INPUT: 
%    s: [scalar] specific humidity at the antenna (in kg/kg)
%    r: [scalar] temperature lapse rate (in kelvin per meter, K/m)
%    opt: [struct] options
%    opt.thin: [scalar, boolean] assume thin layer? defaults to false
%    opt.bennet: [struct] see get_bending_bennet.m
%    opt.* (see get_atm_interf_gen.m for other options)
%
% EXAMPLE:
%    e = 45;  % degrees
%    H = 10;  % meters
%    [dt, ~, ~, Ht] = get_atm_interf_met (e, H)

    if (nargin < 3),  P = [];  end
    if (nargin < 4),  T = [];  end
    if (nargin < 5),  s = [];  end
    if (nargin < 6),  r = [];  end
    if (nargin < 7),  opt = [];  end
    tmp = get_meteo_const();
    if isempty(P),  P = tmp.std_pressure;     end
    if isempty(T),  T = tmp.std_temperature;  end
    if isempty(r),  r = tmp.std_lapse_rate;   end
    if isempty(s),  s = 0;  end
    opt_default = struct('thin',false, 'bennet',[]);
    opt = structmergenonempty(opt_default, opt);

    N = get_refractivity (H, P, T, s, r, opt.thin);
    [de, der] = get_bending_bennet (e, P, T, opt.bennet);

    if (nargout < 4)
        [dt, da, dg] = get_atm_interf_gen (e, H, N, de, der, opt);
    else
        [dt, da, dg, Ht, Ha, Hg] = get_atm_interf_gen (e, H, N, de, der, opt);
    end
end

%%
function N = get_refractivity (H, P, T, s, r, thin)
    if thin
        N = 1e-6*calculate_refractivity (P, T, s);
        %N = calculate_refractivity (P, T, s);  % WRONG!
        return;
    end
    
    Pa = P;  Ta = T;  sa = s;
    dh = -H;  % surface is below antenna
    [Ps, Ts] = reduce_pressure (Pa, Ta, r, dh);  ss = sa;

    Na = 1e-6*calculate_refractivity (Pa, Ta, sa);
    Ns = 1e-6*calculate_refractivity (Ps, Ts, ss);
    N = logavg(Ns, Na);
    %N = (Ns + Na) ./ 2;  % arithmetic mean
end
