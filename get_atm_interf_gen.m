function [dt, da, dg, Ht, Ha, Hg, N, de, der] = get_atm_interf_gen (e, H, N, de, der, opt)
% GET_ATM_INTERF_GEN:  Closed-form interferometric atmospheric delay, generic.
% 
% SYNTAX:
%    [dt, da, dg, Ht, Ha, Hg] = get_atm_interf_gen (e, H, N, de);
%
% INPUT: 
%    e: [vector] satellite elevation angle (in degrees)
%    H: [vector] reflector depth or antenna height above reflector (in meters)
%    N: [vector] layer average refractivity (N=n-1, NOT 1e6*(n-1), unitless)
%    de: [vector] satellite elevation angle bending (in degrees)
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
%    der: [vector] rate of change of elevation bending w.r.t. elevation angle (in degrees per degree)
%    opt: [struct] closed formula options
%    opt.H_approximate [scalar, boolean] obtain altimetry correction approximately? (defaults to false)
%    opt.H_hybrid [scalar, boolean] obtain altimetry correction via numerical derivative of analytical delay? (defaults to false)
%    opt.der_numerical [scalar, boolean] obtain elevation bending rate numerically? (defaults to false)
%    opt.numerical_noend [scalar, boolean] discard end points when calculating results numerically? (defaults to true)
% 
% EXAMPLE:
%    e = 45;  % degrees
%    H = 10;  % meters
%    N = 0;  % unitless
%    de = 0;  % degrees
%    [dt, ~, ~, Ht] = get_atm_interf_gen (e, H, N, de)

    if (nargin < 5),  der = [];  end
    if (nargin < 6),  opt = struct();  end
    opt = structmergenonempty(get_opt_default(), opt);

    da = 2*H*N./sind(e+de);
    dg = 2*H.*(sind(e+de)-sind(e));
    dt = da + dg;

    if (nargout < 4),  return;  end
    
    if opt.H_hybrid
        Ha = -get_height_from_delay (da, e, [], opt.numerical_noend);
        Hg = -get_height_from_delay (dg, e, [], opt.numerical_noend);
        Ht = -get_height_from_delay (dt, e, [], opt.numerical_noend);
        if (nargout >= 9),  der = NaN(size(de));  end
        return;
    end

    if isempty(der) || opt.der_numerical
        der = gradient_all(de, e, [], opt.numerical_noend);
        if isscalar(de) || any(der == 0)
            warning('matlab:get_atm_interf_gen:scalar', ...
                'Input dde_de invalid; output Hg unavailable.');
            der(der == 0) = NaN;
        end
    end
    
    if opt.H_approximate
        tmp = deg2rad(de).*tand(e);
        Ha = H.*N.*cscd(e+de).^2.*(1-tmp).*(1+der);
        Hg = -H.*der + H.*tmp;
    else
        Ha = H.*N.*cscd(e+de).^2.*cosd(e+de)./cosd(e).*(1+der);
        Hg = -H.*der + H.*(sind(de).*tand(e)+1-cosd(de)).*(1+der);
    end
    Ht = Ha + Hg;
end

%%
function opt = get_opt_default ()
    opt = struct();
    opt.H_approximate = false;
    opt.H_hybrid = false;
    opt.der_numerical = false;
    opt.numerical_noend = true;
end
