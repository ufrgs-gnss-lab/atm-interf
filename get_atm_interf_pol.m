function [dt, da, dg, Ht, Ha, Hg, N, de, der] = get_atm_interf_pol (e, H, N_coeff, de_coeff, h, opt)
% GET_ATM_INTERF_POL:  Closed-form interferometric atmospheric delay, using polynomial atmospheric model.
% 
% SYNTAX:
%    [dt, da, dg, Ht, Ha, Hg] = get_atm_interf_pol (e, H);
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
%
% OPTIONAL INPUT: 
%    N_coeff: [vector] linear refraction coefficients (see get_atm_pol.m for details)
%    de_coeff: [vector] angular refraction coefficients (see get_atm_pol.m for details)
%    h: [vector] antenna ellipsoidal height or altitude above ellipsoid (in meters)
%    opt: [struct] options (see get_atm_interf_gen.m for details)
% 
% EXAMPLE:
%    e = 45;  % degrees
%    H = 10;  % meters
%    [dt, ~, ~, Ht] = get_atm_interf_pol (e, H)

% TODO: support time-space coefficient models.

    if (nargin < 3),  N_coeff = [];  end
    if (nargin < 4),  de_coeff = [];  end
    if (nargin < 5),  h = [];  end
    if (nargin < 6),  opt = [];  end

    [N, de, der] = get_atm_pol (e, H, N_coeff, de_coeff, h);
    
    if (nargout < 4)
        [dt, da, dg] = get_atm_interf_gen (e, H, N, de, der, opt);
    else
        [dt, da, dg, Ht, Ha, Hg] = get_atm_interf_gen (e, H, N, de, der, opt);
    end
end

