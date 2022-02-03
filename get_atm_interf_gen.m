function [dt, da, dg, Ht, Ha, Hg] = get_atm_interf_gen (e, H, N, de)
% GET_ATM_INTERF_GEN:  Closed-form interferometric atmospheric delay, generic.
% 
% SYNTAX:
%    [dt, da, dg, Ht, Ha, Hg] = get_atm_interf_gen (e, H, N, de);
%
% INPUT: 
%    e: [vector] satellite elevation angle (in degrees)
%    H: [vector] reflector depth or antenna height above reflector (in meters)
%    N: [vector] layer average refractivity (N=n-1, unitless)
%    de: [vector] satellite elevation angle bending (in degrees)
% 
% OUTPUT: 
%    dt: [vector] total delay (in meters)
%    da: [vector] along-path delay (in meters)
%    dg: [vector] geometric delay (in meters)
%    Ht: [vector] total atmospheric altimetry correction (in meters)
%    Ha: [vector] along-path atmospheric altimetry correction (in meters)
%    Hg: [vector] geometric atmospheric altimetry correction (in meters)
% 
% EXAMPLE:
%    e = 45;  % degrees
%    H = 10;  % meters
%    N = 0;  % unitless
%    de = 0;  % degrees
%    [dt, ~, ~, Ht] = get_atm_interf_gen (e, H, N, de)

    da = 2*H*N./sind(e+de);
    dg = 2*H.*(sind(e+de)-sind(e));
    %dg = 2*H*N.*(sind(e+de)-sind(e));  % WRONG!
    dt = da + dg;

    if (nargout < 4),  return;  end
    %TODO: implement closed-form heights.
    Ht = -get_height_from_delay (dt, e);
    Ha = -get_height_from_delay (da, e);
    Hg = -get_height_from_delay (dg, e);
end

