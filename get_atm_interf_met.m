function [dt, da, dg, Ht, Ha, Hg, N, de] = get_atm_interf_met (e, H, P, T, s, r, thin)
% GET_ATM_INTERF_MET:  Closed-form interferometric atmospheric delay, using in-situ meteorological data.
% 
% SYNTAX:
%    [dt, da, dg, Ht, Ha, Hg, N, de] = get_atm_interf_met (e, H, P, T, s);
%
% INPUT: 
%    e: [vector] satellite elevation angle (in degrees)
%    H: [vector] reflector depth or antenna height above reflector (in meters)
%    P: [scalar] pressure at the antenna (in pascals)
%    T: [scalar] temperature at the antenna (in Kelvin)
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
%
% OPTIONAL INPUT: 
%    s: [scalar] specific humidity at the antenna (in kg/kg)
%    r: [scalar] temperature lapse rate (in K/m, kelvin per meter)
%    thin: [scalar] assume thin layer? (Boolean)
%
% EXAMPLE:
%    e = 45;  % degrees
%    H = 10;  % meters
%    [dt, ~, ~, Ht] = get_atm_interf_met (e, H)

    if (nargin < 3),  P = [];  end
    if (nargin < 4),  T = [];  end
    if (nargin < 5),  s = [];  end
    if (nargin < 6),  r = [];  end
    if (nargin < 7),  thin = [];  end
    tmp = get_meteo_const();
    if isempty(P),  P = tmp.std_pressure;     end
    if isempty(T),  T = tmp.std_temperature;  end
    if isempty(r),  r = tmp.std_lapse_rate;   end
    if isempty(s),  s = 0;  end
    if isempty(thin),  thin = false;  end

    if thin
        Nl = 1e-6*calculate_refractivity (P, T, s);
        %Nl = calculate_refractivity (P, T, s);  % WRONG!
    else
        Pa = P;  Ta = T;  sa = s;
        dh = -H;  % surface is below antenna
        [Ps, Ts] = reduce_pressure (Pa, Ta, r, dh);  ss = sa;

        Na = 1e-6*calculate_refractivity (Pa, Ta, sa);
        Ns = 1e-6*calculate_refractivity (Ps, Ts, ss);
        Nl = exp( (log(Ns) + log(Na)) ./ 2 );  % logarithmic mean
        %Nl = (Ns + Na) ./ 2;  % arithmetic mean
    end
        
    de = get_bending_bennet (e, P, T);

    if (nargout < 4)
        [dt, da, dg] = get_atm_interf_aux (e, H, N, de);
    else
        [dt, da, dg, Ht, Ha, Hg] = get_atm_interf_aux (e, H, N, de);
    end
end

