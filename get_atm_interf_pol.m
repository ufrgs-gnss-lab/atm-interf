function [dt, da, dg, Ht, Ha, Hg, N, de] = get_atm_interf_pol (e, H, N_coeff, de_coeff, h)
% GET_ATM_INTERF_POL:  Closed-form interferometric atmospheric delay, using polynomial atmospheric model.
% 
% SYNTAX:
%    [dt, da, dg, Ht, Ha, Hg] = get_atm_interf_pol (e, H, N, de);
%
% INPUT: 
%    e: [vector] satellite elevation angle (in degrees)
%    H: [vector] reflector depth or antenna height above reflector (in meters)
%    N_coeff: [vector] linear refraction coefficients (mean refractivity and refractivity laspse rate: [N0 dN_dh])
%    de_coeff: [vector] angular refraction coefficients (see code for details)
%    h: [vector] antenna altitude above ellipsoid or ellipsoidal height (in meters)
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
% EXAMPLE:
%    e = 45;  % degrees
%    H = 10;  % meters
%    [dt, ~, ~, Ht] = get_atm_interf_pol (e, H)

% TODO: support time-space coefficient models.

    if (nargin < 3),  N_coeff = [];  end
    if (nargin < 4),  de_coeff = [];  end
    if (nargin < 5),  h = [];  end
    [N_coeff_default, de_coeff_default] = get_coeff_default ();
    if isempty(h),         h = 0;  end
    if isempty(de_coeff),  de_coeff = de_coeff_default;  end
    if isempty(N_coeff),   N_coeff = N_coeff_default;  end

    N = get_refractivity (N_coeff, h, H);
    de = get_elev_bending (de_coeff, e);
    
    if (nargout < 4)
        [dt, da, dg] = get_atm_interf_aux (e, H, N, de);
    else
        [dt, da, dg, Ht, Ha, Hg] = get_atm_interf_aux (e, H, N, de);
    end
end

%%
function [N_coeff, de_coeff] = get_coeff_default ()
    % credit: T. Nikolaidou and F. Geremia-Nievinski (unpublished)
    N_coeff  = [0.00025617 -2.44E-08];
    de_coeff = [5.56947472121108, 1.88401692297586, 1.55363613681730e-05];
end

%%
function Nl = get_refractivity (N_coeff, h, H)
    % refractivity at the antenna:
    Na = get_refractivity_aux (N_coeff, h);

    % refractivity at the surface:
    Ns = get_refractivity_aux (N_coeff, h - H);

    % refractivity at the layer centroid:
    %Nl = (Ns + Na) ./ 2;  % arithmetic mean
    Nl = exp( (log(Ns) + log(Na)) ./ 2 );  % logarithmic mean
end

%%
function Np = get_refractivity_aux (N_coeff, h)
    N0 = N_coeff(1);  % refractivity at ellipsoid

    if (numel(N_coeff)==2)
        dN_dh = N_coeff(2);  % refractivity lapse rate
        DN = h*dN_dh;  % reduction from ellipsoid to point altitude
    else
        dlogN_dh = N_coeff(3);  % log-refractivity lapse rate
        DN = exp(h*dlogN_dh);  % reduction from ellipsoid to point altitude
    end
    
    Np = N0 + DN;  % refractivity at point of interest
end

%%
function de = get_elev_bending (de_coeff, e)
    a = de_coeff(1);
    b = de_coeff(2);
    c = de_coeff(3);

    temp1 = cotd(e + a./(b + e));
    temp2 = cotd(90+ a./(b +90));
    de = c*(1-temp1./temp2);
end

