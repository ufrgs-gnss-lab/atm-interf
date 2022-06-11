function [N, de, der] = get_atm_pol (e, H, N_coeff, de_coeff, h)
    if (nargin < 3),  N_coeff = [];  end
    if (nargin < 4),  de_coeff = [];  end
    if (nargin < 5),  h = [];  end    
    [N_coeff_default, de_coeff_default] = get_atm_pol_coeff_default ();
    if isempty(N_coeff),   N_coeff = N_coeff_default;  end
    if isempty(de_coeff),  de_coeff = de_coeff_default;  end
    if isempty(h),         h = 0;  end

    N = get_refractivity_pol (N_coeff, h, H);
    de = get_elev_bending_pol (de_coeff, e);
    der = get_elev_bending_rate (de_coeff, e);
end

%%
function [N_coeff, de_coeff] = get_atm_pol_coeff_default ()
    % credit: T. Nikolaidou and F. Geremia-Nievinski (unpublished)
    N_coeff  = [0.00025617 -2.44E-08];  % [N0 dN_dh]: mean refractivity (N=n-1) and refractivity laspse rate
    de_coeff = [5.56947472121108, 1.88401692297586, 1.55363613681730e-05]; % [a b c] in de=c*cotd(e+a/(b+e))/cotd(90+a/(b+90))
end

%%
function Nl = get_refractivity_pol (N_coeff, h, H)
    % refractivity at the antenna:
    Na = get_refractivity_pol_aux (N_coeff, h);

    % refractivity at the surface:
    Ns = get_refractivity_pol_aux (N_coeff, h - H);

    % refractivity at the layer centroid:
    %Nl = (Ns + Na) ./ 2;  % arithmetic mean
    Nl = exp( (log(Ns) + log(Na)) ./ 2 );  % logarithmic mean
end

%%
function Np = get_refractivity_pol_aux (N_coeff, h)
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
function de = get_elev_bending_pol (de_coeff, e)
    a = de_coeff(1);
    b = de_coeff(2);
    c = de_coeff(3);

    num = cotd(e + a./(b + e));
    den = cotd(90+ a./(b +90));
    de = c*(1-num./den);
end

%%
function der = get_elev_bending_rate (de_coeff, e)
    a = de_coeff(1);
    b = de_coeff(2);
    c = de_coeff(3);

    num_arg = e + a./(b + e);
    den_arg = 90+ a./(b +90);
    %num = cotd(num_arg);
    den = cotd(den_arg);
    % de = c*(1-num./den);
    dnum_de = -(pi/180)*cscd(num_arg).^2.*(1 + a.*(-1./(b+e).^2));
    %dden_de = 0;
    %dde_de = c*(0 + dnum_de./den + num.*(-1).*(1./den.^2).*dden_de);
    der = c*(-1)*dnum_de./den;  % dde/de
end
