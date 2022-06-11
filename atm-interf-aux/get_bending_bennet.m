function [de, der] = get_bending_bennet (e, P, T, opt)
% GET_BENDING_BENNET Elevation angle bending, via formulation of Bennet (1982).
% 
% INPUT:
%   e: [matrix] elevation angle [degrees] - vacuum or unrefracted
% 
% OUTPUT:
%   de: [matrix] elevation angle bending (refracted minus vacuum) (in degrees)
%   der: [matrix] rate of change of elevation bending w.r.t. elevation angle (in degrees per degree)
% 
% OPTIONAL INPUT:
%   P: [matrix] pressure (in pascal)
%   T: [matrix] temperature (in kelvin)
%   opt: [struct] options
%   opt.tol: [matrix] elevation angle tolerance for inversion (in degrees); 
%     defaults to machine precision for input "e"; set to inf() to disable.
%   opt.form: [char] formulation - 'bowditch' or 'bennet'; defaults to 'bennet'
%   opt.ignore_PT: [scalar, boolean] ignore input P and T? defaults to false
% 
% Bennett, G. (1982). "The Calculation of Astronomical Refraction in Marine
% Navigation". Journal of Navigation, 35(2), 255-259.
% doi:10.1017/S0373463300022037

    if (nargin < 2),  P = [];  end
    if (nargin < 3),  T = [];  end
    if (nargin < 4),  opt = [];  end
    tol_default = eps(e);
    opt_default = struct('tol',tol_default, 'form','bennet', 'ignore_PT',false);
    opt = structmergenonempty(opt_default, opt);

    if isinf(opt.tol)
        [de, der] = get_bending_bennet_aux (e, P, T, opt.form, opt.ignore_PT);
        return;
    end
    
    de0 = 0;
    %de0 = get_bending_saemundsson (e);  % EXPERIMENTAL
    i = 0;
    while true  % (do while)
        h = e + de0;
        %h = e - de0;  % WRONG!
        [de, der] = get_bending_bennet_aux (h, P, T, opt.form, opt.ignore_PT);
        delta = abs(de-de0);
        %if (max(delta) < opt.tol),  return;  end  % WRONG! tol might not be scalar
        if all(delta < opt.tol),  return;  end
        de0 = de;
        i = i + 1;
        %disp([i max(delta)])  % DEBUG
    end
end

%%
function [de, der] = get_bending_bennet_aux (e, P, T, form, ignore_PT)
    % Bennet (1982), p.257, formula G; 
    % "Altitude ranges are expressed in degrees."
    % "unsigned mean refraction (R_M) in minutes of arc"
    % "where h [here, e] is the observed altitude in degrees"
    % Bennet's formula uses the apparent or refracted elevation angle, e'=e+de.
    %Rm = cot( h_rad + (7.31 / (h_rad + 4.4)) );  % WRONG!
    arg = e + (7.31./(e + 4.4));
    arg(arg > 90) = 90;  % avoid negative values near zenith
    Rm = cotd(arg);  % Bennet (1982), p.257, formula G;
    darg_de = 1 - 7.31./(e + 4.4).^2;
    Rmr = -(pi/180)*cscd(arg).^2.*darg_de;  % dRm/de
    
    if isempty(P) || isempty(T) || ignore_PT
        de = Rm./60;
        der = Rmr./60;
        return
    end

    PaToMbar = @(p) p/100;
    %PaToMbar = @(p) p*100;  % WRONG!
    KelvinToCelsius = @(k) k - 273.15;
    %KelvinToCelsius = @(k) k + 273.15;  % WRONG!

    switch lower(char(form))
    case {'bennet','newer'}
        % "For T = 10°C and P = 1010 mb" (...) "If it is accepted that
        % corrections to refraction for non-standard conditions should not
        % be degraded by using formulae that produce significant
        % differences from the generally accepted standard of Garfinkel,
        % then more accurate formulae should be used provided they do not
        % over-complicate the calculation. Rather than use the formulae
        % given by Bowditch or similar forms (see formulae E), it is
        % suggested that the following be used:"
        P = PaToMbar(P);
        T = KelvinToCelsius(T);
        num = (P - 80)./930;
        den = 1 + 8e-5.*(Rm + 39).*(T - 10);
        k = num./den;
        R = k.*Rm;  % Bennet (1982), p.258
        dkr = (num./den.^2).*(-1).*(0+8e-5.*(1+0).*(T-10));  % dk/de
        dRr = dkr.*Rm + k.*Rmr;  % dR/de
        de = R./60;
        der = dRr./60;
    case {'bowditch','older'}
        % "Tables of corrections to altitudes for non-standard conditions
        % are to be found in Bowditch [1977] (...) where t is in degrees F
        % and p in inches of mercury"
        MbarToInchesMercury = @(m) m*0.029530;
        CelsiusToFahrenheit = @(c) 32 + c*9/5;
        p = MbarToInchesMercury(PaToMbar(P));
        t = CelsiusToFahrenheit(KelvinToCelsius(T));
        R = (510 ./ (460 + t)) * (p./29.83) .*Rm;  % Bowditch apud. Bennet (1982), p.258
        de = R./60;
        der = NaN(size(de));
    end
end

%%
function de = get_bending_saemundsson (e)
    % Sæmundsson, Þorsteinn (1986). "Astronomical Refraction". Sky and Telescope. 72: 70.
    % input is vacuum or unrefracted elevation angle.
    R = 1.02*cotd(e + 10.3./(e + 5.11));
    de = R./60;
end
