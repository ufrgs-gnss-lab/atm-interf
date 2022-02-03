function de = get_bending_bennet (e, P, T)
% GET_BENDING_BENNET Elevation angle bending, via formulation of Bennet (1982).
% 
% INPUT:
%   e: elevation angle [degress]
%   P: pressure [Pascal] (optional)
%   T: temperature [Kelvin] (optional)
% 
% OUTPUT:
%   de: elevantion angle bending (refracted minus vacuum) [degrees]
% 
% Bennett, G. (1982). "The Calculation of Astronomical Refraction in Marine Navigation". Journal of Navigation, 35(2), 255-259. doi:10.1017/S0373463300022037

    if (nargin < 3),  P = [];  end
    if (nargin < 4),  T = [];  end

    Rm = cotd( e + (7.31 ./ (e + 4.4)) );  % Bennet (1982), p.257, formula G.
    %Rm = cot( e_rad+ (7.31 / (e_rad+4.4)) ); %WRONG!!

    if isempty(P) || isempty(T)
        de = Rm./60; %arcmin -> deg
        return
    end

    PaToMbar = @(p) p*100;
    MbarToInchesMercury = @(m) m*0.029530;
    P_InHg = MbarToInchesMercury(PaToMbar(P));

    KelvinToCelsius = @(k) k + 273.15;
    CelsiusToFahrenheit = @(c) c*9/5+32;
    T_Fahrn = CelsiusToFahrenheit(KelvinToCelsius(T));

    %R = (510 ./ (460 + T_Fahrn)) * (P_InHg./29.83) .*Rm;  % Bowditch apud. Bennet (1982), p.258
    R = (P_InHg - 80)./930).*(1./(1+8e-5.*(Rm+39).*(T_Fahrn-10))).*Rm;  % Bennet (1982), p.258
    de = R./60;
end

