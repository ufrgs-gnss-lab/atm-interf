function [N, N_hydro, N_nonhydro, N_dry, N_wet] = calculate_refractivity (P, T, s, ref)
% CALCULATE_REFRACTIVITY: Return refractivity given meteorological parameters.
% 
% INPUT:
%   P: pressure, Pa
%   T: temperature, Kelvin
% 
% OUTPUT:
%   N: total (micro-)refracivity, N=1e-6*(n-1), unitless and scaled by one millionth.
%   N_hydro: hydrostatic component
%   N_nonhydro: non-hydrostatic component
%   N_dry: dry component
%   N_wet: wet component
% 
% OPTIONAL INPUT:
%   s: specific humidity, kg/kg (therefore unitless)
%   ref: reference for refractivity coefficients (string)
%      'thayer1974', 'rueguer best average', 'rueguer best available', 'iugg1963', 'bevis1994'

    if (nargin < 3) || isempty(s),  s = 0;  end
    if (nargin < 4),  ref = [];  end
    [k_1, k_2, k_3] = get_refractivity_coeff (ref);
    % (References cited below are given at the bottom of this file.)

    const = get_meteo_const();
    M_dry = const.M_dry;
    M_wet = const.M_wet;
    R_dry = const.R_dry;
    R_wet = const.R_wet;

    % Compressibilities of (dry, moist) air.
    % Langley (1996, p. 118) says: "For typical conditions in 
    % the earth's atmosphere, Z_d and Z_w depart from unity 
    % by less that 1 part in 10^3." Therefore we assume that:
    inv_Z_d = 1;  % Inverse compressibility of dry air
    inv_Z_w = 1;  % Inverse compressibility of moist air

    %%%%%
    P_w = convert_humidity (P, T, s, 'specific humidity', ...
        'partial pressure of water vapor');
    P_d = P - P_w;  % Partial pressure of dry gases

    %%%%%
    % Langley (1996, p. 118)
    N_dry = k_1 * (P_d ./ T) * inv_Z_d; % dry component
    N_wet = ( k_2 * (P_w ./ T) + k_3 * (P_w ./ (T.^2)) ) * inv_Z_w; % wet comp.
    N_total_drywet = N_dry + N_wet;

    %%%%%
    % (Davis et al., 1985, eq. A6):
    k_2_prime = k_2 - (M_wet / M_dry) * k_1;

    % We need the density of the mixed gas to calculate hydrostatic 
    % refractivity:
    density = calculate_density_mixed_gas (T, P, P_w);
    
    % (Davis et al., 1985, first term in eq. A7):
    N_hydro = k_1 * R_dry * density;
    %N_hydro = k_1 .* (P ./ T);  % WRONG!!!
    % (Davis et al., 1985, second and third terms in eq. A7; also eq. A15):
    N_nonhydro = k_2_prime .* (P_w ./ T) + k_3 .* (P_w ./ T.^2);

    N_total_hydrononhydro = N_hydro + N_nonhydro;

    %%%%%
    %myassert(N_total_davis, N_total_thayer, -sqrt(N_total_thayer))  % DEBUG
    N = N_total_drywet;  % Refractivity
    
    %%%%%        
    % References
    % 
    % Langley, R.B. (1996) Propagation of the GPS signals.
    % In Kleusberg; Teunissen (Eds.) GPS for Geodesy 
    % (Lecture Notes in Earth Sciences). Springer. Chap. 3, 
    % pp.103-140.
    % 
    % Matveev, L.T. (1967) Fundamentals of General Meteorology,
    % Physics of the Atmosphere (in Russian, translated by IPST).
    % Dept. of Commerce and NSF, 699 pp.
    %
    % Rocken, C.; Sokolovskiy, S.; Johnson, J.; Hunt, D. (2001) 
    % Improved Mapping of Tropospheric Delays. Journal of 
    % Atmospheric and Oceanic Technology, Vol. 18. July. 
    % pp. 1205-1213.
    %
    % Thayer, G.D. (1974) An improved equation for the radio
    % refractive index of air. Radio Science, Vol. 9, No. 10,
    % pp. 803-807.
    %
    % AMS Glossary
    % <http://amsglossary.allenpress.com/>
    %
    % Davis et al. (1985) Geodesy by radio interferometry: Effects of atmospheric modeling errors on estimates of baseline length. Radio Science 20, 1593-1607.
end

% refractivity coefficients.
function [k_1, k_2, k_3] = get_refractivity_coeff (ref)
    if (nargin < 1) || isempty(ref),  ref = 'thayer1974';  end
    switch lower(ref)
    case 'thayer1974'
        k_1 = 77.60;    % +/- 0.014 [K/mbar]
        k_2 = 64.8;     % +/- 0.08  [K/mbar]
        k_3 = 3.776e5;  % +/- 0.004e5  [K^2/mbar]
    case {'rueguer best average','rueguer2002'}
        k_1 = 77.6890;
        k_2 = 71.2952;
        k_3 = 375463;
    case {'rueguer best available'}
        k_1 = 77.695;
        k_2 = 71.97;
        k_3 = 375406;
    case 'iugg1963'
        k_1 = 77.624;
        k_2 = 64.700;
        k_3 = 371897;
    case 'bevis1994'
        k_1 = 77.60;
        k_2 = 70.4;
        k_3 = 3.739e5;
    otherwise
        error('MATLB:get_refractivity_coeff:badRef', ...
            'Unknown refractivity coefficients reference "%s".', ...
            char(coeff_ref));
    end
    % convert mbar (equivalent to hPa) to Pa, in refractivity constants:
    k_1 = k_1 / 100;
    k_2 = k_2 / 100;
    k_3 = k_3 / 100;
end

%!shared
%! const = get_meteo_const();

%!test
%! % The results below I obtained before modifying calculate_refractivity() 
%! % on Mar 6, 2008. At that time, that routine had been used for two to 
%! % three years.
%! N_hydro = 267.82;
%! N_nonhydro = 776.24;
%! N_dry = 241.03;
%! N_wet = 803.02;
%! tol_abs = 1e-2;
%! 
%! [ignore, N_hydro2, N_nonhydro2, N_dry2, N_wet2] = calculate_refractivity (...
%!     const.std_temperature, const.std_pressure, 1/10);
%! myassert(N_hydro2, N_hydro, -tol_abs);
%! myassert(N_nonhydro2, N_nonhydro, -tol_abs);
%! myassert(N_dry2, N_dry, -tol_abs);
%! myassert(N_wet2, N_wet, -tol_abs);

%!test
%! [ignore, N_hydro, N_nonhydro, N_dry, N_wet] = calculate_refractivity (...
%!     const.std_temperature, const.std_pressure, 1/10);
%! N_total_drywet = N_dry + N_wet;
%! N_total_hydrononhydro = N_hydro + N_nonhydro;
%! tol_abs = eps(N_total_drywet);
%! myassert(N_total_drywet, N_total_hydrononhydro, -tol_abs);

%!test
%! % zero humidity, zero N_wet and N_nonhydro.
%! [ignore, N_hydro, N_nonhydro, N_dry, N_wet] = calculate_refractivity (...
%!     const.std_temperature, const.std_pressure, 0);
%! tol_abs = eps;
%! myassert(N_wet, 0, -tol_abs);
%! myassert(N_nonhydro, 0, -tol_abs);

% any test more meaningful?

