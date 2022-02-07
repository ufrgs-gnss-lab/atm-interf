% T: temperature, in Kelvin
% P: pressure, Pa
% e: partial pressure of water vapor, Pa
function [density1, density2] = calculate_density_mixed_gas (T, P, e)
    const = get_meteo_const();
    c = const.M_wet / const.M_dry;  % auxiliary constant.
    
    % It can be obtained in at least two equivalent ways. 
     
    % Following Saastamoinen (1972, p. 248), 
    %      "The density of the mixture is, of course, equal to 
    %      [the sum of the densities each of the dry-air and water-vapor 
    %      components]":
    % Seeber (2003, p.56) (his P' is our P_d):
    P_d = P - e;  % Partial pressure of dry gases
    density_dry = P_d ./ (const.R_dry .* T);
    density_wet = e   ./ (const.R_wet .* T);
    density1 = density_dry + density_wet;

    if (nargout < 2),  return;  end

    % Following "Gas constant" in AMS Glossary, 
    %     "For moist air, the variable percentage of water vapor is taken 
    %     into account by retaining the gas constant for dry air while using 
    %     the virtual temperature in place of the temperature." 
    T_v = convert_humidity (P, T, e, ...
        'partial pressure of water vapor', 'virtual temperature');
    density2 = (P ./ T_v) ./ const.R_dry;  % state eq. of mixed gases

    % (For references, please see calculate_refractivity.)
end

%!test
%! const = get_meteo_const();
%! T = const.std_temperature;
%! P = const.std_pressure;
%! e = P/10;
%! 
%! [density1, density2] = calculate_density_mixed_gas (T, P, e);
%! % The two formulations above are equivalent:
%! %density1 - density2  % DEBUG
%! 
%! myassert(density1, density2, -eps(density1))

%!test
%! const = get_meteo_const();
%! T = const.std_temperature;
%! P = const.std_pressure;
%! e = P/2;
%! 
%! [density1, density2] = calculate_density_mixed_gas (T, P, e);
%! % The two formulations above are equivalent:
%! %density1 - density2  % DEBUG
%! 
%! myassert(density1, density2, -eps(density1))

