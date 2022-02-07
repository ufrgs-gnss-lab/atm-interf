function [P, T] = reduce_pressure (P0, T0, lapse_rate, height_geopot_rel)
    if isempty(lapse_rate),  lapse_rate = 6.5/1e3;  end
    n = length(P0);
    if isscalar(lapse_rate),  lapse_rate = repmat(lapse_rate, n, 1);  end
    myassert(length(T0), n)
    myassert(length(lapse_rate), n)
    myassert(length(height_geopot_rel), n)

    const = get_meteo_const();

    T = T0 + lapse_rate .* height_geopot_rel;
    temp = (T ./ T0) .^ (-const.g_c ./ (const.R_dry .* lapse_rate));
    P = P0 .* temp;
end

%!shared
%! % U.S. Standard Atmosphere
%! Z = ((80:-5:20)').*1e3;
%! P = [...
%!     8.8627e-3  % 80
%!     2.0679e-2  % 75
%!     4.6342e-2  % 70
%!     9.9220e-2  % 65
%!     2.0314e-1  % 60
%!     3.9969e-1  % 55
%!     7.5944e-1  % 50
%!     1.4313e+0  % 45
%!     2.7752e+0  % 40
%!     5.5892e+0  % 35
%!     1.1718e+1  % 30
%!     2.5110e+1  % 25
%!     5.4748e+1  % 20
%! ] .* 100;  % convert from mbar (=hPa) to Pa
%! T = [...
%!     196.650  % 80
%!     206.650  % 75
%!     217.450  % 70
%!     231.450  % 65
%!     245.450  % 60
%!     259.450  % 55
%!     270.650  % 50
%!     265.050  % 45
%!     251.050  % 40
%!     237.050  % 35
%!     226.650  % 30
%!     221.650  % 25
%!     216.650  % 20
%! ];
%! Z = flipud(Z);
%! P = flipud(P);
%! T = flipud(T);
%! %[Z, P, T]  % DEBUG

%!test
%! % Single point.
%! P0 = P(1);
%! T0 = T(1);
%! Z0 = Z(1);
%! P_true = P(2);
%! Z = Z(2);
%! T = T(2);
%! 
%! height_geopot_rel = Z - Z0;
%! lapse_rate = (T - T0) ./ height_geopot_rel;
%! 
%! P = reduce_pressure (P0, T0, lapse_rate, height_geopot_rel);
%! P2 = integrate_pressure (Z, T, Z0, T0, P0);
%! 
%! %[P_true, P, P2, 100*(P-P_true)./P_true, 100*(P-P2)./P_true]
%! myassert(P, P_true, 1/100)

%!test
%! % Several, independant, points.
%! P0 = P(1:end-1);
%! T0 = T(1:end-1);
%! Z0 = Z(1:end-1);
%! P_true = P(2:end);
%! Z = Z(2:end);
%! T = T(2:end);
%! 
%! height_geopot_rel = Z - Z0;
%! lapse_rate = (T - T0) ./ height_geopot_rel;
%! 
%! P = reduce_pressure (P0, T0, lapse_rate, height_geopot_rel);
%! P2 = zeros(size(lapse_rate));
%! for i=1:length(lapse_rate)
%!     P2(i) = integrate_pressure (Z(i), T(i), Z0(i), T0(i), P0(i));
%! end
%! 
%! %[P_true, P, P2, 100*(P-P_true)./P_true, 100*(P-P2)./P_true]
%! myassert(P, P_true, 1/100)

