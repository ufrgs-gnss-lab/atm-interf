function [dt, da, dg, Ht, Ha, Hg, N, de] = get_atm_interf_gpt (e, H, pos, date, thin, ver)
% GET_ATM_INTERF_GPT:  Closed-form interferometric atmospheric delay, using GPT atmospheric model.
% 
% SYNTAX:
%    [dt, da, dg, Ht, Ha, Hg, N, de] = get_atm_interf_gpt (e, H, pos, date);
%
% INPUT: 
%    e: [vector] satellite elevation angle (in degrees)
%    H: [vector] reflector depth or antenna height above reflector (in meters)
%    pos: [vector] position, in geodetic coordinates (latitude, longitude, altitude; in degrees, degrees, meters)
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
%    date: [vector] date (year month day, in this order)
%    thin: [scalar] assume thin layer? (Boolean)
%    ver: [scalar] GPT version ('1', '2', '2w')
% 
% EXAMPLE:
%    e = 45;  % degrees
%    H = 10;  % meters
%    [dt, ~, ~, Ht] = get_atm_interf_pos (e, H)

    if (nargin < 3),  pos = [];  end
    if (nargin < 4),  date = [];  end
    if (nargin < 5),  thin = [];  end
    if (nargin < 6),  ver = [];  end
    if isempty(thin),  thin = false;  end
    if isempty(ver),  ver = 1;  end

    if isempty(pos)
        P = [];  T = [];  s = [];
    elseif thin
        [P, T, s] = get_met (pos, date, ver);
    else
        [Pa, Ta, sa] = get_met (pos, date, ver);
        [Ps, Ts, ss] = get_met (pos, date, ver, H);
        P = exp((log(Pa)+log(Ps))./2);
        %P = (Pa+Ps)./2;
        T = (Ta+Ts)./2;
        s = (sa+ss)./2;
        thin = false;  % don't redo altitude reduction.
    end    

    if (nargout < 4)
        [dt, da, dg] = get_atm_interf_met (e, H, P, T, s, [], thin);
    else
        [dt, da, dg, Ht, Ha, Hg] = get_atm_interf_met (e, H, P, T, s, [], thin);
    end
end

%%
function [P, T, s] = get_met (pos, date, ver, H)
    if (nargin < 4),  H = [];  end

    if isempty(date)
        temporal = false;
        date = [2000 01 01];
        mjd = 51544;
    else
        temporal = true;
        date = rowvec(date);
        try
            jd = mjuliandate(date);
        catch err
            if strcmp(err.identifier, 'MATLAB:ErrorRecovery:UnlicensedFunction')
                jd = mydatemjd (date);
            else
                rethrow(err)
            end
        end
    end

    pos = rowvec(pos);
    lat = pos(:,1);
    lon = pos(:,2);
    alt = pos(:,3);
    if ~isempty(H),  alt = alt - H;  end

    [P, T, s] = get_met_aux (mjd, lat, lon, alt, ver);
end

function [P, T, s] = get_met_aux (mjd, lat, lon, alt, ver)
    switch lower(ver)
    case {1,'1','v1'}
        [p, t] = gpt (mjd, deg2rad(lat), deg2rad(lon), alt);  e = 0;
    case {2,'2','v2'}
        [p, t, ~, e] = gpt2 (mjd, deg2rad(lat), deg2rad(lon), alt, temporal);
    case {'2w','v2w'}
        [p, t, ~, e] = gpt2_1w (mjd, deg2rad(lat), deg2rad(lon), alt, numel(alt), temporal);
    otherwise
        error('matlab:get_atm_interf_gpt:badVer', 'Unknown version "%s".', char(ver));
    end
    P = p*100;  % from hPa to pascal
    T = t - 273.15;  % from Celsius to kelvin
    e = e*100;  % from hPa to pascal
    if (e == 0)
        s = 0;
    else
        s = convert_humidity (P, T, e, 'partial pressure', 'specific humidity');
    end
end

