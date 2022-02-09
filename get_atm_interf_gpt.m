function [dt, da, dg, Ht, Ha, Hg, N, de, der] = get_atm_interf_gpt (e, H, pos, date, opt)
% GET_ATM_INTERF_GPT:  Closed-form interferometric atmospheric delay, using GPT atmospheric model.
% 
% SYNTAX:
%    [dt, da, dg, Ht, Ha, Hg, N, de, der] = get_atm_interf_gpt (e, H, pos, date);
%
% INPUT: 
%    e: [vector] satellite elevation angle (in degrees)
%    H: [vector] reflector depth or antenna height above reflector (in meters)
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
%    der: [vector] rate of change of elevation bending w.r.t. elevation angle (in degrees per degree)
%
% OPTIONAL INPUT: 
%    pos: [vector] position, in geodetic coordinates (latitude, longitude, altitude; in degrees, degrees, meters)
%    date: [vector] date (year month day, in this order)
%    opt: [struct] closed formula options
%    opt.ver: [char or scalar] GPT version ('1' or 1, '2' or 1, '2w')
%    opt.thin: [scalar, boolean] assume thin layer? defaults to false
%    opt.* (see get_atm_interf_met.m for other options)
% 
% EXAMPLE:
%    e = 45;  % degrees
%    H = 10;  % meters
%    [dt, ~, ~, Ht] = get_atm_interf_gpt (e, H)

    if (nargin < 3),  pos = [];  end
    if (nargin < 4),  date = [];  end
    if (nargin < 5),  opt = [];  end
    opt_default = struct('ver',1, 'thin',false);
    opt = structmergenonempty(opt_default, opt);
    
    [P, T, s, opt] = get_met (H, pos, date, opt);

    if (nargout < 4)
        [dt, da, dg] = get_atm_interf_met (e, H, P, T, s, [], opt);
    else
        [dt, da, dg, Ht, Ha, Hg, N, de, der] = get_atm_interf_met (e, H, P, T, s, [], opt);
    end
end

%%
function [P, T, s, opt] = get_met (H, pos, date, opt)
    if isempty(pos)
        P = [];  T = [];  s = [];
    elseif opt.thin
        [P, T, s] = get_met_aux (pos, date, opt.ver);
    else
        [Pa, Ta, sa] = get_met_aux (pos, date, opt.ver);
        [Ps, Ts, ss] = get_met_aux (pos, date, opt.ver, H);
        P = logavg(Pa, Ps);
        %P = (Pa+Ps)./2;
        T = (Ta+Ts)./2;
        s = (sa+ss)./2;
        opt.thin = false;  % don't redo altitude reduction in get_atm_interf_met
    end
end

%%
function [P, T, s] = get_met_aux (pos, date, ver, H)
    if (nargin < 4),  H = [];  end

    if isempty(pos),  pos = [0 0 0];  end
    pos = rowvec(pos);
    lat = pos(:,1);
    lon = pos(:,2);
    alt = pos(:,3);
    if ~isempty(H),  alt = alt - H;  end
    
    if isempty(date)
        temporal = false;
        mjd = 51544;  % date = [2000 01 01];
    else
        temporal = true;
        date = rowvec(date);
        try
            mjd = mjuliandate(date);
        catch err
            if strcmpi(err.identifier, 'MATLAB:UndefinedFunction') ...
            || strcmpi(err.identifier, 'MATLAB:ErrorRecovery:UnlicensedFunction')              
                mjd = mydatemjd (date);
            else
                rethrow(err);
            end
        end
    end

    [P, T, s] = get_met_aux2 (temporal, mjd, lat, lon, alt, ver);
end

%%
function [P, T, s] = get_met_aux2 (temporal, mjd, lat, lon, alt, ver)
    switch lower(ver)
    case {1,'1'}
        [p, t] = gpt (mjd, deg2rad(lat), deg2rad(lon), alt);  e = 0;
    case {2,'2'}
        file = 'gpt2_5.grd';
        url = 'https://vmf.geo.tuwien.ac.at/codes/';
        check_file(file, url);
        [p, t, ~, e] = gpt2 (mjd, deg2rad(lat), deg2rad(lon), alt, temporal);
    case {'2w','2w5'}
        file = 'gpt2_5w.grd';
        url = 'https://vmf.geo.tuwien.ac.at/codes/';
        check_file(file, url);
        [p, t, ~, e] = gpt2_5w (mjd, deg2rad(lat), deg2rad(lon), alt, numel(alt), temporal);
    otherwise
        error('matlab:get_atm_interf_gpt:badVer', 'Unknown version "%s".', char(ver));
    end
    P = p*100;  % from hPa to pascal
    T = t + 273.15;  % from Celsius to kelvin
    e = e*100;  % from hPa to pascal
    if (e == 0)
        s = 0;
    else
        s = convert_humidity (P, T, e, 'partial pressure', 'specific humidity');
    end
end

%%
function check_file (file, url)
    if ~exist(file, 'file')
        error('matlab:get_atm_interf_gpt:noFile', ...
            'Missing file "%s"; please download from <a href="%s">%s</a>', file, url, url);
    end
end
