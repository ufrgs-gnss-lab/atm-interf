function [p,T,dT,e,ah,aw,undu] = gpt2_artem (dmjd,dlat,dlon,hell,it)

% This subroutine determines pressure, temperature, temperature lapse rate, 
% water vapour pressure, hydrostatic and wet mapping function coefficients 
% ah and aw, and geoid undulation for specific sites near the Earth 
% surface. It is based on a 5 x 5 degree external grid file ('gpt2_5.grd') 
% with mean values as well as sine and cosine amplitudes for the annual and 
% semiannual variation of the coefficients.
%
% c Reference:
% c Lagler, K., Schindelegger, M., Böhm, J., Krásná, H., Nilsson, T. (2013), 
% c GPT2: Empirical slant delay model for radio space geodetic techniques, 
% c Geophys. Res. Lett., Vol. 40, pp. 1069–1073, DOI: 10.1002/grl.50288.
%
% input parameters:
%
% dmjd:  modified Julian date (scalar, only one epoch per call is possible)
% dlat:  ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
% dlon:  longitude in radians [-pi:pi] or [0:2pi] (vector)
% hell:  ellipsoidal height in m (vector)
% nstat: number of stations in dlat, dlon, and hell
%        maximum possible: not relevant for Matlab version
% it:    case 1: no time variation but static quantities
%        case 0: with time variation (annual and semiannual terms)
% 
% output parameters:
%
% p:    pressure in hPa (vector of length nstat) 
% T:    temperature in degrees Celsius (vector of length nstat)
% dT:   temperature lapse rate in degrees per km (vector of length nstat) 
% e:    water vapour pressure in hPa (vector of length nstat)
% ah:   hydrostatic mapping function coefficient at zero height (VMF1) 
%       (vector of length nstat)
% aw:   wet mapping function coefficient (VMF1) (vector of length nstat)
% undu: geoid undulation in m (vector of length nstat)
%
% The hydrostatic mapping function coefficients have to be used with the
% height dependent Vienna Mapping Function 1 (vmf_ht.f) because the
% coefficients refer to zero height.
%
% Example 1 (Vienna, 2 August 2012, with time variation):
%
% dmjd = 56141.d0
% dlat(1) = 48.20d0*pi/180.d0
% dlon(1) = 16.37d0*pi/180.d0
% hell(1) = 156.d0
% nstat = 1
% it = 0
%
% output:
% p = 1002.56 hPa
% T = 22.12 deg Celsius
% dT = -6.53 deg / km
% e = 15.63 hPa
% ah = 0.0012647
% aw = 0.0005726
% undu = 44.06 m
%
% Example 2 (Vienna, 2 August 2012, without time variation, i.e. constant values):
%
% dmjd = 56141.d0
% dlat(1) = 48.20d0*pi/180.d0
% dlon(1) = 16.37d0*pi/180.d0
% hell(1) = 156.d0
% nstat = 1
% it = 1
%
% output:
% p = 1003.49 hPa
% T = 11.95 deg Celsius
% dT = -5.47 deg / km
% e = 9.58 hPa
% ah = 0.0012395
% aw = 0.0005560
% undu = 44.06 m
%
% Klemens Lagler, 2 August 2012
% Johannes Boehm, 6 August 2012, revision
% Klemens Lagler, 21 August 2012, epoch change to January 1 2000
% Johannes Boehm, 23 August 2012, adding possibility to determine constant field
% Johannes Boehm, 27 December 2012, reference added
% Johannes Boehm, 10 January 2013, correction for dlat = -90 degrees
%                                  (problem found by Changyong He)
% Johannes Boehm, 21 May 2013, bug with dmjd removed (input parameter dmjd was replaced
%                 unintentionally; problem found by Dennis Ferguson)
% Artem Karaush (Siberian Research Institute of Metrology, Novosibirsk, Russia) 17 September 2015
%                 Function GPT2 is now more effective and "Matlab-like style": 
%                 It avoids using loops where possible and avoids unneccecary memory allocations, and  
%                 nstat is removed from input parameters, it is calculated
%                 from length of input vectors.
% ---

nstat = length( dlat );

% change the reference epoch to January 1 2000
dmjd1 = dmjd-51544.5;

% mean gravity in m/s**2
gm = 9.80665;
% molar mass of dry air in kg/mol
dMtr = 28.965*10^-3;
% universal gas constant in J/K/mol
Rg = 8.3143;

% factors for amplitudes
if (it==1) % then  constant parameters
    cosfy = 0;
    coshy = 0;
    sinfy = 0;
    sinhy = 0;
else 
    cosfy = cos(dmjd1/365.25*2*pi);
    coshy = cos(dmjd1/365.25*4*pi);
    sinfy = sin(dmjd1/365.25*2*pi);
    sinhy = sin(dmjd1/365.25*4*pi);
end

% read gridfile
fid = fopen( 'gpt2_5.grd','r' );

C = textscan( fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', ...
    'HeaderLines', 1, ...
    'CollectOutput', true );

fclose (fid);

pgrid  = C{1}(:, 3:7);         % pressure in Pascal
Tgrid  = C{1}(:, 8:12);        % temperature in Kelvin
Qgrid  = C{1}(:, 13:17)/1000;  % specific humidity in kg/kg
dTgrid = C{1}(:, 18:22)/1000;  % temperature lapse rate in Kelvin/m
u      = C{1}(:, 23);          % geoid undulation in m
Hs     = C{1}(:, 24);          % orthometric grid height in m
ahgrid = C{1}(:, 25:29)/1000;  % hydrostatic mapping function coefficient, dimensionless
awgrid = C{1}(:, 30:34)/1000;  % wet mapping function coefficient, dimensionless

p    = zeros( nstat, 1 );      % pressure in hPa
T    = zeros( nstat, 1 );      % temperature in degrees Celsius
dT   = zeros( nstat, 1 );      % temperature lapse rate in degrees per km
e    = zeros( nstat, 1 );      % water vapour pressure in hPa
ah   = zeros( nstat, 1 );      % hydrostatic mapping function coefficient at zero height (VMF1) 
aw   = zeros( nstat, 1 );      % wet mapping function coefficient (VMF1)
undu = zeros( nstat, 1 );      % geoid undulation in m


% loop over stations
for k = 1:nstat
    
    % only positive longitude in degrees
    if dlon(k) < 0
        plon = (dlon(k) + 2*pi)*180/pi;
    else
        plon = dlon(k)*180/pi;
    end
    % transform to polar distance in degrees
    ppod = (-dlat(k) + pi/2)*180/pi; 

    % find the index (line in the grid file) of the nearest point
    ipod = floor((ppod+5)/5); 
    ilon = floor((plon+5)/5);
    
    % normalized (to one) differences, can be positive or negative
    diffpod = (ppod - (ipod*5 - 2.5))/5;
    difflon = (plon - (ilon*5 - 2.5))/5;
    % added by HCY
    if ipod == 37
        ipod = 36;
    end

    % get the number of the corresponding line
    indx(1) = (ipod - 1)*72 + ilon;
    
    % near the poles: nearest neighbour interpolation, otherwise: bilinear
    bilinear = 0;
    if ppod > 2.5 && ppod < 177.5 
           bilinear = 1;          
    end          
    
    % case of nearest neighbourhood
    if bilinear == 0

        ix = indx(1);
        
        % transforming ellipsoidial height to orthometric height
        undu(k) = u(ix);
        hgt = hell(k)-undu(k);
            
        % pressure, temperature at the heigtht of the grid
        T0 = Tgrid(ix,1) + ...
             Tgrid(ix,2)*cosfy + Tgrid(ix,3)*sinfy + ...
             Tgrid(ix,4)*coshy + Tgrid(ix,5)*sinhy;
        p0 = pgrid(ix,1) + ...
             pgrid(ix,2)*cosfy + pgrid(ix,3)*sinfy+ ...
             pgrid(ix,4)*coshy + pgrid(ix,5)*sinhy;
         
        % specific humidity
        Q = Qgrid(ix,1) + ...
            Qgrid(ix,2)*cosfy + Qgrid(ix,3)*sinfy+ ...
            Qgrid(ix,4)*coshy + Qgrid(ix,5)*sinhy;
            
        % lapse rate of the temperature
        dT(k) = dTgrid(ix,1) + ...
                dTgrid(ix,2)*cosfy + dTgrid(ix,3)*sinfy+ ...
                dTgrid(ix,4)*coshy + dTgrid(ix,5)*sinhy; 

        % station height - grid height
        redh = hgt - Hs(ix);

        % temperature at station height in Celsius
        T(k) = T0 + dT(k)*redh - 273.15;
        
        % temperature lapse rate in degrees / km
        dT(k) = dT(k)*1000;

        % virtual temperature in Kelvin
        Tv = T0*(1+0.6077*Q);
        
        c = gm*dMtr/(Rg*Tv);
        
        % pressure in hPa
        p(k) = (p0*exp(-c*redh))/100;
        
        % water vapour pressure in hPa
        e(k) = (Q*p(k))/(0.622+0.378*Q);
            
        % hydrostatic coefficient ah 
        ah(k) = ahgrid(ix,1) + ...
                ahgrid(ix,2)*cosfy + ahgrid(ix,3)*sinfy+ ...
                ahgrid(ix,4)*coshy + ahgrid(ix,5)*sinhy;
            
        % wet coefficient aw
        aw(k) = awgrid(ix,1) + ...
                awgrid(ix,2)*cosfy + awgrid(ix,3)*sinfy + ...
                awgrid(ix,4)*coshy + awgrid(ix,5)*sinhy;           
                    
     else % bilinear interpolation
        
        ipod1 = ipod + sign(diffpod);
        ilon1 = ilon + sign(difflon);
        if ilon1 == 73
            ilon1 = 1;
        end
        if ilon1 == 0
            ilon1 = 72;
        end
        
        % get the number of the line
        indx(2) = (ipod1 - 1)*72 + ilon;  % along same longitude
        indx(3) = (ipod  - 1)*72 + ilon1; % along same polar distance
        indx(4) = (ipod1 - 1)*72 + ilon1; % diagonal
        
        % transforming ellipsoidial height to orthometric height:
        % Hortho = -N + Hell
        undul = u(indx);
        hgt = hell(k)-undul;

        % pressure, temperature at the heigtht of the grid
        T0 = Tgrid(indx,1) + ...
             Tgrid(indx,2)*cosfy + Tgrid(indx,3)*sinfy + ...
             Tgrid(indx,4)*coshy + Tgrid(indx,5)*sinhy;
        p0 = pgrid(indx,1) + ...
             pgrid(indx,2)*cosfy + pgrid(indx,3)*sinfy + ...
             pgrid(indx,4)*coshy + pgrid(indx,5)*sinhy;

        % humidity 
        Ql = Qgrid(indx,1) + ...
             Qgrid(indx,2)*cosfy + Qgrid(indx,3)*sinfy + ...
             Qgrid(indx,4)*coshy + Qgrid(indx,5)*sinhy;

        % reduction = stationheight - gridheight
        Hs1 = Hs(indx);
        redh = hgt - Hs1;

        % lapse rate of the temperature in degree / m
        dTl = dTgrid(indx,1) + ...
              dTgrid(indx,2)*cosfy + dTgrid(indx,3)*sinfy + ...
              dTgrid(indx,4)*coshy + dTgrid(indx,5)*sinhy; 

        % temperature reduction to station height
        Tl = T0 + dTl.*redh - 273.15;

        % virtual temperature
        Tv = T0.*(1+0.6077*Ql);  
        c = gm*dMtr./(Rg*Tv);

        % pressure in hPa
        pl = (p0.*exp(-c.*redh))/100;

        % hydrostatic coefficient ah
        ahl = ahgrid(indx,1) + ...
              ahgrid(indx,2)*cosfy + ahgrid(indx,3)*sinfy + ...
              ahgrid(indx,4)*coshy + ahgrid(indx,5)*sinhy;

        % wet coefficient aw
        awl = awgrid(indx,1) + ...
              awgrid(indx,2)*cosfy + awgrid(indx,3)*sinfy + ...
              awgrid(indx,4)*coshy + awgrid(indx,5)*sinhy;
            

        dnpod1 = abs(diffpod); % distance nearer point
        dnpod2 = 1 - dnpod1;   % distance to distant point
        dnlon1 = abs(difflon);
        dnlon2 = 1 - dnlon1;
        
        % pressure
        R1 = dnpod2*pl(1)+dnpod1*pl(2);
        R2 = dnpod2*pl(3)+dnpod1*pl(4);
        p(k) = dnlon2*R1+dnlon1*R2;
            
        % temperature
        R1 = dnpod2*Tl(1)+dnpod1*Tl(2);
        R2 = dnpod2*Tl(3)+dnpod1*Tl(4);
        T(k) = dnlon2*R1+dnlon1*R2;
        
        % temperature in degree per km
        R1 = dnpod2*dTl(1)+dnpod1*dTl(2);
        R2 = dnpod2*dTl(3)+dnpod1*dTl(4);
        dT(k) = (dnlon2*R1+dnlon1*R2)*1000;
            
        % humidity
        R1 = dnpod2*Ql(1)+dnpod1*Ql(2);
        R2 = dnpod2*Ql(3)+dnpod1*Ql(4);
        Q = dnlon2*R1+dnlon1*R2;
        e(k) = (Q*p(k))/(0.622+0.378*Q);
            
        % hydrostatic
        R1 = dnpod2*ahl(1)+dnpod1*ahl(2);
        R2 = dnpod2*ahl(3)+dnpod1*ahl(4);
        ah(k) = dnlon2*R1+dnlon1*R2;
           
        % wet
        R1 = dnpod2*awl(1)+dnpod1*awl(2);
        R2 = dnpod2*awl(3)+dnpod1*awl(4);
        aw(k) = dnlon2*R1+dnlon1*R2;
        
        % undulation
        R1 = dnpod2*undul(1)+dnpod1*undul(2);
        R2 = dnpod2*undul(3)+dnpod1*undul(4);
        undu(k) = dnlon2*R1+dnlon1*R2;
                    
    end 
end

  
