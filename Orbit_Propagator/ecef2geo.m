function GEO = ecef2geo(posvececef, degreeflag) 
%------------------------------------------------------------------------------- 
% function [latitude, longitude, altitude] = ecef2geo(posvececef, degreeflag) 
% 
% references: 
% 
% bowring b., 1976. transformation from spatial to geographical coordinates, 
% survey review 23, pp. 323-327. 
% 
% bowring, b.r. (1985) the accuracy of geodetic latitude and height equations, 
% survey review, vol. 28, no. 218, pp. 202-206. 
% 
% inputs description 
% 
% posvececef earth-centered earth-fixed position vector (m). 
% degreeflag specifies whether output latitude and longitude are in 
% degrees. 
% 
% outputs description 
% 
% latitude geodetic latitude, positive north (deg or rad). 
% longitude geodetic longitude, measured east from greenwich meridian (deg 
% or rad). 
% altitude geodetic altitude above wgs-84 ellipsoid (m). 
%------------------------------------------------------------------------------- 
% wgs-84 defining parameters. 
%------------------------------------------------------------------------------- 
a = 6378137.0; 
f = 1.0 / 298.257223563;

%------------------------------------------------------------------------------- 
% wgs-84 derived parameters. 
%------------------------------------------------------------------------------- 
one_f = 1.0 - f;

b = a * one_f; % semi-minor axis 
e2 = f * (2.0 - f); % first eccentricity squared 
epsilon = e2 / (1.0 - e2); % second eccentricity squared 
b_a = one_f;

%------------------------------------------------------------------------------- 
% extract ecef components from input position vector. 
%------------------------------------------------------------------------------- 
x = posvececef(:, 1); 
y = posvececef(:, 2); 
z = posvececef(:, 3);

%------------------------------------------------------------------------------- 
% initialize outputs. 
%------------------------------------------------------------------------------- 
latitude = zeros(size(x)); 
longitude = latitude; 
altitude = latitude;

%------------------------------------------------------------------------------- 
% quick check for all components zero. 
%------------------------------------------------------------------------------- 
ii0 = (x == 0 & y == 0 & z == 0);

if any(ii0), 
latitude(ii0) = 0; 
longitude(ii0) = 0; 
altitude(ii0) = 0; 
end

%------------------------------------------------------------------------------- 
% quick calculations at poles. 
%------------------------------------------------------------------------------- 
ii1 = (x == 0 & y == 0 & z ~= 0);

if any(ii1), 
latitude(ii1) = sign(z(ii1)) * pi / 2; 
longitude(ii1) = 0; 
altitude(ii1) = abs(z(ii1)) - b; 
end

%------------------------------------------------------------------------------- 
% quick calculations at equator. 
%------------------------------------------------------------------------------- 
ii2 = (~ii0 & ~ii1 & z == 0.0);

if any(ii2), 
longitude(ii2) = atan2(y(ii2), x(ii2)); 
latitude(ii2) = 0; 
p = sqrt(x(ii2).^2 + y(ii2).^2); 
altitude(ii2) = p - a; 
end

%------------------------------------------------------------------------------- 
% main algorithm. in bowring (1985), u is the parametric latitude. it is crucial 
% to maintain the appropriate signs for the sin(u) and sin(lat) in the equations 
% below. 
%------------------------------------------------------------------------------- 
ii = ~ii0 & ~ii1 & ~ii2;

if any(ii),

p2 = x(ii).^2 + y(ii).^2; 
r2 = p2 + z(ii).^2; 
p = sqrt(p2); 
r = sqrt(r2);

%------------------------------------------------------------------------------- 
% equation (17) from bowring (1985), shown to improve numerical accuracy in lat 
%------------------------------------------------------------------------------- 
tanu = b_a * (z(ii) ./ p) .* (1 + epsilon * b ./ r); 
tan2u = tanu .* tanu;

%------------------------------------------------------------------------------- 
% avoid trigonometric functions for determining cos3u and sin3u 
%------------------------------------------------------------------------------- 
cos2u = 1.0 ./ (1.0 + tan2u); 
cosu = sqrt(cos2u); 
cos3u = cos2u .* cosu;

sinu = tanu .* cosu; 
sin2u = 1.0 - cos2u; 
sin3u = sin2u .* sinu;

%------------------------------------------------------------------------------- 
% equation (18) from bowring (1985) 
%------------------------------------------------------------------------------- 
tanlat = (z(ii) + epsilon * b * sin3u) ./ (p - e2 * a * cos3u);

tan2lat = tanlat .* tanlat; 
cos2lat = 1.0 ./ (1.0 + tan2lat); 
sin2lat = 1.0 - cos2lat;

coslat = sqrt(cos2lat); 
sinlat = tanlat .* coslat;

longitude(ii) = atan2(y(ii), x(ii)); 
latitude(ii) = atan(tanlat);

%------------------------------------------------------------------------------- 
% equation (7) from bowring (1985), shown to be numerically superior to other 
% height equations. note that equation (7) from bowring (1985) writes the last 
% term as a^2 / nu, but this reduces to a * sqrt(1 - e^2 * sin(lat)^2), because 
% nu = a / sqrt(1 - e^2 * sin(lat)^2). 
%------------------------------------------------------------------------------- 
altitude(ii) = p .* coslat + z(ii) .* sinlat - a * sqrt(1.0 - e2 * sin2lat);

end

% longitude = unwrap(longitude);

%------------------------------------------------------------------------------- 
% convert outputs if necessary. 
%------------------------------------------------------------------------------- 
if nargin == 2 & degreeflag == 1, 
radtodeg = 180 / pi; 
latitude = latitude * radtodeg; 
longitude = longitude * radtodeg; 
end

GEO = [latitude, longitude, altitude];