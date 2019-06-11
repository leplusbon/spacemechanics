function latlonH = ECEF2latlonH(ECEF)

ea = 6378.1370;
ef = 1 / 298.257223563;
ee = sqrt(ef * (2 - ef));

x = ECEF(1);    y = ECEF(2);    z = ECEF(3);
p = norm([x, y]);

kappaprev = 0;
kappa = 1;
while abs(kappa - kappaprev) > 1e-13
    kappaprev = kappa;
    c = sqrt((p^2 + (1 - ee^2) * z^2 * kappa^2)^3) / ea / ee^2;
    kappa = 1 + (p^2 + (1 - ee^2) * z^2 * kappa^3) / (c - p^2);
end
lat = atan(kappa * z / p);
lon = atan2(y, x);
rhoP = ea / sqrt(1 - ee^2 * sin(lat)^2);
H = p / cos(lat) - rhoP;
latlonH = [lat; lon; H];