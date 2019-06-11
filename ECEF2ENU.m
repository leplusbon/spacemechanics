function [rENU, rElAz] = ECEF2ENU(rECEF, latdeg, londeg, Hm)

lat = latdeg * pi / 180;
lon = londeg * pi / 180;
H = Hm * 0.001;
ea = 6378.1370;
ef = 1 / 298.257223563;
ee = sqrt(ef * (2 - ef));
rhoP = ea / sqrt(1 - ee^2 * sin(lat)^2);
originECEF = [
    (rhoP + H) * cos(lat) * cos(lon);
    (rhoP + H) * cos(lat) * sin(lon);
    (rhoP * (1 - ee^2) + H) * sin(lat)
    ];

drECEF = rECEF - originECEF;
ang1 = pi / 2 - lat;
ang2 = pi / 2 + lon;
R1 = [
    1,  0,          0;
    0,  cos(ang1),  sin(ang1);
    0,  -sin(ang1), cos(ang1)
    ];
R2 = [
    cos(ang2),  sin(ang2),  0;
    -sin(ang2), cos(ang2),  0;
    0,          0,          1
    ];
rENU = R1 * R2 * drECEF;
rElAz = [
    norm(rENU);
    atan2(rENU(3), norm([rENU(1), rENU(2)]));
    atan2(rENU(1), rENU(2))
    ];