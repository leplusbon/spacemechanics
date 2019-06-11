function ECI = ECEF2ECI(ECEF, D)

GMSTdeg = 100.46 + 360.985612288088 * D;
GMSTrad = GMSTdeg * pi / 180;

tf = [
    cos(GMSTrad),   -sin(GMSTrad), 	0;
    sin(GMSTrad),   cos(GMSTrad),   0;
    0,              0,              1
    ];

ECI = tf * ECEF;