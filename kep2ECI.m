function [rv, v] = kep2ECI(kep)

mu = 3.986004418e+5;
p = kep.a * (1 - (kep.e)^2);
nu = kep.u - kep.omega;
r = p / (1 + kep.e * cos(nu));
rPQR = [r * cos(nu); r * sin(nu); 0];
vPQR = sqrt(mu / p) * [-sin(nu); kep.e + cos(nu); 0];

R1 = [
    cos(kep.Omega), -sin(kep.Omega),    0;
    sin(kep.Omega), cos(kep.Omega),     0;
    0,              0,                  1
    ];
R2 = [
    1,  0,          0;
    0,  cos(kep.i), -sin(kep.i);
    0,  sin(kep.i), cos(kep.i)
    ];
R3 = [
	cos(kep.omega), -sin(kep.omega),    0;
    sin(kep.omega), cos(kep.omega),     0;
    0,              0,                  1
    ];

rv = R1 * R2 * R3 * rPQR;
v = R1 * R2 * R3 * vPQR;