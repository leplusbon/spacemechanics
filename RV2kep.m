function kep = RV2kep(r, v)

mu = 3.986004418e+5;
h = cross(r, v);
n = [-h(2); h(1); 0];
ev = cross(v, h) / mu - r / norm(r);
e = norm(ev);
p = norm(h)^2 / mu;

% normalized h, n
hn = h / norm(h);
nn = n / norm(n);

u = atan2(r' * cross(hn, nn), r' * nn);
i = acos(hn(3));
Omega = atan2(nn(2), nn(1));
omega = atan2(ev' * cross(hn, nn), ev' * nn);
a = p / (1 - e^2);

kep.a = a;
kep.e = e;
kep.Omega = Omega;
kep.i = i;
kep.omega = omega;
kep.u = u;