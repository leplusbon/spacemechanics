function kep = gibbs(rv1, rv2, rv3)

r1 = norm(rv1);
r2 = norm(rv2);
r3 = norm(rv3);

N = r3 * cross(rv1, rv2) + r1 * cross(rv2, rv3) + r2 * cross(rv3, rv1);
D = cross(rv1, rv2) + cross(rv2, rv3) + cross(rv3, rv1);
S = rv1 * (r2 - r3) + rv2 * (r3 - r1) + rv3 * (r1 - r2);
eQ = S / norm(S);
eW = N / norm(N);
eP = cross(eQ, eW);
en = [-N(2), N(1), 0];
en = en / norm(en);
p = norm(N) / norm(D);
em = cross(eW, en);

kep = {};
kep.e = norm(S) / norm(D);
kep.a = p / (1 - (kep.e)^2);
kep.i = acos(eW(3));
kep.omega = atan2(-dot(en, eQ), dot(en, eP));
kep.Omega = acos(en(1));
kep.u = atan2(dot(rv2, em), dot(rv2, en));