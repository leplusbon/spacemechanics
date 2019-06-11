%% convertDat
clear; clc; close all;

load('satdata.mat');
satdata = convertDat(satdata);
satdataS1 = convertDat(satdataS1);
mu = 3.986004418e+5;
tm = satdata.D;
tmS1 = satdataS1.D;
sizedat = size(satdata, 1);
sizeS1 = size(satdataS1, 1);
rdat = [satdata.ECI_Xkm,	satdata.ECI_Ykm,	satdata.ECI_Zkm]';
rECEFdat = [satdata.ECEF_Xkm,	satdata.ECEF_Ykm,	satdata.ECEF_Zkm]';
rECEFdatS1 = [satdataS1.ECEF_Xkm, satdataS1.ECEF_Ykm, satdataS1.ECEF_Zkm]';
vdat = [satdata.ECI_V_Xkms,	satdata.ECI_V_Ykms, satdata.ECI_V_Zkms]';
rdatS1 = [satdataS1.ECI_Xkm,	satdataS1.ECI_Ykm,	  satdataS1.ECI_Zkm]';
vdatS1 = [satdataS1.ECI_V_Xkms, satdataS1.ECI_V_Ykms, satdataS1.ECI_V_Zkms]';

%% 1 : r, v -> kep
% 1(a)
tm_1a = satdata.D(1);
ro_1a = [satdata.ECI_Xkm(1);    satdata.ECI_Ykm(1);    satdata.ECI_Zkm(1)];
vo_1a = [satdata.ECI_V_Xkms(1); satdata.ECI_V_Ykms(1); satdata.ECI_V_Zkms(1)];
kep_1a = RV2kep(ro_1a, vo_1a);
meanmot_1a = sqrt(mu / kep_1a.a^3) * 3600 * 24; % rad / day

%1(b)
r_1b = [satdata.ECI_Xkm,	satdata.ECI_Ykm,	satdata.ECI_Zkm]';
v_1b = [satdata.ECI_V_Xkms,	satdata.ECI_V_Ykms, satdata.ECI_V_Zkms]';
a = zeros(size(r_1b, 2), 1);
e = zeros(size(r_1b, 2), 1);
Omega = zeros(size(r_1b, 2), 1);
i = zeros(size(r_1b, 2), 1);
omega = zeros(size(r_1b, 2), 1);
u = zeros(size(r_1b, 2), 1);
for it = 1:size(r_1b, 2)
    r = r_1b(:, it);
    v = v_1b(:, it);
    kep = RV2kep(r, v);
    a(it) = kep.a;
    e(it) = kep.e;
    Omega(it) = kep.Omega;
    i(it) = kep.i;
    omega(it) = kep.omega;
    u(it) = kep.u;
end
tm_1b = satdata.D;
kep_1b = table(a, e, Omega, i, omega, u);
kep_mean_1b.a = mean(kep_1b.a);
kep_mean_1b.e = mean(kep_1b.e);
kep_mean_1b.Omega = wrapToPi(mean(unwrap(kep_1b.Omega)));
kep_mean_1b.i = wrapToPi(mean(unwrap(kep_1b.i)));
kep_mean_1b.omega = wrapToPi(mean(unwrap(kep_1b.omega)));
meanmot_1b = sqrt(mu / kep_mean_1b.a^3) * 3600 * 24; % rad / day


% 1(c)
r_re_1a = zeros(3, sizedat);
r_ECEF_re_1a = zeros(3, sizedat);
latlonH_re_1a = zeros(3, sizedat);
v_re_1a = zeros(3, sizedat);
r_re_1b = zeros(3, sizedat);
r_ECEF_re_1b = zeros(3, sizedat);
latlonH_re_1b = zeros(3, sizedat);
v_re_1b = zeros(3, sizedat);
latlonH_dat = zeros(3, sizedat);

for it = 1:sizedat
    latlonH_dat(:, it) = ECEF2latlonH(rECEFdat(:, it));
end

for it = 1:sizedat
    kep = kep_1a;
    kep.u = kep.u + meanmot_1a * (tm(it) - tm_1a);
    [r_re_1a(:, it), v_re_1a(:, it)] = kep2ECI(kep);
    r_ECEF_re_1a(:, it) = ECI2ECEF(r_re_1a(:, it), tm(it));
    latlonH_re_1a(:, it) = ECEF2latlonH(r_ECEF_re_1a(:, it));
end
med = round(sizedat / 2);
for it = 1:sizedat
    kep = kep_mean_1b;
    kep.u = kep_1b.u(med) + meanmot_1b * (tm(it) - tm_1b(med));
    [r_re_1b(:, it), v_re_1b(:, it)] = kep2ECI(kep);
    r_ECEF_re_1b(:, it) = ECI2ECEF(r_re_1b(:, it), tm(it));
    latlonH_re_1b(:, it) = ECEF2latlonH(r_ECEF_re_1b(:, it));
end
latlonH_dat(1:2, :) = unwrap(latlonH_dat(1:2, :), [], 2);
latlonH_re_1a(1:2, :) = unwrap(latlonH_re_1a(1:2, :), [], 2);
latlonH_re_1b(1:2, :) = unwrap(latlonH_re_1b(1:2, :), [], 2);


%% 1 : groundtrack

tm_x = tm(1):((tm(end) - tm(1)) / 15):tm(end);
lonlat_dat_x = spline(tm, latlonH_dat([2, 1], :), tm_x);
lonlat_re_x_1a = spline(tm, latlonH_re_1a([2, 1], :), tm_x);
lonlat_re_x_1b = spline(tm, latlonH_re_1b([2, 1], :), tm_x);

textcomment = datestr(datenum(2000, 1, 1, 0, 0, 0) + tm_x', 'HH:MM:SS AM');
dx = 3; dy = 0;

figure(100);
% worldmap([-90, 90], [-180, 180]);
% lakes = shaperead('worldlakes','UseGeoCoords',true);
h = zeros(6, 1);
load coastlines;
%geoshow(coastlat, coastlon); 
hold on;
ter = plot([coastlon; NaN; coastlon+360; NaN; coastlon-360], ...
    [coastlat; NaN; coastlat; NaN; coastlat], 'b');
ter.Color = [0.5, 0.5, 0.5];
axis equal;
grid on;
h(1:3) = plot( ...
    latlonH_dat(2, :) * 180 / pi, latlonH_dat(1, :) * 180 / pi, 'k', ...
    latlonH_re_1a(2, :) * 180 / pi, latlonH_re_1a(1, :) * 180 / pi, 'r', ...
    latlonH_re_1b(2, :) * 180 / pi, latlonH_re_1b(1, :) * 180 / pi, 'b' ...
    );
hold on;
hp = plot(lonlat_dat_x(1, :) * 180 / pi, lonlat_dat_x(2, :) * 180 / pi, 'k+');
hq = plot(lonlat_re_x_1a(1, :) * 180 / pi, lonlat_re_x_1a(2, :) * 180 / pi, 'r+');
hr = plot(lonlat_re_x_1b(1, :) * 180 / pi, lonlat_re_x_1b(2, :) * 180 / pi, 'b+');
hp.MarkerSize = 10;
hq.MarkerSize = 10;
hr.MarkerSize = 10;
text(lonlat_dat_x(1, :)' * 180 / pi + dx, lonlat_dat_x(2, :)' * 180 / pi + dy, textcomment);
legend(h(1:3), 'Measurement', 'Regeneration (1a)', 'Regeneration (1b)');
ylim([-90, 90]);
xlim([-260, 160]);
title('SNUGLITE Ground Track, April 11, 2019');
xlabel('Longitude [\circ]');
ylabel('Latitude [\circ]');

%% 1 : 3dplot ECI
figure(101);
hold on;
h = zeros(1, 3);
h(1) = plot3(rdat(1, :), rdat(2, :), rdat(3, :), 'k');
h(2) = plot3(r_re_1a(1, :), r_re_1a(2, :), r_re_1a(3, :), 'r');
h(3) = plot3(r_re_1b(1, :), r_re_1b(2, :), r_re_1b(3, :), 'b');
[ex, ey, ez] = sphere;
ex = ex * 6371; ey = ey * 6371; ez = ez * 6371;
sh = surf(ex, ey, ez);
sh.FaceColor = [1, 1, 1];
sh.EdgeColor = [0.7, 0.7, 0.7];
axis equal;
legend(h(1:3), 'Measurement', 'Regeneration (1a)', 'Regeneration (1b)');
title('SNUGLITE Orbit Visualization, April 11, 2019');
xlabel('X_{ECI} [km]');
ylabel('Y_{ECI} [km]');
zlabel('Z_{ECI} [km]');

%% 1 : 3dplot ECEF
figure(102);
hold on;
h = zeros(1, 3);
h(1) = plot3(rECEFdat(1, :), rECEFdat(2, :), rECEFdat(3, :), 'k');
h(2) = plot3(r_ECEF_re_1a(1, :), r_ECEF_re_1a(2, :), r_ECEF_re_1a(3, :), 'r');
h(3) = plot3(r_ECEF_re_1b(1, :), r_ECEF_re_1b(2, :), r_ECEF_re_1b(3, :), 'b');
he = plot3( ...
    6371 * cosd(coastlon).*cosd(coastlat), ...
    6371 * sind(coastlon).*cosd(coastlat), 6371 * sind(coastlat));
he.Color = [0.5, 0.5, 0.5];
[ex, ey, ez] = sphere;
ex = ex * 6371; ey = ey * 6371; ez = ez * 6371;
sh = surf(ex, ey, ez);
sh.FaceColor = [1, 1, 1];
sh.EdgeColor = [0.8, 0.8, 0.8];
axis equal;
legend(h(1:3), 'Measurement', 'Regeneration (1a)', 'Regeneration (1b)');
title('SNUGLITE Orbit Visualization, April 11, 2019');
xlabel('X_{ECEF} [km]');
ylabel('Y_{ECEF} [km]');
zlabel('Z_{ECEF} [km]');

%% 2
% (a)
rv1 = [satdata.ECI_Xkm(1), satdata.ECI_Ykm(1), satdata.ECI_Zkm(1)];
rv2 = [satdata.ECI_Xkm(2), satdata.ECI_Ykm(2), satdata.ECI_Zkm(2)];
rv3 = [satdata.ECI_Xkm(3), satdata.ECI_Ykm(3), satdata.ECI_Zkm(3)];
kep_2a = gibbs(rv1, rv2, rv3);
tm_2a = satdata.D(2);

% (b)
kep_2b_raw = [];
for starti = 1:3:598
    dat = satdata(starti:(starti+2), :);
    rv = [dat.ECI_Xkm, dat.ECI_Ykm, dat.ECI_Zkm];
    rv1 = rv(1, :);
    rv2 = rv(2, :);
    rv3 = rv(3, :);
    k = gibbs(rv1, rv2, rv3);
    kep_2b_raw = [kep_2b_raw; [k.a, k.e, k.Omega, k.i, k.omega, k.u]];
end
kep_2b_raw(:, 6) = unwrap(kep_2b_raw(:, 6));

tm_2b = satdata.D(2:3:599);
kep_2b.a = mean(kep_2b_raw(:, 1));
kep_2b.e = mean(kep_2b_raw(:, 2));
kep_2b.Omega = mean(unwrap(kep_2b_raw(:, 3)));
kep_2b.i = mean(unwrap(kep_2b_raw(:, 4)));
kep_2b.omega = mean(unwrap(kep_2b_raw(:, 5)));

% (c)
kep_2c_raw = [];
for starti = 1:30:571
    dat = satdata(starti:10:(starti+20), :);
    rv = [dat.ECI_Xkm, dat.ECI_Ykm, dat.ECI_Zkm];
    rv1 = rv(1, :);
    rv2 = rv(2, :);
    rv3 = rv(3, :);
    k = gibbs(rv1, rv2, rv3);
    kep_2c_raw = [kep_2c_raw; [k.a, k.e, k.Omega, k.i, k.omega, k.u]];
end
kep_2c_raw(:, 6) = unwrap(kep_2c_raw(:, 6));

tm_2c = satdata.D(11:30:581);
kep_2c.a = mean(kep_2c_raw(:, 1));
kep_2c.e = mean(kep_2c_raw(:, 2));
kep_2c.Omega = mean(unwrap(kep_2c_raw(:, 3)));
kep_2c.i = mean(unwrap(kep_2c_raw(:, 4)));
kep_2c.omega = mean(unwrap(kep_2c_raw(:, 5)));

%% 2d : regen
tm_re = satdata.D;

% re_2a
meanmot_2a = sqrt(mu / kep_2a.a^3) * 3600 * 24;
r_re_2a = zeros(3, sizedat);
r_ECEF_re_2a = zeros(3, sizedat);
latlonH_re_2a = zeros(3, sizedat);
v_re_2a = zeros(3, sizedat);
for it = 1:sizedat
    kep = kep_2a;
    kep.u = kep_2a.u + meanmot_2a * (tm_re(it) - tm_2a);
    [r_re_2a(:, it), v_re_2a(:, it)] = kep2ECI(kep);
    r_ECEF_re_2a(:, it) = ECI2ECEF(r_re_2a(:, it), tm_re(it));
    latlonH_re_2a(:, it) = ECEF2latlonH(r_ECEF_re_2a(:, it));
end
latlonH_re_2a(1, :) = unwrap(latlonH_re_2a(1, :));
latlonH_re_2a(2, :) = unwrap(latlonH_re_2a(2, :));

% re_2b
meanmot_2b = sqrt(mu / kep_2b.a^3) * 3600 * 24;
med = round(size(tm_2b, 1) / 2);
tm_mean_2b = tm_2b(med);
kep_2b.u = kep_2b_raw(med, 6);

r_re_2b = zeros(3, sizedat);
r_ECEF_re_2b = zeros(3, sizedat);
latlonH_re_2b = zeros(3, sizedat);
v_re_2b = zeros(3, sizedat);
for it = 1:sizedat
    kep = kep_2b;
    kep.u = kep_2b.u + meanmot_2b * (tm_re(it) - tm_mean_2b);
    [r_re_2b(:, it), v_re_2b(:, it)] = kep2ECI(kep);
    r_ECEF_re_2b(:, it) = ECI2ECEF(r_re_2b(:, it), tm_re(it));
    latlonH_re_2b(:, it) = ECEF2latlonH(r_ECEF_re_2b(:, it));
end
latlonH_re_2b(1, :) = unwrap(latlonH_re_2b(1, :));
latlonH_re_2b(2, :) = unwrap(latlonH_re_2b(2, :));

% re_2c
meanmot_2c = sqrt(mu / kep_2c.a^3) * 3600 * 24;
med = round(size(tm_2c, 1) / 2);
tm_mean_2c = tm_2c(med);
kep_2c.u = kep_2c_raw(med, 6);

r_re_2c = zeros(3, sizedat);
r_ECEF_re_2c = zeros(3, sizedat);
latlonH_re_2c = zeros(3, sizedat);
v_re_2c = zeros(3, sizedat);
for it = 1:sizedat
    kep = kep_2c;
    kep.u = kep_2c.u + meanmot_2c * (tm_re(it) - tm_mean_2c);
    [r_re_2c(:, it), v_re_2c(:, it)] = kep2ECI(kep);
    r_ECEF_re_2c(:, it) = ECI2ECEF(r_re_2c(:, it), tm_re(it));
    latlonH_re_2c(:, it) = ECEF2latlonH(r_ECEF_re_2c(:, it));
end
latlonH_re_2c(1, :) = unwrap(latlonH_re_2c(1, :));
latlonH_re_2c(2, :) = unwrap(latlonH_re_2c(2, :));

%% 2d : groundtrack

tm_x = tm(1):((tm(end) - tm(1)) / 15):tm(end);
lonlat_dat_x = spline(tm, latlonH_dat([2, 1], :), tm_x);
lonlat_re_x_2a = spline(tm, latlonH_re_2a([2, 1], :), tm_x);
lonlat_re_x_2b = spline(tm, latlonH_re_2b([2, 1], :), tm_x);
lonlat_re_x_2c = spline(tm, latlonH_re_2c([2, 1], :), tm_x);

textcomment = datestr(datenum(2000, 1, 1, 0, 0, 0) + tm_x', 'HH:MM:SS AM');
dx = 3; dy = 0;

figure(200);
% worldmap([-90, 90], [-180, 180]);
% lakes = shaperead('worldlakes','UseGeoCoords',true);
h = zeros(6, 1);
load coastlines;
%geoshow(coastlat, coastlon); 
hold on;
ter = plot([coastlon; NaN; coastlon+360; NaN; coastlon-360], ...
    [coastlat; NaN; coastlat; NaN; coastlat], 'b');
ter.Color = [0.5, 0.5, 0.5];
axis equal;
grid on;
h(1:4) = plot( ...
    latlonH_dat(2, :) * 180 / pi, latlonH_dat(1, :) * 180 / pi, 'k', ...
    latlonH_re_2a(2, :) * 180 / pi, latlonH_re_2a(1, :) * 180 / pi, 'r', ...
    latlonH_re_2b(2, :) * 180 / pi, latlonH_re_2b(1, :) * 180 / pi, 'b', ...
    latlonH_re_2c(2, :) * 180 / pi, latlonH_re_2c(1, :) * 180 / pi, 'm' ...
    );
hold on;
hp = plot(lonlat_dat_x(1, :) * 180 / pi, lonlat_dat_x(2, :) * 180 / pi, 'k+');
hq = plot(lonlat_re_x_2a(1, :) * 180 / pi, lonlat_re_x_2a(2, :) * 180 / pi, 'r+');
hr = plot(lonlat_re_x_2b(1, :) * 180 / pi, lonlat_re_x_2b(2, :) * 180 / pi, 'b+');
hs = plot(lonlat_re_x_2c(1, :) * 180 / pi, lonlat_re_x_2c(2, :) * 180 / pi, 'm+');
hp.MarkerSize = 10;
hq.MarkerSize = 10;
hr.MarkerSize = 10;
hs.MarkerSize = 10;
text(lonlat_dat_x(1, :)' * 180 / pi + dx, lonlat_dat_x(2, :)' * 180 / pi + dy, textcomment);
legend(h(1:4), 'Measurement', 'Regeneration (2a)', 'Regeneration (2b)', 'Regeneration (2c)');
ylim([-90, 90]);
xlim([-260, 160]);
title('SNUGLITE Ground Track, April 11, 2019');
xlabel('Longitude [\circ]');
ylabel('Latitude [\circ]');


%% 1 : 3dplot ECI
figure(201);
hold on;
h = zeros(1, 4);
h(1) = plot3(rdat(1, :), rdat(2, :), rdat(3, :), 'k');
h(2) = plot3(r_re_2a(1, :), r_re_2a(2, :), r_re_2a(3, :), 'r');
h(3) = plot3(r_re_2b(1, :), r_re_2b(2, :), r_re_2b(3, :), 'b');
h(4) = plot3(r_re_2c(1, :), r_re_2c(2, :), r_re_2c(3, :), 'm');
[ex, ey, ez] = sphere;
ex = ex * 6371; ey = ey * 6371; ez = ez * 6371;
sh = surf(ex, ey, ez);
sh.FaceColor = [1, 1, 1];
sh.EdgeColor = [0.7, 0.7, 0.7];
axis equal;
legend(h(1:4), 'Measurement', 'Regeneration (2a)', 'Regeneration (2b)', 'Regeneration (2c)');
title('SNUGLITE Orbit Visualization, April 11, 2019');
xlabel('X_{ECI} [km]');
ylabel('Y_{ECI} [km]');
zlabel('Z_{ECI} [km]');

%% 1 : 3dplot ECEF
figure(202);
hold on;
h = zeros(1, 4);
h(1) = plot3(rECEFdat(1, :), rECEFdat(2, :), rECEFdat(3, :), 'k');
h(2) = plot3(r_ECEF_re_2a(1, :), r_ECEF_re_2a(2, :), r_ECEF_re_2a(3, :), 'r');
h(3) = plot3(r_ECEF_re_2b(1, :), r_ECEF_re_2b(2, :), r_ECEF_re_2b(3, :), 'b');
h(4) = plot3(r_ECEF_re_2c(1, :), r_ECEF_re_2c(2, :), r_ECEF_re_2c(3, :), 'm');
he = plot3( ...
    6371 * cosd(coastlon).*cosd(coastlat), ...
    6371 * sind(coastlon).*cosd(coastlat), 6371 * sind(coastlat));
he.Color = [0.5, 0.5, 0.5];
[ex, ey, ez] = sphere;
ex = ex * 6371; ey = ey * 6371; ez = ez * 6371;
sh = surf(ex, ey, ez);
sh.FaceColor = [1, 1, 1];
sh.EdgeColor = [0.8, 0.8, 0.8];
axis equal;
legend(h(1:4), 'Measurement', 'Regeneration (2a)', 'Regeneration (2b)', 'Regeneration (2c)');
title('SNUGLITE Orbit Visualization, April 11, 2019');
xlabel('X_{ECEF} [km]');
ylabel('Y_{ECEF} [km]');
zlabel('Z_{ECEF} [km]');


%% 4


tm_4 = (datenum(2019, 4, 12, 0, 51, 29) + (0:(2 * 3600 * 24)) / (3600 * 24))' - datenum(2000, 1, 1, 0, 0, 0);
size_4 = size(tm_4, 1);
r_re_1b_4 = zeros(3, size_4);
r_ECEF_re_1b_4 = zeros(3, size_4);
latlonH_re_1b_4 = zeros(3, size_4);
v_re_1b_4 = zeros(3, size_4);
r_ENU_re_1b_4 = zeros(3, size_4);
rElAz_4 = zeros(3, size_4);

med = round(sizedat / 2);
for it = 1:size_4
    kep = kep_mean_1b;
    kep.u = kep_1b.u(med) + meanmot_1b * (tm_4(it) - tm_1b(med));
    [r_re_1b_4(:, it), v_re_1b_4(:, it)] = kep2ECI(kep);
    r_ECEF_re_1b_4(:, it) = ECI2ECEF(r_re_1b_4(:, it), tm_4(it));
    latlonH_re_1b_4(:, it) = ECEF2latlonH(r_ECEF_re_1b_4(:, it));
    [r_ENU_re_1b_4(:, it), rElAz_4(:, it)] = ECEF2ENU(r_ECEF_re_1b_4(:, it), ...
        37.448347027777770, 126.9524995277778, 277.6);
end

% skyplot
figure(401);
skyplot_4 = polarplot(rElAz_4(3, :), wrapToPi(rElAz_4(2, :)) * 180 / pi, 'k', ...
    rElAz_4(3, 1:240:end), wrapToPi(rElAz_4(2, 1:240:end)) * 180 / pi, 'r+'...
    );
textr = [];
texttheta = [];
textcomment = [];
for it = 1:240:size_4
    if rElAz_4(2, it) > 0
        textr = [textr; rElAz_4(2, it) * 180 / pi];
        texttheta = [texttheta; rElAz_4(3, it)];
        textcomment = [textcomment; datestr(datenum(2000, 1, 1, 0, 0, 0) + tm_4(it), 'mm/dd HH:MM:SS AM')];
    end
end
rlim([0, 90]);
text(texttheta, textr, textcomment);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.RDir = 'reverse';
title('SNUGLITE Skyplot, April 12-14, 2019');
% ax.RAxis.Label.String = 'El [\circ]';
% ax.ThetaAxis.Label.String = 'Az [\circ]';

%%

r_re_1bS1 = zeros(3, sizeS1);
r_ECEF_re_1bS1 = zeros(3, sizeS1);
latlonH_re_1bS1 = zeros(3, sizeS1);
v_re_1bS1 = zeros(3, sizeS1);
latlonH_datS1 = zeros(3, sizeS1);

for it = 1:sizeS1
    latlonH_datS1(:, it) = ECEF2latlonH(rECEFdatS1(:, it));
end

med = round(sizedat / 2);
for it = 1:sizeS1
    kep = kep_mean_1b;
    kep.u = kep_1b.u(med) + meanmot_1b * (tmS1(it) - tm_1b(med));
    [r_re_1bS1(:, it), v_re_1bS1(:, it)] = kep2ECI(kep);
    r_ECEF_re_1bS1(:, it) = ECI2ECEF(r_re_1bS1(:, it), tmS1(it));
    latlonH_re_1bS1(:, it) = ECEF2latlonH(r_ECEF_re_1bS1(:, it));
end
latlonH_dat(1:2, :) = unwrap(latlonH_dat(1:2, :), [], 2);
latlonH_re_1b(1:2, :) = unwrap(latlonH_re_1b(1:2, :), [], 2);

%% 4 : groundtrack

tm_xS1 = tmS1(1):((tmS1(end) - tmS1(1)) / 8):tmS1(end);
lonlat_dat_xS1 = spline(tmS1, latlonH_datS1([2, 1], :), tm_xS1);
lonlat_re_x_1bS1 = spline(tmS1, latlonH_re_1bS1([2, 1], :), tm_xS1);

textcomment = datestr(datenum(2000, 1, 1, 0, 0, 0) + tm_xS1', 'HH:MM:SS AM');
dx = 3; dy = 0;

figure(400);
% worldmap([-90, 90], [-180, 180]);
% lakes = shaperead('worldlakes','UseGeoCoords',true);
h = zeros(2, 1);
load coastlines;
%geoshow(coastlat, coastlon); 
hold on;
ter = plot([coastlon; NaN; coastlon+360; NaN; coastlon-360], ...
    [coastlat; NaN; coastlat; NaN; coastlat], 'b');
ter.Color = [0.5, 0.5, 0.5];
axis equal;
grid on;
h(1:2) = plot( ...
    latlonH_datS1(2, :) * 180 / pi, latlonH_datS1(1, :) * 180 / pi, 'k', ...
    latlonH_re_1bS1(2, :) * 180 / pi, latlonH_re_1bS1(1, :) * 180 / pi, 'b' ...
    );
hold on;
hp = plot(lonlat_dat_xS1(1, :) * 180 / pi, lonlat_dat_xS1(2, :) * 180 / pi, 'k+');
hq = plot(lonlat_re_x_1bS1(1, :) * 180 / pi, lonlat_re_x_1bS1(2, :) * 180 / pi, 'b+');
hp.MarkerSize = 10;
hq.MarkerSize = 10;
text(lonlat_dat_xS1(1, :)' * 180 / pi + dx, lonlat_dat_xS1(2, :)' * 180 / pi + dy, textcomment);
tx2 = text(lonlat_re_x_1bS1(1, :)' * 180 / pi + dx, lonlat_re_x_1bS1(2, :)' * 180 / pi + dy, textcomment);
set(tx2, 'Color', [0, 0, 1]);
legend(h(1:2), 'Measurement', 'Prediction (1b)');
ylim([-40, 90]);
xlim([100, 180]);
title('SNUGLITE Ground Track, April 12, 2019');
xlabel('Longitude [\circ]');
ylabel('Latitude [\circ]');

%% 5

%% 6

J2 = 1082.635854e-06;
r_re_1b_6 = zeros(3, sizedat);
r_ECEF_re_1b_6 = zeros(3, sizedat);
latlonH_re_1b_6 = zeros(3, sizedat);
v_re_1b_6 = zeros(3, sizedat);
med = round(sizedat / 2);
for it = 1:sizedat
    kep = kep_mean_1b;
    p = kep.a * (1 - (kep.e)^2);
    meanmot_6 = meanmot_1b - ...
        0.75 * J2 * (6378.1363 / p)^2 * sqrt(1 - (kep.e)^2) * meanmot_1b ...
        * (3 * (sin(kep.i))^2 - 2);
    kep = kep_mean_1b;
    kep.u = kep_1b.u(med) + meanmot_6 * (tm(it) - tm_1b(med));
    kep.Omega = kep_mean_1b.Omega - 1.5 * J2 * (6378.1363 / p)^2 ...
        * meanmot_1b * cos(kep.i) * (tm(it) - tm(med));
    kep.omega = kep_mean_1b.omega + 1.5 * J2 * (6378.1363 / p)^2 ...
        * meanmot_1b * (2 - 2.5 * (sin(kep.i))^2) * (tm(it) - tm(med));
    [r_re_1b_6(:, it), v_re_1b_6(:, it)] = kep2ECI(kep);
    r_ECEF_re_1b_6(:, it) = ECI2ECEF(r_re_1b_6(:, it), tm(it));
    latlonH_re_1b_6(:, it) = ECEF2latlonH(r_ECEF_re_1b_6(:, it));
end
latlonH_re_1b_6(1:2, :) = unwrap(latlonH_re_1b_6(1:2, :), [], 2);

%% 6 : groundtrack

tm_x = tm(1):((tm(end) - tm(1)) / 15):tm(end);
lonlat_dat_x = spline(tm, latlonH_dat([2, 1], :), tm_x);
lonlat_re_x_1b_6 = spline(tm, latlonH_re_1b_6([2, 1], :), tm_x);

textcomment = datestr(datenum(2000, 1, 1, 0, 0, 0) + tm_x', 'HH:MM:SS AM');
dx = 3; dy = 0;

figure(600);
% worldmap([-90, 90], [-180, 180]);
% lakes = shaperead('worldlakes','UseGeoCoords',true);
h = zeros(6, 1);
load coastlines;
%geoshow(coastlat, coastlon); 
hold on;
ter = plot([coastlon; NaN; coastlon+360; NaN; coastlon-360], ...
    [coastlat; NaN; coastlat; NaN; coastlat], 'b');
ter.Color = [0.5, 0.5, 0.5];
axis equal;
grid on;
h(1:2) = plot( ...
    latlonH_dat(2, :) * 180 / pi, latlonH_dat(1, :) * 180 / pi, 'k', ...
    latlonH_re_1b_6(2, :) * 180 / pi, latlonH_re_1b_6(1, :) * 180 / pi, 'b' ...
    );
hold on;
hp = plot(lonlat_dat_x(1, :) * 180 / pi, lonlat_dat_x(2, :) * 180 / pi, 'k+');
hr = plot(lonlat_re_x_1b_6(1, :) * 180 / pi, lonlat_re_x_1b_6(2, :) * 180 / pi, 'b+');
hp.MarkerSize = 10;
hr.MarkerSize = 10;
text(lonlat_dat_x(1, :)' * 180 / pi + dx, lonlat_dat_x(2, :)' * 180 / pi + dy, textcomment);
legend(h(1:2), 'Measurement', 'Regeneration (6)');
ylim([-90, 90]);
xlim([-260, 160]);
title('SNUGLITE Ground Track, April 11, 2019');
xlabel('Longitude [\circ]');
ylabel('Latitude [\circ]');


%%
% 
% %% 2a : regen
% 
% med = round(size(tm_2a / 2));
% u_re_2a = kep_2a_raw(med, 6) + meanmot_2a * (tm_2a - tm_2a(med));
% 
% tm_re_2a = satdata.D;
% 
% r_ECEF_re_2a = zeros(3, size(tm_re_2a, 1));
% latlonH_re_2a = zeros(3, size(tm_re_2a, 1));
% med = round(sizedat / 2);
% meanmot_2a = sqrt(mu / a^3) * 3600 * 24;
% for it = 1:sizedat
%     kep = kep_2a;
%     kep.u = tm_re_2a;
%     [r_re_1b(:, it), v_re_1b(:, it)] = kep2ECI(kep);
%     r_ECEF_re_1b(:, it) = ECI2ECEF(r_re_1b(:, it), tm(it));
%     latlonH_re_1b(:, it) = ECEF2latlonH(r_ECEF_re_1b(:, it));
% end
% 
% for it = 1:size(tm_re, 1)
%     kk = {};
%     kk.a = a;
%     kk.e = e;
%     kk.Omega = Omega;
%     kk.i = i;
%     kk.omega = omega;
%     kk.u = u_re(it);
%     [ECI_r_rekm, ECI_v_rekms] = kep2ECI(kk);
%     ECEF_r_rekm = ECI2ECEF(ECI_r_rekm, tm_re(it));
%     ECEF_v_rekms = ECI2ECEF(ECI_v_rekms, tm_re(it));
%     ECEF_X_rekm(it) = ECEF_r_rekm(1);
%     ECEF_Y_rekm(it) = ECEF_r_rekm(2);
%     ECEF_Z_rekm(it) = ECEF_r_rekm(3);
%     latlonH_re = ECEF2latlonH(ECEF_r_rekm);
%     lat_redeg(it) = latlonH_re(1) * 180 / pi;
%     lon_redeg(it) = latlonH_re(2) * 180 / pi;
%     H_rekm(it) = latlonH_re(3);
% end
% 
% kep = [];
% for starti = 1:3:118
%     dat = satdata(starti:(starti+2), :);
%     rv = [dat.ECI_Xkm, dat.ECI_Ykm, dat.ECI_Zkm];
%     rv1 = rv(1, :);
%     rv2 = rv(2, :);
%     rv3 = rv(3, :);
%     kS1 = gibbs(rv1, rv2, rv3);
%     kep = [kep; [kS1.a, kS1.e, kS1.Omega, kS1.i, kS1.omega, kS1.u]];
% end
% 
% tmS1 = satdata.D(2:3:119);
% aS1 = mean(kep(:, 1));
% eS1 = mean(kep(:, 2));
% OmegaS1 = mean(kep(:, 3));
% iS1 = mean(kep(:, 4));
% omegaS1 = mean(kep(:, 5));
% uS1 = wrapToPi(kep(:, 6));
% 
% latS1deg = zeros(size(satdata.D));
% lonS1deg = zeros(size(satdata.D));
% HS1deg = zeros(size(satdata.D));
% for it = 1:size(satdata.D, 1)
%     latlonHS1 = ECEF2latlonH([satdata.ECEF_Xkm(it); satdata.ECEF_Ykm(it); satdata.ECEF_Zkm(it)]);
%     latS1deg(it) = latlonHS1(1) * 180 / pi;
%     lonS1deg(it) = latlonHS1(2) * 180 / pi;
%     HS1km(it) = latlonHS1(3);
% end
% 
% %% GIBBS : Plot
% 
% figure(1);
% subplot(3, 1, 1);
% plot(satdataS1.D, satdataS1.ECEF_Xkm, 'k', tm_re, ECEF_X_rekm, 'r');
% grid on;
% ylabel('X [km] (ECEF)');
% legend('Measurement', 'Reconstructed');
% subplot(3, 1, 2);
% plot(satdataS1.D, satdataS1.ECEF_Ykm, 'k', tm_re, ECEF_Y_rekm, 'r');
% grid on;
% ylabel('Y [km] (ECEF)');
% subplot(3, 1, 3);
% plot(satdataS1.D, satdataS1.ECEF_Zkm, 'k', tm_re, ECEF_Z_rekm, 'r');
% grid on;
% ylabel('Z [km] (ECEF)');
% xlabel('Time [days from 0h, 01 Jan 2000, UT1]');
% 
% %% GIBBS : groundTrack
% lonlatS1x = spline(satdataS1.D, [lonS1deg, latS1deg]', (7041.078:0.002:7041.09)');
% lonlat_rex = spline(tm_re, [lon_redeg, lat_redeg]', (7041.078:0.002:7041.09)');
% textcomment = datestr(datenum(2000, 1, 1, 0, 0, 0) + (7041.078:0.002:7041.09)', 'HH:MM:SS AM');
% dx = 3; dy = 0;
% 
% figure(2);
% % worldmap([-90, 90], [-180, 180]);
% % lakes = shaperead('worldlakes','UseGeoCoords',true);
% h = zeros(5, 1);
% load coastlines;
% h(5) = geoshow(coastlat, coastlon); 
% grid on;
% hold on;
% h(1:2) = plot(lonS1deg, latS1deg, 'k', lon_redeg, lat_redeg, 'r');
% h(3:4) = plot(lonlatS1x(1, :), lonlatS1x(2, :), 'k+', ...
%     lonlat_rex(1, :), lonlat_rex(2, :), 'r+');
% text(lonlatS1x(1, :)' + dx, lonlatS1x(2, :)' + dy, textcomment);
% legend(h(1:2), 'Measurement', 'Reconstructed');
% ylim([0, 90]);
% xlim([120, 180]);
% title('SNUGLITE Ground Track, April 12, 2019');
% xlabel('Longitude [\circ]');
% ylabel('Latitude [\circ]');
% 
% 
