%% convertDat
clear; clc; close all;

load('satdata.mat');
satdata = convertDat(satdata);
satdataS1 = convertDat(satdataS1);
mu = 3.986004418e+5;
J2 = 1082.635854e-06;
tm = satdata.D;
tmS1 = satdataS1.D;
sizedat = size(satdata, 1);
sizedatS1 = size(satdataS1, 1);
rdat = [satdata.ECI_Xkm,	satdata.ECI_Ykm,	satdata.ECI_Zkm]';
rECEFdat = [satdata.ECEF_Xkm,	satdata.ECEF_Ykm,	satdata.ECEF_Zkm]';
rECEFdatS1 = [satdataS1.ECEF_Xkm, satdataS1.ECEF_Ykm, satdataS1.ECEF_Zkm]';
vdat = [satdata.ECI_V_Xkms,	satdata.ECI_V_Ykms, satdata.ECI_V_Zkms]';
rdatS1 = [satdataS1.ECI_Xkm,	satdataS1.ECI_Ykm,	  satdataS1.ECI_Zkm]';
vdatS1 = [satdataS1.ECI_V_Xkms, satdataS1.ECI_V_Ykms, satdataS1.ECI_V_Zkms]';

kep_raw = [];
for starti = 1:3:598
    dat = satdata(starti:(starti+2), :);
    rv = [dat.ECI_Xkm, dat.ECI_Ykm, dat.ECI_Zkm];
    rv1 = rv(1, :);
    rv2 = rv(2, :);
    rv3 = rv(3, :);
    k = gibbs(rv1, rv2, rv3);
    kep_raw = [kep_raw; [k.a, k.e, k.Omega, k.i, k.omega, k.u]];
end
kep_raw(:, 6) = unwrap(kep_raw(:, 6));

tm_7 = satdata.D(2:3:599);
kep_mean.a = mean(kep_raw(:, 1));
kep_mean.e = mean(kep_raw(:, 2));
kep_mean.Omega = mean(kep_raw(:, 3));
kep_mean.i = mean(kep_raw(:, 4));
kep_mean.omega = mean(unwrap(kep_raw(:, 5)));
kep_mean.u = mean(kep_raw(:, 6));

med = round(sizedat / 2);
meanmot_0 = sqrt(mu / (kep_mean.a)^3) * 3600 * 24;
%meanmot_0 = 93.94103218;
%kep_mean.a = nthroot(mu / (meanmot_0 / 3600 / 24)^2, 3);
% at = [tm_7 - tm(med), ones(size(tm_7))];
% reab = (at' * at) \ (at' * kep_raw(:, 6));
% meanmot = reab(1);
% p = kep_mean.a * (1 - (kep_mean.e)^2);
% meanmot_0 = meanmot / ...
%         (1 - 0.75 * J2 * (6378.1363 / p)^2 * sqrt(1 - (kep_mean.e)^2) ...
%         * (3 * (sin(kep_mean.i))^2 - 2));
% 
% maxiter = 100;
% Re = 6378.1363;
% temp = 1 - (kep_mean.e)^2;
% meanmot_0 = meanmot;
% for it = 1:maxiter
%     a2 = (mu / (meanmot_0 / 3600 / 24)^2)^(2 / 3);
%     meanmot_0 = meanmot / ( ...
%         1 - 0.75 * J2 * Re^2 / (a2 * temp^2) * sqrt(temp) * (3 * sin(kep_mean.i)^2 - 2)...
%         );
% end
% kep_mean.a = sqrt(a2);

r_re = zeros(3, sizedat);
r_ECEF_re = zeros(3, sizedat);
r_ENU_re = zeros(3, sizedat);
rElAz_re = zeros(3, sizedat);
latlonH_re = zeros(3, sizedat);
v_re = zeros(3, sizedat);
latlonH_dat = zeros(3, sizedat);
for it = 1:sizedat
    latlonH_dat(:, it) = ECEF2latlonH(rECEFdat(:, it));
end
for it = 1:sizedat
    kep = kep_mean;
    p = kep_mean.a * (1 - (kep_mean.e)^2);
    meanmot = meanmot_0 - ...
        0.75 * J2 * (6378.1363 / p)^2 * sqrt(1 - (kep_mean.e)^2) * meanmot_0 ...
        * (3 * (sin(kep_mean.i))^2 - 2);
    kep.u = kep_mean.u + meanmot * (tm(it) - tm(med));
    kep.Omega = kep_mean.Omega - 1.5 * J2 * (6378.1363 / p)^2 ...
        * meanmot_0 * cos(kep_mean.i) * (tm(it) - tm(med));
    kep.omega = kep_mean.omega + 1.5 * J2 * (6378.1363 / p)^2 ...
        * meanmot_0 * (2 - 2.5 * (sin(kep_mean.i))^2) * (tm(it) - tm(med));
    [r_re(:, it), v_re(:, it)] = kep2ECI(kep);
    r_ECEF_re(:, it) = ECI2ECEF(r_re(:, it), tm(it));
    latlonH_re(:, it) = ECEF2latlonH(r_ECEF_re(:, it));
    [r_ENU_re(:, it), rElAz_re(:, it)] = ECEF2ENU(r_ECEF_re(:, it), ...
        37.448347027777770, 126.9524995277778, 277.6);
end
latlonH_dat(1:2, :) = latlonH_dat(1:2, 1) + unwrap(latlonH_dat(1:2, :) - latlonH_dat(1:2, 1), [], 2);
latlonH_re(1:2, :) = latlonH_dat(1:2, 1) + unwrap(latlonH_re(1:2, :) - latlonH_dat(1:2, 1), [], 2);

%% 7 : groundtrack
tm_x = tm(1):((tm(end) - tm(1)) / 15):tm(end);
lonlat_dat_x = spline(tm, latlonH_dat([2, 1], :), tm_x);
lonlat_re_x = spline(tm, latlonH_re([2, 1], :), tm_x);

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
h(1:2) = plot( ...
    latlonH_dat(2, :) * 180 / pi, latlonH_dat(1, :) * 180 / pi, 'k', ...
    latlonH_re(2, :) * 180 / pi, latlonH_re(1, :) * 180 / pi, 'r' ...
    );
hold on;
hp = plot(lonlat_dat_x(1, :) * 180 / pi, lonlat_dat_x(2, :) * 180 / pi, 'k+');
hq = plot(lonlat_re_x(1, :) * 180 / pi, lonlat_re_x(2, :) * 180 / pi, 'r+');
hp.MarkerSize = 10;
hq.MarkerSize = 10;
text(lonlat_dat_x(1, :)' * 180 / pi + dx, lonlat_dat_x(2, :)' * 180 / pi + dy, textcomment);
legend(h(1:2), 'Measurement', 'Regeneration (7)');
ylim([-90, 90]);
xlim([-260, 160]);
title('SNUGLITE Ground Track, April 11, 2019');
xlabel('Longitude [\circ]');
ylabel('Latitude [\circ]');


%% verify
%J2 = 0;
J2 = 1082.635854e-06;
r_reS1 = zeros(3, sizedatS1);
r_ECEF_reS1 = zeros(3, sizedatS1);
r_ENU_reS1 = zeros(3, sizedatS1);
rElAz_reS1 = zeros(3, sizedatS1);
r_ENU_datS1 = zeros(3, sizedatS1);
rElAz_datS1 = zeros(3, sizedatS1);
latlonH_reS1 = zeros(3, sizedatS1);
v_reS1 = zeros(3, sizedatS1);
latlonH_datS1 = zeros(3, sizedatS1);
for it = 1:sizedatS1
    latlonH_datS1(:, it) = ECEF2latlonH(rECEFdatS1(:, it));
end
medS1 = round(sizedatS1 / 2);
for it = 1:sizedatS1
    kep = kep_mean;
    p = kep_mean.a * (1 - (kep_mean.e)^2);
    meanmot = meanmot_0 - ...
        0.75 * J2 * (6378.1363 / p)^2 * sqrt(1 - (kep_mean.e)^2) * meanmot_0 ...
        * (3 * (sin(kep_mean.i))^2 - 2);
    kep.u = kep_mean.u + meanmot * (tmS1(it) - tm(med));
    kep.Omega = kep_mean.Omega - 1.5 * J2 * (6378.1363 / p)^2 ...
        * meanmot_0 * cos(kep_mean.i) * (tmS1(it) - tm(med));
    kep.omega = kep_mean.omega + 1.5 * J2 * (6378.1363 / p)^2 ...
        * meanmot_0 * (2 - 2.5 * (sin(kep_mean.i))^2) * (tmS1(it) - tm(med));
    [r_reS1(:, it), v_reS1(:, it)] = kep2ECI(kep);
    r_ECEF_reS1(:, it) = ECI2ECEF(r_reS1(:, it), tmS1(it));
    latlonH_reS1(:, it) = ECEF2latlonH(r_ECEF_reS1(:, it));
    [r_ENU_reS1(:, it), rElAz_reS1(:, it)] = ECEF2ENU(r_ECEF_reS1(:, it), ...
        37.448347027777770, 126.9524995277778, 277.6);
    [r_ENU_datS1(:, it), rElAz_datS1(:, it)] = ECEF2ENU(rECEFdatS1(:, it), ...
        37.448347027777770, 126.9524995277778, 277.6);
end
latlonH_datS1(1:2, :) = unwrap(wrapTo2Pi(latlonH_datS1(1:2, :)), [], 2);
latlonH_reS1(1:2, :) = unwrap(wrapTo2Pi(latlonH_reS1(1:2, :)), [], 2);

%% 7 : groundtrack
tm_xS1 = tmS1(1):((tmS1(end) - tmS1(1)) / 8):tmS1(end);
lonlat_dat_xS1 = spline(tmS1, latlonH_datS1([2, 1], :), tm_xS1);
lonlat_re_xS1 = spline(tmS1, latlonH_reS1([2, 1], :), tm_xS1);

textcomment = datestr(datenum(2000, 1, 1, 0, 0, 0) + tm_xS1', 'HH:MM:SS AM');
dx = 3; dy = 0;

figure(201);
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
    latlonH_datS1(2, :) * 180 / pi, latlonH_datS1(1, :) * 180 / pi, 'k', ...
    latlonH_reS1(2, :) * 180 / pi, latlonH_reS1(1, :) * 180 / pi, 'r' ...
    );
hold on;
hp = plot(lonlat_dat_xS1(1, :) * 180 / pi, lonlat_dat_xS1(2, :) * 180 / pi, 'k+');
hq = plot(lonlat_re_xS1(1, :) * 180 / pi, lonlat_re_xS1(2, :) * 180 / pi, 'r+');
hp.MarkerSize = 10;
hq.MarkerSize = 10;
text(lonlat_dat_xS1(1, :)' * 180 / pi + dx, lonlat_dat_xS1(2, :)' * 180 / pi + dy, textcomment);
legend(h(1:2), 'Measurement', 'Prediction (7)');
ylim([-10, 90]);
xlim([110, 190]);
title('SNUGLITE Ground Track, April 11, 2019');
xlabel('Longitude [\circ]');
ylabel('Latitude [\circ]');

%% skyplot
figure(401);
skyplot_7 = zeros(4, 1);
skyplot_7(1:4) = polarplot(rElAz_datS1(3, :), wrapToPi(rElAz_datS1(2, :)) * 180 / pi, 'k', ...
    rElAz_datS1(3, 1:10:end), wrapToPi(rElAz_datS1(2, 1:10:end)) * 180 / pi, 'k+', ...
    rElAz_reS1(3, :), wrapToPi(rElAz_reS1(2, :)) * 180 / pi, 'r', ...
    rElAz_reS1(3, 1:10:end), wrapToPi(rElAz_reS1(2, 1:10:end)) * 180 / pi, 'r+' ...
    );
textr = [];
texttheta = [];
textcomment = [];
for it = 1:10:sizedatS1
    if rElAz_reS1(2, it) > 0
        textr = [textr; rElAz_reS1(2, it) * 180 / pi];
        texttheta = [texttheta; rElAz_reS1(3, it)];
        textcomment = [textcomment; datestr(datenum(2000, 1, 1, 0, 0, 0) + tmS1(it), 'DD HH:MM:SS AM')];
    end
end
rlim([0, 90]);
text(texttheta, textr, textcomment);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.RDir = 'reverse';
title('SNUGLITE Skyplot, April 11-12, 2019');
legend(skyplot_7([1, 3]), 'Measurement', 'Prediction (7)');
