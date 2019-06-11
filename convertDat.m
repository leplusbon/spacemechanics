function satdata = convertDat(satdata)

satdata.YY = satdata.YY + 2000;

satdata.D = datenum( ...
    satdata.YY, satdata.MM, satdata.DD, ...
    satdata.HH, satdata.MM1, satdata.SS) - datenum(2000, 1, 1, 0, 0, 0) + datenum(0, 0, 0, 0, 0, 5);

satdata.GMSTdeg = 100.46 + 360.985612288088 * satdata.D;

for i = 1:size(satdata, 1)
    ECEF_rkm = [satdata.ECEF_Xkm(i); satdata.ECEF_Ykm(i); satdata.ECEF_Zkm(i)];
    ECEF_Vkms = [satdata.ECEF_V_Xkms(i); satdata.ECEF_V_Ykms(i); satdata.ECEF_V_Zkms(i)];
    GMSTrad = satdata.GMSTdeg(i) * pi / 180;
    tf = [
        cos(GMSTrad),   -sin(GMSTrad),  0;
        sin(GMSTrad),   cos(GMSTrad),   0;
        0,              0,              1
        ]; % ECEF->ECI, clockwise rotation
    ECI_rkm = tf * ECEF_rkm;
    ECI_Vkms = tf * ECEF_Vkms;
    ECI_rkm = ECEF2ECI(ECEF_rkm, satdata.D(i));
    %ECI_Vkms = ECEF2ECI(ECEF_Vkms, satdata.D(i));
    ECI_Vkms = ECEF2ECI(ECEF_Vkms, satdata.D(i))+(cross(7.292115857*10^(-5)*[0 0 1],ECI_rkm.').');
    satdata.ECI_Xkm(i) = ECI_rkm(1);
    satdata.ECI_Ykm(i) = ECI_rkm(2);
    satdata.ECI_Zkm(i) = ECI_rkm(3);
    satdata.ECI_V_Xkms(i) = ECI_Vkms(1);
    satdata.ECI_V_Ykms(i) = ECI_Vkms(2);
    satdata.ECI_V_Zkms(i) = ECI_Vkms(3);
end

