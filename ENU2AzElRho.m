function AzElRho = ENU2AzElRho(enu)

Rho = norm(enu);
Az = atan2(enu(1), enu(2));
El = atan2(enu(3), sqrt(enu(1)^2 + enu(2)^2));
AzElRho = [Az; El; Rho];
if enu(3) < 0
    AzElRho = [NaN; NaN; NaN];
end