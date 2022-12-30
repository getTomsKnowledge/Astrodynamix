%{

    z1 = 450;
    z2 = 1800;


%}


function [deltai_1, deltai_2, deltaV1, deltaV2, deltaVT] = get_inclination_split(deltai, r1_m, r2_m, z1_m, z2_m)

r1_m = consts.R_earth_km*1e3 + z1_m;
r2_m = consts.R_earth_km*1e3 + z2_m;
a12_m = (r1_m + r2_m)/2;
mu_m = consts.mu_earth_km*1e9;

v_cs1 = sqrt(mu_m/r1_m);
v_cs2 = sqrt(mu_m/r2_m);
v_pt = sqrt(2*mu_m*(1/r1_m - 1/(2*a12_m)));
v_at = sqrt(2*mu_m*(1/r2_m - 1/(2*a12_m)));

syms alpha;

deltai = 15;

dV_tot = @(alpha) ...
    sqrt(v_pt^2 + v_cs1^2 - 2*v_pt*v_cs1*cosd(alpha)) ...
    + sqrt(v_at^2 + v_cs2^2 - 2*v_at*v_cs2*cosd(deltai-alpha));

options = optimset('display', 'off');

deltai_1 = fmincon(dV_tot, 0);
deltai_2 = deltai - deltai_1;

deltaV1 = sqrt(v_pt^2 + v_cs1^2 - 2*v_pt*v_cs1*cosd(deltai_1));
deltaV2 = sqrt(v_at^2 + v_cs2^2 - 2*v_at*v_cs2*cosd(deltai_2));
deltaVT = deltaV1 + deltaV2;



% Problem 2:
fprintf('  Inclination Split\n\n');
fprintf('     The inclination split is as follows:\n');
fprintf('        alpha = %0.3f deg.', deltai_1);
fprintf('        beta = %0.3f deg.\n\n', deltai_2);

fprintf('        The split maneuever dVs are:\n');
fprintf('          deltaV_1 = %0.3f km/s\n', deltaV1);
fprintf('          deltaV_2 = %0.3f km/s\n', deltaV2);
fprintf('          totalDV  = %0.3f km/s\n\n', deltaVT);

end