%% Determines position of GPS satellites from TLE / RINEX data

% NOTE:  ANGULAR INPUT IN RADIANS!!

%{
    % TEST DATA:

    % Obtained from RINEX elements:
    % line,pos:  2,4; 1,3; 2,2; 4,3; 1,4; 3,3; 4,4; 4,1;
    %            5,1; 1,2; 4,2; 3,4; 3,2; 2,3; 2,1


    % Sample RINEX tle:
    t = 5000;
    tGPSw = 500000;

    rGPS = get_gps_satellite_positions(...
    5.15366711998*1e3, 4.010524357057*1e-9, 8.152227383107*1e-3,...
    6.99666045202*1e-1, 9.29852083934*1e-1, 9.00856367346*1e-1,...
    -7.769966892113*1e-9, 9.729982837983*1e-1, 4.014452936740*1e-10,...
    -2.562500000000*1e1, 1.663750000000*1e2, 5.774199962616*1e-8,...
    -1.415610313416*1e-7, 1.143291592598*1e-5, -1.413747668266*1e-6,...
    t, tGPSw)
%}

function gpsOut = get_gps_satellite_positions(...
    aRoot, nDelta, e, o, MNought, ONought, ODot, iNought, iDot, ...
    Crs, Crc, Cis, Cic, Cus, Cuc, t, tGPSw)

    consts = OrbitConstants();
    earthRate = consts.oRateEarth_radsec;
    mu = consts.mu_earth_km*1e-9;

    % 1.) Calculate (corrected) mean motion:
    meanMotion = sqrt(mu) / (aRoot^3) + nDelta;

    % 2.) Mean anomaly:
    Mk = MNought + meanMotion*t;

    % 3.) Alternate mean anomaly:
    syms E;
    y(E) = E - e*sin(E) - Mk;
    Ek = solve(y == 0);

    % 4.) True anomaly:
    nu = atan(sqrt(1 - e^2)*sin(Ek), cos(Ek) - e);

    % 5.) Arg. of latitude, uncorrected:
    uPrimek = nu + o;

    % 6.) Arg. of latitude, corrected:
    uk = uPrimek + Cus*sin(2*uPrimek) + Cuc*cos(2*uPrimek);

    % 7.) Radius:
    rk = (aRoot^2)*(1 - e*cos(Ek));

    % 8.) Inclination:
    ik = iNought + iDot*t + Cis*sin(2*uPrimek) + Cic*cos(2*uPrimek);

    % 9.) Lon. of Ascending Node:
    Ok = ONought + ODot*t - earthRate*tGPSw;

    % 10.) In-plane position:
    xp = rk*cos(uk);
    yp = rk*sin(uk);

    % 11.) ECEF coordinates:
    X = xp*cos(Ok) - yp*cos(ik)*sin(Ok);
    Y = xp*sin(Ok) + yp*cos(ik)*cos(Ok);
    Z = yp*sin(ik);

    gpsOut = [X,Y,Z];


end