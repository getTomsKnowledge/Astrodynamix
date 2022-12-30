%% get_topo_los - Transforms from AZ EL Range to ENZ to TCEF
% inputs:  Azimuth, Elevation
% outputs: rhoHat = (x, y, z)_topo

function [rhoHat] = get_topo_los( ...
    AZ, EL, siteLLA)

    % AZELRANGE to ENZ:

    % Set of input angles:
    sphereAngles = [AZ, EL];

    % Simplify matrix calcs (prettify for eye):
    cAZ = cosd(sphereAngles(1));
    cEL = cosd(sphereAngles(2));
    sAZ = sind(sphereAngles(1));
    sEL = sind(sphereAngles(2));

    % Make transformation:
    rhoENZ = [cEL*sAZ;
               cEL*cAZ;
               sEL];

    % TIME:

    % Get local sidereal time:
    stOut = siderealtime(siteLLA.epoch, siteLLA); % st = [thetaLST thetaGST] deg
    thLST = stOut(1); % in degrees
    phi   = siteLLA.latitude_deg; % in degrees

    % ENZ - TCEF TRANSFORMATION:

    % Rotation vectors:
    Rz = [cosd(-90 - thLST)   sind(-90 - thLST)     0;
          -sind(-90 - thLST)    cosd(-90 - thLST)     0; 
               0                    0                1];

    Rx = [     1                    0            0;
               0              cosd(phi - 90)    sind(phi - 90);
               0              -sind(phi - 90)  cosd(phi - 90)];

    % TCE unit vector (Line-of-Sight vector)
    rhoHat = Rz * Rx * rhoENZ;

end