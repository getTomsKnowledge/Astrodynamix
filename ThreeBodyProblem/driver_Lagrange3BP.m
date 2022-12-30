%% Lagrange Points and Zero-Velocity Contours

%{

This driver and set of methods will allow the user
to calculate various quantities associated with the
three-body problem of celestial mechanics.  Quantities
include:

- Lagrange equillibrium points in the restricted orbital plane
- Surfaces of zero relative veolicity
- Stability of Lagrange points
- Trajectories for a test mass in the restricted 3BP
- Elliptical orbits for three massive bodies in the general problem

%}

%{

References:

[1] Schaub, Hanspeter; Junkins, John L.  _Analytical Mechanics of Space
Systems_ .  4th ed.  2018.  Reston, VA: American Institute of Aeronautics
and Astronautics, Inc., 2018. pp. 577-620.
[2] Kasdin, N. Jeremy; Paley, Derek A. _Engineering Dynamics: A Comprehensive Introduction_. Princeton,
NJ:  Princeton University Press, 2011.  pp. 324 - 328.
[3] Bate, Roger R.; Mueller, Donald D.; White, Jerry E.; Saylor, William W.
_Fundamentals of Astrodynamics_.  2nd ed. Garden City, NY: Dover Publications, 2020.  pp. 4-10.
[4] _An Introduction to the Mathematics and Methods of Astrodynamics_. Rev. ed. Reston, VA: American Institute of Aeronautics and
Astronautics, 1999.  pp. 365-408.

%}

%% Setup (Formatting, Global Constants)

run('TidyWorkspace.m');
%run('CheckSNaG.m');
printing = PrintFormatting();
heading = HeadingGenerator(...
    'EML1',...
    'Lagrange 3BP Study',...
    'Calculates Lagrange points, zero-velocity contours, etc...');

heading.print_heading();
consts = OrbitConstants();

%% Principles

fprintf('\n\nPRINCIPLES:\n\n');
fprintf("   'Prop. 73-74: A body inside a homogeneous sphere, each of whose points\n");
fprintf("    is the center of a centripetal inverse-square law, is attracted by a force\n");
fprintf("    directed to the center of the sphere, proportional to the distance to the\n");
fprintf("    center; while a body outside the sphere is attracted by a force inversely\n");
fprintf("    proportional to the square of distance to the center.'\n        - Principia, Sir Isaac Newton, 1687\n");

%% Frame Rotation Rate

omegaSq = (consts.mu_earth_km + consts.mu_moon_km) ...
    / (consts.r_earthMoon_km^3);
revRate_sec = (2*pi) / sqrt(omegaSq);
revRate_days = revRate_sec/(24*60*60);

fprintf('\n\n     Instantaneous Rotation-Rate of Earth-Moon Frame: \n');
fprintf('        omega = %0.5d rad/s\n', sqrt(omegaSq));
fprintf('        1 revolution = %0.3f days\n\n', revRate_days);

%% Barycentric Coordinates

rho_m = consts.r_earthMoon_km*1e3;
muEarth_ms = consts.mu_earth_km*1e9;
muMoon_ms = consts.mu_moon_km*1e9;

RG_Earth_km = ...
    (muMoon_ms/(muMoon_ms + muEarth_ms))...
    *rho_m;

xi_Moon  = rho_m - RG_Earth_km;
xi_Earth = xi_Moon - rho_m;

fprintf('\n     In Earth-Moon barycentric frame:\n');
fprintf('        rEarth = %0.2f km\n', xi_Earth*1e-3);
fprintf('        rMoon  = %0.0f km\n\n', xi_Moon*1e-3);


%% Lagrange Points

%{
% [2] Kasdin & Paley:  restricted 3BP planar solution, (x,y)_B

% Rename constants for brevity:
x1 = xi_Moon;
mu1 = consts.mu_moon_km*1e9;
x2 = -xi_Earth;
mu2 = consts.mu_earth_km*1e9;

% Lagrange Quintic Equation:
a5 = omegaSq;
a4 = 2*omegaSq*(x2 - x1);
a3 = omegaSq*(x2^2 - 4*x1*x2 + x1^2);
a2 = 2*omegaSq*((x1^2)*x2 - x1*(x2^2)) - (mu1 + mu2);
a1 = omegaSq*(x1^2)*(x2^2) - 2*(mu1*x2 - mu2*x1);
a0 = -(mu1*(x2^2) + mu2*(x1^2));

quintCoeff = [a5, a4, a3, a2, a1, a0];

quinticSolutions = roots(quintCoeff);
index = 1;
solution = 0;
for counter = 1:4
    solution = quinticSolutions(counter);
    if (isreal(solution))
        lagrangePoints(index) = solution;
        index = index + 1;
    end
end

lagrangePoints


% Alternate solution:
M1 = mu2/consts.G_Nmkg;
M2 = mu1/consts.G_Nmkg;
R = rho_m;
syms r;
f(r) = (M1 / ((R - r)^2)) - M2/(r^2) - (((M1 * R) / (M1 + M2)) - r) * ((M1 + M2)/(R^3));

altQuintic = double(solve(f(r) == 0));

%}

% Earth-Moon Bipolar Coordinates, rho1/rho2

% [4] p. 380, eqns. 8.37 - 8.39

% Initial Conditions:
xi1 = xi_Earth;
mu1 = muEarth_ms;
xi2 = xi_Moon;
mu2 = muMoon_ms;

% Collinear Solutions:

% L1:

% Quintic Coefficients
a5Battin = -(mu1 + mu2);
a4Battin =  (3*mu1 + 2*mu2)*rho_m;
a3Battin = -(3*mu1 + mu2)*(rho_m^2);
a2Battin =  mu2*(rho_m^3);
a1Battin = -2*mu2*(rho_m^4);
a0Battin =  mu2*(rho_m^5);

battinCoefficients = ...
    [a5Battin, a4Battin, a3Battin, a2Battin, a1Battin, a0Battin];

rho2Roots_L1 = roots(battinCoefficients);

for counter = 1:5
    currentRoot = rho2Roots_L1(counter);
    if (isreal(currentRoot))
        rho2_L1 = currentRoot;
    end
end

rho1_L1 = rho_m - rho2_L1;
L1 = [rho1_L1 + xi1, 0];

fprintf('\n\n     L1 distance from Earth along Earth-Moon axis:\n');
fprintf('        R_L1 = %0.0f km', rho1_L1*1e-3);

% Get L2:

clear a5Battin a4Battin a3Battin ...
    a2Battin a1Battin a0Battin ...
    battinCoefficients;

% Quintic Coefficients
a5Battin = -(mu1 + mu2);
a4Battin = -(3*mu1 + 2*mu2)*rho_m;
a3Battin = -(3*mu1 + mu2)*(rho_m^2);
a2Battin =  mu2*(rho_m^3);
a1Battin =  2*mu2*(rho_m^4);
a0Battin =  mu2*(rho_m^5);

battinCoefficients = ...
    [a5Battin, a4Battin, a3Battin, a2Battin, a1Battin, a0Battin];

rho2Roots_L2 = roots(battinCoefficients);

for counter = 1:5
    currentRoot = rho2Roots_L2(counter);
    if (isreal(currentRoot))
        rho2_L2 = currentRoot;
    end
end

rho1_L2 = rho_m + rho2_L2;
L2 = [rho1_L2 + xi1, 0];

fprintf('\n\n     L2 distance from Earth along Earth-Moon axis:\n');
fprintf('        R_L2 = %0.0f km', rho1_L2*1e-3);

% Get L3:

clear a5Battin a4Battin a3Battin ...
    a2Battin a1Battin a0Battin ...
    battinCoefficients;

% Quintic Coefficients
a5Battin =   mu1 + mu2;
a4Battin =  (2*mu1 + 3*mu2)*rho_m;
a3Battin =  (mu1 + 3*mu2)*(rho_m^2);
a2Battin =  -mu1*(rho_m^3);
a1Battin =  -2*mu1*(rho_m^4);
a0Battin =  -mu1*(rho_m^5);

battinCoefficients = ...
    [a5Battin, a4Battin, a3Battin, a2Battin, a1Battin, a0Battin];

rho1Roots_L3 = roots(battinCoefficients);

for counter = 1:5
    currentRoot = rho1Roots_L3(counter);
    if (isreal(currentRoot))
        rho1_L3 = currentRoot;
    end
end

rho2_L3 = rho_m + rho1_L3;
L3 = [-rho1_L3 + xi1, 0];

fprintf('\n\n     L3 distance from Earth along Earth-Moon axis:\n');
fprintf('        R_L3 = %0.0f km\n\n', -rho1_L3*1e-3);

% Equilateral Triangle Solutions:
alpha = pi/3;
xiEquilat_m  = rho_m*cos(alpha);
etaEquilat_m = rho_m*sin(alpha);

% L4:
fprintf('     L4 Coordinates (relative to earth):\n');
fprintf('        L4: <%0.0f, %0.0f> km\n\n', ...
    xiEquilat_m*1e-3, etaEquilat_m*1e-3);

L4 = [xiEquilat_m + xi1, etaEquilat_m];

% L5:
fprintf('     L5 Coordinates (relative to earth):\n');
fprintf('        L5: <%0.0f, %0.0f> km\n\n', ...
    xiEquilat_m*1e-3, -etaEquilat_m*1e-3);

L5 = [xiEquilat_m + xi1, -etaEquilat_m];

%{

syms xi eta zeta C;

rho1Sq(xi, eta, zeta) = (xi - xi1_Earth)^2 + eta^2 + zeta^2;
rho2Sq(xi, eta, zeta) = (xi - xi2_Moon)^2 + eta^2 + zeta^2;


vRelSq(xi, eta, C) = ...
    omegaSq*(eta^2 + eta^2) ...
    + 2*(mu1/sqrt(rho1Sq) + mu2/sqrt(rho2Sq)) - C;

% Jacobi Integral:
J(xi, eta, zeta) =...
    (mu1/2)*( rho1Sq/(rho^3) + 2/sqrt(rho1Sq)) ...
    + (mu2/2)*( rho2Sq/(rho^3) + 2/sqrt(rho2Sq));
%}

%% L4/L5 Stability Criterion:

% Equilateral Lagrange Points:

L4Check = mu1/mu2 + mu2/mu1;
L4StabilityCriterion = 25;
if (L4Check >= L4StabilityCriterion)
    fprintf('\n\n     L4, L5 meet stability per:\n');
    fprintf('        m1/m2 + m2/m1 = %0.2f >= 25\n\n', L4Check);
end


%% Loci of Zero Relative Velocity

lociNumber = 25;
cStarVector = zeros(lociNumber,1);
cStar0 = 3*(rho_m^2)*omegaSq;
cStarVector = cStar0*[...
    1.000001, 1.000002, 1.00004, 1.00005, 1.0001, ...
    1.001, 1.002, 1.003, 1.005, 1.0075 ...
    1.0076, 1.0077, 1.0079, 1.008, 1.085, ...
    1.01, 1.05, 1.1, 1.2, 1.5, ...
    1.75, 2, 2.5, 3, 4]'; % all values in m/s
dimension = 2;
timeStep = 5000;
xiEtaOutput = zeros(dimension*timeStep, lociNumber*dimension);

thisX2 = 0;
thisX1 = 0;
maxX1 = 0;
minX1 = 0;
maxX2 = 0;
minX2 = 0;
firstSol = 0;
secondSol = 0;

for outerLoopCount = 1:lociNumber


    
    % 1.  Set C*, A, B:
    cStar = cStarVector(outerLoopCount);
    B = 3 + (rho_m/mu2)*(cStar - 3*(rho_m^2)*omegaSq);

    if (B < 0)
        continue;
    end

    % Solve for range of x2 values:

    x_Coefficients = [1, 0, -B, 2];
    x_Solutions = roots(x_Coefficients);
    
    index = 1;
    for xCount = 1:3
        xSol = x_Solutions(xCount);
        if (isreal(xSol) && (xSol > 0))
            realX(index) = xSol;
            index = index + 1;
        end
    end

    [minX2, maxX2] = bounds(realX, 'all');

    if (minX2 == maxX2)
        continue;
    end

    % Solve for two values of x1:
    for innerLoopTwo = 1:timeStep

        deltaX2 = (maxX2 - minX2)/timeStep;
        x2_Values = minX2:deltaX2:maxX2;
        currentX2 = x2_Values(innerLoopTwo);

        A = ((rho_m*cStar)/mu1) ...
            - (mu2/mu1)*(currentX2^2 + 2/currentX2);

        if (A < 3)
            continue;
        end

        x1_Coefficients = [1, 0, -A, 2];
    
        % 2.  Get x1, x2:
        x1_Solutions = roots(x1_Coefficients);

        index = 1;
        for x1Count = 1:3
            x1Sol = x1_Solutions(x1Count);
            if (isreal(x1Sol) && (x1Sol > 0))
                realX1(index) = x1Sol;
                index = index + 1;
            end
        end

        [minX1, maxX1] = bounds(realX1, 'all');

        xPairs(innerLoopTwo, 1:4) = [minX1, currentX2, maxX1, currentX2];

    end

    % 3.  Transform to bipolar coords. r1,r2:
    rhoVectorTop = [xPairs(:,1), xPairs(:,2)];
    rhoVectorBottom = [xPairs(:,3), xPairs(:,4)];
    rhoVector = rho_m*[rhoVectorTop; rhoVectorBottom];

    % 4.  Get alpha:
    badEtaCount = 0;
    badXiCount = 0;
    for rhoLoop = 1:length(rhoVector)
        rho1Now = rhoVector(rhoLoop, 1);
        rho2Now = rhoVector(rhoLoop, 2);
        alpha = real(acosd( ( (rho1Now^2 + rho_m^2) - rho2Now^2 ) / (2*rho1Now*rho_m)));
        beta  = real(acosd( ( (rho2Now^2 + rho_m^2) - rho1Now^2 ) / (2*rho2Now*rho_m)));
        % 5. Transform to barycentric rotating frame:
        etaNow1 = rho1Now*sind(alpha);
        etaNow2 = rho2Now*sind(beta);
        if (etaNow1 ~= etaNow2)
            badEtaCount = badEtaCount + 1;
            badSin(badEtaCount, 1) = alpha;
            badSin(badEtaCount, 2) = beta;
        end
        xiNow1  = ( etaNow1 / tand(alpha) ) - abs(xi1);
        xiNow2  = abs(xi2) - ( etaNow2 / tand(beta) );
        if (xiNow1 ~= xiNow2)
            badXiCount = badXiCount + 1;
            badTan(badXiCount, 1) = alpha;
            badTan(badXiCount, 2) = beta;
        end
        xiEtaOutput(rhoLoop, (2*outerLoopCount - 1):(2*outerLoopCount)) = [xiNow1, etaNow1];
    end

end

badEtaCount
badXiCount

% Plot 

lagrangePointsSet = [L1; L2; L3; L4; L5];

C1 = sprintf('%0.3d', (1/cStar0)*cStarVector(1));
C2 = sprintf('%0.3d', (1/cStar0)*cStarVector(2));
C3 = sprintf('%0.3d', (1/cStar0)*cStarVector(3));
C4 = sprintf('%0.3d', (1/cStar0)*cStarVector(4));
C5 = sprintf('%0.3d', (1/cStar0)*cStarVector(5));
C6 = sprintf('%0.3d', (1/cStar0)*cStarVector(6));
C7 = sprintf('%0.3d ', (1/cStar0)*cStarVector(7));
C8 = sprintf('%0.3d', (1/cStar0)*cStarVector(8));
C9 = sprintf('%0.3d', (1/cStar0)*cStarVector(9));
C10 = sprintf('%0.3d', (1/cStar0)*cStarVector(10));
C11 = sprintf('%0.3d', (1/cStar0)*cStarVector(11));
C12 = sprintf('%0.3d', (1/cStar0)*cStarVector(12));
C13 = sprintf('%0.3d', (1/cStar0)*cStarVector(13));
C14 = sprintf('%0.3d', (1/cStar0)*cStarVector(14));
C15 = sprintf('%0.3d', (1/cStar0)*cStarVector(15));
C16 = sprintf('%0.3d', (1/cStar0)*cStarVector(16));
C17 = sprintf('%0.3d', (1/cStar0)*cStarVector(17));
C18 = sprintf('%0.3d', (1/cStar0)*cStarVector(18));
C19 = sprintf('%0.3d', (1/cStar0)*cStarVector(19));
C20 = sprintf('%0.3d', (1/cStar0)*cStarVector(20));
C21 = sprintf('%0.3d', (1/cStar0)*cStarVector(21));
C22 = sprintf('%0.3d', (1/cStar0)*cStarVector(22));
C23 = sprintf('%0.3d', (1/cStar0)*cStarVector(23));
C24 = sprintf('%0.3d', (1/cStar0)*cStarVector(24));
C25 = sprintf('%0.3d', (1/cStar0)*cStarVector(25));

earthCircle = circle(xi1, 0, consts.R_earth_km*1e3);
moonCircle = circle(xi2, 0, consts.R_moon_km*1e3);

figure('Name', 'Zero-Velocity Surfaces for Restricted 3BP');
scatter(lagrangePointsSet(:,1), lagrangePointsSet(:,2), 'x');
hold on;
circle(xi1, 0, consts.R_earth_km*1e3);
hold on;
circle(xi2, 0, consts.R_moon_km*1e3);
hold on;
for plotCount = 1:lociNumber
    xiIndex  = 2*plotCount - 1;
    etaIndex = 2*plotCount;
    scatter(xiEtaOutput(:, xiIndex), xiEtaOutput(:, etaIndex), 0.5,'*', 'b');
    scatter(xiEtaOutput(:, xiIndex), -xiEtaOutput(:, etaIndex), 0.5, '*', 'b');
    hold on;
end
axis equal;
xlabel('$\xi$ (m)', 'Interpreter', 'latex');
ylabel('$\eta$ (m)', 'Interpreter', 'latex');
title({'Zero-Velocity Curves for the Earth-Moon System', "Work by Tom West, based on Battin's treatment, pp. 376-380"})
legend(...
    'Lagrange Points (L1 - L5)', 'Earth', 'Moon',...
    'Surfaces of Zero Relative Velocity', ...
    'Interpreter', 'latex');

%% Lagrange Point Stability

% Plot 1 month of deviation, starting with small kick 
% -->  1 m/s along Earth-Moon axis:  x0 = [0, 0, 1, 0]
daySpan   = 24*60*60; % one day in seconds
monthSpan = 30*daySpan; % one month in seconds
weekSpan  = 7*daySpan; % one week in seconds
twoWeeks  = 2*weekSpan; 
threeDays = 3*daySpan;
fourDays = 4*daySpan;
fiveDays = 5*daySpan;
fiveHours = 5*60*60;
sixHours  = 6*60*60;
tenHours  = 10*60*60;
twelveHours = 12*60*60;
timeSpan = [0 monthSpan];
xi_L1 = rho1_L1 + xi1;
omega = sqrt(omegaSq);
X0 = [xi_L1 - 3715080, 0, 10500000, -4000, 1820, 0]'; % Trying to jump back to earth :-)
[tL1, disturbedL1] = ode45(@propagate_lagrange_stability, timeSpan, X0, [], xi1, xi2, mu1, mu2, omega);

% Collinear Lagrange Points:
run('solve_eml1_halo.m');

lZeros = zeros(length(lagrangePointsSet), 1);
zZeros = zeros(length(xiEtaOutput), 1);

figure('Name', '5 m/s Impulse at L1 over 12 Hours');
%scatter(lagrangePointsSet(:,1)*1e-3, lagrangePointsSet(:,2)*1e-3, 'x');
%hold on;
%circle(xi1*1e-3, 0, consts.R_earth_km);
%hold on;
%circle(xi2*1e-3, 0, consts.R_moon_km);
%hold on;
scatter3(lagrangePointsSet(:,1), lagrangePointsSet(:,2), lZeros(:), 'x');
hold on;
%circle(xi1, 0, consts.R_earth_km*1e3);
%hold on;
%circle(xi2, 0, consts.R_moon_km*1e3);
%hold on;
for plotCount = 1:lociNumber
    xiIndex  = 2*plotCount - 1;
    etaIndex = 2*plotCount;
    scatter3(xiEtaOutput(:, xiIndex), xiEtaOutput(:, etaIndex), zZeros(:), 0.5,'*', 'b');
    scatter3(xiEtaOutput(:, xiIndex), -xiEtaOutput(:, etaIndex), zZeros(:), 0.5, '*', 'b');
    hold on;
end
%plot3(disturbedL1(:, 1), disturbedL1(:, 2), disturbedL1(:, 3), 'r');
plot3(rHalo(:,1), rHalo(:,2), rHalo(:,3), 'r');
hold on;
title({'Reference Orbit - Halo at 25,000 km altitude over EML1', ...
    "Work by Tom West, based on Richardson's 'Analytic Construction of Periodic Orbits about the Collinear Points'"});
xlabel('$\xi$ (m)', 'Interpreter', 'latex');
ylabel('$\eta$ (m)', 'Interpreter', 'latex');
legend(...
    'Lagrange Points (L1 - L5)', 'Earth', 'Moon',...
    'Surfaces of Zero Relative Velocity', ...
    'Interpreter', 'latex');
axis equal;
zlim([-100000000, 100000000]);

%% Trajectories in Bipolar Coordinates

%% FUNCTIONS:
function h = circle(x,y,r)
    hold on;
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit, 'k');
end

function xDot = propagate_lagrange_stability(t, xIn, xi1In, xi2In, mu1In, mu2In, o)

    % Radial vectors:
    rho1 = sqrt((xIn(1) - xi1In).^2 + xIn(2).^2 + xIn(3).^2);
    rho2 = sqrt((xIn(1) - xi2In).^2 + xIn(2).^2 + xIn(3).^2);

    % Forcing matrix:
    f0Coeff = (mu1In/(rho1^3) + mu2In/(rho2^3));
    F0 = [2  0  0;
          0 -1  0;
          0  0 -1];
    F0 = f0Coeff*F0;

    % Frame angular velocity:
    O  = [0 -o  0;
          o  0  0;
          0  0  0];

    % State-Space Transition Matrix:
    M = [zeros(3), eye(3);
        (F0 - O*O),  -2*O];

    xDot = M * xIn;

end