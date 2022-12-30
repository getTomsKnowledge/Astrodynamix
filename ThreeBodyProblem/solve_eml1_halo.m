%% Solve for Earth-Moon L1 Halo Orbit

%run('TidyWorkspace.m');
%run('CheckSNaG.m');
%heading = HeadingGenerator('EML1', 'Halo Solution',"Richardson's 'Analytic Construction...' (1979)");
%heading.print_heading();
%printing = PrintFormatting();
%consts = OrbitConstants();

%{
% Sun test values:
mu2 = consts.mu_sun_km*1e9;
mu1 = consts.mu_earth_km*1e9;
mu  = mu1 / (mu1 + mu2);
r1 = (consts.r_earthSun_km - consts.L1_earthSun_km)*1e3;
n1 = 1.99099*1e-7;
a1 = consts.r_earthSun_km*1e3;
gLSun = r1/a1;


sunC2 = (1/(gLSun^3)) * (  mu + (((1 - mu)*(gLSun^3))/((1 - gLSun)^3))  );
sunLambda = roots([1, 0, (c2 - 2), 0, -((c2 - 1)*(1 + 2*c2))])

%}


%% Set first Legendre polynomial coefficient, c2:

% Normalized parameters for dimensionless analysis:
muEarth = consts.mu_earth_km*1e9;
muMoon = consts.mu_moon_km*1e9;
mu = muMoon / (muMoon + muEarth);
sma = consts.r_earthMoon_km*1e3;
r1 = (consts.r_earthMoon_km - consts.L1_earthMoon_km)*1e3;
sma = consts.r_earthMoon_km*1e3;
n1 = sqrt((muMoon + muEarth)/(sma^3));
gL1 = r1 / sma;


% set c2:
c2 = l1_cn_generator(2, gL1, mu);

%% Get allowed orbital frequency, lambda:

lambdaCoefficients = [1, 0, (c2 - 2), 0, -( (c2 - 1)*(1 + 2*c2) )];
lambdaSolutions = roots(lambdaCoefficients);

for n = 1:4
    if (isreal(lambdaSolutions(n)))
        lambda = abs(lambdaSolutions(n));
    end
end

%% Calculate constants:

% Legendre coefficients:
c3 = l1_cn_generator(3, gL1, mu);
c4 = l1_cn_generator(4, gL1, mu);

% wave number?
k = ((lambda^2) + 1 + 2*c2)/(2*lambda);

% d:
d1 = ((3*(lambda^2))/k)*(k*(6*(lambda^2) - 1) - 2*lambda);
d2 = ((8*(lambda^2))/k)*(k*(11*(lambda^2) - 1) - 2*lambda);

% d2:
d21 = -(c3 / (2*(lambda^2)) );

% a2:
a21 = (3*c3*(k^2 - 2))/(4*(1 + 2*c2));
a22 = (3*c3)/(4*(1 + 2*c2));
a23 = -((3*c3*lambda)/(4*k*d1))*(3*(k^3)*lambda - 6*k*(k - lambda) + 4);
a24 = -((3*c3*lambda)/(4*k*d1))*(2 + 3*k*lambda);

% b2:
b21 = -((3*c3*lambda)/(2*d1))*(3*k*lambda - 4);
b22 = (3*c3*lambda)/d1;

% b3:
b31 = (3/(8*d2))* ...
    (...
        ( 8*lambda* ...
            ( 3*c3*(k*b21 - 2*a23) - c4*(2 + 3*(k^2))) ...
        )...
        + (9*(lambda^2) + 1 + 2*c2)* ...
            (4*c3*(k*a23 - b21) + k*c4*(4 + (k^2))...
        ) ...
     );
b32 = (1/d2)...
        *(...
            9*lambda*( c3*(k*b22 + d21 - 2*a24) -c4)...
            + (3/8)*(9*(lambda^2) + 1 + 2*c2)*( 4*c3*(k*a24 - b22) + k*c4)  ...
         );

% d3:
d31 = (3 / (64*(lambda^2))) * (4*c3*a24 + c4);
d32 = (3/(64 * (lambda^2)))*(4*c3*(a23 - d21) + c4*(4 + (k^2)));

% a3:
a31 = -((9*lambda)/(4*d2))*(4*c3*(k*a23 - b21) + k*c4*(4 + k^2)) ...
    + ((9*(lambda^2) + 1 - c2)/(2*d2))*(3*c3*(2*a23 - k*b21) + c4*(2 + 3*k^2));
a32 = -(1/d2)*...
    (...
        ((9*lambda)/4)*(4*c3*(k*a24 - b22) + k*c4) + ...
        (3/2)*(9*(lambda^2) + 1 - c2)*(c3*(k*b22 + d21 - 2*a24) - c4) ...
     );

% a:
a1 = -(3/2)*c3*(2*a21 + a23 + 5*d21) ...
    - (3/8)*c4*(12 - k^2);
a2 = (3/2)*c3*(a24 - 2*a22) + (9/8)*c4;

%% Solve for frequency correction, s1, s2:
sCoeff = (1 / (2*lambda*(lambda*(1 + (k^2)) - 2*k)));
s1 = sCoeff * ...
    (...
        (3/2)*c3*( 2*a21*((k^2) - 2) - a23*((k^2) + 2) - 2*k*b21 ) - ...
        (3/8)*c4*(3*(k^4) + 8*(k^2) + 8)...
    );
s2 = sCoeff * ...
    (...
        (3/2)*c3*( 2*a22*((k^2) - 2) + a24*((k^2) + 2) + 2*k*b22 + 5*d21 ) + ...
        (3/8)*c4*(12 - (k^2)) ...
    );

%% Solve amplitude-constraint relationship:

l1 = a1 + 2*(lambda^2)*s1;
l2 = a2 + 2*(lambda^2)*s2;



%% Get minimum x amplitude:

Delta = lambda^2 - c2;

axBar = sqrt(abs(Delta / l1));
axMin_m = r1 * axBar;
axMin = axMin_m/r1;

%% Get desired z amplitude:

badInput = 1;

% Get good input:
while (badInput)
    catchAz_m = input('Please enter initial z amplitude in meters:\n');
    catchAz = catchAz_m / r1;
    catchAx = sqrt(-(l2*(catchAz^2) + Delta)/l1);
    if (catchAx < axBar)
        fprintf('\nSorry, chief!  Bad input.  Ax > Ax_min = %0.1f km\n', axMin_m*1e-3);
        Az = 0;
        Ax = axBar;
    else
        badInput = 0;
        Az = catchAz;
        Ax = catchAx;
    end
end

% How angelic is it?
haloTilt = atand(Az/Ax);
fprintf('\n     Your halo is tilted %0.5f degrees about eTheta.\n', haloTilt);

%% Solve frequency correction term:
o2 = s1*(Ax^2) + s2*(Az^2);
omega = 1 + o2;

%% solution in synodic x,y,z coordinates:

% Get solution class:
enn = input('Please enter the desired solution class (enn = 1 for Class I, enn = 3 for Class II):\n');

% Phase variables:
phi0_deg = input('Enter initial phase in degrees on 0 to 360:\n');
t0 = input('Enter initial time in seconds:\n');
phi0 = deg2rad(phi0_deg);
psi0 = phi0 + enn*(pi/2);
tau10 = lambda*t0 + phi0;
delSwitch = 2 - enn;

% Synodic coordinate set:
x0_nonDim = a21*(Ax^2) ...
    + a22*(Az^2) ...
    - Ax*cos(tau10) ...
    + (a23*(Ax^2) - a24*(Az^2))*cos(2*tau10) ...
    + (a31*(Ax^3) - a32*Ax*(Az^2))*cos(3*tau10);
y0_nonDim = k*Ax*sin(tau10) ...
    + (b21*(Ax^2) - b22*(Az^2))*sin(2*tau10)...
    + (b31*(Ax^3) - b32*Ax*(Az^2))*sin(3*tau10);
z0_nonDim = delSwitch* ...
    ( ...
        Az*cos(tau10) ...
        + d21*Ax*Az*(cos(2*tau10) - 3) ...
        + ( d32*Az*(Ax^2) - d31*(Az^3) )*cos(3*tau10) ...
    );

% Set initial velocity:
xDot0_nonDim = 0;
yDot0_nonDim = (lambda*omega*n1)*( k*Ax + 2*(b21*(Ax^2) - b22*(Az^2)) + 3*(b31*(Ax^3) - b32*Ax*(Az^2)) );
zDot0_nonDim = 0;

x0_km = x0_nonDim*r1*1e-3;
y0_km = y0_nonDim*r1*1e-3;
z0_km = z0_nonDim*r1*1e-3;
yDot0_kms = yDot0_nonDim*r1*1e-3;

%% Propagate solution:

tauMax = 100;
tauSteps = 5000;
dTau = tauMax / tauSteps;
tauRange = [0:dTau:tauMax];
tau1Range = zeros(tauMax, 1);
rHalo = zeros(tauMax, 3);
for f = 1:tauSteps
    tau1Now =  lambda*tauRange(f) + phi0;
    tau1Range(f) = tau1Now;
    x = a21*(Ax^2) ...
        + a22*(Az^2) ...
        - Ax*cos(tau1Now) ...
        + (a23*(Ax^2) - a24*(Az^2))*cos(2*tau1Now) ...
        + (a31*(Ax^3) - a32*Ax*(Az^2))*cos(3*tau1Now);
    y = k*Ax*sin(tau1Now) ...
        + (b21*(Ax^2) - b22*(Az^2))*sin(2*tau1Now)...
        + (b31*(Ax^3) - b32*Ax*(Az^2))*sin(3*tau1Now);
    z = delSwitch* ...
        ( ...
            Az*cos(tau1Now) ...
            + d21*Ax*Az*(cos(2*tau1Now) - 3) ...
            + ( d32*Az*(Ax^2) - d31*(Az^3) )*cos(3*tau1Now) ...
        );
    rHalo(f, 1:3) = [x, y, z];
end

zAltitude = sprintf('Altitude above L1 on 3BP plane, Az = %0.2f km', catchAz_m*1e-3);
figure('Name', 'Look at this Graph');
plot3(rHalo(:,1), rHalo(:,2), rHalo(:,3));
title('Reference Orbit - Halo at L1', "Work by Tom West, based on Richardson's 'Analytic Construction of Periodic Orbits about the Collinear Points'");
xlabel('x (-)');
ylabel('y (-)');
zlabel('z (-)');

rHalo = rHalo * (consts.r1_moonToL1_km*1e3);
rHalo(:,1) = rHalo(:,1) + consts.xi_L1_km*1e3;

%{
x0 = ;
y0 = ;
z0 = ;
%}

%{
Ax > 0
Az >= 0 --> set Az to zero, solve for Ax min
%}

%% FUNCTIONS:

function cOut = l1_cn_generator(nIn, gIn, muIn)

    cOut = (1/(gIn^3))* ...
        (muIn +    ((-1)^nIn)*(((1 - muIn)*(gIn^(nIn+1)))/((1 - gIn)^(nIn+1)))  );
end