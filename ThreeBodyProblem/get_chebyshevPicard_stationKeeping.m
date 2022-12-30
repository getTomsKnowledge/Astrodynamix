%% L1 Station Keeping - Chebyshev-Picard Iteration

%{

Parameters:
tau
rCM
omega
|r1| -> distance betwen Earth and L1

%}

run('TidyWorkspace.m');
%run('CheckSNaG.m');
heading = HeadingGenerator('EML1', 'Station Keeping', 'MCPI Solutions About L1');
printing = PrintFormatting();
consts = OrbitConstants();

% Constant parameters:
mu1 = consts.mu_earth_km*1e9;
mu2 = consts.mu_moon_km*1e9;
r12 = consts.r_earthMoon_km*1e3;
r1  = consts.L1_earthMoon_km*1e3;
r2  = r12 - r1;
rCM = r1 - consts.rCM_earthMoon_km*1e3; % location of L1 wrt CoM
mu  = mu2 / (mu1 + mu2);
omegaSq = (mu1 + mu2) / (r12^3);

x1 = r1 / r1;
x2 = r2 / r1;
x12 = r12 / r1;
xCM = rCM / r1;
x0 = rx/r1;
y0 = ry/r1;
z0 = rz/r1;
X0 = [x0, y0, z0, vx0, vy0, vz0]';
revNum = 10;
revSpan = revNum * 2 * pi;
tRange = 0:(1/sqrt(omegaSq)):(revSpan/sqrt(omegaSq));
tauSpan = [0 revSpan];

%tau = omega*t0;

[tauOut, Y] = ode45(@get_MCPI_ode, tauSpan, X0, [], x1, x2, x12, mu);

function xDot =  get_MCPI_ode(t, xIn, x1, x2, x12, xCM, muIn)

    rho1 = sqrt((xIn(1) - x1).^2 + xIn(2).^2 + xIn(3).^2);
    rho2 = sqrt((xIn(1) - x2).^2 + xIn(2).^2 + xIn(3).^2);

    % Velocities:    
    xDot(1) = xIn(4);
    xDot(2) = xIn(5);
    xDot(3) = xIn(6);
    
    % EOMs:
    xDot(4) = 2*xIn(5) + (xIn(1)-xCM) ...
        -(x12^3) * ...
        (  ((1 - muIn)/(rho1^3))*(xIn(1) - x1)  +  (muIn/(rho2^3))*(xIn(1) - x2) );
    
    xDot(5) = -2*xIn(4) + xIn(2) ...
        - (x12^3) * ...
        ( ((1 - muIn)/(rho1^3)) + (muIn/(rho2^3))  ) * xIn(2);
    
    xDot(6) = -(x12^3) * ...
        (  ((1 - muIn)/(rho1^3)) + (muIn/(rho2^3))  )*xIn(3);


end

function output = nonDim()


end