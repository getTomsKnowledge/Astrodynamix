%% Leapfrog 3BP Propagator

clc;
clear all;
close all;
home;

%% Setup:

% Gravitational parameters:
muSun = 1.327124400189e11; % km^3/s^2
muEarth = 3.9860044188e5; % km^3/s^2
muMoon = 4.90486959e3; % km^3/s^2

%{

Mom's birthday:

My birthday:

Chris's birthday:

Peggy's birthday:

%}

% Initial states:
% rEarth, my birthday = [1.085654661170196E+07 -1.516246665530904E+08 -9.740893969707191E+03]';
r0_Earth = [4.781615230795520E+07 1.393640709905750E+08 -1.437791641300917E+04]';
% vEarth, my birthday = [2.924778860543529E+01  2.018192849080797      1.110912062926461E-03]';
v0_Earth = [-2.871524177256561E+01 9.373338592355815E+00 -1.463442458238262E-05]';
X0_Earth = [r0_Earth; v0_Earth];

% rMoon, my birthday  = [1.059481977955440E+07 -1.513568853278044E+08 -1.365246557713300E+04]';
r0_Moon  = [4.817632258673497E+07 1.392504306482499E+08 1.733991004971415E+04]';
% vMoon, my birthday  = [2.845478370962656E+01  1.335568444795943     -9.403977943926289E-02]';
v0_Moon  = [-2.846702572208396E+01 1.037557201569116E+01 3.302227521669909E-02]';
X0_Moon  = [r0_Moon; v0_Moon];

% rSun, my birthday   = [3.891963415068314E+04  5.550923454368106E+04 -6.493728399841038E+03]';
r0_Sun   = [8.445244453936435E+05 -4.063250008632522E+05 -2.860802401718454E+04]';
% vSun, my birthday   = [8.758026516243883E-03  2.695060513075030E-03 -1.589335583160812E-04]';
v0_Sun   = [1.028055855381472E-02 6.431591197170886E-03 -2.545093399360199E-04]';
X0_Sun   = [r0_Sun; v0_Sun];

yBody2 = X0_Earth(1:3,1)';
yBody3 = X0_Moon(1:3,1)';
yBody1 = X0_Sun(1:3,1)';
vBody2 = X0_Earth(4:6,1)';
vBody3 = X0_Moon(4:6,1)';
vBody1 = X0_Sun(4:6,1)';

% year = 31536000;
year = 31536000*68;
%minute = 60;
minute = 60*60;

%Initialize trajectory matrices for plotting:
earthOut = zeros(year/minute, 3);
earthOut(1,1:3) = X0_Earth(1:3,1)';
moonOut = zeros(year/minute, 3);
moonOut(1,1:3) = X0_Moon(1:3,1)';
sunOut = zeros(year/minute, 3);
sunOut(1,1:3) = X0_Sun(1:3,1)';
bigMu = muSun + muEarth + muMoon;

h = minute;

%% Run:




for n = 0:minute:year

    R_Body2Body3 = yBody3 - yBody2;
    r_Body2Body3 = norm(R_Body2Body3);
    R_Body2Body1  = yBody1 - yBody2;
    r_Body2Body1 = norm(R_Body2Body1);
    R_Body1Body3   = yBody3 - yBody1;
    r_Body1Body3 = norm(R_Body1Body3);

    % Velocity Verlet:
    vBody2_nPlus1Half = vBody2 + (h/2) ...
        .* (  ((muSun/(r_Body2Body1^3))   .*  R_Body2Body1) ...
        +  ( (muMoon/(r_Body2Body3^3))   .*  R_Body2Body3)  );

    vBody3_nPlus1Half  = vBody3  + (h/2) ...
        .* (  ((muSun/(r_Body1Body3^3))    .* -R_Body1Body3) ...
        +  ( (muEarth/(r_Body2Body3^3))  .* -R_Body2Body3)  );

    vBody1_nPlus1Half   = vBody1   + (h/2) ...
        .* (  ((muEarth/(r_Body2Body1^3)) .* -R_Body2Body1) ...
        +  ( (muMoon/(r_Body1Body3^3))     .*  R_Body1Body3)  );

    yBody2_nPlusOne = yBody2  + h .* vBody2_nPlus1Half;
    yBody3_nPlusOne  = yBody3   + h .* vBody3_nPlus1Half;
    yBody1_nPlusOne   = yBody1    + h .* vBody1_nPlus1Half;
    
    % Update interplanetary distances:
    R_Body2Body1   = yBody1_nPlusOne - yBody2_nPlusOne;
    r_Body2Body1   = norm(R_Body2Body1);
    R_Body2Body3  = yBody3_nPlusOne - yBody2_nPlusOne;
    r_Body2Body3  = norm(R_Body2Body3);
    R_Body1Body3    = yBody3_nPlusOne - yBody1_nPlusOne;
    r_Body1Body3    = norm(R_Body1Body3);
        
    vBody2_nPlus3Halves = (h/2) ...
        .* (  ((muSun/(r_Body2Body1^3))   .*  R_Body2Body1) ...
        +  ( (muMoon/(r_Body2Body3^3))   .*  R_Body2Body3)  );
    vBody2_nPlusOne = vBody2_nPlus1Half + vBody2_nPlus3Halves;
    
    vBody3_nPlus3Halves  = (h/2) ...
        .* (  ((muSun/(r_Body1Body3^3))    .* -R_Body1Body3) ...
        +  ( (muEarth/(r_Body2Body3^3))  .* -R_Body2Body3)  );
    vBody3_nPlusOne = vBody3_nPlus1Half + vBody3_nPlus3Halves;
    
    vBody1_nPlus3Halves   = (h/2) ...
        .* (  ((muEarth/(r_Body2Body1^3)) .* -R_Body2Body1) ...
        +  ( (muMoon/(r_Body1Body3^3))     .*  R_Body1Body3)  );    
    vBody1_nPlusOne = vBody1_nPlus1Half + vBody1_nPlus3Halves;

    if (n ~= 0)
        r_baryOld = (1/bigMu) * (muSun .* yBody1 + muEarth .* yBody2 + muMoon .* yBody3);
        r_baryNew = (1/bigMu) * (muSun * yBody1_nPlusOne + muEarth * yBody2_nPlusOne + muMoon * yBody3_nPlusOne);
        vBary = (r_baryNew - r_baryOld)/minute;
        RBody2_curr = yBody2_nPlusOne - r_baryNew;
        RBody3_curr = yBody3_nPlusOne - r_baryNew;
        RBody1_curr = yBody1_nPlusOne - r_baryNew;
        VBody2_curr = vBody2_nPlusOne - vBary;
        VBody3_curr = vBody3_nPlusOne - vBary;
        VBody1_curr = vBody1_nPlusOne - vBary;

        yBody2 = RBody2_curr;
        vBody2 = VBody2_curr;
        yBody3  = RBody3_curr;
        vBody3  = VBody3_curr;
        yBody1   = RBody1_curr;
        vBody1   = VBody1_curr;
 
    end

    if (n == 0)
        for m = 1:3
            earthOut(2,m) = yBody2(1,m);
            moonOut(2,m)  = yBody3(1,m);
            sunOut(2,m)   = yBody1(1,m);
        end
    else
        for o = 1:3
            earthOut( ((n/minute) + 2), o) = yBody2(1,o);
            moonOut( ((n/minute) + 2), o)  = yBody3(1,o);
            sunOut( ((n/minute) + 2), o)   = yBody1(1,o);
        end
    end  
end

figure('Name', '3BP Propagator Solution (By: TVW)');
plot3(earthOut(:,1), earthOut(:,2), earthOut(:,3), 'g');
hold on;
plot3(moonOut(:,1), moonOut(:,2), moonOut(:,3), 'b');
hold on;
plot3(sunOut(:,1), sunOut(:,2), sunOut(:,3), 'r');
grid on;
title({'3BP Propagator Solution', 'per Kasdin & Paley (2011), App. C, p. 656 "velocity Verlet" algorithm'});
xlabel('x');
ylabel('y');
zlabel('z');

%% FUNCTIONS:

%{
% Gravitational parameters:
muBody1 = 2.01e-4; % km^3/s^2
% 1.327124400189e11
muBody2 = 1.9e-4; % km^3/s^2
% 3.9860044188e5
muBody3 = 1.07e-4; % km^3/s^2
% 4.90486959e3

% Initial states:
% r0_Body2 = [1.085654661170196E+07 -1.516246665530904E+08 -9.740893969707191E+03]';
r0_Body2 = [1 1 1]';
%v0_Body2 = [2.924778860543529E+01  2.018192849080797      1.110912062926461E-03]';
v0_Body2 = [0  0  -0.01]';
X0_Body2 = [r0_Body2; v0_Body2];

%r0_Body3  = [1.059481977955440E+07 -1.513568853278044E+08 -1.365246557713300E+04]';
r0_Body3  = [1 0.1 -0.9]';
%v0_Body3  = [2.845478370962656E+01  1.335568444795943     -9.403977943926289E-02]';
v0_Body3  = [-0.01  0  0]';
X0_Body3  = [r0_Body3; v0_Body3];

%r0_Body1   = [3.891963415068314E+04  5.550923454368106E+04 -6.493728399841038E+03]';
r0_Body1   = [-1.1  -0.61  1.01]';
%v0_Body1   = [8.758026516243883E-03  2.695060513075030E-03 -1.589335583160812E-04]';
v0_Body1   = [0 0.01 -0.001]';
X0_Body1   = [r0_Body1; v0_Body1];

yBody2 = X0_Body2(1:3,1)';
yBody3 = X0_Body3(1:3,1)';
yBody1 = X0_Body1(1:3,1)';
vBody2 = X0_Body2(4:6,1)';
vBody3 = X0_Body3(4:6,1)';
vBody1 = X0_Body1(4:6,1)';

% year = 31536000;
year = 40000;
%minute = 60;
minute = 1;

%Initialize trajectory matrices for plotting:
earthOut = zeros(year/minute, 3);
earthOut(1,1:3) = X0_Body2(1:3,1)';
moonOut = zeros(year/minute, 3);
moonOut(1,1:3) = X0_Body3(1:3,1)';
sunOut = zeros(year/minute, 3);
sunOut(1,1:3) = X0_Body1(1:3,1)';
bigMu = muBody1 + muBody2 + muBody3;

h = minute;
%}

