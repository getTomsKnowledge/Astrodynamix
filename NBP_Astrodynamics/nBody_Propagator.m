%% n-Body Propagator

clear all;
close all;
clc;
home;

lineBreak = '=============================================================================';
programName = 'N-Body Propagator (A Leapfrog-Algorithm Implementation)';
author = 'Tom West';
date = '05/26/2022';


fprintf('\n%s\n%s\n\n', lineBreak, lineBreak);
fprintf('     Program: %s\n     Author: %s\n     Date: %s\n\n', programName, author, date);
fprintf('%s\n%s\n\n', lineBreak, lineBreak);
%% Solar System Data:

% Table of gravitational parameters:
muSun       = 1.327124400189e11;
muMercury   = 2.20329e4;
muVenus     = 3.248599e5;
muEarth     = 3.9860044188e5;
muMoon      = 4.90486959e3;
muMars      = 4.2828372e4;
muCeres     = 6.26325e1;
muJupiter   = 1.266865349e8;
muSaturn    = 3.79311879e7;
muUranus    = 5.7939399e6;
muNeptune   = 6.8365299e6;
muPluto     = 8.719e2;
muEris      = 1.1089e3;

% J2000 data (PV):
% Position:
RSun_J2000      = [-1.067706805381631E+06 -4.182752718185146E+05 3.086181725478008E+04]; 
RMercury_J2000  = [-2.052943316392625E+07 -6.733155053453046E+07 -3.648992526181437E+06];
RVenus_J2000    = [-1.085242008576727E+08 -5.303290245135058E+06 6.166496117013918E+06];
REarth_J2000    = [-2.756674048064499E+07 1.442790215211286E+08 3.025066782881320E+04]; 
RMoon_J2000     = [-2.785834886487951E+07 1.440040417795086E+08 6.652186445663124E+04]; 
RMars_J2000     = [2.069805421179929E+08 -2.425286841879086E+06 -5.125451306189817E+06];
RCeres_J2000    = [-3.570100661713269E+08 1.185847406508015E+08 6.929550952765751E+07];
RJupiter_J2000  = [5.974999178522581E+08 4.391864046755430E+08 -1.519599985574219E+07];
RSaturn_J2000   = [9.573176521108806E+08 9.824380076870195E+08 -5.518211788151336E+07];
RUranus_J2000   = [2.157907331373991E+09 -2.055043522898880E+09 -3.559460196112704E+07];
RNeptune_J2000  = [2.513978490096488E+09 -3.739132780869018E+09 1.906330132972622E+07];
RPluto_J2000    = [-1.478398655393869E+09 -4.182993264891145E+09 8.752463989487143E+08];
REris_J2000     = [1.322247531122232E+10 4.602025019881923E+09 -3.903660777776308E+09];
% Velocity:
VSun_J2000      = [9.312571926508239E-03 -1.282475570795343E-02 -1.633507186347347E-04];
VMercury_J2000  = [3.700430442865286E+01 -1.117724068322721E+01 -4.307791469481385E+00];
VVenus_J2000    = [1.391218600360602E+00 -3.515311993219235E+01 -5.602056889533600E-01];
VEarth_J2000    = [-2.978494749858966E+01 -5.482119695038260E+00 1.843295966752478E-05];
VMoon_J2000     = [-2.914141610973326E+01 -6.213103677848060E+00 -1.148803176307656E-02];
VMars_J2000     = [1.171984953605987E+00 2.628323979722416E+01 5.221336703898150E-01];
VCeres_J2000    = [-6.196623898412994E+00 -1.834193788479624E+01 5.778897659018130E-01];
VJupiter_J2000  = [-7.900547720232828E+00 1.114339277066948E+01 1.307023308633424E-01];
VSaturn_J2000   = [-7.421900386834246E+00 6.723930997204315E+00 1.775749426204376E-01];
VUranus_J2000   = [4.646586369459156E+00 4.614774391558801E+00 -4.308124107669631E-02];
VNeptune_J2000  = [4.474858465459663E+00 3.063881605796575E+00 -1.659044011083001E-01];
VPluto_J2000    = [5.269124016493589E+00 -2.669250607326493E+00 -1.250716402199096E+00];
VEris_J2000     = [-3.431877929177132E-01 1.676377476595594E+00 1.504390854069972E+00];

J2000 = cell(13,1);
J2000{1}.body  = {'Sun', muSun, RSun_J2000, VSun_J2000};
J2000{2}.body  = {'Mercury', muMercury, RMercury_J2000, VMercury_J2000};
J2000{3}.body  = {'Venus', muVenus, RVenus_J2000, VVenus_J2000};
J2000{4}.body  = {'Earth', muEarth, REarth_J2000, VEarth_J2000};
J2000{5}.body  = {'Moon', muMoon, RMoon_J2000, VMoon_J2000};
J2000{6}.body  = {'Mars', muMars, RMars_J2000, VMars_J2000};
J2000{7}.body  = {'Ceres', muCeres, RCeres_J2000, VCeres_J2000};
J2000{8}.body  = {'Jupiter', muJupiter, RJupiter_J2000, VJupiter_J2000};
J2000{9}.body  = {'Saturn', muSaturn, RSaturn_J2000, VSaturn_J2000};
J2000{10}.body = {'Uranus', muUranus, RUranus_J2000, VUranus_J2000};
J2000{11}.body = {'Neptune', muNeptune, RNeptune_J2000, VNeptune_J2000};
J2000{12}.body = {'Pluto', muPluto, RPluto_J2000, VPluto_J2000};
J2000{13}.body = {'Eris', muEris, REris_J2000, VEris_J2000};

J2000_muTable = table( ...
    {cell2mat(J2000{1}.body(1,2))}, ...
    {cell2mat(J2000{2}.body(1,2))}, ...
    {cell2mat(J2000{3}.body(1,2))}, ...
    {cell2mat(J2000{4}.body(1,2))}, ...
    {cell2mat(J2000{5}.body(1,2))}, ...
    {cell2mat(J2000{6}.body(1,2))}, ...
    {cell2mat(J2000{7}.body(1,2))}, ...
    {cell2mat(J2000{8}.body(1,2))}, ...
    {cell2mat(J2000{9}.body(1,2))}, ...
    {cell2mat(J2000{10}.body(1,2))}, ...
    {cell2mat(J2000{11}.body(1,2))}, ...
    {cell2mat(J2000{12}.body(1,2))}, ...
    {cell2mat(J2000{13}.body(1,2))}, ...
    'VariableNames', ...
    {cell2mat(J2000{1}.body(1,1)) ...
    cell2mat(J2000{2}.body(1,1)) ...
    cell2mat(J2000{3}.body(1,1)) ...
    cell2mat(J2000{4}.body(1,1)) ...
    cell2mat(J2000{5}.body(1,1)) ...
    cell2mat(J2000{6}.body(1,1)) ...
    cell2mat(J2000{7}.body(1,1)) ...
    cell2mat(J2000{8}.body(1,1)) ...
    cell2mat(J2000{9}.body(1,1)) ...
    cell2mat(J2000{10}.body(1,1)) ...
    cell2mat(J2000{11}.body(1,1)) ...
    cell2mat(J2000{12}.body(1,1)) ...
    cell2mat(J2000{13}.body(1,1))}, ... 
    'RowName', ...
    {'GM (km^3/s^2)'});

disp(J2000_muTable);

J2000_pvTable = table( ...
    {cell2mat(J2000{1}.body(1,3)); cell2mat(J2000{1}.body(1,4))}, ...
    {cell2mat(J2000{2}.body(1,3)); cell2mat(J2000{2}.body(1,4))}, ...
    {cell2mat(J2000{3}.body(1,3)); cell2mat(J2000{3}.body(1,4))}, ...
    {cell2mat(J2000{4}.body(1,3)); cell2mat(J2000{4}.body(1,4))}, ...
    {cell2mat(J2000{5}.body(1,3)); cell2mat(J2000{5}.body(1,4))}, ...
    {cell2mat(J2000{6}.body(1,3)); cell2mat(J2000{6}.body(1,4))}, ...
    {cell2mat(J2000{7}.body(1,3)); cell2mat(J2000{7}.body(1,4))}, ...
    {cell2mat(J2000{8}.body(1,3)); cell2mat(J2000{8}.body(1,4))}, ...
    {cell2mat(J2000{9}.body(1,3)); cell2mat(J2000{9}.body(1,4))}, ...
    {cell2mat(J2000{10}.body(1,3)); cell2mat(J2000{10}.body(1,4))}, ...
    {cell2mat(J2000{11}.body(1,3)); cell2mat(J2000{11}.body(1,4))}, ...
    {cell2mat(J2000{12}.body(1,3)); cell2mat(J2000{12}.body(1,4))}, ...
    {cell2mat(J2000{13}.body(1,3)); cell2mat(J2000{13}.body(1,4))}, ...
    'VariableNames', ...
    {cell2mat(J2000{1}.body(1,1)) ...
    cell2mat(J2000{2}.body(1,1)) ...
    cell2mat(J2000{3}.body(1,1)) ...
    cell2mat(J2000{4}.body(1,1)) ...
    cell2mat(J2000{5}.body(1,1)) ...
    cell2mat(J2000{6}.body(1,1)) ...
    cell2mat(J2000{7}.body(1,1)) ...
    cell2mat(J2000{8}.body(1,1)) ...
    cell2mat(J2000{9}.body(1,1)) ...
    cell2mat(J2000{10}.body(1,1)) ...
    cell2mat(J2000{11}.body(1,1)) ...
    cell2mat(J2000{12}.body(1,1)) ...
    cell2mat(J2000{13}.body(1,1))}, ... 
    'RowName', ...
    {'Position (km)', 'Velocity (km/s)'});

disp(J2000_pvTable);
%% Gather input:


n = input('Enter number of bodies:\n');
runTime = input('\nEnter propagation time (in seconds, multiple of 60):\n');
stepSize = input('\nEnter step size (in seconds):\n');
entries = runTime/stepSize;
bodyData = cell(n, 1);
ephemerides = cell(n+1,1);
trajectoryData = cell(n, 1);
angularMomentumData = cell(2,1);
angMom_curr = [0 0 0];

h = stepSize;

for index = 1:n
    bodyData{index}.bodyName = input('\n\nEnter body name:\n', 's');
    bodyData{index}.muBody   = input('\nEnter gravitational parameter:\n');
    bodyData{index}.RBody    = input('\nEnter initial barycentric position:\n');
    bodyData{index}.VBody    = input('\nEnter initial barycentric velocity:\n');
    bodyData{index}.output   = zeros(entries, 3);
    bodyData{index}.f_n        = 0;
end % end for

for manyIndices = 1:n
    ephemerides{manyIndices}.R_ij = [0 0 0];
    ephemerides{manyIndices}.r_ij = [0 0 0];
end

ephemerides{(n+1)}.barycenterOLD = [0 0 0]; % barycenter vector relative to initial barycenter
ephemerides{(n+1)}.barycenterNEW = [0 0 0];
ephemerides{(n+1)}.VBarycenter = [0 0 0];

for indices = 1:n
    trajectoryData{indices}.nthTrajectory = zeros((runTime/stepSize), 3);
end

angularMomentumData{1} = zeros(runTime/stepSize, 3);
angularMomentumData{2} = zeros(runTime/stepSize, 1);

bigMu = 0;

for s = 1:n
    bigMu = bigMu + bodyData{s}.muBody;
end

for indicesAnywhere = 1:n
    ephemerides{(n+1)}.barycenterOLD = ephemerides{(n+1)}.barycenterOLD + (bodyData{indicesAnywhere}.muBody * bodyData{indicesAnywhere}.RBody);
end
ephemerides{(n+1)}.barycenterOLD = (1 / bigMu) * ephemerides{(n+1)}.barycenterOLD;

%% Processing:

for index = 1:stepSize:(runTime + 1)
%{
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
   %} 
    if (index ~= 1)
        % Get radial vectors, distances between attracting bodies:
        for t = 1:n
            bodyData{t}.f_n = 0;
        end
        
        for i = 1:n
            for j = 1:n
                if (strcmp(bodyData{j}.bodyName,bodyData{i}.bodyName) == 0) % 1 true, 0 false
                    ephemerides{i}.R_ij = bodyData{j}.RBody - bodyData{i}.RBody;  % radial vector R_ij
                    ephemerides{i}.r_ij = norm(ephemerides{i}.R_ij); % radial distance r_ij
                    bodyData{i}.f_n     = bodyData{i}.f_n + (bodyData{j}.muBody/(ephemerides{i}.r_ij^3)) * ephemerides{i}.R_ij;
                end % end if
            end % end for
        end % end for
        
        % current position relative to old barycenter, acceleration
        % applied, v_n now known:
        
        for k = 1:n   
            % Velocity Verlet
            bodyData{k}.VBody = bodyData{k}.VBody + (h/2) * bodyData{k}.f_n; % v_n+1/2
            % v_n+1/2 now known, moving on to position:
            bodyData{k}.RBody = bodyData{k}.RBody + h * bodyData{k}.VBody; % y_n+1
            % new position relative to old barycenter known
        end % end for

        % Get new barycenter position, velocity (accurate because vBary 
        % is unaccelerated and thus a straight line):
        ephemerides{(n+1)}.barycenterNEW = [0 0 0];
        for l = 1:n      
            ephemerides{(n+1)}.barycenterNEW = ephemerides{(n+1)}.barycenterNEW + bodyData{l}.muBody * bodyData{l}.RBody;
        end % end for
        
        ephemerides{(n+1)}.barycenterNEW = (1 / bigMu) * ephemerides{(n+1)}.barycenterNEW;
        ephemerides{(n+1)}.VBarycenter = (  ephemerides{(n+1)}.barycenterNEW - ephemerides{(n+1)}.barycenterOLD ) / h;

        
        % New coordinate frame established around new barycenter.
        % Get states of bodies in this new frame:
        for indicial = 1:n
            bodyData{indicial}.RBody = bodyData{indicial}.RBody - ephemerides{(n+1)}.barycenterNEW;
        end
        
        % Prepare for next step in trajectory:
        for w = 1:n
            bodyData{w}.f_n = 0;
        end % end for
        
        for i = 1:n
            for j = 1:n
                if (strcmp(bodyData{j}.bodyName, bodyData{i}.bodyName) == 0)
                    ephemerides{i}.R_ij = bodyData{j}.RBody - bodyData{i}.RBody;  % radial vector R_ij
                    ephemerides{i}.r_ij = norm(ephemerides{i}.R_ij); % radial distance r_ij
                    bodyData{i}.f_n     = bodyData{i}.f_n + (bodyData{j}.muBody/(ephemerides{i}.r_ij^3)) * ephemerides{i}.R_ij;
                end % end if
            end % end for
        end % end for

        % n+1th velocity, no change in barycenter:
        
        for m = 1:n
            bodyData{m}.VBody = bodyData{m}.VBody + (h/2) * bodyData{m}.f_n;
        end % end for
        
        ephemerides{(n+1)}.barycenterOLD = ephemerides{(n+1)}.barycenterNEW;
        ephemerides{(n+1)}.barycenterNEW = [0 0 0];        
        for indicesEverywhere = 1:n
            ephemerides{(n+1)}.barycenterNEW = ephemerides{(n+1)}.barycenterNEW + bodyData{indicesEverywhere}.muBody * bodyData{indicesEverywhere}.RBody;
        end
        ephemerides{(n+1)}.barycenterNEW = (1 / bigMu) * ephemerides{(n+1)}.barycenterNEW;

        % Update barycentric velocity of bodies
        if (index ~= 1)
            for o = 1:n
                bodyData{o}.VBody = bodyData{o}.VBody - ephemerides{(n+1)}.VBarycenter;
            end  % end for
        end % end if

        ephemerides{(n+1)}.barycenterOLD = ephemerides{(n+1)}.barycenterNEW;
    end % end if
        
    % Get trajectory points:
    if (index == 1)
        for p = 1:n
            trajectoryData{p}.nthTrajectory(index,1) = bodyData{p}.RBody(1,1);
            trajectoryData{p}.nthTrajectory(index,2) = bodyData{p}.RBody(1,2);
            trajectoryData{p}.nthTrajectory(index,3) = bodyData{p}.RBody(1,3);
        end % end for
        for p = 1:n
            angMom_curr = angMom_curr + cross(bodyData{p}.RBody, bodyData{p}.VBody);
        end % end for
        angularMomentumData{1}(index, :) = angMom_curr;
        angularMomentumData{2} = norm(angMom_curr);
    else
        for q = 1:n
            trajectoryData{q}.nthTrajectory((index - 1)/stepSize + 1,1) = bodyData{q}.RBody(1,1);
            trajectoryData{q}.nthTrajectory((index - 1)/stepSize + 1,2) = bodyData{q}.RBody(1,2);
            trajectoryData{q}.nthTrajectory((index - 1)/stepSize + 1,3) = bodyData{q}.RBody(1,3);            
        end % end for
        angMom_curr = [0 0 0];
        for q = 1:n
            angMom_curr = angMom_curr + cross(bodyData{q}.RBody, bodyData{q}.VBody);
        end % end for
        angularMomentumData{1}((index-1)/stepSize + 1, :) = angMom_curr;
        angularMomentumData{2} = norm(angMom_curr);
    end % end if   
    
    
    
end % End processing loop

%% Plotting:
Names = cell(1, n); %For legend

for finalIndex = 1:n
    Names{1, finalIndex} = bodyData{finalIndex}.bodyName;
end

figure('Name', 'N-Body Propagator Solution (By: TVW)');
for moreIndices = 1:n
    plot3(trajectoryData{moreIndices}.nthTrajectory(:,1), ...
        trajectoryData{moreIndices}.nthTrajectory(:,2), ...
        trajectoryData{moreIndices}.nthTrajectory(:,3));
    %{
    text( ...
        trajectoryData{moreIndices}.nthTrajectory(1,1), ...
        trajectoryData{moreIndices}.nthTrajectory(1,1), ...
        trajectoryData{moreIndices}.nthTrajectory(1,1), ...
        bodyData{moreIndices}.bodyName);
    %}
    hold on;
end
grid on;
axis equal;
title({'Barycentric-Frame N-Body Propagator', 'Leapfrog Algorithm - Kasdin & Paley (2011), p. 656'});
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
legend(Names{:});

time = [0:stepSize:runTime];

figure('Name', 'Separation of Total Angular Momentum');
subplot(4,1,1);
plot(time(:), angularMomentumData{2}(:,1));
subplot(4,1,2);
plot(time(:), angularMomentumData{1}(:,1));
subplot(4,1,3);
plot(time(:), angularMomentumData{1}(:,2));
subplot(4,1,4);
plot(time(:), angularMomentumData{1}(:,3));
sgtitle('Separation of Total Angular Momentum');
xlabel('Time (s)');