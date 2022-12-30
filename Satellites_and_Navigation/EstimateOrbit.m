%% ESTIMATION

% Setup workspace, printing, text formatting:
run('TidyWorkspace.m');
run("CheckSNaG.m");
heading = HeadingGenerator(...
    'ENAE441', 'Final Project', 'Pt.2 - Orbit Estimation');
printing = PrintFormatting();

% Get global constants:
consts = OrbitConstants();
chemConsts = ChemicalConstants();

waitMsgESIMATION = "Beginning Orbit Estimation (3.2).  Please wait...";
xWaitESIMATION   = 0;
waitBarESTIMATION = waitbar(xWaitESIMATION, waitMsgESIMATION);


%% DATA / STATION / ORBIT PREPARATION:

xWaitESIMATION = 0.05;
waitMsgESTIMATION = "Loading data, force models, etc...";
waitbar(xWaitESIMATION, waitBarESTIMATION, waitMsgESTIMATION);

% Load 1st & 2nd nights:
night1 = load('opt2satCset4.mat');
night2 = load('opt3satCset4.mat');

% Get desired properties (AZ, EL, Datetimes, Obs#))
ignoredProperties = [...
    "right_ascension_deg",
    "declination_deg",
    "site_latitude_deg",
    "site_longitude_deg",
    "site_altitude_m"];

% Full observation set:
observationSet1 = select_columns(night1.opt2satCset4, ignoredProperties, true);
observationSet1.datetime = datetime(observationSet1.datetime, 'TimeZone','UTC');
observationSet2 = select_columns(night2.opt3satCset4, ignoredProperties, true);
observationSet2.datetime = datetime(observationSet2.datetime, 'TimeZone','UTC');

% Observation subset (N = 25):
obsSubsetDim = size(observationSet1);
earlyObsNight1 = observationSet1(1:4:100, :);
earlyObsNight2 = observationSet2(1:4:100, :);
lateObsNight1 = observationSet1(774:4:874, :);
%od_earlyObsNight1 = determine_orbit();

tInitial_N1 = observationSet1.datetime(1);
tInitial_N2 = observationSet2.datetime(1);
tFinal_N2 = observationSet2.datetime(end);

% Set up observation station object for OAP, Chile:

% Observation site LLA for OAP, Chile:
% Double-precision values (10 sig figs)
chileLLA.latitude_deg  = -30.1428030000;
chileLLA.longitude_deg = -70.6945280000;
chileLLA.altitude_m    =  1500.000000;
chileName = "OAP-Chile";

chileLLA.epoch = earlyObsNight1.datetime(1);

% Station tool from SNaG-app:
oapChile = make_station(chileName, ...
    chileLLA.latitude_deg, ...
    chileLLA.longitude_deg, ...
    chileLLA.altitude_m);


% Set initial orbit determination:
obsVector = table2array(night1.opt2satCset4(:,1));
goodIndex1 = 101; %find(obsVector == 2);
goodIndex2 = 176; %find(obsVector == 48);
goodIndex3 = 251; %find(obsVector == 319);

goodObs1 = night1.opt2satCset4(goodIndex1,:);
goodObs2 = night1.opt2satCset4(goodIndex2,:);
goodObs3 = night1.opt2satCset4(goodIndex3,:);

initialOrbit = orbit_checker( ...
    goodObs1,...
    goodObs2,...
    goodObs3,...
    chileLLA,...
    consts.mu_earth_km);


xWaitESIMATION = 0.1;
waitMsgESTIMATION = "Calculating perturbations.";
waitbar(xWaitESIMATION, waitBarESTIMATION, waitMsgESTIMATION);

%% MODELS:

% Mass:
satMass_kg = 1000; % kg
hispaSat3Wet = 6200; % Launch mass of Amazonas sat
hispaSat3Dry = 2819; % Dry mass of ''
echoStarWet = 2885; % Launch mass of EchoStar II
echoStarDry = 2000; % dry mass of EchoStar II
smallMass = 500; % 0.5 mt
avgMass = 1000; % 1 mt
bigMass = 2000; % 2 mt


% FORCES:
forceModel_twoBody  = constants.force_twobody;

% Gravity perturbations only:
forceModel_2x0      = force_model(2,0,0,0,0,0,satMass_kg);
forceModel_2x2      = force_model(2,2,0,0,0,0,satMass_kg);
forceModel_20x20    = force_model(20,20,0,0,0,0,satMass_kg);
forceModel_40x40    = force_model(40,40,0,0,0,0,satMass_kg);

% Solar Radiation Pressure:

% Reflectivity Coefficients:
lowCRef = 0.5;
avgCRef = 1.2;
highCRef = 1.5;

% Area:
echostarArea = 20;

smallSatArea = 0.5;
medSatArea   = 1;
largeSatArea = 10;

craftRadArea = 20;

% set:  [reflectiveArea, reflectivityCoefficient, s/c mass]
srp_set1 = [1, lowCRef,  satMass_kg];
srp_set2 = [1, avgCRef,  satMass_kg];
srp_set3 = [1, highCRef, satMass_kg];

% Luni-Solar 4BP:
%sunCode = bodycode('Sun');
sunCode = 10;
%moonCode = bodycode('Luna');
moonCode = 301;


%% DETERMINE ORBIT:

waitMsgESTIMATION = 'Determining orbit...';
xWaitESIMATION = 0.3;
waitbar(xWaitESIMATION, waitBarESTIMATION, waitMsgESTIMATION);

initialPVT_N1 = pvt(goodObs2.datetime, initialOrbit.R2, initialOrbit.V2);
estimationModel = force_model(20,20,0,0,1,1,1000,0,0);
earlyOD_N1 = determine_orbit(...
    initialPVT_N1, oapChile, earlyObsNight1, estimationModel);
estimatedPVT_N1 = earlyOD_N1.estimated;
tic
estimatedOrbit = propagate(...
    estimatedPVT_N1, tInitial_N2, tFinal_N2, 30, estimationModel, 0);
estimationTime = toc;
interpolatedNight2 = ephemeris_interp(estimatedOrbit, observationSet2.datetime);

%% DATA PROCESSING:

waitMsgESTIMATION = 'Converting interpolated PVT data to AER...';
aerCounter = 0;
xWaitESIMATION = 0.6;
waitbar(xWaitESIMATION, waitBarESTIMATION, waitMsgESTIMATION);
pause(1);

% GET RMS VALUES:

N = length(night2.opt3satCset4.datetime);
interpolatedStates.position_m = interpolatedNight2.position_m;
interpolatedStates.epoch = interpolatedNight2.epoch;

tic
% Return to AER Frame:
for n = 1:N
    waitMsgESTIMATION = sprintf('AER:  %d/%d', n, N);
    waitbar(xWaitESIMATION, waitBarESTIMATION, waitMsgESTIMATION);
    chileLLA.epoch = interpolatedNight2.epoch(n);
    thisAER = eci_to_azelrn(interpolatedStates.epoch(n), interpolatedStates.position_m(n, 1:3), chileLLA);
    interpolatedNight2.aer(n) = thisAER;
end

waitMsgESTIMATION = 'Calculating residuals/RMS from night 2...';
xWaitESIMATION = 0.3;
waitbar(xWaitESIMATION, waitBarESTIMATION, waitMsgESTIMATION);

for n = 1:N
    waitMsgESTIMATION = sprintf('Residual:  %d/%d', n, N);
    waitbar(xWaitESIMATION, waitBarESTIMATION, waitMsgESTIMATION);
    currentAZobs  = night2.opt3satCset4.azimuth_deg(n);
    currentELobs  = night2.opt3satCset4.elevation_deg(n);
    currentAZpred = interpolatedNight2.aer(n).azimuth_deg;
    currentELpred = interpolatedNight2.aer(n).elevation_deg;
    interpResidualDistribution.datetime(n) = night2.opt3satCset4.datetime(n);
    interpResidualDistribution.azimuth_deg(n) = currentAZobs - currentAZpred;
    interpResidualDistribution.elevation_deg(n) = currentELobs - currentELpred;
    interpResidualDistribution.RMS(n) = ...
        (interpResidualDistribution.azimuth_deg(n))^2 ...
        + (interpResidualDistribution.elevation_deg(n))^2;
end

totalRMS = sqrt((1/N)*sum(interpResidualDistribution.RMS));

processingTime = toc;

%% GRAPHING/VISUALS:

figure('Name', 'Histogram of RMS Values for Night 2 data vs. Night 1 estimate');
histogram(interpResidualDistribution.RMS);
xlabel('Total Residual Error (deg.)');
ylabel('Count (# observations)');
title('Night 2 Observations vs. Predicted Night 2 Values');

figure('Name', 'Azimuth Residuals vs. Time (Night 2 Obs., Night 2 Pred.)');
plot(interpResidualDistribution.azimuth_deg(:), 'Color', [0.4940 0.1840 0.5560]);
xlabel('Time (Julian Date)');
ylabel('Azimuthal Residual (deg.)');
title('Night 2 Observations vs. Predicted Night 2 Values from Night 1 Determination, ephemeris_interp()');

figure('Name', 'Elevation Residuals vs. Time (Night 2 Obs., Night 2 Pred.)');
plot(interpResidualDistribution.elevation_deg(:), 'Color', [0.4660 0.6740 0.1880]);
xlabel('Time (Julian Date)');
ylabel('Azimuthal Residual (deg.)');
title('Night 2 Observations vs. Predicted Night 2 Values from Night 1 Determination, ephemeris_interp()');

figure('Name', 'RMS Values over time');
plot(second(interpResidualDistribution.datetime(:), 'secondofday'), interpResidualDistribution.RMS(:));
xlabel('Time (Julian Date)');
ylabel('Total RMS Error (deg.)');
title('Total RMS Over the Course of Propagated Night 2 Values');


%% Check for maneuver

% Set initial orbit determination:
obsVector_N2 = table2array(night2.opt3satCset4(:,1));
goodIndex1_N2 = find(obsVector_N2 == 613);
goodIndex2_N2 = find(obsVector_N2 == 1091);
goodIndex3_N2 = find(obsVector_N2 == 1186);

goodObs1_N2 = night2.opt3satCset4(goodIndex1_N2,:);
goodObs2_N2 = night2.opt3satCset4(goodIndex2_N2,:);
goodObs3_N2 = night2.opt3satCset4(goodIndex3_N2,:);

initialOrbit_N2 = orbit_checker( ...
    goodObs1_N2,...
    goodObs2_N2,...
    goodObs3_N2,...
    chileLLA,...
    consts.mu_earth_km);

initialPVT_N2 = pvt(goodObs2_N2.datetime, initialOrbit_N2.R2, initialOrbit_N2.V2);
estimationModel_N2 = force_model(20,20,0,0,1,1,1000,0,0);
earlyOD_N2 = determine_orbit(...
    initialPVT_N2, oapChile, earlyObsNight2, estimationModel_N2);
estimatedPVT_N2 = earlyOD_N2.estimated;
estimatedOrbit_N1 = propagate(...
    estimatedPVT_N1, tInitial_N1, tFinal_N2, 30, estimationModel, 0);
tic
estimatedOrbit_N2 = propagate(...
    estimatedPVT_N2, tInitial_N1, tFinal_N2, 30, estimationModel, 0);
estimatedOrbit_N2.epoch(:) = datetime(estimatedOrbit_N2.epoch(:), 'TimeZone', 'UTC');
estimationTime_Conjunction = toc;

%interpolatedNight2 = ephemeris_interp(estimatedOrbit, observationSet2.datetime);

%conjunctionTable = conjunction(...
%    estimatedPVT_N1, estimatedPVT_N2, tInitial_N1, tFinal_N2);

%% PRINT TO CONSOLE:
fprintf('\n\n     Time to Estimation Convergence:  tConverge = %0.2f sec \n\n', estimationTime);
fprintf('     Time for AER, Residual Calculation: tProcessing = %0.2f sec\n\n', processingTime);
fprintf('    The final RMS value is:  %0.8f degrees\n\n', totalRMS);