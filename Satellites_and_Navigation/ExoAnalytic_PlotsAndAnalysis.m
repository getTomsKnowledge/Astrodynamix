%% INITIAL DATA ANALYSIS:

run('TidyWorkspace.m');
run('CheckSNaG.m');
printing = PrintFormatting();
heading = HeadingGenerator('ENAE441', 'Orbit Determination Project', 'Determine GEO orbit based on observation data');
heading.print_heading;

consts = OrbitConstants();

%% USER INPUT:

% Select from multiple group data sets for comparison:
chosenSet = input("Type '12' for group 12, '1' for Group 1, '11' for Group 11.");
%chosenSet = str2num(groupInput);
switch chosenSet
    case 1
        selectedGroupN1 = 'opt2satAset1.mat';
        selectedGroupN2 = 'opt2satAset1.mat';
    case 11
        selectedGroupN1 = 'opt2satCset3.mat';
        selectedGroupN2 = 'opt2satCset3.mat';
    case 12
        selectedGroupN1 = 'opt2satCset4.mat';
        selectedGroupN2 = 'opt3satCset4.mat';
    otherwise
        fprintf('INVALID ENTRY!');
        return;
end

night1Data = load(selectedGroupN1);
night2Data = load(selectedGroupN2);

%% DATA PREPARATION:

% Convert datetimes to Julian Date:
t0_NightOne = juliandate(night1Data.opt2satCset4.datetime(1));
times_NightOne = juliandate(night1Data.opt2satCset4.datetime(:)) - t0_NightOne;
t0_NightTwo = juliandate(night2Data.opt3satCset4.datetime(1));
times_NightTwo = juliandate(night2Data.opt3satCset4.datetime(:)) - t0_NightTwo;

% Get AZ, EL Data:
obsAZ_NightOne = night1Data.opt2satCset4.azimuth_deg;
obsEL_NightOne = night1Data.opt2satCset4.elevation_deg;
obsAZ_NightTwo = night2Data.opt3satCset4.azimuth_deg;
obsEL_NightTwo = night2Data.opt3satCset4.elevation_deg;

%% GRAPHING:

% Build Figures

% Night One:
figure('Name', 'Night One: October 29th, 2020 UTC - OAP Chile, 299*E Cluster');
subplot(2,1,1);
scatter(times_NightOne(:), obsAZ_NightOne(:), 5, [0.4940 0.1840 0.5560], 'Marker', '*');
xlabel('Times [Julian Date Since Epoch]');
ylabel('Azimuth, AZ (deg.)');
title({'299*E Cluster Sat C Azimuth vs. Time (10/29/2020, OAP Chile)', 'Courtesy: ExoAnalytic Solutions, 11/10/2020'});
subplot(2,1,2);
scatter(times_NightOne(:), obsEL_NightOne(:), 5, [0.4660 0.6740 0.1880], 'Marker', '*');
xlabel('Times [Julian Date Since Epoch]');
ylabel('Elevation, EL (deg.)');
title({'299*E Cluster Sat C Elevation vs. Time (10/29/2020, OAP Chile)','Courtesy: ExoAnalytic Solutions, 11/10/2020'});

% Night Two:
figure('Name', 'Night One: October 29th, 2020 UTC - OAP Chile, 299*E Cluster');
subplot(2,1,1);
scatter(times_NightTwo(:), obsAZ_NightTwo(:), 5, [0.4940 0.1840 0.5560], 'Marker', '*');
xlabel('Times [Julian Date Since Epoch]');
ylabel('Azimuth, AZ (deg.)');
title({'299*E Cluster Sat C Azimuth vs. Time (10/30/2020, OAP Chile)', 'Courtesy: ExoAnalytic Solutions, 11/10/2020'});
subplot(2,1,2);
scatter(times_NightTwo(:), obsEL_NightTwo(:), 5, [0.4660 0.6740 0.1880], 'Marker', '*');
xlabel('Times [Julian Date Since Epoch]');
ylabel('Elevation, EL (deg.)');
title({'299*E Cluster Sat C Elevation vs. Time (10/30/2020, OAP Chile)','Courtesy: ExoAnalytic Solutions, 11/10/2020'})
