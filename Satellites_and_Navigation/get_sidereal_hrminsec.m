%% Sidereal Time Converter - Tom West, 2022

% datetimeIn - Local datetime in UTC for comparison
% lstDegIn   - The local sidereal time in degrees
% timeOut    - The converted lst in hrs/mins/secs

% Prints steps of computation

% Demo:  
% nowUTC = nowutc(); now_lst = siderealtime(nowUTC, stddef.umd_lla);
% timeOut = get_sidereal_hrminsec(nowUTC, now_lst) 


function timeOut_lst = get_sidereal_hrminsec(datetimeIn, lstDegIn_deg)
        
    % Convert to hrs mins sec with datetime():
    timeHMS = [deg2rad(lstDegIn_deg(1)) deg2rad(lstDegIn_deg(2))];
    timeHMS = [(timeHMS(1) / (2*pi))  (timeHMS(2) / (2*pi))];
    t = timeHMS(1);
    hrs = floor(t * 24);
    t = (t*24 - hrs)/24;
    mins = floor((t*24) * 60);
    t = (t*24*60 - mins)/(24*60);
    secs = t*24*60*60;
    timeOut_lst = sprintf('%d hrs %d mins %0.4f sec [sidereal]', hrs, mins, secs);
    
    % Output results:
    fprintf('     The local sidereal time at %s UTC is:  \n       %s\n', datetimeIn, timeOut_lst);
    fprintf('       (-OR- %0.3f deg)\n\n', lstDegIn_deg(1));

end