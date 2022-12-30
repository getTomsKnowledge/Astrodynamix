function results = orbit_checker(obs1, obs2, obs3, siteLLA, mu_m)

    muBody2 = mu_m*1e9;

    siteLLA.epoch = obs1.datetime;
    site1 = eci(siteLLA);
    site1.latitude_deg = siteLLA.latitude_deg;
    site1.longitude_deg = siteLLA.longitude_deg;
    site1.altitude_m = siteLLA.altitude_m;

    siteLLA.epoch = obs2.datetime;
    site2 = eci(siteLLA);
    site2.latitude_deg = siteLLA.latitude_deg;
    site2.longitude_deg = siteLLA.longitude_deg;
    site2.altitude_m = siteLLA.altitude_m;

    siteLLA.epoch = obs3.datetime;
    site3 = eci(siteLLA);
    site3.latitude_deg = siteLLA.latitude_deg;
    site3.longitude_deg = siteLLA.longitude_deg;
    site3.altitude_m = siteLLA.altitude_m;
    
    % Get line of sight vectors:
    los1 = get_topo_los(...
        obs1.azimuth_deg, obs1.elevation_deg, site1);
    los2 = get_topo_los(...
        obs2.azimuth_deg, obs2.elevation_deg, site2);
    los3 = get_topo_los(...
        obs3.azimuth_deg, obs3.elevation_deg, site3);

    %{
    Apparent Right Ascension / Declination technique:
    [los1, los2, los3] = get_los( ...
        obs1.right_ascension_deg, obs1.declination_deg, ...
        obs2.right_ascension_deg, obs2.declination_deg, ...
        obs3.right_ascension_deg, obs3.declination_deg);
    %}


    t1 = second(site1.epoch, 'secondofday');
    t2 = second(site2.epoch, 'secondofday');
    t3 = second(site3.epoch, 'secondofday');

    [R1Out, R2Out, R3Out] = gauss_solver(...
        los1, los2, los3, ...
        site1.position_m, site2.position_m, site3.position_m, ...
        t1, t2, t3);

    [V1, V2, V3] = gibbs_solver(...
        R1Out, R2Out, R3Out, muBody2);

    dayRange = [0 24*60*60];

    state1 = [R1Out;V1];
    [tOut1, xOut1] = ode45(@propagate_2BP, dayRange, state1, [], muBody2);

    state2 = [R2Out;V2];
    [tOut2, xOut2] = ode45(@propagate_2BP, dayRange, state2, [], muBody2);

    state3 = [R3Out;V3];
    [tOut3, xOut3] = ode45(@propagate_2BP, dayRange, state3, [], muBody2);

    % Process results:
    result1.position_m = R1Out;
    result1.epoch = site1.epoch;
    results.R1 = R1Out;
    results.V1 = V1;
    result2.position_m = R2Out;
    result2.epoch = site2.epoch;
    results.R2 = R2Out;
    results.V2 = V2;
    result3.position_m = R3Out;
    result3.epoch = site3.epoch;
    results.R3 = R3Out;
    results.V3 = V3;
    
    result1.aer = aer(result1, site1);
    result2.aer = aer(result2, site2);
    result3.aer = aer(result3, site3);
    
    az1res = obs1.azimuth_deg - result1.aer.azimuth_deg;
    el1res = obs1.elevation_deg - result1.aer.elevation_deg;
    
    az2res = obs2.azimuth_deg - result2.aer.azimuth_deg;
    el2res = obs2.elevation_deg - result2.aer.elevation_deg;
    
    az3res = obs3.azimuth_deg - result3.aer.azimuth_deg;
    el3res = obs3.elevation_deg - result3.aer.elevation_deg;
    
    RMS = (az1res)^2 + (el1res)^2; 
    RMS = (az2res)^2 + (el2res)^2;
    RMS = (az3res)^2 + (el3res)^2;
    results.RMS = sqrt( RMS/3 );

    rmsResultsString = sprintf('RMS = 0.4f deg. for <%d, %d, %d>', ...
        results.RMS, obs1.observation_number, obs2.observation_number, obs3.observation_number);
    results.r1 = norm(R1Out);
    results.r2 = norm(R2Out);
    results.r3 = norm(R3Out);

    r1Label = sprintf('r1 = %0.3f (km)', results.r1*1e-3);
    r2Label = sprintf('r2 = %0.3f (km)', results.r2*1e-3);
    r3Label = sprintf('r3 = %0.3f (km)', results.r3*1e-3);

    %{
    figure('Name', 'Orbit Check Plot');
    plot3(xOut1(:,1), xOut1(:,2), xOut1(:,3));
    hold on;
    plot3(xOut2(:,1), xOut2(:,2), xOut2(:,3));
    plot3(xOut3(:,1), xOut3(:,2), xOut3(:,3));
    title('Orbit Check');
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title(rmsResultsString);
    legend(r1Label, r2Label, r3Label);
    %}

end