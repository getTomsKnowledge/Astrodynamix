%% Determines location of receiver via GPS clock biases, pseudoranges, positions

%{
    % TEST DATA:
    % ECEF satellite positions:
    r1_m = [791.015754, -26491.796133, 195.542196]'*1e3;
    r2_m = [-10034.534619, -12199.565845, 21805.910224]'*1e3;
    r3_m = [6740.994371, -23259.598032, 10425.385636]'*1e3;
    r4_m = [-8056.794719, -22862.235722, 10725.417482]'*1e3;

    % Satellite pseudoranges:
    R1_m = 21954447.398;
    R2_m = 22050447.992;
    R3_m = 20355865.297;
    R4_m = 21249501.531;

    % Satellite clock biases:
    cd1_m = -87.480116*1e3;
    cd2_m =  288.118412*1e3;
    cd3_m = -93.540637*1e3;
    cd4_m = -269.563631*1e3;

%}

function gpsOutput = get_gps_solution(...
    r1_m, r2_m, r3_m, r4_m,...
    R1_m, R2_m, R3_m, R4_m,...
    cd1_m, cd2_m, cd3_m, cd4_m,...
    satName1, satName2, satName3, satName4);

    thisConsts. = OrbitConstants();

    % Vectorize inputs:
    r = [r1_m, r2_m, r3_m, r4_m];    
    R = [R1_m, R2_m, R3_m, R4_m]';    
    satClock = [cd1_m, cd2_m, cd3_m, cd4_m]';
    rho = zeros(4,1);
    for j = 1:4
        rho(j) = norm(r(:,j));
    end


    % Estimation:
    
    % 1. Observer initial position guess:
    xGuess = zeros(4,1); % origin always works
    
    % Set up residual vector:
    epsilon = zeros(4,1);
    
    % Loop control variables:
    diff = 1;
    tol = 1*1e-10;
    chiSqPrevious = 100;
    chiSquared = 10;
    deltaChi = 100;
    
    while (deltaChi > tol)
    
        % 2. Compute residual:
        for n = 1:4
            epsilon(n) = ...
                R(n) - sqrt( ...
                (r(1,n) - xGuess(1))^2 + (r(2,n) - xGuess(2))^2 + (r(3,n) - xGuess(3))^2 ...
                ) + satClock(n) - xGuess(4);
        end
    
        % 3. Check residual:
        % 4. Compute Jacobian:
        A = [(xGuess(1) - r(1,1))/rho(1), (xGuess(2) - r(2,1))/rho(1), (xGuess(3) - r(3,1))/rho(1) 1;
             (xGuess(1) - r(1,2))/rho(2), (xGuess(2) - r(2,2))/rho(2), (xGuess(3) - r(3,2))/rho(2) 1;
             (xGuess(1) - r(1,3))/rho(3), (xGuess(2) - r(2,3))/rho(3), (xGuess(3) - r(3,3))/rho(3) 1;
             (xGuess(1) - r(1,4))/rho(4), (xGuess(2) - r(2,4))/rho(4), (xGuess(3) - r(3,4))/rho(4) 1];
    
        % 5. Newton's method:
        xGuess = xGuess + inv(A)*epsilon;
    
        % 6. Check for convergence:
        chiSqPrevious = chiSquared;
        chiSquared = sum(norm(epsilon(:)))/4;
        prevDiff = diff;
        diff = abs(chiSquared - chiSqPrevious);
        deltaChi = abs(diff - prevDiff);
    
    end
    
    % Assign output struct fields:
    gpsOutput.x_m           = xGuess(1);
    gpsOutput.y_m           = xGuess(2);
    gpsOutput.z_m           = xGuess(3);
    gpsOutput.biasDist_m    = xGuess(4);
    gpsOutput.clockBias_sec = xGuess(4)/thisConsts.c;
    gpsOutput.chiSquared    = chiSquared;
    gpsOutput.tol           = tol;

    % Problem 3
    fprintf('  GPS Results\n\n');
    fprintf('     For a tolerance of %4.0d, the location is:\n', tol);
    fprintf('        ECEF Position:  <%0.3f, %0.3f, %0.3f> [km]\n\n', ...
        xGuess(1)*1e-3, xGuess(2)*1e-3, xGuess(3)*1e-3);
    fprintf('     and a clock bias (in both distance and time):\n');
    fprintf('        Bias Distance: %0.3f km\n', xGuess(4)*1e-3);
    fprintf('        Time Offset:   %0.3f millisec\n\n', (xGuess(4)/thisConsts.c)*1e3);
    fprintf('     Minimized chi-squared:\n');
    fprintf('        chiSq = %6.6d sq. m', chiSquared)
    fprintf('     Satellites Considered:\n');
    fprintf('        %s , %s , %s , %s\n\n', ...
        satName1, satName2, satName3, satName4);

end