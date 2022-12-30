%% Batch estimation
function initialState_nonlinear = get_batch_nonlinear_least_squares(datasetName, obsDim, stateDim, estimate)

    %% PREPARE DATA
    dataset = load(datasetName);
    N = length(dataset.meas);

    %% Unweighted Nonlinear Least Squares
    
    % In book nomenclature,
    % Dimensions are as follows:
    % yObs,i: 2 x 1,  yBigObs:            (2*11) x 1
    % yPred,i: 2 x 1, yBigPred:           (2*11) x 1
    %                 epsilon (residual): (2*11) x 1 
        
    % 1. Set up initial Newton-Raphson input, x_0:
    x_n = estimate;
    
    % Transform to new (radial) observation space:
    rho_0 = sqrt(x_n(1)^2 + x_n(2)^2); % 2-D Range
    rhoHat_0 = x_n(1:2) / rho_0; % 2-D unit vector
    rhoDot_0 = dot(x_n(3:4), rhoHat_0); % 2-D Range Rate
    rho_n = zeros(N,1);
    
    % Compute initial RMS value of chi-squared:
    lastRMS = 10;
    deltaRMS = 1;
    
    % Let's get loopy!
    tolerance = 1*1e-13;
    nthState = [0, 0, 0, 0]; % 2-D state vector
    counter = 0;
    
    while (deltaRMS > tolerance)
    
        % Loop counter:
        counter = counter + 1;
    
        % Transform prediction to states:
        statePred = bigA*x_n;
    
        for n = 1:N
            
            % Get nth coordinate:
            x_i = statePred(2*n-1);
            y_i = statePred(2*n);
    
            % Compute nth predicted range:
            rho_n(n,1) = ...
                sqrt(x_i^2 + y_i^2);
    
            % 4. Compute nth residual:
            res_n(n,1) = a9p2.meas(n) - rho_n(n);
    
            % 5. Get nth Jacobian (little model):
            c = cos(meanMotion*a9p1.times(n));
            s = sin(meanMotion*a9p1.times(n));
            drhodx =  (x_i*(4 - 3*c) + y_i*6*(s - meanMotion*a9p1.times(n))) ...
                / rho_n(n);
            drhody = (y_i) ...
                / rho_n(n);
            drhodxdot = (x_i*(s/meanMotion) + y_i*(-2/meanMotion)*(1 - c)) ...
                / rho_n(n);
            drhodydot = (x_i*(2/meanMotion)*(1 - c) + y_i*((4*s/meanMotion) - 3*a9p1.times(n))) ...
                / rho_n(n);
            bigH(n, 1:4) = [drhodx, drhody,  drhodxdot, drhodydot];
    
        end
    
        % Prepare output state vector:
        nthState = x_n;
    
        % 6. Solve for dX by dX = inv(A^T*A)*A^T*nEpsilon
        P = bigH'*bigH;
        dx = inv(P)*bigH'*res_n;
    
        % 7. Update Gauss-Newton Method, (n+1)X = nX -(-dX)
        x_n = x_n + dx;
    
        % 8. Check for convergence w/ eqn. 9.38.  If not, return to 3).  
    
        % Get current RMS:
        nthChiSquared = 0;
        nthChiSquared = sum(res_n.^2);
        thisRMS = sqrt( nthChiSquared / N);
    
        % Find difference in RMS values:
        deltaRMS = abs(thisRMS - lastRMS);
    
        % update:
        lastRMS = thisRMS;
    
    
    end
    
    
    %% PRINT RESULTS TO CONSOLE:
    
    fprintf('Problem 1 - Unweighted Linear Least-Squares Estimation\n\n');
    fprintf('     The predicted initial state is:\n');
    fprintf('        x_0    = %0.2f   m\n', estimate(1,1));
    fprintf('        y_0    = %0.2f   m\n', estimate(2,1));
    fprintf('        xDot_0 =   %0.6f m/s\n', estimate(3,1));
    fprintf('        yDot_0 =  %0.5f  m/s\n\n\n', estimate(4,1));
    
    fprintf('%s\n\n', printing.lineBreak);
    
    fprintf('  Problem 2 - Nonlinear Least-Squares Solution\n\n');
    fprintf('     The predicted initial state after %d iterations is:\n\n', ...
        counter);
    fprintf('        x_0    = %0.3f m\n', nthState(1));
    fprintf('        y_0    = %0.3f m\n', nthState(2));
    fprintf('        xdot_0 = %0.7f m/s\n', nthState(3));
    fprintf('        ydot_0 = %0.6f m/s\n\n', nthState(4));
    fprintf('             with an rms range chi-squared of %0.4f m \n\n', lastRMS);
    
end
