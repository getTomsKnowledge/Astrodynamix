%% Sequential Least-Squares Estimation - Kalman Filtering

%{

% TEST DATA:

% Given parameters:
meanMotion = 1*1e-3; % rad/s
timeStep = 4*60;
% Initial position estimate from ELMS announcement:
x0 = [133.2657, 135.92806, 0.06794, -0.26640]';

%}

function estimatedStates = get_sequential_kalman_least_squares(datasetName, obsDim, stateDim, estimate,  qVals, rVals)

    % SET UP DATA:
    dataset = load(datasetName); 

    % State Covariance, Q
    Q = zeros(stateDim);
    for n = 1:stateDim
        Q(n,n) = qVals(n);
    end
    
    % Observation Covariance, R
    R = zeros(obsDim, stateDim)
    for n = 1:obsDim
        R(n,n) = rVals(n);
    end

    N = length(dataset.meas / obsDim);

    % Measurement components:
    times = dataset.times;
    xi = zeros(N,1);
    xii = zeros(N,1);
    rhoi = zeros(N,1);
    for n = 1:N
        xi(n) = dataset.meas(2*n - 1);
        xii(n) = dataset.meas(2*n);
        rhoi(n) = dataset.meas(n);
    end
    
    % HCW curve preparation:
    phiNow = zeros(stateDim);
    propTime = N*timeStep;
    X_hcwCurve = zeros(propTime, stateDim);
    for t = 0:propTime % 100-minute propagation
        phiNow = clohessyWiltshire_transition_matrix(...
            meanMotion, t);
        phiXY = [phiNow(1,1) phiNow(1,2), phiNow(1,4), phiNow(1,5);
            phiNow(2,1) phiNow(2,2), phiNow(2,4), phiNow(2,5);
            phiNow(4,1) phiNow(4,2), phiNow(4,4), phiNow(4,5);
            phiNow(5,1) phiNow(5,2), phiNow(5,4), phiNow(5,5)];
        X_hcwCurve(t+1, :) = transpose(phiXY*x0);
    end
    
    l = 1;
    hcwDiscrete = zeros(N,stateDim);
    for t = 0:timeStep:((N-1)*timeStep)
        phiNow = clohessyWiltshire_transition_matrix(...
            meanMotion, t);
        phiXY = [phiNow(1,1) phiNow(1,2), phiNow(1,4), phiNow(1,5);
            phiNow(2,1) phiNow(2,2), phiNow(2,4), phiNow(2,5);
            phiNow(4,1) phiNow(4,2), phiNow(4,4), phiNow(4,5);
            phiNow(5,1) phiNow(5,2), phiNow(5,4), phiNow(5,5)];
        hcwDiscrete(l, :) = phiXY*x0;
        l = l + 1;
    end
    
    
    % KALMAN FILTER:
    
    % Kalman Filter Preparation:
    
    % Get little state-transition matrix for 4-min propagation, Phi_k:
    bigPhi = clohessyWiltshire_transition_matrix(...
        meanMotion, timeStep);
    kPhi = [bigPhi(1,1) bigPhi(1,2), bigPhi(1,4), bigPhi(1,5);
            bigPhi(2,1) bigPhi(2,2), bigPhi(2,4), bigPhi(2,5);
            bigPhi(4,1) bigPhi(4,2), bigPhi(4,4), bigPhi(4,5);
            bigPhi(5,1) bigPhi(5,2), bigPhi(5,4), bigPhi(5,5)];
    
    % Vector to store Kalman-filter estimated states for graphing:
    X_kalmEstimate = zeros(N,stateDim);
    K_kalmGain = zeros(N,1);
    % Setup:
    % Set initial estimate:
    xHat = x0;
    % Set initial covariance matrix:
    a = 125;
    chiSq = 1000;
    %{
    da = 1;
    diffPrev = 1000;
    tol = 1;
    diff = 1000;
    delta = 1000;
    %}
    P = a*[   R,     zeros(obsDim);
        zeros(obsDim),  eye(obsDim)];
    % Set state-observation transform:
    H = [eye(obsDim), zeros(obsDim)]; % observation and state are direct; no dotted vals
    
    %{
    while (abs(delta) > tol)
    %}
        for n = 1:N
        
            % 1)  PROPAGATE:
            
            % Estimate current state based on previous estimate:
            x_pred = kPhi*xHat;
            % Estimate current covariance based on previous estimate: 
            P_pred = kPhi*P*kPhi';
        
            % 2) UPDATE:
            
            % Compute Kalman gain:
            K_prev = P_pred*H'*inv(H*P_pred*H' + R);
            K_kalmGain(n) = sqrt(norm(K_prev(:,1))^2 + norm(K_prev(:,2))^2);
            % Get nth observation:
            z = dataset.meas((2*n-1) : (2*n));
            % Compute state estimate:
            xHat = x_pred + K_prev*(z - H*x_pred);
            X_kalmEstimate(n,:) = xHat';
            res(n, 1:2) = (z - H*x_pred)';
            % Compute state-obs covariance:
            P = (eye(4) - K_prev*H)*P_pred;
            minP = kPhi*P*kPhi';
        
        end
    
    
        chiSqPrev = chiSq;
        chiSq = norm(res(:,1).^2 + res(:,2).^2) / N;
        %{
        diffPrev = diff;
        diff = chiSq - chiSqPrev;
    
        if (diff < 0)
            a = a + da;
        else
            a = a - da;
            da = da/10;
        end
    
        delta = diff - diffPrev;
    end
        %}
    
    % Prepare plot:
    
    resNorm(:,1) = res(:,1).^2 + res(:,2).^2;
    resNorm(:) = sqrt(resNorm(:));
    
    
    %% PRINT TO CONSOLE
    
    figure('Name', 'Problem 1: Kalman Filter - Tom West, UID 117659399');
    subplot(3,1,1);
    plot(xi(:), xii(:), 'Color', 'r', 'Marker', '+', 'LineStyle', 'none');
    title(...
        {'Kalman-Filtered Estimates of Position for Relative Orbital Motion', printing.myName},...
        'Interpreter', 'latex');
    xlabel('X-Axis Value (m)', 'Interpreter', 'latex');
    ylabel('Y-Axis Value (m)', 'Interpreter', 'latex');
    hold on;
    plot(X_kalmEstimate(:,1), X_kalmEstimate(:,2), 'Color', 'g', 'Marker', 'diamond', 'LineStyle','none');
    plot(X_hcwCurve(:,1), X_hcwCurve(:,2), 'Color', 'b');
    plot(hcwDiscrete(:,1), hcwDiscrete(:,2), 'Color', 'b', 'Marker', 'o', 'LineStyle', 'none');
    legend(...
        'Observed Position, $z_{obs}$ [m]', ...
        'Kalman-Estimated Position, $x_{est}$ [m]', ...
        'H-C-W Solution', 'Interpreter', 'latex');
    subplot(3,1,2);
    plot(0:240:(propTime - timeStep), resNorm(:));
    title('Residual (error) vs. Time', 'Interpreter', 'latex');
    xlabel('Time, t (s)', 'Interpreter', 'latex');
    ylabel('$\varepsilon$ (m)', 'Interpreter', 'latex');
    subplot(3,1,3);
    plot(0:240:(propTime-timeStep), K_kalmGain(:));
    title('Kalman Gain vs. Time', 'Interpreter', 'latex');
    xlabel('Time, t (s)', 'Interpreter', 'latex');
    ylabel('$\left\Vert K_{Kalman} \right\Vert$ = $\sqrt{ (K_x)^2 + (K_y)^2}$', 'Interpreter', 'latex');
    sgtitle({'Problem 1 - Observation Data vs. HCW Fit w/ P0 Tuning Parameter a = 125', printing.myName},...
        'Interpreter', 'latex');
    
    
    %{
    subplot(3,1,2);
    plot(times(:), rhoi(:), 'Color', 'b', 'Marker', 'x', 'LineStyle', 'none');
    xlabel('time, t (s)');
    ylabel('Range (m)');
    %}
    
    
    
    % Problem 1
    fprintf('  Problem 1 - Kalman Filter\n\n');
    fprintf('     Tuning factor, a, for initial covariance, P0 = a*f(R):\n');
    fprintf('            a = %0.1f\n\n', a);
    fprintf('     Minimized chi-squared value:\n');
    fprintf('        chiSq = %0.3f sq. m\n\n', chiSq);
    
    fprintf('%s\n\n', printing.lineBreak);

end