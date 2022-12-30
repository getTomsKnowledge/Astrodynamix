%% Batch Least Squares - Orbit Determination

function initialState_linear = get_batch_linear_least_squares(datasetName, obsDim, stateDim)
    
    % 1. PREPARE BATCH:

    % Get data from given .mat file:
    dataset = load(datasetName);
    obsSpace = zeros(length(dataset.observations)/obsDim, 1); % Column vector of observation times
    N = length(dataset.observations); % Total measurements
    
    % Set up independent variable vector, s:
    for n = 1:N
        obsSpace(n) = dataset.times(n);
    end
    
    % 2. APPLY LEAST-SQUARES METHOD
    
    % a) Prepare observation space, yObs:
    yObs = dataset.meas;
    
    % b) Prepare big model matrix, A, using eqn. 5.91:
    bigA = zeros(obsDim*N,stateDim);
    
    % STATE TRANSITION MATRIX:

    % i.f.f. obsDim < 6
    for n = 1:N
        bigPhi = clohessyWiltshire_transition_matrix(meanMotion, obsSpace(n));
        bigA((2*n-2), 1:4)   = [bigPhi(1,1), 0, bigPhi(1,4), bigPhi(1,5)];
        bigA((2*n-1),   1:4) = [bigPhi(2,1), 1, bigPhi(2,4), bigPhi(2,5)];
    end

    % OBSERVATION SPACE TRANSFORMATION:

    % Observation space transformation is identity:
    % H = [eye(2) zeros(2)];

    % COVARIANCE:

    % c) Get covariance matrix, P:
    P1 = bigA'*bigA;

    % NORMAL EQUATION:
    
    % d) Apply Normal Equation, eqn. 9.16:
    initialState_linear = inv(P1)*bigA'*yObs;

end