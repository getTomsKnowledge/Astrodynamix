    %% Problem 2 - Conjunction/Collisions
function collisionData = did_they_collide()
    % Problem values:
    cosmos1408Debris_ID = 50979;
    meteor_ID = 11251;
    TCA = datetime('2022-11-18 01:34:31.386000');
    tStart = TCA - minutes(20);
    tEnd = TCA + minutes(20);
    %stepSize = seconds(30);
    MIN_RNG = 164;
    
    cosmosElements = spacetrack_orbit(cosmos1408Debris_ID, TCA);
    meteorElements = spacetrack_orbit(meteor_ID, TCA);
    
    cosmosConsequences = ...
        conjunction(cosmosElements, meteorElements, tStart, tEnd);
    
    % Problem 2
    fprintf('  Problem 2 - Conjunction/Collision\n\n');
    fprintf('     According to conjunction.m:\n\n');
    fprintf('     The closest approach occurs at:\n');
    fprintf('        TCA_SNaG-app = %s\n', cosmosConsequences.epoch);
    fprintf('     at a distance of:\n')
    fprintf('        rho_SNaG-app = %0.3f m\n\n', cosmosConsequences.distance_m);
    fprintf('     This agrees moderately well with the Space-Track report:\n');
    fprintf('        TCA_spaceTrack = %s\n', TCA);
    fprintf('        rho_spaceTrack = %0.3f m\n\n', MIN_RNG);
    
    fprintf('%s\n\n', printing.lineBreak);
end