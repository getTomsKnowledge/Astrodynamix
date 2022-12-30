% Clohessy-Wiltshire relations for relative motion in
% circular, Restricted 2BP orbit

function bigPhi = get_phi_hcw(nIn, tIn)

    consts = OrbitConstants();
    
    % Parameters:
    c = cos(nIn*tIn);
    s = sin(nIn*tIn);
    
    % Quadrants (position, velocity):
    rrPhi = [      (4-3*c)               0             0;
                6*(s-nIn*tIn)            1             0;
                      0                  0             c];
    rvPhi = [      s/nIn          (2/nIn)*(1-c)        0;
              -(2/nIn)*(1-c)   (((4*s)/nIn)-3*tIn)     0
                      0                  0           s/nIn];
    vrPhi = [      3*nIn*s               0             0;
                 -6*nIn*(1-c)            0             0
                      0                  0          -nIn*s];
    vvPhi = [       c                   2*s            0;
                  -2*s                (4*c-3)          0;
                    0                    0             c];
    
    % Concatenate quadrants:
    bigPhi = [rrPhi, rvPhi;
              vrPhi, vvPhi];

end