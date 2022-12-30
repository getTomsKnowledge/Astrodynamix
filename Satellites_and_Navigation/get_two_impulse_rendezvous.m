
%{

TEST DATA:
% Problem values:
tof = 2400; % seconds
z = 900; % altitude in km
r = 900 + consts.R_earth_km;

drTwoImp = [200, 150, -150]'; %km
drdotTwoImp = [-0.1, 0.050, 0.08]'; %km/s

%}
function [drdotDelta1, drdotDelta2, drdotDeltaTot] = get_two_impulse_rendezvous(tof, zTwoImpulse_m, drTwoImp, drdotTwoImp)

rTwoImpulse_m = zTwoImpulse_m + consts.R_earth_km*1e3;
nTwoImpulse = sqrt( (consts.mu_earth_km*1e9) / (rTwoImpulse_m^3) );
cTwo = cos(nTwoImpulse*tof);
sTwo = sin(nTwoImpulse*tof);


rrPhi = [    (4-3*cTwo)                             0                   0;
            6*(sTwo-nTwoImpulse*tof)                1                   0;
                0                                   0                 cTwo];
rvPhi = [      sTwo/nTwoImpulse          (2/nTwoImpulse)*(1-cTwo)       0;
          -(2/nTwoImpulse)*(1-cTwo)   (((4*sTwo)/nTwoImpulse)-3*tof)    0;
                0                                    0           sTwo/nTwoImpulse];
vrPhi = [     3*nTwoImpulse*sTwo                     0                  0;
            -6*nTwoImpulse*(1-cTwo)                  0                  0;
                0                                    0        -nTwoImpulse*sTwo];
vvPhi = [       cTwo                              2*sTwo                0;
              -2*sTwo                           (4*cTwo-3)              0;
                0                                    0               cTwo];

bigPhi = [rrPhi, rvPhi;
          vrPhi, vvPhi];

% Compute the separate delta v's:

% Position match:
drdotDelta1 = -inv(rvPhi)*rrPhi*drTwoImp - drdotTwoImp;
% Velocity match @ target:
drdotDelta2 = -(vrPhi - vvPhi*inv(rvPhi)*rrPhi) * drdotTwoImp;

drdotDeltaTot = norm(drdotDelta1) + norm(drdotDelta2);


end