function [sigmaa]=sigmaa(cryst,reso,rms)
%

fsol=0.95;
bsol=300;

%

rcryst = 1./cryst;   % // orthonormal
c2f = [rcryst(1),0 ; 0,rcryst(2)];

grid = reso/3;
gridXY = [ cryst(1)/ceil(cryst(1)/grid), cryst(2)/ceil(cryst(2)/grid) ];

% // reciprocal space lattice
[H,K] = meshgrid ([0:floor(gridXY(2)/2-0.5),-floor(gridXY(2)/2):-1] , ...
                  [0:floor(gridXY(1)/2-0.5),-floor(gridXY(1)/2):-1] );


% // calc real space density
S2 = (H*rcryst(1)).^2 +(K*rcryst(2)).^2
sigmaa = sqrt ( (1-fsol*exp(-bsol*S2))*exp(-8*pi*pi/3*rms*rms*S2) );
