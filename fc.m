function [f_c,resomask]=fc(atoms,cryst,reso,B)



rcryst = 1./cryst;   % // orthonormal
c2f = [rcryst(1),0 ; 0,rcryst(2)]; % diagonal linear algebra operation

grid = reso/3;
gridXY = [ ceil(cryst(1)/grid), ceil(cryst(2)/grid) ];

% // reciprocal space lattice
[H,K] = meshgrid ([0:floor(gridXY(1)/2-0.5),-floor(gridXY(1)/2):-1] , ...
                  [0:floor(gridXY(2)/2-0.5),-floor(gridXY(2)/2):-1] );

% // scattering -- call everything carbon
C_scat = [0.215600, 2.310000, 1.020000, 1.588600, 0.865000, 20.843899, 10.207500, 0.568700, 51.651199];

% // calc real space density
S2 = (H*rcryst(1)).^2 +(K*rcryst(2)).^2;
resomask = (S2 <= 1/(reso^2));
f_0 = C_scat(1) + ...
      C_scat(2) * exp(-C_scat(6).*S2/4) + ...
      C_scat(3) * exp(-C_scat(7).*S2/4) + ...
      C_scat(4) * exp(-C_scat(8).*S2/4) + ...
      C_scat(5) * exp(-C_scat(9).*S2/4);
f_c = zeros( size(H) );
Bscale  = f_0 .* exp( -B/4*S2 );

natoms = size(atoms,1);
for ii=1:natoms
	fracX = mod( c2f*atoms(ii,:)' , 1 );
	dotProd = -(fracX(2)*H + fracX(1)*K);
	f_c = f_c + Bscale .* cos(2*pi*dotProd) + i * (Bscale .* sin(2*pi*dotProd));
end

