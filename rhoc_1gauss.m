function [rho_c,atommask]=rhoc(atoms,cryst,reso,B)
%

rcryst = 1./cryst;   % // orthonormal
c2f = [rcryst(1),0 ; 0,rcryst(2)];

grid = reso/3;
gridXY = [ cryst(1)/ceil(cryst(1)/grid), cryst(2)/ceil(cryst(2)/grid) ];

% // real space lattice
[Y,X] = meshgrid ([0:gridXY(1):cryst(1)-gridXY(1)] , ...
                  [0:gridXY(2):cryst(2)-gridXY(2)] );

% // reciprocal space lattice
% //[H,K] = meshgrid ([0:floor(gridXY(2)/2-0.5),-floor(gridXY(2)/2):-1] , ...
% //                  [0:floor(gridXY(1)/2-0.5),-floor(gridXY(1)/2):-1] );

% // scattering
% //k = (pi/reso);
k = min( 0.5/(0.6+0.006*B).^2 , (pi/reso) );

% // calc real space density
% //S2 = (H*rcryst(1)).^2 +(K*rcryst(2)).^2;
% //resomask = (S2 <= 1/(reso^2));
% //f_0 = C_scat(1) + ...
% //      C_scat(2) * exp(-C_scat(6).*S2/4) + ...
% //      C_scat(3) * exp(-C_scat(7).*S2/4) + ...
% //      C_scat(4) * exp(-C_scat(8).*S2/4) + ...
% //      C_scat(5) * exp(-C_scat(9).*S2/4);
% //f_c = zeros( size(H) );
% //p_c = zeros( size(H) );
% //Bscale  = f_0 .* exp( -B/4*S2 );

rho_c = zeros( size(X) );
atommask = false( size(X) );
natoms = size(atoms,1);

for ii=1:natoms
	d2 = ( mod(X-atoms(ii,1)+cryst(1)/2,cryst(1))-cryst(1)/2 ).^2 + ...
	     ( mod(Y-atoms(ii,2)+cryst(2)/2,cryst(2))-cryst(2)/2 ).^2;

	rho_c = rho_c + exp(-d2*k);
	atommask = atommask | d2<(3.2*3.2);
end

