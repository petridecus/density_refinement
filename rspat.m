function [p_c,atommask, epsilon]=rspat(atoms,cryst,r_min,r_max,reso,B)
%

rcryst = 1./cryst;   % // orthonormal
c2f = [rcryst(1),0 ; 0,rcryst(2)];
ATOM_MASK = 3.2;

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
k = min( 0.25/(0.6+0.006*B).^2 , (pi/reso).^2 );   % k_convolution = 0.5*k

p_c = zeros( size(X) );
inv_p_mask = ones( size(X) );
natoms = size(atoms,1);

d2 = ( mod(X+cryst(1)/2,cryst(1))-cryst(1)/2 ).^2 + ...
     ( mod(Y+cryst(2)/2,cryst(2))-cryst(2)/2 ).^2;
if (r_min <= 0) 
	epsilon = (1./(1+exp( d2 - r_max*r_max)));
else
	epsilon = (1./(1+exp( r_min*r_min - d2 ))) .* (1./(1+exp( d2 - r_max*r_max)));
end


%for ii=1:natoms
%for jj=1:natoms
%	x_ij = atoms(ii,:) - atoms(jj,:);
%	d2 = ( mod(X-x_ij(1)+cryst(1)/2,cryst(1))-cryst(1)/2 ).^2 + ...
%	     ( mod(Y-x_ij(2)+cryst(2)/2,cryst(2))-cryst(2)/2 ).^2;
%	inv_msk = 1./(1+exp( d2 - (ATOM_MASK)*(ATOM_MASK)  ));
%
%	p_c = p_c + exp(-d2*k);
%	inv_p_mask = inv_p_mask .* (1 - inv_msk);
%end
%end

atommask = epsilon.*(1-inv_p_mask);
