function [rho_c,atommask]=rhoc(atoms,cryst,reso,B)
% rho_c    - output intensities, 2D array (will be 3D for our purposes)
% resomask - every reflection has an associated resolution, boolean data
%            use if true

% atoms    - XYZ cartestion coordinate of atom, 
% cryst    - crystal size, only 2D here
% reso     - resolution for this atom and we return resomask value that 
%            is based on this compared to some target
% B        - atomic B factor. will eventually be fitted, we can ignore for the
%            time being

% ex crystal is dimensions [4, 5]
% reso is .01 (1/100)
grid = reso/3; % grid = .00333333 (1/300)

% what does ceil do?
ceil1 = ceil(cryst(1)/grid); % 4 * 300 = 1200
ceil2 = ceil(cryst(2)/grid); % 5 * 300 = 1500

gridXY = [ cryst(1)/ceil1, cryst(2)/ceil2 ]; % 4/1200, 5/1500 = .00333, .00333
% ^ are the above lines meant to round the step size in some way

% // real space lattice
[Y,X] = meshgrid (0:gridXY(1):cryst(1)-gridXY(1) , ...
                  0:gridXY(2):cryst(2)-gridXY(2));
% ^ iterating thru crystal size with step size effectively equal to reso/3

% Y - Should have 1500 rows and 1200 columns 
% -----------------------------------------------
%            0  : [0, .0033, .0066, ... , 3.9967]
%  .0033    (1) : [0, .0033, .0066, ... , 3.9967]
%              ...
% 4.9967 (1499) : [0, .0033, .0066, ... , 3.9967]

% X - Should have 1200 rows and 1500 columns 
% -----------------------------------------------
%            0  : [0, .0033, .0066, ... , 4.9967]
%  .0033    (1) : [0, .0033, .0066, ... , 4.9967]
%              ...
% 3.9967 (1199) : [0, .0033, .0066, ... , 4.9967]

% // reciprocal space lattice
% //[H,K] = meshgrid ([0:floor(gridXY(2)/2-0.5),-floor(gridXY(2)/2):-1] , ...
% //                  [0:floor(gridXY(1)/2-0.5),-floor(gridXY(1)/2):-1] );

% // scattering -- call everything carbon
C_scat = [0.215600, 2.310000, 1.020000, 1.588600, 0.865000, 20.843899, 10.207500, 0.568700, 51.651199];

% // calc real space density
% //S2 = (H*rcryst(1)).^2 +(K*rcryst(2)).^2;
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
    % atom at [1, 4]

    % -------------------------------------------------------------------
    d2X = mod(X-atoms(ii,1)+cryst(1)/2,cryst(1)) - cryst(1)/2;
    % For every value in matrix X:
    % - subtract 1
    % - add (4 / 2)
    % - take that value % 4
    % - subtract 2 again

    % ((1.0033 - 1 + 2) % 4) - 2 = (2.0033 % 4) - 2 = 2.0033 - 2 = .0033
    % ((3.5 - 1 + 2) % 4) - 2 = (4.5 % 4) - 2 = .5 - 2 = -1.5
    % ((.6 - 1 + 2) % 4) - 2 = (1.6 % 4) - 2 = 1.6 - 2 = -.4

    % d2X is an array of dimensions X, with values that are closer to 0
    % the closer the value in X is to the xcoord of atoms(ii)
    % -------------------------------------------------------------------

    d2Y = mod(Y-atoms(ii,2)+cryst(2)/2,cryst(2)) - cryst(2) / 2;
    
    d2 = ( d2X ).^2 + ...
	     ( d2Y ).^2;

	rho_c = rho_c + C_scat(1) * sqrt(pi*4/B) * exp(-d2*pi*pi*4/B) + ...
                    C_scat(2) * sqrt(pi*4/(B+C_scat(6))) * exp(-d2*pi*pi*4/(B+C_scat(6)) ) + ...
                    C_scat(3) * sqrt(pi*4/(B+C_scat(7))) * exp(-d2*pi*pi*4/(B+C_scat(7)) ) + ...
                    C_scat(4) * sqrt(pi*4/(B+C_scat(8))) * exp(-d2*pi*pi*4/(B+C_scat(8)) ) + ...
                    C_scat(5) * sqrt(pi*4/(B+C_scat(9))) * exp(-d2*pi*pi*4/(B+C_scat(9)) ) ;
	atommask = atommask | d2<(3.2*3.2);
end

