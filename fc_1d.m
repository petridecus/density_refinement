function [rho_c,x]=fc(atoms,cryst,reso,B)
%

rcryst = 1./cryst;   % // orthonormal
c2f = [rcryst(1)];

grid = reso/3;
gridXY = [ ceil(cryst/grid)];

% // reciprocal space lattice
H = [0:floor(gridXY/2-0.5),-floor(gridXY/2):-1] ;
x = mod( [0:cryst/gridXY:cryst-cryst/gridXY]+cryst/2, cryst) - cryst/2;

% // scattering -- call everything carbon
C_scat = [0.215600, 2.310000, 1.020000, 1.588600, 0.865000, 20.843899, 10.207500, 0.568700, 51.651199];
O_scat = [0.250800, 3.048500, 2.286800, 1.546300, 0.867000, 13.277100, 5.701100, 0.323900, 32.908897];

scat=C_scat;

S2 = (H*rcryst).^2;
resomask = (S2 <= 1/(reso^2));
f_0 = scat(1) + ...
      scat(2) * exp(-scat(6).*S2/4) + ...
      scat(3) * exp(-scat(7).*S2/4) + ...
      scat(4) * exp(-scat(8).*S2/4) + ...
      scat(5) * exp(-scat(9).*S2/4);
f_c = zeros( size(H) );
Bscale  = f_0 .* exp( -B/4*S2 );

natoms = length(atoms);
for ii=1:natoms
	fracX = mod( c2f*atoms(ii)' , 1 );
	dotProd = -(fracX*H);
	f_c = f_c + Bscale .* cos(2*pi*dotProd) + i * (Bscale .* sin(2*pi*dotProd));
end

rho_c = real(ifft(f_c.*resomask));
