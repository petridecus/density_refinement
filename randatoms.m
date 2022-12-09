function atoms=randatoms(N,cryst)
%

bondlength=1.5;

atoms = zeros(N,2);

atoms(1,:) = rand(1,2);
theta = rand()*360;

for ii=2:N
	atoms(ii,:) = atoms(ii-1,:) + bondlength*[cosd(theta),sind(theta)];
	theta = theta + (-1)^ii*(30 + rand()*100);
end

% // recenter
shift = -mean(atoms,1) + cryst/2;
for ii=1:N
	atoms(ii,:) = atoms(ii,:) + shift;
end
