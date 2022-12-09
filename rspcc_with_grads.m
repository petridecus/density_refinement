function [pcc,grads,  dCOdx,dOdx,dO2dx,dCdx,dC2dx,dVdx]=rscc(p_o,atoms,cryst,r_min,r_max,reso,B)
%

rcryst = 1./cryst;   % // orthonormal
c2f = [rcryst(1),0 ; 0,rcryst(2)];
ATOM_MASK = 3.2;

grid = reso/3;
gridXY = [ cryst(1)/ceil(cryst(1)/grid), cryst(2)/ceil(cryst(2)/grid) ];

% // real space lattice
[Y,X] = meshgrid ([0:gridXY(1):cryst(1)-gridXY(1)] , ...
                  [0:gridXY(2):cryst(2)-gridXY(2)] );

% // scattering
% //k = (pi/reso);
k = min( 0.25/(0.6+0.006*B).^2 , (pi/reso)^2 );   % k_convolution = 0.5*k

p_c = zeros( size(X) );
inv_p_mask = ones( size(X) );
natoms = size(atoms,1);


% // 0 initialize
for ii=1:natoms
	%
	rho_dx_mask{ii} = zeros( size(X) );
	rho_dy_mask{ii} = zeros( size(X) );
	rho_dx_atm{ii} = zeros( size(X) );
	rho_dy_atm{ii} = zeros( size(X) );
end

% // 0.5 epsilon mask
d2 = ( mod(X+cryst(1)/2,cryst(1))-cryst(1)/2 ).^2 + ...
     ( mod(Y+cryst(2)/2,cryst(2))-cryst(2)/2 ).^2;
if (r_min <= 0) 
	epsilon = (1./(1+exp( d2 - r_max*r_max)));
else
	epsilon = (1./(1+exp( r_min*r_min - d2 ))) .* (1./(1+exp( d2 - r_max*r_max)));
end


% // 1 compute p_c
for ii=1:natoms
	for jj=1:natoms
		x_ij = atoms(ii,:) - atoms(jj,:);
		dX = mod(x_ij(1)-X+cryst(1)/2,cryst(1))-cryst(1)/2;
		dY = mod(x_ij(2)-Y+cryst(2)/2,cryst(2))-cryst(2)/2;
		d2 = ( dX.^2 + dY.^2 );

		% // 2 masks ... radial mask from r_min to r_max
		% //         ... mask around each gaussian for segmentation
		% // Both smoothly decay
		% // epsilon = (1./(1+exp( r_min*r_min - dot(x_ij,x_ij) ))) .* ...
		% //           (1./(1+exp( dot(x_ij,x_ij) - r_max*r_max)));
		inv_msk = 1./(1+exp( d2 - (ATOM_MASK*ATOM_MASK)  ));
		atm = exp(-d2*k);

		% // aggregate
		p_c = p_c + atm;
		inv_p_mask = inv_p_mask .* (1 - inv_msk);

		if (ii ~= jj)
			% // derivative caching
			rho_dx_mask{ii} = rho_dx_mask{ii} + ( -inv_msk .* dX );
			rho_dy_mask{ii} = rho_dy_mask{ii} + ( -inv_msk .* dY );
			rho_dx_atm{ii} = rho_dx_atm{ii} + ( atm .* dX );
			rho_dy_atm{ii} = rho_dy_atm{ii} + ( atm .* dY );
	
			% // as j moves it also changes this vector
			rho_dx_mask{jj} = rho_dx_mask{jj} + ( -inv_msk .* -dX );
			rho_dy_mask{jj} = rho_dy_mask{jj} + ( -inv_msk .* -dY );
			rho_dx_atm{jj} = rho_dx_atm{jj} + ( atm .* -dX );
			rho_dy_atm{jj} = rho_dy_atm{jj} + ( atm .* -dY );
		end
	end
end

% // 2 COMPUTE SUMMARY STATISTICS
p_mask = epsilon.*(1-inv_p_mask);
sumCO = sum(sum( p_o.*p_c.*p_mask ));
sumO  = sum(sum( p_o.*p_mask ));
sumO2 = sum(sum( p_o.*p_o.*p_mask ));
sumC  = sum(sum( p_c.*p_mask ));
sumC2 = sum(sum( p_c.*p_c.*p_mask ));
sumV  = sum(sum( p_mask ));

varC = (sumC2 - sumC*sumC / sumV );
varO = (sumO2 - sumO*sumO / sumV );
if (varC == 0 || varO == 0)
	pcc = 0;
else
	pcc = (sumCO - sumC*sumO/sumV) / sqrt( varC*varO );
end

% // 3 CALCULATE PER-ATOM DERIVATIVES
grads = zeros(natoms,2);
for ii=1:natoms
	delx_mask = 2.0*epsilon.*inv_p_mask.*rho_dx_mask{ii};
	dVdx  = sum(sum( delx_mask ));
	dOdx  = sum(sum( delx_mask.*p_o ));
	dO2dx = sum(sum( delx_mask.*p_o.*p_o ));
	% // delx_rhoc = -2.0*k*rho_dx_atm{ii};
	delx_eps = -2.0*k*epsilon.*rho_dx_atm{ii};
	dCdx  = sum(sum( delx_eps ));
	dCOdx = sum(sum( delx_eps.*p_o ));
	dC2dx = sum(sum( 2.0*delx_eps.*p_c ));

	dely_mask = 2.0*epsilon.*inv_p_mask.*rho_dy_mask{ii};
	dVdy  = sum(sum( dely_mask ));
	dOdy  = sum(sum( dely_mask.*p_o ));
	dO2dy = sum(sum( dely_mask.*p_o.*p_o ));
	% // dely_rhoc = -2.0*k*rho_dy_atm{ii};
	dely_eps = -2.0*k*epsilon.*rho_dy_atm{ii};
	dCdy  = sum(sum( dely_eps ));
	dCOdy = sum(sum( dely_eps.*p_o ));
	dC2dy = sum(sum( 2.0*dely_eps.*p_c ));

	% // finally compute dCC/dx_ij
	f = ( sumCO - sumC*sumO / sumV );
	g = sqrt( varC*varO );

	% // dfdx = dCOdx - 1/(sumV*sumV) * ( dOdx*sumC*sumV - sumO*sumC*dVdx);
	% // dgdx = 0.5 * ( ...
	% // 		sqrt(varO)/sqrt(varC) * ( dC2dx + ( sumC*sumC*dVdx/(sumV*sumV) ) )  + ...
	% // 		sqrt(varC)/sqrt(varO) * ( dO2dx - ( 1/(sumV*sumV) * ( 2*sumV*sumO*dOdx - sumO*sumO*dVdx ) ) ) );
	dfdx = dCOdx - 1/(sumV*sumV) * ( dOdx*sumC*sumV + dCdx*sumO*sumV - sumO*sumC*dVdx);
	dgdx = 0.5 * ( ...
			sqrt(varO)/sqrt(varC) * ( dC2dx - ( 1/(sumV*sumV) * ( 2*sumV*sumC*dCdx - sumC*sumC*dVdx ) ) )  + ...
			sqrt(varC)/sqrt(varO) * ( dO2dx - ( 1/(sumV*sumV) * ( 2*sumV*sumO*dOdx - sumO*sumO*dVdx ) ) ) );
	dfdy = dCOdy - 1/(sumV*sumV) * ( dOdy*sumC*sumV + dCdy*sumO*sumV - sumO*sumC*dVdy);
	dgdy = 0.5 * ( ...
			sqrt(varO)/sqrt(varC) * ( dC2dy - ( 1/(sumV*sumV) * ( 2*sumV*sumC*dCdy - sumC*sumC*dVdy ) ) )  + ...
			sqrt(varC)/sqrt(varO) * ( dO2dy - ( 1/(sumV*sumV) * ( 2*sumV*sumO*dOdy - sumO*sumO*dVdy ) ) ) );

	grads(ii,:) = [ (g*dfdx - f*dgdx)/(g*g) , (g*dfdy - f*dgdy)/(g*g) ];
end


