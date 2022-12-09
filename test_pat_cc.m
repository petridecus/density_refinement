cryst = [100.0,100.0];
natoms = 1500;
B = 25;
reso = 2.0;
patt_mask = [5.0,20.0]; % // ignored for now

moving_atom = 6
H = 0.1;

pcc_a = 0;
pcc_n = 0;
dx_a = 0;
dx_n = 0;
dy_a = 0;
dy_n = 0;

for ii=1:1
	atoms1 = randatoms(natoms,cryst,1.5);
	atoms2 = atoms1 + 4.0 * rand(size(atoms1)) - 2.0;
	
	% // patterson 1
	p_1 = rspat(atoms1,cryst,patt_mask(1),patt_mask(2),reso,B);

	% // patterson 2
	[p_2, p_mask,epsilon] = rspat(atoms2,cryst,patt_mask(1),patt_mask(2),reso,B);

	pcc_n(ii) = patterson_correl_masked( p_1, p_2, p_mask);

	% // numeric derivs
	atoms2_xm1 = atoms2;
	atoms2_xm1(moving_atom,:) = atoms2_xm1(moving_atom,:) - [H,0];
	atoms2_xp1 = atoms2;
	atoms2_xp1(moving_atom,:) = atoms2_xp1(moving_atom,:) + [H,0];
	[p_2_xm1, p_mask_xm1] = rspat(atoms2_xm1,cryst,patt_mask(1),patt_mask(2),reso,B);
	[p_2_xp1, p_mask_xp1] = rspat(atoms2_xp1,cryst,patt_mask(1),patt_mask(2),reso,B);

	[pcc_xm1, sumCO_xm1,sumO_xm1,sumO2_xm1,sumC_xm1,sumC2_xm1,sumV_xm1] = ...
		patterson_correl_masked( p_1, p_2_xm1, p_mask_xm1);
	[pcc_xp1, sumCO_xp1,sumO_xp1,sumO2_xp1,sumC_xp1,sumC2_xp1,sumV_xp1] = ...
		patterson_correl_masked( p_1, p_2_xp1, p_mask_xp1);

	ddx_n = ([sumCO_xp1,sumO_xp1,sumO2_xp1,sumC_xp1,sumC2_xp1,sumV_xp1] - ...
		[ sumCO_xm1,sumO_xm1,sumO2_xm1,sumC_xm1,sumC2_xm1,sumV_xm1]) / (2*H);

	atoms2_ym1 = atoms2;
	atoms2_ym1(moving_atom,:) = atoms2_ym1(moving_atom,:) - [0,H];
	atoms2_yp1 = atoms2;
	atoms2_yp1(moving_atom,:) = atoms2_yp1(moving_atom,:) + [0,H];
	[p_2_ym1, p_mask_ym1] = rspat(atoms2_ym1,cryst,patt_mask(1),patt_mask(2),reso,B);
	[p_2_yp1, p_mask_yp1] = rspat(atoms2_yp1,cryst,patt_mask(1),patt_mask(2),reso,B);
	pcc_ym1 = patterson_correl_masked( p_1, p_2_ym1, p_mask_ym1);
	pcc_yp1 = patterson_correl_masked( p_1, p_2_yp1, p_mask_yp1);
	dpcc_dx_numeric = [ (pcc_xp1 - pcc_xm1) / (2*H) , (pcc_yp1 - pcc_ym1) / (2*H) ];
	dx_n(ii) = dpcc_dx_numeric(1);
	dy_n(ii) = dpcc_dx_numeric(2);

	% // analytic derivs
	[pcc_a(ii), grads, dCOdx,dOdx,dO2dx,dCdx,dC2dx,dVdx] = ...
		rspcc_with_grads(p_1,atoms2,cryst,patt_mask(1),patt_mask(2),reso,B);

	ddx_a=[dCOdx,dOdx,dO2dx,dCdx,dC2dx,dVdx];

	dpcc_dx_analytic = grads(moving_atom,:);
	dx_a(ii) = dpcc_dx_analytic(1);
	dy_a(ii) = dpcc_dx_analytic(2);
end

corrcoef(pcc_n,pcc_a)
corrcoef(dx_n,dx_a)
corrcoef(dy_n,dy_a)
