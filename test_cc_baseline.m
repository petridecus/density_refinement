x = [0:.1:14];
[X,Y] = meshgrid(x,x);

ds = [0:.1:7,-7:.1:-.1];
[xd,yd] = meshgrid(ds);
DS = sqrt(xd.^2+yd.^2);
DS_xm1 = sqrt((xd-1).^2+yd.^2);
DS_xp1 = sqrt((xd+1).^2+yd.^2);

tol = 1e-5;
MASKWIDTH = 2.5;
	
for i=1:16
	X_curr = [4.5+i/10,9.3];
	currpos = (abs(X-X_curr(1))<tol).*(abs(Y-X_curr(2))<tol);
	
	% // true atom locations
	X1 = [5.4,9.3];
	X2 = [5.3,10.7];
	X3 = [3.8,8.8];
	
	STDEV = .2 +3*rand();
	HEIGHT = 1; % + 4*rand();
	stdevs(i) = STDEV;
	
	% // observed density
	d1 = ((X-X1(1)).^2+(Y-X1(2)).^2);
	d2 = ((X-X2(1)).^2+(Y-X2(2)).^2);
	d3 = ((X-X3(1)).^2+(Y-X3(2)).^2);
	fo = HEIGHT*exp(-d1/(STDEV*STDEV)) + ...
	     HEIGHT*exp(-d2/(STDEV*STDEV)) + ...
	     HEIGHT*exp(-d3/(STDEV*STDEV)) ;
	
	
	atmmask = abs(DS)<MASKWIDTH;
	atm   = exp(-DS.^2/(STDEV*STDEV)).*atmmask;
	% //atmm1 = exp(-DS_xm1.^2/(STDEV*STDEV)).*atmmask;
	% //atmp1 = exp(-DS_xp1.^2/(STDEV*STDEV)).*atmmask;
	
	% // put a spike at each atom location ...
	impulse_m1 = (abs(X-(X_curr(1)-0.1))<tol).*(abs(Y-X_curr(2))<tol)  + ...
	             (abs(X-X2(1))<tol).*(abs(Y-X2(2))<tol) + ...
	             (abs(X-X3(1))<tol).*(abs(Y-X3(2))<tol);
	impulse    = (abs(X-X_curr(1)    )<tol).*(abs(Y-X_curr(2))<tol)  + ...
	             (abs(X-X2(1))<tol).*(abs(Y-X2(2))<tol) + ...
	             (abs(X-X3(1))<tol).*(abs(Y-X3(2))<tol);
	impulse_p1 = (abs(X-(X_curr(1)+0.1))<tol).*(abs(Y-X_curr(2))<tol)  + ...
	             (abs(X-X2(1))<tol).*(abs(Y-X2(2))<tol) + ...
	             (abs(X-X3(1))<tol).*(abs(Y-X3(2))<tol);
	
	% // ... and convolute to compute fc
	fc_m1 = real(ifft2( fft2(atm).*fft2(impulse_m1) ));
	fc    = real(ifft2( fft2(atm).*fft2(impulse) ));
	fc_p1 = real(ifft2( fft2(atm).*fft2(impulse_p1) ));
	
	% // ... also mask ea. fc
	resmask    = abs(fc)>tol;
	resmask_m1 = abs(fc_m1)>tol;
	resmask_p1 = abs(fc_p1)>tol;
	
	% // Compute the true CC at each position ...
	true_cc_p1 = corrcoef( fc_p1(resmask_p1==1) , fo(resmask_p1==1) );
	true_cc    = corrcoef(    fc(resmask==1)    , fo(resmask==1) );
	true_cc_m1 = corrcoef( fc_m1(resmask_m1==1) , fo(resmask_m1==1) );
	
	% // ... as well as the (numeric) derivative
	true_dcc_dx(i) = ( true_cc_p1 - true_cc_m1 ) / 0.2;
	
	% // stats over mask
	V = sum(sum(resmask));
	C = sum(fc(resmask));
	O = sum(fo(resmask));
	C2 = sum(fc(resmask).^2);
	O2 = sum(fo(resmask).^2);
	CO = sum( fo(resmask==1).*fc(resmask==1) );
	
	% // more stats over mask
	% // convolutions	
	ddx = zeros(141,141); 
	 ddx(1,140) = 0; ddx(1,141) = -1;    ddx(1,2) = 1; ddx(1,3) = 0; 
	conv_atm_map  = real( ifft2(fft2(atm).*fft2(fo)) );
	 ddx_atm_map   = real( ifft2(fft2(ddx).*fft2(conv_atm_map)) ) / 0.2;
	conv_map_mask  = real( ifft2(fft2(atmmask).*fft2(fo)) );
	 ddx_map_mask   = real( ifft2(fft2(ddx).*fft2(conv_map_mask)) ) / 0.2;
	conv_map2_mask  = real( ifft2(fft2(atmmask).*fft2(fo.^2)) );  
	 ddx_map2_mask   = real( ifft2(fft2(ddx).*fft2(conv_map2_mask)) ) / 0.2;
	

	% // stats over _shifted_ resmask
	Vp1 = sum(sum(resmask_p1));
	Vm1 = sum(sum(resmask_m1));
	Op1 = sum(fo(resmask_p1));
	Om1 = sum(fo(resmask_m1));
	O2p1 = sum(fo(resmask_p1).^2);
	O2m1 = sum(fo(resmask_m1).^2);
	Cp1 = sum(fc_p1(resmask_p1));
	Cm1 = sum(fc_m1(resmask_m1));
	C2p1 = sum(fc_p1(resmask_p1).^2);
	C2m1 = sum(fc_m1(resmask_m1).^2);
	COp1 = sum( fo(resmask_p1).*fc_p1(resmask_p1) );
	COm1 = sum( fo(resmask_m1).*fc_m1(resmask_m1) );
	
	% // numeric deriv of stats over resmask
	delO  = ( Op1 - Om1 )   / 0.2;
	delO2 = ( O2p1 - O2m1 ) / 0.2;
	delC  = ( Cp1 - Cm1 )   / 0.2;
	delC2 = ( C2p1 - C2m1 ) / 0.2;
	delV =  ( Vp1 - Vm1 )   / 0.2;
	delCO = ( COp1 - COm1 ) / 0.2;
	
	% // true_delO(i) = delO;
	% // pred_delO(i) = ddx_map_mask(impulse1==1);
	
	% // true_delO2(i) = delO2;
	% // pred_delO2(i) = ddx_map2_mask(impulse1==1);
	
	% // true_delCO(i) = delCO;
	% // pred_delCO(i) = ddx_atm_map(impulse1==1);
	
	varC = (V*C2 - C*C);
	varO = (V*O2 - O*O);

	f       =  ( V*CO - C*O ) / (V);
	g       =  ( sqrt(varO*varC) ) / (V);

	% //f_prime =  ddx_atm_map(impulse1==1)*V + C*ddx_map_mask(impulse1==1);
	% //g_prime =  0.5/(g) * sqrt(varO) * (V*ddx_map2_mask(impulse1==1) - 2*C*ddx_map_mask(impulse1==1));
	
	f_prime = (V*delCO + delV*CO - C*delO) / (V);
	g_prime = ( (.5/sqrt(varC) * ( delV*C2 + V*delC2 )) * sqrt(varO) + ... 
	            (.5/sqrt(varO) * ( delV*O2 + V*delO2 - 2*O*delO )) * sqrt(varC) ) / (V);
	          
	f_p1 = (Vp1*COp1 - Cp1*Op1) / (V);
	f_m1 = (Vm1*COm1 - Cm1*Om1) / (V);
	true_fp(i) = (f_p1 - f_m1)/.2;
	pred_fp(i) = f_prime;
	
	varCp1 = (Vp1*C2p1 - Cp1*Cp1);
	varOp1 = (Vp1*O2p1 - Op1*Op1);
	varCm1 = (Vm1*C2m1 - Cm1*Cm1);
	varOm1 = (Vm1*O2m1 - Om1*Om1);
	g_p1 = ( sqrt(varCp1 * varOp1))  / (V);
	g_m1 = ( sqrt(varCm1 * varOm1)) / (V);
	true_gp(i) = (g_p1 - g_m1)/.2;
	pred_gp(i) = g_prime;
	
	true_ccs1(i) =  f/g;
	true_ccs2(i) = true_cc;
	
	true_ccs_p1_1(i) =  f_p1/g_p1;
	true_ccs_p1_2(i) = true_cc_p1;
	
	true_ccs_m1_1(i) =  f_m1/g_m1;
	true_ccs_m1_2(i) = true_cc_m1;
	
	% //pred_dcc_dx1(i) =  (V * ddx_atm_map(currpos==1) ) / D;
	pred_dcc_dx1(i) =  (f_prime*g - g_prime*f)/(g*g);
	pred_dcc_dx2(i) =  (true_fp(i)*g - true_gp(i)*f)/(g*g);
end 

cor(true_ccs1, true_ccs2)
polyfit(true_ccs1, true_ccs2,1)
cor(true_ccs_p1_1, true_ccs_p1_2)
polyfit(true_ccs_p1_1, true_ccs_p1_2,1)
cor(true_ccs_m1_1, true_ccs_m1_2)
polyfit(true_ccs_m1_1, true_ccs_m1_2,1)
disp('---')
cor(true_dcc_dx,pred_dcc_dx1)
cor(true_dcc_dx,pred_dcc_dx2)
polyfit(true_dcc_dx,pred_dcc_dx1,1)
polyfit(true_dcc_dx,pred_dcc_dx2,1)
disp('---')
cor(true_fp,pred_fp)
cor(true_gp,pred_gp)
polyfit(true_fp,pred_fp,1)
polyfit(true_gp,pred_gp,1)
