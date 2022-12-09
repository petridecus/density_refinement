STEPSIZE = 1.0; % // must evenly divide 8

x = [0:STEPSIZE:16];
[Y,X] = meshgrid(x,x);

ds = [0:STEPSIZE:8,-8:STEPSIZE:-STEPSIZE];
[yd,xd] = meshgrid(ds);

DS = sqrt(xd.^2+yd.^2);

tol = 1e-5;

MASKWIDTH = 4.0
SIGMOID_STEEPNESS  = 3;
STDEV     =  1

for i=1:1000
	offset = 4*rand(1,2)-2;
	off_dist(i) = norm(offset);
	% // true atom locations
	X1 = [8,8]; %//+ rand(1,2);
%//X1 = round(X1/STEPSIZE)*STEPSIZE;
	theta1 = rand()*360;
	thetas(i) = theta1;
	theta2 = mod( theta1+120 , 360 );
	X2 = X1 + 1.4*[cosd(theta1),sind(theta1)]; 
	X3 = X1 + 1.4*[cosd(theta2),sind(theta2)]; 
	theta3 = mod( theta2-180 + 120 , 360);
	theta4 = mod( theta2-180 + 240 , 360);
	X4 = X3 + 1.5*[cosd(theta3),sind(theta3)]; 
	X5 = X3 + 1.5*[cosd(theta4),sind(theta4)];
%//X2 = round(X2/STEPSIZE)*STEPSIZE;
%//X3 = round(X3/STEPSIZE)*STEPSIZE;

	X_curr = X1+offset;
%//X_curr = round(X_curr/STEPSIZE)*STEPSIZE;
	
	HEIGHT = 1; % // + 4*rand();
	
	% // observed density
	d1 = (X-(X1(1))).^2 + (Y-(X1(2))).^2;
	d2 = (X-(X2(1))).^2 + (Y-(X2(2))).^2;
	d3 = (X-(X3(1))).^2 + (Y-(X3(2))).^2;
	d4 = (X-(X4(1))).^2 + (Y-(X4(2))).^2;
	d5 = (X-(X5(1))).^2 + (Y-(X5(2))).^2;
	fo = HEIGHT*exp(-d1/(STDEV*STDEV)) + ...
	     HEIGHT*exp(-d2/(STDEV*STDEV)) + ...
	     HEIGHT*exp(-d3/(STDEV*STDEV)) + ...
	     HEIGHT*exp(-d4/(STDEV*STDEV)) + ...
	     HEIGHT*exp(-d5/(STDEV*STDEV)) ;
	
	
	% //atmmask = abs(DS)<MASKWIDTH;
	atmmask = 1./(1 + exp(-SIGMOID_STEEPNESS * (MASKWIDTH-DS)));
	atm   = exp(-DS.^2/(STDEV)^2).*atmmask;
	
	% // put a spike at each atom location ...
	XClow = floor( X_curr / STEPSIZE + 1 );
	X2low = floor( X2 / STEPSIZE + 1 );
	X3low = floor( X3 / STEPSIZE + 1 );
	X4low = floor( X4 / STEPSIZE + 1 );
	X5low = floor( X5 / STEPSIZE + 1 );
	XCcoeff =  X_curr / STEPSIZE + 1 - XClow;
	X2coeff =  X2 / STEPSIZE + 1 - X2low;
	X3coeff =  X3 / STEPSIZE + 1 - X3low;
	X4coeff =  X4 / STEPSIZE + 1 - X4low;
	X5coeff =  X5 / STEPSIZE + 1 - X5low;
	
	currpos = zeros( size(X) );
		currpos( XClow(1)   , XClow(2) )   = (1-XCcoeff(1))*(1-XCcoeff(2));
		currpos( XClow(1)+1 , XClow(2) )   = (XCcoeff(1))*(1-XCcoeff(2));
		currpos( XClow(1)   , XClow(2)+1 ) = (1-XCcoeff(1))*(XCcoeff(2));
		currpos( XClow(1)+1 , XClow(2)+1 ) = (XCcoeff(1))*(XCcoeff(2));

	impulse = zeros( size(X) );
		impulse( X2low(1)   , X2low(2) )   = impulse( X2low(1)   , X2low(2) )  + (1-X2coeff(1))*(1-X2coeff(2));
		impulse( X2low(1)+1 , X2low(2) )   = impulse( X2low(1)+1 , X2low(2) )  + (X2coeff(1))*(1-X2coeff(2));
		impulse( X2low(1)   , X2low(2)+1 ) = impulse( X2low(1)   , X2low(2)+1 ) + (1-X2coeff(1))*(X2coeff(2));
		impulse( X2low(1)+1 , X2low(2)+1 ) = impulse( X2low(1)+1 , X2low(2)+1 ) + (X2coeff(1))*(X2coeff(2));
	% //---
		impulse( X3low(1)   , X3low(2) )   = impulse( X3low(1)   , X3low(2) ) + (1-X3coeff(1))*(1-X3coeff(2));
		impulse( X3low(1)+1 , X3low(2) )   = impulse( X3low(1)+1 , X3low(2) ) + (X3coeff(1))*(1-X3coeff(2));
		impulse( X3low(1)   , X3low(2)+1 ) = impulse( X3low(1)   , X3low(2)+1 ) + (1-X3coeff(1))*(X3coeff(2));
		impulse( X3low(1)+1 , X3low(2)+1 ) = impulse( X3low(1)+1 , X3low(2)+1 ) + (X3coeff(1))*(X3coeff(2));
	% //---
		impulse( X4low(1)   , X4low(2) )   = impulse( X4low(1)   , X4low(2) ) + (1-X4coeff(1))*(1-X4coeff(2));
		impulse( X4low(1)+1 , X4low(2) )   = impulse( X4low(1)+1 , X4low(2) ) + (  X4coeff(1))*(1-X4coeff(2));
		impulse( X4low(1)   , X4low(2)+1 ) = impulse( X4low(1)   , X4low(2)+1 ) + (1-X4coeff(1))*(  X4coeff(2));
		impulse( X4low(1)+1 , X4low(2)+1 ) = impulse( X4low(1)+1 , X4low(2)+1 ) + (  X4coeff(1))*(  X4coeff(2));
	% //---
		impulse( X5low(1)   , X5low(2) )   = impulse( X5low(1)   , X5low(2) ) + (1-X5coeff(1))*(1-X5coeff(2));
		impulse( X5low(1)+1 , X5low(2) )   = impulse( X5low(1)+1 , X5low(2) ) + (  X5coeff(1))*(1-X5coeff(2));
		impulse( X5low(1)   , X5low(2)+1 ) = impulse( X5low(1)   , X5low(2)+1 ) + (1-X5coeff(1))*(  X5coeff(2));
		impulse( X5low(1)+1 , X5low(2)+1 ) = impulse( X5low(1)+1 , X5low(2)+1 ) + (  X5coeff(1))*(  X5coeff(2));
	% //---
	impulse_m1 = impulse;
	impulse_p1 = impulse;
		impulse( XClow(1)   , XClow(2) )   = impulse( XClow(1)   , XClow(2) ) + (1-XCcoeff(1))*(1-XCcoeff(2));
		impulse( XClow(1)+1 , XClow(2) )   = impulse( XClow(1)+1 , XClow(2) ) + (XCcoeff(1))*(1-XCcoeff(2));
		impulse( XClow(1)   , XClow(2)+1 ) = impulse( XClow(1)   , XClow(2)+1 ) + (1-XCcoeff(1))*(XCcoeff(2));
		impulse( XClow(1)+1 , XClow(2)+1 ) = impulse( XClow(1)+1 , XClow(2)+1 ) + (XCcoeff(1))*(XCcoeff(2));
	% //---
		impulse_m1( XClow(1)-1 , XClow(2) )   = impulse_m1( XClow(1)-1 , XClow(2) ) + (1-XCcoeff(1))*(1-XCcoeff(2));
		impulse_m1( XClow(1)   , XClow(2) )   = impulse_m1( XClow(1)   , XClow(2) ) + (XCcoeff(1))*(1-XCcoeff(2));
		impulse_m1( XClow(1)-1 , XClow(2)+1 ) = impulse_m1( XClow(1)-1 , XClow(2)+1 ) + (1-XCcoeff(1))*(XCcoeff(2));
		impulse_m1( XClow(1)   , XClow(2)+1 ) = impulse_m1( XClow(1)   , XClow(2)+1 ) + (XCcoeff(1))*(XCcoeff(2));
	% //---
		impulse_p1( XClow(1)+1 , XClow(2) )   = impulse_p1( XClow(1)+1 , XClow(2) ) + (1-XCcoeff(1))*(1-XCcoeff(2));
		impulse_p1( XClow(1)+2 , XClow(2) )   = impulse_p1( XClow(1)+2 , XClow(2) ) + (XCcoeff(1))*(1-XCcoeff(2));
		impulse_p1( XClow(1)+1 , XClow(2)+1 ) = impulse_p1( XClow(1)+1 , XClow(2)+1 ) + (1-XCcoeff(1))*(XCcoeff(2));
		impulse_p1( XClow(1)+2 , XClow(2)+1 ) = impulse_p1( XClow(1)+2 , XClow(2)+1 ) + (XCcoeff(1))*(XCcoeff(2));
	
	% // ... and convolute to compute fc
	fc_m1 = real(ifft2( fft2(atm).*fft2(impulse_m1) ));
	fc    = real(ifft2( fft2(atm).*fft2(impulse) ));
	fc_p1 = real(ifft2( fft2(atm).*fft2(impulse_p1) ));

	d1 = (X-(X_curr(1))).^2 + (Y-(X_curr(2))).^2;
	fc_true = HEIGHT*exp(-d1/(STDEV*STDEV)) + ...
	     HEIGHT*exp(-d2/(STDEV*STDEV)) + ...
	     HEIGHT*exp(-d3/(STDEV*STDEV)) + ...
	     HEIGHT*exp(-d4/(STDEV*STDEV)) + ...
	     HEIGHT*exp(-d5/(STDEV*STDEV)) ;

	% // ... also mask ea. fc
	resmask    = abs(fc)>tol;
	resmask_m1 = abs(fc_m1)>tol;
	resmask_p1 = abs(fc_p1)>tol;

	fcerrs(i) = cor (fc_true(resmask), fc(resmask));


	% // Compute the true CC at each position ...
	true_cc_p1 = corrcoef( fc_p1(resmask_p1==1) , fo(resmask_p1==1) );
	true_cc    = corrcoef(    fc(resmask==1)    , fo(resmask==1) );
	true_cc_m1 = corrcoef( fc_m1(resmask_m1==1) , fo(resmask_m1==1) );
	
	% // ... as well as the (numeric) derivative
	true_dcc_dx(i) = ( true_cc_p1 - true_cc_m1 ) / (2*STEPSIZE);
	
	% // stats over mask
	V = sum(sum(resmask));
	C = sum(fc(resmask));
	O = sum(fo(resmask));
	C2 = sum(fc(resmask).^2);
	O2 = sum(fo(resmask).^2);
	CO = sum( fo(resmask).*fc(resmask) );
	
	% // more stats over mask
	% // convolutions	
	ddx = zeros( size(X) ); 
	 ddx(length(X),length(X)) = 0; ddx(2,length(X)) = 0; 
	 ddx(length(X),1)         = 1; ddx(2,1)         = -1; 
	 ddx(length(X),2)         = 0; ddx(2,2)         = 0; 
	conv_atm_map  = real( ifft2(fft2(atm).*fft2(fo)) );
	 ddx_atm_map   = real( ifft2(fft2(ddx).*fft2(conv_atm_map)) ) / (2*STEPSIZE);
	conv_atm_fc  = real( ifft2(fft2(atm).*fft2(fc)) );
	 ddx_atm_fc   = real( ifft2(fft2(ddx).*fft2(conv_atm_fc)) ) / (2*STEPSIZE);


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
	delO  = ( Op1  - Om1  )    / (2*STEPSIZE);
	delO2 = ( O2p1 - O2m1 )  / (2*STEPSIZE);
	delC  = ( Cp1  - Cm1  )    / (2*STEPSIZE);
	delC2 = ( C2p1 - C2m1 )  / (2*STEPSIZE);
	delV =  ( Vp1  - Vm1  )    / (2*STEPSIZE);
	delCO = ( COp1 - COm1 )  / (2*STEPSIZE);

	if (delC > 1.0) 
		return
	end

	true_delC(i) = delC;

	varC = (C2 - C*C/V);
	varO = (O2 - O*O/V);

	f  =  ( V*CO - C*O ) / (V);
	g  =  ( sqrt(varO*varC) );

	f_prime1 =  sum(sum( ddx_atm_map .* currpos ) );
	g_prime1 =  0.0;
	
	f_prime2 =  sum(sum( ddx_atm_map .* currpos ) );
	g_prime2 = ( sqrt(varO)/sqrt(varC) * sum(sum( ddx_atm_fc .* currpos ))  ) ;
	
	f_prime3 =  sum(sum( ddx_atm_map .* currpos ) );
	g_prime3 = ( .5/sqrt(varC) * ( V*delC2 ) * sqrt(varO) ) / V;

	true_ccs1(i) =  f/g;
	true_ccs2(i) = true_cc;
	
	% //pred_dcc_dx1(i) =  (V * ddx_atm_map(currpos==1) ) / D;
	pred_dcc_dx1(i) =  (f_prime1*g - g_prime1*f)/(g*g);
	pred_dcc_dx2(i) =  (f_prime2*g - g_prime2*f)/(g*g);
	pred_dcc_dx3(i) =  (f_prime3*g - g_prime3*f)/(g*g);
end 

cor(true_ccs1, true_ccs2)

disp('---')
cor(true_dcc_dx,pred_dcc_dx1)
polyfit(true_dcc_dx,pred_dcc_dx1,1)
disp('---')
cor(true_dcc_dx,pred_dcc_dx2)
polyfit(true_dcc_dx,pred_dcc_dx2,1)
disp('---')
cor(true_dcc_dx,pred_dcc_dx3)
polyfit(true_dcc_dx,pred_dcc_dx3,1)
