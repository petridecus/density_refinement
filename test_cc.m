STEPSIZE = 1; % // must evenly divide 8

x = [0:STEPSIZE:16];
[Y,X] = meshgrid(x,x);

ds = [0:STEPSIZE:8,-8:STEPSIZE:-STEPSIZE];
[yd,xd] = meshgrid(ds);

% //params
H = 0.1
RESO = 5
ATOM_MASK = 6
k = (pi/RESO)
C = (k/pi)^1.5
A = 8.0
tol = 0.000001

for i=1:100
	offset = 4*rand(1,2)-2;
	off_dist(i) = norm(offset);

	% // true atom locations
	X1 = [8,8] + rand(1,2);
	theta1 = rand()*360;
	thetas(i) = theta1;
	theta2 = mod( theta1+120 , 360 );
	X2 = X1 + 1.4*[cosd(theta1),sind(theta1)]; 
	X3 = X1 + 1.4*[cosd(theta2),sind(theta2)]; 
	theta3 = mod( theta2-180 + 120 , 360);
	theta4 = mod( theta2-180 + 240 , 360);
	X4 = X3 + 1.5*[cosd(theta3),sind(theta3)]; 
	X5 = X3 + 1.5*[cosd(theta4),sind(theta4)];
	X_curr = X1+offset;
	
	% // observed density
	d1 = (X-(X1(1))).^2 + (Y-(X1(2))).^2;
	d2 = (X-(X2(1))).^2 + (Y-(X2(2))).^2;
	d3 = (X-(X3(1))).^2 + (Y-(X3(2))).^2;
	d4 = (X-(X4(1))).^2 + (Y-(X4(2))).^2;
	d5 = (X-(X5(1))).^2 + (Y-(X5(2))).^2;
	fo = 5*randn(size(X));%//C*A*exp(-k*d1) + C*A*exp(-k*d2) + C*A*exp(-k*d3) + C*A*exp(-k*d4) + C*A*exp(-k*d5) + .1*randn(size(X));
	
	dcurr = (X-(X_curr(1))).^2 + (Y-(X_curr(2))).^2;
	dcurr_p1 = (X-(X_curr(1)-H)).^2 + (Y-(X_curr(2))).^2;
	dcurr_m1 = (X-(X_curr(1)+H)).^2 + (Y-(X_curr(2))).^2;
	fc = C*A*exp(-k*dcurr) + C*A*exp(-k*d2) + C*A*exp(-k*d3) + C*A*exp(-k*d4) + C*A*exp(-k*d5);
	fc_p1 = C*A*exp(-k*dcurr_p1) + C*A*exp(-k*d2) + C*A*exp(-k*d3) + C*A*exp(-k*d4) + C*A*exp(-k*d5);
	fc_m1 = C*A*exp(-k*dcurr_m1) + C*A*exp(-k*d2) + C*A*exp(-k*d3) + C*A*exp(-k*d4) + C*A*exp(-k*d5);

	atm = ( C*A*exp(-k*dcurr) );
	ddx = ( 2*(X-X_curr(1))*(-k) );
	d_fx(i) = sum(sum( ddx .* atm .* fo ) );
	d_fc(i) = sum(sum( 2 * ddx .* atm .* fc ) );

	% // ... also mask ea. fc
	resmask_2345 = ( 1 - 1./(1+exp( ( d2 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) ) .* ...
	               ( 1 - 1./(1+exp( ( d3 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) ) .* ...
	               ( 1 - 1./(1+exp( ( d4 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) ) .* ...
	               ( 1 - 1./(1+exp( ( d5 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) );
	inv_atmmask    = 1./(1+exp( ( dcurr ) - (ATOM_MASK-1)*(ATOM_MASK-1) ));
	resmask    = 1 - (( 1 - inv_atmmask ) .* resmask_2345 ); %//abs(fc)>tol;
	resmask_m1 = 1 - (( 1 - 1./(1+exp( ( dcurr_m1 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) ) .* resmask_2345 ); %//abs(fc)>tol;
	resmask_p1 = 1 - (( 1 - 1./(1+exp( ( dcurr_p1 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) ) .* resmask_2345 ); %//abs(fc)>tol;

	d_fo(i) = sum(sum( 2*(X_curr(1)-X) .* inv_atmmask .* (1-resmask) .* fo ) );
	d_fo2(i) = sum(sum( 2*(X_curr(1)-X) .* inv_atmmask .* (1-resmask) .* fo .* fo ) );

	% // Compute the true CC at each position ...
	V_p1 = sum(sum(resmask_p1));
	C_p1 = sum(sum(fc_p1.*resmask_p1));
	O_p1 = sum(sum(fo.*resmask_p1));
	C2_p1 = sum(sum( (fc_p1.^2).*resmask_p1));
	O2_p1 = sum(sum( (fo.^2).*resmask_p1));
	OC_p1 = sum(sum( fo.*fc_p1 .*resmask_p1));
	V_m1 = sum(sum(resmask_m1));
	C_m1 = sum(sum(fc_m1.*resmask_m1));
	O_m1 = sum(sum(fo.*resmask_m1));
	C2_m1 = sum(sum( (fc_m1.^2).*resmask_m1));
	O2_m1 = sum(sum( (fo.^2).*resmask_m1));
	OC_m1 = sum(sum( fo.*fc_m1 .*resmask_m1));

	varO_p1 = ( O2_p1 - O_p1.*O_p1 / V_p1 );
	varC_p1 = ( C2_p1 - C_p1.*C_p1 / V_p1 );
	true_cc_p1 = ( OC_p1 - O_p1.*C_p1 / V_p1 ) / sqrt( varC_p1 * varO_p1 );
	varO_m1 = ( O2_m1 - O_m1.*O_m1 / V_m1 );
	varC_m1 = ( C2_m1 - C_m1.*C_m1 / V_m1 );
	true_cc_m1 = ( OC_m1 - O_m1.*C_m1 / V_m1 ) / sqrt( varC_m1 * varO_m1 );

	% // ... as well as the (numeric) derivative
	true_dcc_dx(i) = ( true_cc_p1 - true_cc_m1 ) / (2*H);

	% // pred_cc
	%pred_dcc_dx1(i) = d_fx;

	V_ = sum(sum(resmask));
	C_ = sum(sum(fc.*resmask));
	O_ = sum(sum(fo.*resmask));
	C2_ = sum(sum( (fc.^2).*resmask));
	O2_ = sum(sum( (fo.^2).*resmask));
	OC_ = sum(sum( fo.*fc .*resmask));

	true_delOC(i) = (OC_p1 - OC_m1) / (2*H);
	true_delC2(i) = (C2_p1 - C2_m1) / (2*H);
	true_delO2(i) = (O2_p1 - O2_m1) / (2*H);
	true_delO(i) = (O_p1 - O_m1) / (2*H);
	true_delC(i) = (C_p1 - C_m1) / (2*H);

	f    = ( OC_ - O_.*C_ / V_ );
	varO = ( O2_ - O_.*O_ / V_ );
	varC = ( C2_ - C_.*C_ / V_ );
	g = sqrt( varC * varO );

	dfdx1 =  d_fx(i) - d_fo(i)*C_/V_;
	dgdx1 = 0.5 * ( sqrt(varO)/sqrt(varC) * d_fc(i)  +  sqrt(varC)/sqrt(varO) * (d_fo2(i) + 2*(O_/V_)*d_fo(i))  ) ;

	dfdx2 = d_fx(i);
	dgdx2 = dgdx1;

	pred_dcc_dx1(i) =  (dfdx1*g - dgdx1*f)/(g*g);
	pred_dcc_dx2(i) =  (dfdx2*g - dgdx2*f)/(g*g);
end 

cor(true_dcc_dx,pred_dcc_dx1)
polyfit(true_dcc_dx,pred_dcc_dx1,1)
disp('---')
cor(true_dcc_dx,pred_dcc_dx2)
polyfit(true_dcc_dx,pred_dcc_dx2,1)
disp('***')
cor(d_fo,true_delO)
polyfit(d_fx,true_delOC,1)
disp('---')
cor(d_fx,true_delOC)
polyfit(d_fx,true_delOC,1)
disp('---')
cor(d_fc,true_delC2)
polyfit(d_fc,true_delC2,1)
