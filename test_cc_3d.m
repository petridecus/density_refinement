[Y,X,Z] = meshgrid(1:10,1:10,1:9);

% //params
H = 0.0005
RESO = 5
ATOM_MASK = 5
k = (pi/RESO)^2
C = (k/pi)^1.5
A = 8.0
tol = 0.00001

for i=1:1
 	X_curr = X_1;
	X2 = X_2;
	X3 = X_3;
	X4 = X_4;
	X5 = X_4;

	% // observed density
	fo = rhoObs ; %//5*randn(size(X)); %//C*A*exp(-k*d1) + C*A*exp(-k*d2) + C*A*exp(-k*d3) + C*A*exp(-k*d4) + C*A*exp(-k*d5) + .1*randn(size(X));

	% // calc density
	d2 = 4 * ( (X-(X2(1))).^2 + (Y-(X2(2))).^2 + (Z-(X2(3))).^2 );
	d3 = 4 * ( (X-(X3(1))).^2 + (Y-(X3(2))).^2 + (Z-(X3(3))).^2 );
	d4 = 4 * ( (X-(X4(1))).^2 + (Y-(X4(2))).^2 + (Z-(X4(3))).^2 );
	d5 = 4 * ( (X-(X5(1))).^2 + (Y-(X5(2))).^2 + (Z-(X5(3))).^2 );
	dcurr    = 4 * ( (X-(X_curr(1))).^2   + (Y-(X_curr(2))).^2 + (Z-(X_curr(3))).^2 );
	dcurr_p1 = 4 * ( (X-(X_curr(1)-H)).^2 + (Y-(X_curr(2))).^2 + (Z-(X_curr(3))).^2 );
	dcurr_m1 = 4 * ( (X-(X_curr(1)+H)).^2 + (Y-(X_curr(2))).^2 + (Z-(X_curr(3))).^2 );

	fc = C*A*exp(-k*dcurr) + C*A*exp(-k*d2) + C*A*exp(-k*d3) + C*A*exp(-k*d4) + C*A*exp(-k*d5);
	fc_p1 = C*A*exp(-k*dcurr_p1) + C*A*exp(-k*d2) + C*A*exp(-k*d3) + C*A*exp(-k*d4) + C*A*exp(-k*d5);
	fc_m1 = C*A*exp(-k*dcurr_m1) + C*A*exp(-k*d2) + C*A*exp(-k*d3) + C*A*exp(-k*d4) + C*A*exp(-k*d5);

	atm = ( C*A*exp(-k*dcurr) );
	ddx = ( 2*2*(X-X_curr(1))*(-k) );
	d_fx(i) = sum(sum(sum( ddx .* atm .* fo )));
	d_fc(i) = sum(sum(sum( 2 * ddx .* atm .* fc )));

	% // ... also mask ea. fc
	resmask_2345 = ( 1 - 1./(1+exp( ( d2 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) ) .* ...
	               ( 1 - 1./(1+exp( ( d3 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) ) .* ...
	               ( 1 - 1./(1+exp( ( d4 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) ) .* ...
	               ( 1 - 1./(1+exp( ( d5 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) );
	inv_atmmask    = 1./(1+exp( ( dcurr ) - (ATOM_MASK-1)*(ATOM_MASK-1) ));
	resmask    = 1 - (( 1 - inv_atmmask ) .* resmask_2345 ); %//abs(fc)>tol;
	resmask_m1 = 1 - (( 1 - 1./(1+exp( ( dcurr_m1 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) ) .* resmask_2345 ); %//abs(fc)>tol;
	resmask_p1 = 1 - (( 1 - 1./(1+exp( ( dcurr_p1 ) - (ATOM_MASK-1)*(ATOM_MASK-1) )) ) .* resmask_2345 ); %//abs(fc)>tol;

	d_fo(i) = sum(sum(sum( 2*2*(X_curr(1)-X) .* inv_atmmask .* (1-resmask) .* fo )));
	d_fo2(i) = sum(sum(sum( 2*2*(X_curr(1)-X) .* inv_atmmask .* (1-resmask) .* fo .* fo )));

	% // Compute the true CC at each position ...
	V_p1  = sum(sum(sum( resmask_p1 )));
	C_p1  = sum(sum(sum( fc_p1.*resmask_p1 )));
	O_p1  = sum(sum(sum( fo.*resmask_p1 )));
	C2_p1 = sum(sum(sum( (fc_p1.^2).*resmask_p1 )));
	O2_p1 = sum(sum(sum( (fo.^2).*resmask_p1 )));
	OC_p1 = sum(sum(sum( fo.*fc_p1.*resmask_p1 )));
	V_m1  = sum(sum(sum( resmask_m1 )));
	C_m1  = sum(sum(sum( fc_m1.*resmask_m1 )));
	O_m1  = sum(sum(sum( fo.*resmask_m1 )));
	C2_m1 = sum(sum(sum( (fc_m1.^2).*resmask_m1 )));
	O2_m1 = sum(sum(sum( (fo.^2).*resmask_m1 )));
	OC_m1 = sum(sum(sum( fo.*fc_m1 .*resmask_m1 )));

	varO_p1    = ( O2_p1 - O_p1.*O_p1 / V_p1 );
	varC_p1    = ( C2_p1 - C_p1.*C_p1 / V_p1 );
	true_cc_p1 = ( OC_p1 - O_p1.*C_p1 / V_p1 ) / sqrt( varC_p1 * varO_p1 );
	varO_m1    = ( O2_m1 - O_m1.*O_m1 / V_m1 );
	varC_m1    = ( C2_m1 - C_m1.*C_m1 / V_m1 );
	true_cc_m1 = ( OC_m1 - O_m1.*C_m1 / V_m1 ) / sqrt( varC_m1 * varO_m1 );

	% // ... as well as the (numeric) derivative
	true_dcc_dx(i) = ( true_cc_p1 - true_cc_m1 ) / (2*2*H);

	% // pred_cc
	% //pred_dcc_dx1(i) = d_fx;

	V_  = sum(sum(sum( resmask )));
	C_  = sum(sum(sum( fc.*resmask )));
	O_  = sum(sum(sum( fo.*resmask )));
	C2_ = sum(sum(sum( (fc.^2).*resmask )));
	O2_ = sum(sum(sum( (fo.^2).*resmask )));
	OC_ = sum(sum(sum( fo.*fc .*resmask )));

	true_delOC(i) = (OC_p1 - OC_m1) / (2*2*H);
	true_delC2(i) = (C2_p1 - C2_m1) / (2*2*H);
	true_delO2(i) = (O2_p1 - O2_m1) / (2*2*H);
	true_delO(i) = (O_p1 - O_m1) / (2*2*H);
	true_delC(i) = (C_p1 - C_m1) / (2*2*H);
	true_delV(i) = (V_p1 - V_m1) / (2*2*H);

	d_V(i) = sum(sum(sum( 2*2*(X_curr(1)-X) .* inv_atmmask .* (1-resmask) )));

	f    = ( OC_ - O_.*C_ / V_ );
	varO = ( O2_ - O_.*O_ / V_ );
	varC = ( C2_ - C_.*C_ / V_ );
	g = sqrt( varC * varO );

	true_delF(i) = (( OC_p1 - O_p1.*C_p1 / V_p1 ) - ( OC_m1 - O_m1.*C_m1 / V_m1 )) / (2*2*H);
	true_delG(i) = (sqrt( varC_p1 * varO_p1 ) - sqrt( varC_m1 * varO_m1 )) / (2*2*H);

d_fo=-35.5239;
d_fx=-0.2365;
d_fc=0.19059;
d_fo2=-156.019;
d_V=-8.31459;

O_ = 188.665;
C_ = 4.9608;
O2_ = 754.053;
C2_ = 1.59895;
CO_ = 20.2422;
V_ = 67.7138;

varO = ( O2_ - O_.*O_ / V_ );
varC = ( C2_ - C_.*C_ / V_ );

	pred_delF(i) = d_fx(i) - (1/V_^2) * ( d_fo(i)*C_*V_ - O_*C_*d_V(i));
	pred_delG(i) = 0.5 * ( ...
	      sqrt(varO)/sqrt(varC) * ( d_fc(i) + ( C_*C_*d_V(i)/V_^2 ) )  +  ...
	      sqrt(varC)/sqrt(varO) * ( d_fo2(i) - ( (1/V_^2) * ( 2*V_*O_*d_fo(i) - O_*O_*d_V(i) ) ) ) );
	dfdx1 = pred_delF(i);
	dgdx1 = pred_delG(i);

	pred_dcc_dx1(i) =  (dfdx1*g - dgdx1*f)/(g*g);
end 

pred_dcc_dx1
true_dcc_dx

disp('***')

pred_delF
pred_delG

disp('***')

O_
C_
O2_
C2_
CO_
V_

disp('***')

d_fx
d_fc
d_fo
d_fo2
d_V