STEPSIZE = 1.0; % // must evenly divide 8

x = [0:STEPSIZE:16];
[Y,X] = meshgrid(x,x);

ds = [0:STEPSIZE:8,-8:STEPSIZE:-STEPSIZE];
[yd,xd] = meshgrid(ds);

DS = sqrt(xd.^2+yd.^2);

tol = 1e-5;

MASKWIDTH = 3.0
SIGMOID_STEEPNESS  = 1;
STDEV     =  1

X1 = [8,8]; %//+ rand(1,2);
theta1 = rand()*360;
theta2 = mod( theta1+120 , 360 );
X2 = X1 + 1.4*[cosd(theta1),sind(theta1)]; 
X3 = X1 + 1.4*[cosd(theta2),sind(theta2)]; 
theta3 = mod( theta2-180 + 120 , 360);
theta4 = mod( theta2-180 + 240 , 360);
X4 = X3 + 1.5*[cosd(theta3),sind(theta3)]; 
X5 = X3 + 1.5*[cosd(theta4),sind(theta4)];

% // observed density
d1 = (X-(X1(1))).^2 + (Y-(X1(2))).^2;
d2 = (X-(X2(1))).^2 + (Y-(X2(2))).^2;
d3 = (X-(X3(1))).^2 + (Y-(X3(2))).^2;
d4 = (X-(X4(1))).^2 + (Y-(X4(2))).^2;
d5 = (X-(X5(1))).^2 + (Y-(X5(2))).^2;
fo = exp(-d1/(STDEV*STDEV)) + exp(-d2/(STDEV*STDEV)) + ...
     exp(-d3/(STDEV*STDEV)) + exp(-d4/(STDEV*STDEV)) + ...
     exp(-d5/(STDEV*STDEV)) ;

% // 1 = sharp mask
% // 2 = smooth mask
atmmask1 = abs(DS)<(MASKWIDTH);
atmmask2 = 1./(1 + exp(-SIGMOID_STEEPNESS * (MASKWIDTH-DS)));
atm1   = exp(-DS.^2/(STDEV)^2).*atmmask1;
atm2   = exp(-DS.^2/(STDEV)^2).*atmmask2;

X2low = floor( X2 / STEPSIZE + 1 );
X3low = floor( X3 / STEPSIZE + 1 );
X4low = floor( X4 / STEPSIZE + 1 );
X5low = floor( X5 / STEPSIZE + 1 );
X2coeff =  X2 / STEPSIZE + 1 - X2low;
X3coeff =  X3 / STEPSIZE + 1 - X3low;
X4coeff =  X4 / STEPSIZE + 1 - X4low;
X5coeff =  X5 / STEPSIZE + 1 - X5low;

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


for noffset = 1:201
	offset = 0.02*(noffset-101);
	X_dx = X1+[offset , 0];
	X_dy = X1+[0 , offset];

	% // put a spike at each atom location ...
	Xdxlow = floor( X_dx / STEPSIZE + 1 );
	Xdylow = floor( X_dy / STEPSIZE + 1 );
	Xdxcoeff =  X_dx / STEPSIZE + 1 - Xdxlow;
	Xdycoeff =  X_dy / STEPSIZE + 1 - Xdylow;

	impulse_dx = impulse;
	impulse_dy = impulse;
		impulse_dx( Xdxlow(1)   , Xdxlow(2) )   = impulse_dx( Xdxlow(1)   , Xdxlow(2) ) + (1-Xdxcoeff(1))*(1-Xdxcoeff(2));
		impulse_dx( Xdxlow(1)+1 , Xdxlow(2) )   = impulse_dx( Xdxlow(1)+1 , Xdxlow(2) ) + (Xdxcoeff(1))*(1-Xdxcoeff(2));
		impulse_dx( Xdxlow(1)   , Xdxlow(2)+1 ) = impulse_dx( Xdxlow(1)   , Xdxlow(2)+1 ) + (1-Xdxcoeff(1))*(Xdxcoeff(2));
		impulse_dx( Xdxlow(1)+1 , Xdxlow(2)+1 ) = impulse_dx( Xdxlow(1)+1 , Xdxlow(2)+1 ) + (Xdxcoeff(1))*(Xdxcoeff(2));
	% //---
		impulse_dy( Xdylow(1)   , Xdylow(2) )   = impulse_dy( Xdylow(1)   , Xdylow(2) ) + (1-Xdycoeff(1))*(1-Xdycoeff(2));
		impulse_dy( Xdylow(1)+1 , Xdylow(2) )   = impulse_dy( Xdylow(1)+1 , Xdylow(2) ) + (Xdycoeff(1))*(1-Xdycoeff(2));
		impulse_dy( Xdylow(1)   , Xdylow(2)+1 ) = impulse_dy( Xdylow(1)   , Xdylow(2)+1 ) + (1-Xdycoeff(1))*(Xdycoeff(2));
		impulse_dy( Xdylow(1)+1 , Xdylow(2)+1 ) = impulse_dy( Xdylow(1)+1 , Xdylow(2)+1 ) + (Xdycoeff(1))*(Xdycoeff(2));

	% // ... and convolute to compute fc
	fc_dx = real(ifft2( fft2(atm1).*fft2(impulse_dx) ));
	fc_dy = real(ifft2( fft2(atm1).*fft2(impulse_dy) ));

	d1x = (X-(X_dx(1))).^2 + (Y-(X_dx(2))).^2;
	d1y = (X-(X_dy(1))).^2 + (Y-(X_dy(2))).^2;
	d2 = (X-(X2(1))).^2 + (Y-(X2(2))).^2;
	d3 = (X-(X3(1))).^2 + (Y-(X3(2))).^2;
	d4 = (X-(X4(1))).^2 + (Y-(X4(2))).^2;
	d5 = (X-(X5(1))).^2 + (Y-(X5(2))).^2;
	fc_dx = exp(-d1x/(STDEV*STDEV)) + exp(-d2/(STDEV*STDEV)) + ...
	     exp(-d3/(STDEV*STDEV)) + exp(-d4/(STDEV*STDEV)) + ...
	     exp(-d5/(STDEV*STDEV)) ;
	fc_dy = exp(-d1y/(STDEV*STDEV)) + exp(-d2/(STDEV*STDEV)) + ...
	     exp(-d3/(STDEV*STDEV)) + exp(-d4/(STDEV*STDEV)) + ...
	     exp(-d5/(STDEV*STDEV)) ;

	% // ... also mask ea. fc
	resmask_dx_smooth = 1./(1 + exp(-SIGMOID_STEEPNESS * (fc_dx)));
	resmask_dx_sharp = abs(fc_dx)>tol;
	resmask_dy_smooth = 1./(1 + exp(-SIGMOID_STEEPNESS * (fc_dy)));
	resmask_dy_sharp = abs(fc_dx)>tol;

	% // wted CC
	Cdx_sm(noffset) = sum(sum( fc_dx .* resmask_dx_smooth ));
	Cdx_sh(noffset) = sum(sum( fc_dx .* resmask_dx_sharp ));
	Cdy_sm(noffset) = sum(sum( fc_dy .* resmask_dy_smooth ));
	Cdy_sh(noffset) = sum(sum( fc_dy .* resmask_dy_sharp ));
	Odx_sm(noffset) = sum(sum( fo .* resmask_dx_smooth ));
	Odx_sh(noffset) = sum(sum( fo .* resmask_dx_sharp ));
	Ody_sm(noffset) = sum(sum( fo .* resmask_dy_smooth ));
	Ody_sh(noffset) = sum(sum( fo .* resmask_dy_sharp ));
	Cdx2_sm(noffset) = sum(sum( fc_dx .* fc_dx .* resmask_dx_smooth ));
	Cdx2_sh(noffset) = sum(sum( fc_dx .* fc_dx .* resmask_dx_sharp ));
	Cdy2_sm(noffset) = sum(sum( fc_dy .* fc_dy .* resmask_dy_smooth ));
	Cdy2_sh(noffset) = sum(sum( fc_dy .* fc_dy .* resmask_dy_sharp ));
	Odx2_sm(noffset) = sum(sum( fo .* fo  .* resmask_dx_smooth ));
	Odx2_sh(noffset) = sum(sum( fo .* fo  .* resmask_dx_sharp ));
	Ody2_sm(noffset) = sum(sum( fo .* fo  .* resmask_dy_smooth ));
	Ody2_sh(noffset) = sum(sum( fo .* fo  .* resmask_dy_sharp ));
	Ndx_sm(noffset) = sum(sum( resmask_dx_smooth ));
	Ndx_sh(noffset) = sum(sum( resmask_dx_sharp ));
	Ndy_sm(noffset) = sum(sum( resmask_dy_smooth ));
	Ndy_sh(noffset) = sum(sum( resmask_dy_sharp ));
	OCdx_sm(noffset) = sum(sum( fc_dx .* fo  .* resmask_dx_smooth ));
	OCdx_sh(noffset) = sum(sum( fc_dx .* fo  .* resmask_dx_sharp ));
	OCdy_sm(noffset) = sum(sum( fc_dy .* fo  .* resmask_dy_smooth ));
	OCdy_sh(noffset) = sum(sum( fc_dy .* fo  .* resmask_dy_sharp ));

	cc_dx_sm(noffset) = (Ndx_sm(noffset)*OCdx_sm(noffset) - Odx_sm(noffset) .* Cdx_sm(noffset)) / ...
	                               sqrt( ( Ndx_sm(noffset)*Cdx2_sm(noffset) - Cdx_sm(noffset)*Cdx_sm(noffset) )* ...
	                                     ( Ndx_sm(noffset)*Odx2_sm(noffset) - Odx_sm(noffset)*Odx_sm(noffset) ) );
	cc_dx_sh(noffset) = (Ndx_sh(noffset)*OCdx_sh(noffset) - Odx_sh(noffset) .* Cdx_sh(noffset)) / ...
	                               sqrt( ( Ndx_sh(noffset)*Cdx2_sh(noffset) - Cdx_sh(noffset)*Cdx_sh(noffset) )* ...
	                                     ( Ndx_sh(noffset)*Odx2_sh(noffset) - Odx_sh(noffset)*Odx_sh(noffset) ) );
	cc_dy_sm(noffset) = (Ndy_sm(noffset)*OCdy_sm(noffset) - Ody_sm(noffset) .* Cdy_sm(noffset)) / ...
	                               sqrt( ( Ndy_sm(noffset)*Cdy2_sm(noffset) - Cdy_sm(noffset)*Cdy_sm(noffset) )* ...
	                                     ( Ndy_sm(noffset)*Ody2_sm(noffset) - Ody_sm(noffset)*Ody_sm(noffset) ) );
	cc_dy_sh(noffset) = (Ndy_sh(noffset)*OCdy_sh(noffset) - Ody_sh(noffset) .* Cdy_sh(noffset)) / ...
	                               sqrt( ( Ndy_sh(noffset)*Cdy2_sh(noffset) - Cdy_sh(noffset)*Cdy_sh(noffset) )* ...
	                                     ( Ndy_sh(noffset)*Ody2_sh(noffset) - Ody_sh(noffset)*Ody_sh(noffset) ) );
	offsets(noffset) = offset;
end 

 plot( offsets, cc_dx_sh,'r', offsets, cc_dx_sm,'b', offsets, cc_dy_sh,'g', offsets, cc_dy_sm,'k')
