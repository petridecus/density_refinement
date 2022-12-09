STEPSIZE = 0.5; % // must evenly divide 2

ds = [0:STEPSIZE:8,-6:STEPSIZE:-STEPSIZE];
[Xd,Yd,Zd] = ndgrid(ds);

% //params
% //RESO = 5

RS=[];
effRs=[];
CC_ala=[];
CC_leu=[];

for RESO=2:.5:15
	
	k = (pi/RESO)^2;
	C = (k/pi)^1.5;
	A = 8.0;
	
	x_N  = [-0.66,	0.86,	0.96];
	x_CA = [0,	0,	0];
	x_C  = [-0.45,	-1.45,	0.16];
	x_O  = [-0.94,	-2.06,	-0.8];
	x_CB = [1.52,	0.09,	0.16];
	x_CG = [2.24,	1.27,	-0.49];
	x_CD1 = [3.69,	1.26,	-0.08];
	x_CD2 = [2.11,	1.18,	-2.01];
	
	% // observed density
	d_N  = (Xd-(x_N(1))).^2  + (Yd-(x_N(2))).^2  + (Zd-(x_N(3))).^2;
	d_CA = (Xd-(x_CA(1))).^2 + (Yd-(x_CA(2))).^2 + (Zd-(x_CA(3))).^2;
	d_C  = (Xd-(x_C(1))).^2  + (Yd-(x_C(2))).^2  + (Zd-(x_C(3))).^2;
	d_O  = (Xd-(x_O(1))).^2  + (Yd-(x_O(2))).^2  + (Zd-(x_O(3))).^2;
	d_CB = (Xd-(x_CB(1))).^2 + (Yd-(x_CB(2))).^2 + (Zd-(x_CB(3))).^2;
	d_CG = (Xd-(x_CG(1))).^2 + (Yd-(x_CG(2))).^2 + (Zd-(x_CG(3))).^2;
	d_CD1 = (Xd-(x_CD1(1))).^2 + (Yd-(x_CD1(2))).^2 + (Zd-(x_CD1(3))).^2;
	d_CD2 = (Xd-(x_CD2(1))).^2 + (Yd-(x_CD2(2))).^2 + (Zd-(x_CD2(3))).^2;

	rho_ala = C*A*exp(-k*d_N) + C*A*exp(-k*d_CA) + C*A*exp(-k*d_C) + C*A*exp(-k*d_O) + C*A*exp(-k*d_CB);
	rho_leu = C*A*exp(-k*d_N) + C*A*exp(-k*d_CA) + C*A*exp(-k*d_C) + C*A*exp(-k*d_O) + C*A*exp(-k*d_CB) + C*A*exp(-k*d_CG) + C*A*exp(-k*d_CD1) + C*A*exp(-k*d_CD2);
	
	bestCC=0;
	bestK2 = -1;
	%//for effR=4:.2:12
	effR = 0.8*RESO+2.4;
	k2 = (pi/effR)^2;
	C2 = (k2/pi)^1.5;
	rho1 = C2*A*exp(-k2*d_CA);
	thisCC_ala = cor( reshape(rho_ala,1,[]), reshape(rho1,1,[]));
	thisCC_leu = cor( reshape(rho_leu,1,[]), reshape(rho1,1,[]));
	%//end

	%//[RESO,pi/sqrt(bestK2),bestCC];
	RS = [RS;RESO];
	effRs = [effRs;pi/sqrt(bestK2)];
	CC_ala = [CC_ala;thisCC_ala];
	CC_leu = [CC_leu;thisCC_leu];
end