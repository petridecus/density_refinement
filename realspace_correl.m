function [rscc]=realspace_correl(rho_1,rho_2,mask)
%

c = corrcoef( rho_1(mask),rho_2(mask) );
rscc = c(1,2);