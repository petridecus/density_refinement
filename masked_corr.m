function pcc=masked_corr(p1,p2,p_mask)
%

sumCO = sum(sum( p1.*p2.*(p_mask) ));
sumO  = sum(sum( p1.*(p_mask) ));
sumO2 = sum(sum( p1.*p1.*(p_mask) ));
sumC  = sum(sum( p2.*(p_mask) ));
sumC2 = sum(sum( p2.*p2.*(p_mask) ));
sumV  = sum(sum( p_mask ));

varC = (sumC2 - sumC*sumC / sumV );
varO = (sumO2 - sumO*sumO / sumV );
pcc = (sumCO - sumC*sumO/sumV) / sqrt( varC*varO );
