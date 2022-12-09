CA1 = [ 40.439  46.023  48.742  ];
C1 = [  39.361  45.009  48.382  ] - CA1;
N1 = [  41.793  45.458  48.727 ] - CA1;
CA2 = [ 38.536  43.316  46.944 ] ;
C2 = [  38.270  42.206  47.973  ] - CA2;
N2 = [  39.578  44.214  47.353 ] - CA2;


C1 = [-1.036 , -0.257 , 1.11 ]
N1 = [-0.423 , 0.874 , -1.081 ]
C2 = [-0.99 , -1.056 , -0.521 ]
N2 = [ 1.017 , 0.269 , -1.006 ]

X1 = (C1+N1) / norm(C1+N1)
X2 = (C2+N2) / norm(C2+N2)

beta1 = acos( X1(3) );
beta2 = acos( X2(3) );
gamma1 = atan2( X1(1) , X1(2) );
gamma2 = atan2( X2(1) , X2(2) );

sin(beta1)
cos(beta1)
sin(beta2)
cos(beta2)

sin(gamma1)
cos(gamma1)
sin(gamma2)
cos(gamma2)

O1 = [  38.362  44.875  49.092 ];
O2 = [  37.129  41.742  48.147 ];

R1gb = [  cos(gamma1) , cos(beta1)*sin(gamma1) , sin(beta1)*sin(gamma1) ; ...
		 -sin(gamma1) , cos(beta1)*cos(gamma1) , sin(beta1)*cos(gamma1) ; ...
		  0           , -sin(beta1)            , cos(beta1) ]
R2gb = [  cos(gamma2) , cos(beta2)*sin(gamma2) , sin(beta2)*sin(gamma2) ; ...
		 -sin(gamma2) , cos(beta2)*cos(gamma2) , sin(beta2)*cos(gamma2) ; ...
		  0           , -sin(beta2)            , cos(beta2) ]

Rgb_N1 = R1gb^-1*N1';
Rgb_N2 = R2gb^-1*N2';

alpha1 = atan2( Rgb_N1(1) , Rgb_N1(2) )
alpha2 = atan2( Rgb_N2(1) , Rgb_N2(2) )

R1 =  [ - sin(alpha1)*cos(beta1)*sin(gamma1) + cos(alpha1)*cos(gamma1) , cos(alpha1)*cos(beta1)*sin(gamma1) + sin(alpha1)*cos(gamma1) , sin(beta1)*sin(gamma1) ; ...
		- sin(alpha1)*cos(beta1)*cos(gamma1) - cos(alpha1)*sin(gamma1) , cos(alpha1)*cos(beta1)*cos(gamma1) - sin(alpha1)*sin(gamma1) , sin(beta1)*cos(gamma1) ; ...
		  sin(alpha1)*sin(beta1) , -cos(alpha1)*sin(beta1) , cos(beta1)  ]

R2 =  [ - sin(alpha2)*cos(beta2)*sin(gamma2) + cos(alpha2)*cos(gamma2) , cos(alpha2)*cos(beta2)*sin(gamma2) + sin(alpha2)*cos(gamma2) , sin(beta2)*sin(gamma2) ; ...
		- sin(alpha2)*cos(beta2)*cos(gamma2) - cos(alpha2)*sin(gamma2) , cos(alpha2)*cos(beta2)*cos(gamma2) - sin(alpha2)*sin(gamma2) , sin(beta2)*cos(gamma2) ; ...
		  sin(alpha2)*sin(beta2) , -cos(alpha2)*sin(beta2) , cos(beta2)  ]

R = R2*(R1^-1)

err1 = norm( O2' - (R*(O1-CA1)'+CA2'))
err2 = norm(R*C1' - C2' )
err3 = norm( R*N1' - N2' )