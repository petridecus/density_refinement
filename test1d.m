SS = 1;
X=[0:SS:8.0-SS];

for i=1:401
	x1 = 0.01*(i-1)+2.0;
	x2 = 5.07;
	impulse = zeros(size(X));

	X1l = floor( x1 / SS + 1 );
	X1c = x1 / SS + 1 - X1l;
	X2l = floor( x2 / SS + 1 );
	X2c = x2 / SS + 1 - X2l;
	
	impulse(X1l) = (1-X1c);
	impulse(X1l+1) = (X1c);
	impulse(X2l) = impulse(X2l) + (1-X2c);
	impulse(X2l+1) = impulse(X2l+1) + (X2c);
	
	atm = exp(-0.5*min(X,8.0-X).^2/(.4)^2);

	fc1 = real(ifft( fft(atm) .* fft(impulse) ));
	fc2 = exp(-0.5*(X-x1).^2/(.4)^2) + exp(-0.5*(X-x2).^2/(.4)^2);
	fc3 = real(ifft( fft(atm) .* fft(impulse.*impulse) ));

	a(i) = sum(fc1.^2);
	b(i) = sum(fc2.^2);
	c(i) = sum(fc3);
end