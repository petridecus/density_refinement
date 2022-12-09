coords = [[2.0, 1.1]; ...
          [4.2, 3.7]; ...
          [2.2, 1.3]];

[rho, mask] = rhoc(coords, [ 5 , 5 ], .01, rand());

imagesc(rho);