function W = hann2_circular_shrink(nx, nz, gamma, dx, dz)

x = ((1:nx) - (nx+1)/2) * dx;
z = ((1:nz) - (nz+1)/2) * dz;
[X, Z] = meshgrid(x, z);
r = sqrt(X.^2 + Z.^2);

Rx = max(abs(x));
Rz = max(abs(z));
R = gamma * min(Rx, Rz);

W = zeros(nz, nx);
mask = (r <= R);
W(mask) = 0.5 * (1 + cos(pi * r(mask) / R));

end