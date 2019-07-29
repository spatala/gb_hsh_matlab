# Code Checks

## RHS of first equation on page 20.

1. `clebschgordan(j1,m1,j2,m2,j,m)`
	+ Code to compute Clebsh-Gordan coefficients (from Oliver)
2. `compute_complex_Yl(th,phi,l)`
	+ Compute the vector of spherical harmonics for a given theta(`th`) and phi (`phi`) for a fixed `l`. The output is a column vector with `m = -l,...0,...,m`.
3. `sph_harm_Mab_formula(ab_ord, l, th, phi)`
	+ Compute the sum for a fixed `e=l` in the RHS of the formula. The normalization by `PI_e` is done in the function.
4. 