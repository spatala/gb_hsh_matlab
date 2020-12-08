# Code Descriptions

## Folders
1. `compute_mvec`
2. `mbp_indices`
3. `misc_funcs`
4. `param_conversions`
5. `rotation_Umat`
6. `symmetry_generators`

## `compute_mvec`
1. `mbp_basis`: Returns the basis functions `$M^{a,b}$` for the grain boundary space using the MBP parameterization.
2. `convert_gbrots`: Convert (g1,g2) rotation matrices (of GBs) into quaternions and Cayley-Klein parameters
3. `compute_Mvec`: Computes MBP basis functions for given set of (a, b, gamma, alpha, beta) for a given order (a+b \leq Nmax) of symmetrized MBP basis function.
4. `compute_Mvec_full`: Computes MBP basis functions for given set of (a, b, gamma, alpha, beta) for a given order (a+b \leq Nmax) of symmetrized MBP basis function.
5. `calc_Mvec_symm`: % Computes symmetrized MBP basis functions for given set of(a, b, gamma, alpha, beta) for a given order (a+b \leq Nmax) of symmetrized MBP basis function.


## symmetry_generators
1. `get_symm_rots.m`
2. `get_symmgen_angs.m`: Function to set the generators (as angles - (w, th, ph)) for each point group (only proper rotation elements are set.)
3. `get_symmgen_mats.m`: Function to set the generators (as matrices) for each point group (only proper rotation elements are set.)

## mbp_indices
1. `mbp_ab.m`: Get indices for MBP functions with fixed a and b.
2. `mbp_inds_array.m`: Get indices for MBP function with values of (a,b) given by `symm_orders`.

## misc_funcs
1. `clean.m` : Attempts to minimize the number of nonzero coefficients of the eigenvectors in `S` while maintaining orthonormality.
2. `clebsch_gordan.m`
3. `spnull.m`

## param_conversions
1. `check_conv.m`
2. `mbp_to_rots.m`: Convert GB parameters from Misorientation-Boundary plane (MBP) to SO3xSO3 parametrization (rots).
3. `rotmat_to_angs.m`
4. `rots_to_angs.m`
5. `rots_to_mbp.m`


## rotation_Umat

This folder contains functions to construct SO(3) irreducible representatives (matrices).

### Functions
1. `angles_to_CK.m`:returns `U^{1/2}` for the rotation given by a rotation angle and axis. Rows and cols ordered by increasing `m'` and `m`.
2. `CK_to_angles.m`: returns the rotation angle and axis for the given `U^{1/2}` with rows and cols ordered by increasing m' and m.
3. `rotation_mat.m`: returns the rotation matrix `U` of the order `j`, with a rotation angle `w` and rotation axis given by polar and azimuthal angles `th` and `ph`.
4. `wigner_little_d.m`: returns the Wigner `d` matrix for the given order and rotation angle, with rows and columns labeled in increasing values of `m`. The second argument is an estimate of the magnitude of numerical error.

### Functions To be modified or removed:
1. **To be modified**
	+ `rotation_wo_svd.m`: Function to compute elements of SO(3) Irreducible representative  `U^{a_val}_{alp_val, al_val}` for a rotation given by (q,Q) quaternion.
2. **To be removed**
	+ `rotation.m`: Constructs irreducible representation of SO(3) of dimension (N+1) for the rotation specified by axis `ax` and angle `an`.





# Functions List

## Util_functions

### To Document

1. `spnull`
2. `sp_null`
3. `rotation`

### Doucmented

1. `get_symmgen_mats`
2. `get_symm_rots`
3. `check_conv`
4. `mbp_inds_ab_array`
5. `mbp_ab`
6. `mbp_basis`
7. `convert_gbrots`
8. `compute_Mvec`
9. `clean`
10. `clebsch_gordan`
11. `sp_orth`
12. `mbp_to_rots`
13. `rotmat_to_angs`
14. `rots_to_angs`
15. `rots_to_mbp`
16. `angles_to_CK`
17. `CK_to_angles`
18. `rotation_mat`
19. `rotation_wo_svd`
20. `wigner_little_d`
21. `get_symmgen_angs`

 












