# Function Descriptions

## `covert_gbrots`

Convert `(g1,g2)` rotation matrices (of GBs) into 
quaternions and Cayley-Klein parameters

- Input:
  + rot_mats: Array of size `3X6XN`
  		`N` GB parameters in `SO3XSO3` parameterization

- Output:
  + `s`: struct of N-arrays (8 elements)
  	- ((`s.a1`, `s.lb1`, `s.q1`, `s.Q1`), (`s.a2`, `s.lb2`, `s.q2`, `s.Q2`))
  	- `s.a`  :   `q_0 - i*q_3`
  	- `s.lb` : `-q_2 - i*q_1`
  	- `s.q`  : `q_0`
  	- `s.Q`  : `(q_1, q_2, q_3)`

- Notes:
	+ This struct is used to compute MBP basis functions

## `compute_Mvec`

Computes MBP basis functions for given set of (`a`, `b`, `gamma`, `alpha`, `beta`)
for a given order (`a+b \leq Nmax`) of symmetrized MBP basis function.

- Input:
  + s: struct of N(=nrots)-arrays (8 elements)
    - ((`s.a1`, `s.lb1`, `s.q1`, `s.Q1`), (`s.a2`, `s.lb2`, `s.q2`, `s.Q2`))
  + nrots: Number of GB parameters
  + num_cols: Number of columns in the given order (a+b \leq Nmax) of 
      symmetrized MBP basis function
  + vec_inds: The non-zero `indices` in the symmetrized MBP basis 
      functions.

- Output:
  + Mvec: array of size `nrots x num_cols`
    - The rows coresspond to different rotations.
    - The columns correspond to all possible combinations of 
      `a`, `b`, `gamma`, `alpha`, `beta`.
    - Only the non-zero indices (stored in `vec_inds`) are calculated.

- Notes:
  + The idea is to compute MBP functions only for those combinations of 
      (`a`, `b`, `gamma`, `alpha`, `beta`) that are necessary for the 
      symmetrized MBP basis functions.
  + For any given order of symmetrized MBP basis functions 
      (`a+b \leq Nmax`), we can determine all the non-zero indices 
      (corresponding to unique a, b, gamma, alpha, beta). 
      This is stored in the array `vec_inds`.


## `mbp_basis`

Returns the basis functions $M^{a,b}$ for the grain boundary space
using the MBP parameterization.

- Input:
  + a, b: Indices for $M^{a,b}$ function.
  + mbp_angs: (omega_m, theta_m, phi_m, omega_b, phi_b)

- Output:
  + M: MBP basis function with order (a,b). 
	- Rows ordered by (\gamma, \alpha, \beta) in lexicographic order, i.e., starting with negative values, ending with positive values).

- Notes:
  + Follows Equation () of the manuscript (finalize eqn. number after publishing).


## `calc_Mvec_full`


Returns the full and symmetrized basis functions for the grain boundary 
space using the SO3 X SO3 parameterization.

- Input:
  + top_dir: Directory to find the symmetrized basis functions
  + pt_grp: The crystallographic point group ('Oh')
  + Nmax: all basis functions such that a+b \leq Nmax
  + coeffs_typ: 'aPLUSb_max'
  + null_typ: 'zero' or 'const'
    - 'const': Crystal and GES and constant no-boundary condition
    - 'null':  Crystal and GES and null no-boundary condition
  + rots: gA and gB, orientations of the two crystals defining the GB

- Output:
  + Mvec: Full basis vectors
  + SMvec: Symmetrized basis vectors



## `calc_Mvec_symm`

Computes symmetrized MBP basis functions for given set of
(a, b, gamma, alpha, beta) for a given order (a+b \leq Nmax)
of symmetrized MBP basis function.

- Input:
  + top_dir: Directory to find the symmetrized basis functions
  + pt_grp: The crystallographic point group ('Oh')
  + Nmax: all basis functions such that a+b \leq Nmax
  + coeffs_typ: 'aPLUSb_max'
  + symm_typ: string ('cryst' or 'cryst_ges' or 'const' or 'zero')
    - specify which level of symmetry to use to compute basis functions
    - 'cryst': Only crystal point group symmetries
    - 'cryst_ges': Crystal and GES
    - 'const': Crystal and GES and constant no-boundary condition
    - 'null':  Crystal and GES and null no-boundary condition
  + s_angs: struct of N(=nrots)-arrays (8 elements)
    - ((s.a1, s.lb1, s.q1, s.Q1), (s.a2, s.lb2, s.q2, s.Q2))

- Output:
  + SMvec: symmetrized basis functions
    - The rows coresspond to different rotations.
    - The columns correspond to the number of symmetrized basis functions.


