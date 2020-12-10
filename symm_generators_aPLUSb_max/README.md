# Function Descriptions

## `generate_gb_cryst_symm`

- Inputs:
  + pt_grp	: Point group symmetry of the underlying crystal
  + Nmax	: Corresponds to max(a,b) for the basis functions in M^{a,b}
  + TOL		: threshold for a pivot to be considered significant in sp_null code


## `generate_cryst_symm_ab`

Basis functions symmetrized with crystal symmetries

## `generate_ges_mat`

Generate the matrix operator for grain exchange symmetry

## `combine_cryst_ges`

Generate symmetrized basis functions with both crystal point-group and grain-exchange symmetries applied

## `generate_gb_null`

Generate the matrix operator for Null-boundary singularity (i.e. the functions take the value **zero** at zero misorientation)

## `combine_cryst_ges_gbnull`

Generate symmetrized basis functions with both crystal point-group and grain-exchange symmetries applied


## `generate_yp_left_ab`

The matrix equivalent for `(g1, g2) -> (Ypi*g1, Ypi*g2)`
`M_ab(Ypi*g1, Ypi*g2)` = `M_ab(g1,g2) * ypileft_mat`

- Input:
  + a,b: Integers. Order for M_{a,b} function.

- Output:
  + ypileft_mat:
      Matrix operation with size `(Na+1)(Nb+1)(Nc+1)`
      `Na = 2a`; `Nb = 2b`; `Nc = 2*min(a,b)`

- Notes:
  + For a column corresponding to `gamma1, alpha1, beta1`,
  the row corresponding to `-gamma1, alpha1, beta1` contains
  the value `(-1)^(a+b)`


## `generate_gb_null`


## `save_symmvec_MabInds`

Save, "tot_Uprops" and "ind_ranges", that store the indices corresponding
to non-zero components of the symmetrized basis functions.

- Input:
  + top_dir: Path to identify data_files that contains symmetrized basis
  functions.
  + pt_grp: Point group symmetry of the underlying crystal.
  + Nmax: Order of the basis functions (a,b) s.t. max(a+b) \leq Nmax.

- Output (stored in .mat file):
  + tot_Uprops: Indices corresponding to the non-zero components of the
  basis-vector. This array contains all the non-zero components of all 
  the symmetrized basis-vector. They are vertically stacked. The 
  index-range for each basis-vector is provided in ind_ranges variable.
  The columns are 1(#),3(a),4(b),5(gamma),6(alpha),7(beta) of tot_inds.

  + ind_ranges: The number of rows are the number of symmetrized basis
  functions. Each row corresponds to the "start" and "stop" indices in 
  the tot_Uprops array containing the non-zero components.
  For example, if the third-row in "ind_ranges" is [296, 544]. Then the
  rows 296:544 contain the (a,b,gamma,alpha,beta) indices of non-zero
  components for the "third" basis vector.


### `save_symmvec_MabInds_cryst`

Indices of basis functions symmetrized with just the crystal point-group symmetries.

### `save_symmvec_MabInds_cryst_ges`

Indices of basis functions symmetrized with just the crystal point-group symmetries and grain-exchange symmetry.

### `save_symmvec_MabInds_gbnull_const`

Indices of basis functions symmetrized with just the crystal point-group, grain-exchange symmetries, and a **Constant** constraint for zero misorientation.

### `save_symmvec_MabInds_gbnull_zero`

Indices of basis functions symmetrized with just the crystal point-group, grain-exchange symmetries, and a **Zero** constraint for zero misorientation.