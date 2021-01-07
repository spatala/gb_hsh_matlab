## `get_top_dir.m`

Provides the directory path for `gb_hsh_matlab`.

## `check_symms`

This code will call the `symm_checks` and `symm_checks_nullgb_constraint` for various symmetries, for a given `pt_grp`, `Nmax` and `coeffs_typ`

## `symm_checks`

Steps:
   1) Load a random grain boundary (g1, g2).
   2) Create all the symmetrically equivalent grain boundaries
   3) Compute MBP basis functions for all the symmetrically equivalent
      grain boundaries.
   4) Check that all the basis functions return the same vectors.

Input:
   1) pt_grp
   2) Nmax
   3) coeffs_typ: string
       'nmax' or 'aPLUSb_max'
   4) symm_typ: string
       'cryst', 'cryst_ges', 'zero', 'const'


## `symm_checks_nullgb_constraint`

1) Computes the full MBP basis function for a random GB parameter with
identity misorientation
2) Computes the norm of the "constrained" basis-function
3) If constraint is:
  +  'zero': Checks that the norm is equal to zero, i.e. f(g,g) = 0
  +  'const': Checks that the basis functions have the same value, i.e.
  f(g,g) = C

