# Function Descriptions

## `sp_orth`


Finds the **left** and **right** null spaces of the sparse matrix `A`.

- Inputs:
  + `A`  : a sparse matrix with `m` rows and `n` columns.
  + `opt`: `0` or `1` specifies that the left or right null 
            space is desired.
  + `TOL`: threshold for a pivot to be considered significant

- Outputs:
  + `N`: the left or right null space of the matrix A depending on the
        value of opt. The sparsity is increased at the expense of the 
        columns not being orthonormal.

- Notes:
  + The algorithm is not well-established in the literature, 
      but performs the decomposition `A = L U Q`, 
	- where `L` and `Q` are invertible matrices, 
	- `U` and `A` are the same size, and 
	- `U` has the block structure `U = [U_{11}, 0; 0, 0]`. 
  + If `U` is rank `r`, then since `A Q^{-1} = L U` the `n - r` columns 
      of `Q^{-1}` corresponding to pivotless columns of `U` are a basis 
      for the right null space. 
  + Similarly, since `L^{-1} A = U Q` the `m - r` rows of `L^{-1}` 
      corresponding to pivotless rows of `U` are a basis for the 
      left null space.


## `sp_null`

Finds an orthogonal basis for the column space of `A`, but is 
actually just a wrapper for `qr`.

- Inputs:
  + `A` - matrix specifying the column space.

- Outputs:
  `Q` - matrix whose columns form an orthonormal column space for `A`.


## `clean`

clean attempts to construct an orthonormal column space for `A` with the 
minimum number of nonzero entries. Assumes that `A` has full column rank.

- Inputs:
  + `A`  : matrix specifying the column space.
  + `TOL`: threshold below which an entry is considered insignificant.

- Outputs:
  + `A`: matrix whose columns form an orthonormal column space for `A`.

## `clebsch_gordan`

Computes all of the relevant Clebsch-Gordan coefficients for 
the specified `j1`, `j2`, `j` and `m`.

- Inputs:
  + `j1`, `j2`, `j`, `m`: indices

- Outputs:
  + `C`  : Clebsch-Gordan coefficients
  + `m1` : indices of the corresponding m1
  + `m2` : indices of the corresponding m2

- Notes:
  + Follows the approach of W. Straub in viXra:1403.0263.
  + Output is ordered in increasing values of m1.

