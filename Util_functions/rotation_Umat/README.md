# Function Descriptions

## angles_to_CK


Computes `U^{1/2}` for the rotation given by a rotation angle (`w`) and axis (`th` and `ph`). 
Rows and cols ordered by increasing `m'` and `m`.
 
- Input
   + `w`, `th`, `ph`:
       The rotation angle, polar and azimuthal angles of the rotation
       matrix.

- Output
   + `CK`: Cayley-Klein Parameters

## CK_to_angles

Computes the rotation angle (`w`) and axis (`th`, `ph`) for the 
given `U^{1/2}` with rows and cols ordered by increasing `m'` and `m`.

- Input
   + `CK`: Cayley-Klein Parameters

- Output
   + `w`, `th`, `ph`:
       The rotation angle, polar and azimuthal angles of the rotation
       matrix.

## wigner_little_d

Computes the Wigner `d` matrix for the given order and rotation angle, with rows and colums labelled in increasing values of m.

- Input
   + `j`: order of the matrix
   + `theta`: rotation angle


- Output
   + `d`: Wigner (small) d-matrix
   + `err`: estimate of the magnitude of numerical error.

- Notes  
  + Follows the approach of X. M. Feng et al in 10.1103/PhysRevE.92.043307.

## rotation_wo_svd

Computes elements of SO(3) Irreducible representative 
`U^{a_val}_{alp_val, al_val}` for rotations given by (q,Q) quaterion.

Optimized for array of rotations as input (suppose N rotations).

- Input:
	+ `a`: `q-1i*Q(3)`
	+ `lb`: `-Q(2)-1i*Q(1)`
	+ `q`: q0 component of the rotation quaternion. Array of size N X 1. N is the number rotations.
	+ `Q`: [q1,q2,q3] components array of the rotation quaternion. Array of size N X 3.
	+ `a_val`: The order of the rotation matrix
	+ `alp_val`: Corresponds to row-index. 
	+ `al_val`:  Corresponds to column-index.

- Output:
	+ `U_val`: Value of U in the alp_val row and al_val column. 
		Size N X 1, where N is the number of rotations.

- Notes:
	+ `U^{a_val}_{alp_val, al_val}`: a_val is the order of the rotation matrix.
	+ The "-" sign is added because "rotation.m" had the convention of decreasing order of alp_val and al_val.


## rotation_mat

Computes the `(2j + 1)`-dimensional irreducible representation (irrep) of `SU(2)`, indicated by `U^{a}_{\alpha' \alpha}`

- Input Parameters
  + `j`: Umatrix of order j.
  + rot_angs: `3 X 1` array. `w` = rot_angs(1); `th` = rot_angs(2); `ph` = rot_angs(3);
	+ `w`: rotation angle;
    + `th`: polar angle of rotation axis
    + `ph`: azimuthal angle of rotation axis

- Output
  + `[U]`: rotation matrix `U` of the order `j`. Size: `(2*j+1, 2*j+1)`

- Notes:
  + Follows Eq. 6 on page 81 of D. A. Varshalovich et al, Quantum Theory of Angular Momentum, 1988.

