# Function Descriptions

## `rots_to_mbp`

Converts GB parameters from `SO3xSO3` parameterization (rots) to 
Misorientation-Boundary plane (MBP).

- Input
	+ `rots`: 1 X 3 X 6 matrix
        - First three columns are the first orientation matrix (g1) and
        - the second three are the second orientaiton matrix (g2)

- Output
	+ `mbp`: 3 X 4 matrix
        - First three columns are the misorientation matrix and 
        - the last column is the boundary-plane.

## `rotmat_to_angs`

Convert rotation maxtrix to angles

- Input:
  + g: $3 \times 3$ rotation matrix.

- Output:
  + rot_angs: $1 \times 3$ rotation angles
	- w: Omega (rotation angle)
	- th: polar angle (theta), and 
	- ph: azimuthal angle (phi)

## `mbp_to_rots`

Convert GB parameters from Misorientation-Boundary plane (MBP) to
SO3xSO3 parametersation (rots).

- Input
  + mbp: 3 X 4 matrix
        First three columns are the misorientation matrix and the
        last column is the boundary-plane.

- Output
  + rots: 3 X 6 matrix
        First three columns are the first orientation matrix (g1)
        and the second three are the second orientaiton matrix (g2)



## `rots_to_angs`

Convert `SO(3) X SO(3)` parameterization into five angles

- Input
  + `g1`, `g2`: `3 X 3` rotation matrices.

- Output
  + mbp_angs: `1 X 5` matrix. The columns are as follows:
	- `omega_m`: rotation angle
	- `theta_m`: polar angle
	- `phi_m`  : azimuthal angle
	- `omega_b`: rotation angle
	- `phi_b`  : azimuthal angle

- Notes
	+ The first three angles give the misorientation matrix. The last two give the rotation matrix to bring the boundary-plane normal parallel to z-axis. The rotation is along a vector that is in the X-Y plane.


## `check_conv`


Function checks that the angles (w, th, ph) match with the rotation
matrix (r1)

- Input:
  + r1: Rotation matrix
  + rot_angs: array of rotation angles (w, th, ph)

- Output: None
  + Prints an error message if the check fails.

