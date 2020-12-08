# Function Descriptions

## `get_symmgen_angs`

Function to set the symmetry generators for each point group.
The symmetries are for grain boundaries in the `SO(3) \times SO(3)`
parameterization, i.e. `(O_a, O_b) ~ (O_a S_i, O_b S_j)`

- Input:
  + `pt_grp`:	(string) 
  			Point-group of the underlying crystal

- Output:
  + `ga_s`, `gb_s`: (Cell arrays)
   		grain boundary symmetry generators.
  + `num_gen`: (integer)
   		Number of generators for the point-group.
  + `Laue`: (string)
   		`0` for non-Laue and `1` for Laue group.

- Notes:
  + Only proper rotation elements are set.
  + The symmetry operation is provided using three angles
  	- `w` : rotation angle
  	- `th`: polar angle of the axis
  	- `ph`: azimuthal angle of the axis


## `get_symm_rots`

Computes all the symmetrically equivalent GB parameters.
That is, for a given (g1, g2) all the equivalent (`g1*S`, `g2*`S) are
computed.

- Inputs:
   + g1, g2: `3 \times 3` rotation matrices of the 
       `SO(3) \times SO(3)` parameterization
   + pt_grp: string
       Crystal point group to apply symmetries
   + data_fname: string
       Location of the symmetry elements (matrices)
   + opt:
       1: crystal symmetries only (if Laue (M,n) = (M,-n))
       2: crystal symmetries + Grain Exchange symmetry

- Output:
  + symm_rots: `3,6,(1,2,4)*nrot_csymm`
      - All the symmetrically equivalent (g1, g2) pairs. The total number
        of equivalent parameters depends on opt (1 or 2).
      - `opt == 1`, has only crystal point group symmetries.
      - `opt == 2`, has both crystal + GES.

