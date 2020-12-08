# Function Descriptions

## `mbp_inds_ab_array`

Get indices for MBP function with values 
of (a,b) given by symm_orders.

- Input
  + `symm_orders`: n X 2 array.
        The possible values of (a,b) are provided in the symm_orders
        array.

- Output
  + `mbp`: N X 7 matrix.
        The columns are as follows:
        1) Index number for fixed MBP function.
        2) Index number for fixed a, b.
        3) Value of `a`
        4) Value of `b`
        5) Value of `gamma` - range [-c, c], where $c = min(a,b)$
        6) Value of `alpha` - range [-a, a]
        7) Value of `beta`  - range [-b, b]


## `mbp_ab`

Get indices for MBP functions with fixed a and b.

- Input
  + a, b: Integers
        order of MBP function (`\tensor[_{\gamma}]{M}{^a_{\alpha}^b_{\beta}}`)

- Output
  + mbp: N X 6 matrix. The columns are as follows:
    - Index number for fixed a, b.
    - Value of 'a'
    - Value of 'b'
    - Value of 'gamma' - range [-c, c], where $c = min(a,b)$
    - Value of 'alpha' - range [-a, a]
    - Value of 'beta'  - range [-b, b]
