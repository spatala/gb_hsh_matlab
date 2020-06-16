function Mvec = compute_Mvec(s, nrots, num_cols, vec_inds)
% 
% Computes MBP basis functions for given set of (a, b, gamma, alpha, beta)
% for a given order (a+b \leq Nmax) of symmetrized MBP basis function.
% 
% - Input:
%   + s: struct of N(=nrots)-arrays (8 elements)
%     - ((s.a1, s.lb1, s.q1, s.Q1), (s.a2, s.lb2, s.q2, s.Q2))
%   + nrots: Number of GB parameters
%   + num_cols: Number of columns in the given order (a+b \leq Nmax) of 
%       symmetrized MBP basis function
%   + vec_inds: The non-zero indices in the symmetrized MBP basis 
%       functions.
% 
% - Output:
%   + Mvec: array of size nrots x num_cols
%     - The rows coresspond to different rotations.
%     - The columns correspond to all possible combinations of 
%       a, b, gamma, alpha, beta.
%     - Only the non-zero indices (stored in vec_inds) are calculated.
% 
% - Notes:
%   + The idea is to compute MBP functions only for those combinations of 
%       (a, b, gamma, alpha, beta) that are necessary for the 
%       symmetrized MBP basis functions.
%   + For any given order of symmetrized MBP basis functions 
%       (a+b \leq Nmax), we can determine all the non-zero indices 
%       (corresponding to unique a, b, gamma, alpha, beta). 
%       This is stored in the array vec_inds.
% 

a1_rots  = s.a1;
a2_rots  = s.a2;
lb1_rots = s.lb1;
lb2_rots = s.lb2;
q1_rots  = s.q1;
q2_rots  = s.q2;
Q1_rots  = s.Q1;
Q2_rots  = s.Q2;


a_val   = vec_inds(:,3);
alp_val = vec_inds(:,6);
al_val  = vec_inds(:,5);
b_val   =  vec_inds(:,4);
bep_val =  vec_inds(:,7);
be_val  = -vec_inds(:,5);


Mvec = sparse(nrots,num_cols);

for ct2 = 1:size(vec_inds,1)
    a = a_val(ct2); b = b_val(ct2);
    U1 = rotation_wo_svd(a1_rots,lb1_rots, ...
        q1_rots, Q1_rots, a, alp_val(ct2), al_val(ct2));
    U2 = rotation_wo_svd(a2_rots,lb2_rots, ...
        q2_rots, Q2_rots, b, bep_val(ct2), be_val(ct2));
    PI_ab = realsqrt((2*a+1)*(2*b+1));
    % M = PI_ab*M/(2*pi*pi);
    Mvec(:,vec_inds(ct2,1)) = PI_ab*(U1.*U2)/(sqrt(2*(pi^3)));
end
end