function Mrot = calc_rot_mat_Rab(Rab_rot, a,b)
% function [Mrot, inds, Rab_inds] = calc_rot_mat_Rab(Rab, a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Get the value of MBP functions for a given $(a,b)$ and R^(ab)
%%%%%%% matrix
%%%%%%%
%%%%%%% Input
%%%%%%% Rab_rot: (2*a+1)*(2*b+1) square matrix
%%%%%%%         Is given by kron(Ua(g1), Ub(g2)
%%%%%%%
%%%%%%% a, b: Integers
%%%%%%%         order of MBP function (gamma_MBP^{a,b}_{alpha, beta})
%%%%%%%
%%%%%%%
%%%%%%% Output
%%%%%%% Mrot: (2*a+1)*(2*b+1)*(2c+1) square matrix.
%%%%%%%         Rotation matrix for MBP^{a,b}.
%%%%%%%         This is well-defined only for right multiplications.
%%%%%%%         Left multiplication works for a few cases. (Ypi, Zth etc.)
%%%%%%%

c = min(a, b); Na = 2*a; Nb = 2*b; Nc = 2*min(a,b);
tf1 = (Na+1)*(Nb+1);
R_col_inds = zeros(tf1,Nc+1);
for ct1=1:Nc+1
    g1 = ct1-(c+1);
    ai = g1 + a + 1; bi = -g1 + b + 1;
    ind_row = (ai-1)*(2*b+1)+bi;
    R_col_inds(:,ct1) = (ind_row-1)*(tf1)+1 : ind_row*tf1;
end
tcols = R_col_inds(:); Mrot = Rab_rot(tcols, tcols);
end