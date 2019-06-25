function M = calc_Mfunc(a,b, rots)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Get the value of MBP functions for a given $(a,b)$.
%%%%%%%
%%%%%%% Input
%%%%%%% a, b: Integers
%%%%%%%         order of MBP function (gamma_MBP^{a,b}_{alpha, beta})
%%%%%%%
%%%%%%% rots: 3 X 6 matrix
%%%%%%%         First three columns are the first orientation matrix (g1)
%%%%%%%         and the second three are the second orientaiton matrix (g2)
%%%%%%%
%%%%%%% Output
%%%%%%% M: 1 X N row vector.
%%%%%%%         The value of the function MBP function 
%%%%%%%         (gamma_MBP^{a,b}_{alpha, beta})
%%%%%%%         for different values of a, b, gamma, alpha, beta
%%%%%%%         are calculated for the GB parameter (given by rot)
%%%%%%%
r1=rots(:,1:3);r2=rots(:,4:6);

c = min(a, b);
ax_ang1 = vrrotmat2vec(r1');
ax_ang2 = vrrotmat2vec(r2');

inds_mat = mbp_ab(a,b);
num_inds = size(inds_mat,1);


Na = 2*a; Nb = 2*b;
Ua_r1 = rotation(ax_ang1(1:3), ax_ang1(4), Na);
Ub_r2 = rotation(ax_ang2(1:3), ax_ang2(4), Nb);

M = zeros(1, num_inds);
for ct1=1:num_inds
    g1 = inds_mat(ct1, 4);
    a1 = inds_mat(ct1, 5);
    b1 = inds_mat(ct1, 6);
    g1_inda =  g1 + c + 1;
    g1_indb = -g1 + c + 1;
    a1_ind = a1 + a + 1;
    b1_ind = b1 + b + 1;
    M(ct1) = Ua_r1(a1_ind, g1_inda)*Ub_r2(b1_ind, g1_indb);
end

% if a == b
%     PI_ab = sqrt(2*a+1);
% else
%     if (a > b)
%         PI_ab = sqrt(prod(2*b+1:2*a+1));
%     else
%         PI_ab = sqrt(prod(2*a+1:2*b+1));
%     end
% end
% 
% M = PI_ab*M/(2*pi*pi);
end