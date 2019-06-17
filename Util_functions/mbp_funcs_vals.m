function M_vals = mbp_funcs_vals(r1, r2, N)
% a = 1; b = 1;
% tot_inds = mbp_inds(N);
tot_num = 0;
for ct1=0:N
    for ct2=0:N
        ct3 = min(ct1, ct2);
        tot_num = tot_num + (2*ct1+1)*(2*ct2+1)*(2*ct3+1);
    end
end
tot_inds = zeros(tot_num,7);
M_vals = zeros(tot_num,1);
ind_start = 1;
for ct1=0:N
    for ct2=0:N
        M = calc_Mfunc(ct1,ct2,r1,r2);
        inds_mat = mbp_ab(ct1,ct2);
        ind_stop = ind_start+size(inds_mat,1)-1;
        tot_inds(ind_start:ind_stop,1) = ind_start:ind_stop;
        tot_inds(ind_start:ind_stop,2:7) = inds_mat;
        M_vals(ind_start:ind_stop,1) = M';
        ind_start = ind_stop+1;
    end
end
end

function M = calc_Mfunc(a,b,r1,r2)
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

if a == b
    PI_ab = sqrt(2*a+1);
else
    if (a > b)
        PI_ab = sqrt(prod(2*b+1:2*a+1));
    else
        PI_ab = sqrt(prod(2*a+1:2*b+1));
    end
end

M = PI_ab*M/(2*pi*pi);
end

% M = zeros(10, num_inds);
%
% for ct2 = 1:10
%     th = 2*pi*rand();
%     gz = vrrotvec2mat([0,0,1,th]);
%     r1 = gz*r1; ax_ang1 = vrrotmat2vec(r1');
%     r2 = gz*r2; ax_ang2 = vrrotmat2vec(r2');
%
%     Na = 2*a_val;
%     Ua_r1 = rotation(ax_ang1(1:3), ax_ang1(4), Na);
%
%     Nb = 2*b_val;
%     Ub_r2 = rotation(ax_ang2(1:3), ax_ang2(4), Nb);
%
%
%
%
%
%     for ct1=1:num_inds
%         g1 = inds_mat(ct1, 4);
%         a1 = inds_mat(ct1, 5);
%         b1 = inds_mat(ct1, 6);
%         g1_inda =  g1 + c_val + 1;
%         g1_indb = -g1 + c_val + 1;
%         a1_ind = a1 + a_val + 1;
%         b1_ind = b1 + b_val + 1;
%         M(ct2, ct1) = Ua_r1(a1_ind, g1_inda)*Ub_r2(b1_ind, g1_indb);
%     end
% end


% U_r1_mats = cell(2,1); U_r2_mats = cell(2,1);
%
% N = 2*0;
% U_r1_mats{1} = rotation(ax_ang1(1:3), ax_ang1(4), N);
% U_r2_mats{1} = rotation(ax_ang2(1:3), ax_ang2(4), N);
%
%
% N = 2*1;
% U_r1_mats{2} = rotation(ax_ang1(1:3), ax_ang1(4), N);
% U_r2_mats{2} = rotation(ax_ang2(1:3), ax_ang2(4), N);




%%%%% To check rotation g_z
% for ct1=1:100
% th = 2*pi*rand();
% % th = 0;
%
% aval = 1; bval = 1;
% rots = rot_mats(:,:,1);
% r1 = rots(:,1:3);
% r1 = vrrotvec2mat([0,0,1,th])*r1;
% ax_ang1 = vrrotmat2vec(r1');
% U1 = rotation(ax_ang1(1:3), ax_ang1(4), 2*aval);
% r2 = rots(:,4:6);
% r2 = vrrotvec2mat([0,0,1,th])*r2;
% ax_ang2 = vrrotmat2vec(r2');
% U2 = rotation(ax_ang2(1:3), ax_ang2(4), 2*aval);
%
% U1(1,3)*U1(2,1)
% end