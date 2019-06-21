clear all; clc;

fname = get_dir_name();

Nmax = 8;
pt_grp = 'O';
% pt_grp = 'C2';
mat_name = [fname,'/ptgrp_',pt_grp,'/cryst_symm/symm_ab_',...
    pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1, b1);
num_rows = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
% Nmax = max(a1);

tot_inds = mbp_inds_ab_array(symm_orders, num_rows);


ges_mat = sparse(num_rows, num_rows);
% ges_mat = zeros(num_rows, num_rows);

for ct1=1:num_rows
%     ct1
    a1 = tot_inds(ct1,3);
    b1 = tot_inds(ct1,4);
    gamma1 = tot_inds(ct1,5);
    alpha1 = tot_inds(ct1,6);
    beta1 = tot_inds(ct1,7);
    
    a2 = b1;
    b2 = a1;
    gamma2 = gamma1;
    alpha2 = beta1;
    beta2 = alpha1;
    
    ind1 = find((tot_inds(:,3) == a2) & ...
        (tot_inds(:,4) == b2) & ...
        (tot_inds(:,5) == gamma2) & ...
        (tot_inds(:,6) == alpha2) & ...
        (tot_inds(:,7) == beta2));
    
    ges_mat(ct1,ind1) = 1;
    ges_mat(ct1,ct1) = (-1)^(a1+b1);
    
end

[Q,R] = qr(ges_mat);

st1 = 0;
j_inds = zeros(num_rows, 1);
ct2 = 1;
for ct1=1:num_rows
    ct1
    ind1 = find(abs(R(:,ct1)), 1, 'last');
    if ct1 == 1
        st1 = ind1;
    else
        if (ind1 > st1)
            st1 = ind1;
        else
            j_inds(ct2) = ct1;
            ct2 = ct2 + 1;
        end
    end
end
j_inds(ct2:end) =[];
col1 = ges_mat;
col1(:,j_inds) = [];

mat_name = [fname,'/ptgrp_',pt_grp,'/ge_symm/Y_ges_Nmax_',...
    num2str(Nmax),'.mat'];
save(mat_name, 'col1');
% save('Y_ges.mat','col1');

Q_mat = col1*col1';
[v, d] = eig(full(Q_mat));
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
S = orth(v(:,col));
if any(col)
    S = orth(v(:,col));
end
mat_name = [fname,'/ptgrp_',pt_grp,'/ge_symm/Sarr_ges_Nmax_',...
    num2str(Nmax),'.mat'];
save(mat_name,'S');