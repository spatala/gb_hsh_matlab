function [] = generate_gb_const()

top_dir = get_top_dir();
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(util_dir));

pt_grp = 'Oh'; Nmax = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'aPLUSb_max_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_aPLUSb_max_',...
    num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1 = symm_orders(:,1); b1 = symm_orders(:,2);
na1 = 2*a1+1;nb1 = 2*b1+1; ng1 = 2*min(a1,b1)+1; inds1 = na1.*nb1.*ng1;
nsymm = numel(inds1);
row_inds = zeros(nsymm,2);
col_inds = zeros(nsymm,2);
for ct1 = 1:nsymm
    ta1 = symm_orders(ct1,1); tb1 = symm_orders(ct1,2);
    row_inds(ct1,:) = [(ta1-tb1)^2+1,(ta1+tb1+1)^2];
    if ct1 == 1
        col_inds(ct1,:) = [1,sum(inds1(1:ct1))];
    else
        col_inds(ct1,:) = [sum(inds1(1:ct1-1))+1,sum(inds1(1:ct1))];
    end
end

null_mat = sparse(max(max(row_inds)),sum(inds1));
e_eps_inds = zeros(max(max(row_inds)),nsymm*2);
for ct1 = 1:nsymm
    % ct1
    ta1 = symm_orders(ct1,1); tb1 = symm_orders(ct1,2);
    [mat_ab, e_eps_vals] = null_mat_ab(ta1,tb1);
    r1 = row_inds(ct1,1); r2 = row_inds(ct1,2);
    c1 = col_inds(ct1,1); c2 = col_inds(ct1,2);
    null_mat(r1:r2,c1:c2) = mat_ab;
    e_eps_inds(r1:r2,2*ct1-1:2*ct1) = e_eps_vals;
end

rmpath(genpath(util_dir));

end
function [mat_ab, e_eps_vals] = null_mat_ab(a,b)
c = min(a,b);
na = 2*a+1; 
nb = 2*b+1;
ng = 2*c+1;

e_range = abs(a-b):a+b; 
num_e = numel(e_range);
mat_ab = sparse(sum(2*e_range+1),na*nb*ng);

r_ct = 1;
e_eps_vals = zeros(sum(2*e_range+1),2);
for e_ct = 1:num_e
    e = e_range(e_ct);
    [Cz, ~, ~] = clebsch_gordan(a, b, e, 0);
    eps_range = -e:e;
    ne = 2 * e + 1;
    for eps_ct = 1:ne
        eps1 = eps_range(eps_ct);
        [Ce, m1, m2] = clebsch_gordan(a, b, e, -eps1);
        ind = (m1 + a) * nb + m2 + b + 1;
        C_val = Ce*(Cz.');
        C_val = C_val(:);
        ind1 = ind + (0:ng-1)*na*nb;
        ind1 = ind1(:);
        mat_ab(r_ct, ind1) = C_val/(sqrt(2*e+1));
        e_eps_vals(r_ct,:) = [e, eps1];
        r_ct = r_ct + 1;
    end
end
mat_ab = mat_ab*sqrt(2)*sqrt((2*a+1)*(2*b+1))/pi;
end
