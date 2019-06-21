clear all; clc;

curr_pwd = split(pwd,'/');
top_dir = '';
for ct1=1:length(curr_pwd)
    top_dir = strcat(top_dir,curr_pwd{ct1},'/');
    if (strcmp(curr_pwd{ct1},'MATLAB_Codes'))
        break;
    end
end
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(top_dir));

s1 = set_vars();
Nmax = s1.Nmax; pt_grp = s1.pt_grp;

data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/cryst_symm/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ct1=0:Nmax
    for ct2=0:Nmax
        a_val = ct1; Na = 2*a_val;
        b_val = ct2; Nb = 2*b_val;
        c_val = min(ct1, ct2); Nc = 2*c_val;
        [ct1, ct2]
        
        tot_inds = mbp_ab(a_val,b_val);  num_inds = size(tot_inds,1);
        ges_mat = sparse(num_inds, num_inds);
        for ct3=1:num_inds
            a1 = tot_inds(ct3,2); b1 = tot_inds(ct3,3); 
            gamma1 = tot_inds(ct3,4); alpha1 = tot_inds(ct3,5);
            beta1 = tot_inds(ct3,6);
            
            a2 = b1; b2 = a1; gamma2 = gamma1;
            alpha2 = beta1; beta2 = alpha1;
            
            ind1 = find((tot_inds(:,2) == a2) & ...
                (tot_inds(:,3) == b2) & ...
                (tot_inds(:,4) == gamma2) & ...
                (tot_inds(:,5) == alpha2) & ...
                (tot_inds(:,6) == beta2));
            
            ges_mat(ct3,ind1) = 1;
            ges_mat(ct3,ct3) = (-1)^(a1+b1);
        end
        
        [Q,R] = qr(ges_mat);
        st1 = 0; j_inds = zeros(num_inds, 1); ct4 = 1;
        for ct5=1:num_inds
            ind1 = find(abs(R(:,ct5)), 1, 'last');
            if ct5 == 1
                st1 = ind1;
            else
                if (ind1 > st1)
                    st1 = ind1;
                else
                    j_inds(ct4) = ct1;
                    ct4 = ct4 + 1;
                end
            end
        end
        j_inds(ct4:end) =[]; col1 = ges_mat; col1(:,j_inds) = [];
        
        mat_name = [fname,'/ptgrp_',pt_grp,...
            '/ge_symm/Y_ges_a_', num2str(a_val),...
            '_b_',num2str(b_val), num2str(Nmax),'.mat'];
        save(mat_name, 'col1');
        
        Q_mat = col1*col1';
        [v, d] = eig(full(Q_mat));
        col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);

        if any(col)
            S = orth(v(:,col));
        end
        mat_name = [fname,'/ptgrp_',pt_grp,...
            '/ge_symm/Sarr_ges_Nmax_a_', num2str(a_val),...
            '_b_',num2str(b_val), num2str(Nmax),'.mat'];
        save(mat_name,'S');
    end
end