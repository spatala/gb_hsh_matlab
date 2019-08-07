clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curr_pwd = split(pwd,'/');
top_dir = '';
for ct1=1:length(curr_pwd)
    top_dir = strcat(top_dir,curr_pwd{ct1},'/');
    if (strcmp(curr_pwd{ct1},'gb_hsh_matlab'))
        break;
    end
end
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(util_dir));
% addpath(genpath('Spherical-Harmonic-Transform/'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pt_grp = 'C1'; Nmax = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_inds = mbp_inds_ab_array(symm_orders);

a_val = symm_orders(:,1)'; b_val = symm_orders(:,2)';
c_val = min(a_val, b_val);
num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));

null_mat = sparse((2*Nmax+1)^2,num_cols);

for a=a_val
    for b=b_val
        c = min(a,b);
        na = 2*a+1; nb = 2*b+1; ng = 2*c+1;
        alpha = -a:a; beta = -b:b; gamma = -c:c;
        
        cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
        ind_start = find(cond1,1);
        ind_stop  = find(cond1,1,'last');
        
        for e = abs(a - b):(a + b)
            epsilon = -e:e;
            ne = 2 * e + 1;
            
            [C, m1, m2] = clebsch_gordan(a, b, e, 0);
            ind = (m1 + a) * nb + m2 + b + 1;
            Cz = sparse(ind, 1, C, na * nb, 1);
            for p = 1:ne % \epsilon
                M3 = zeros(na * nb, ng);
                
                [C, m1, m2] = clebsch_gordan(a, b, e, -epsilon(p));
                ind = (m1 + a) * nb + m2 + b + 1;
                Ce = sparse(ind, 1, C, na * nb, 1);
                for s = 1:na % \alpha
                    for t = 1:nb % \beta
                        row = (alpha(s) + a) * nb + beta(t) + b + 1;
                        for u = 1:ng % \gamma
                            M3(row, u) = M3(row, u) + ...
                                (Ce(row)/ sqrt(2 * e + 1)) * ...
                                Cz((gamma(u) + a) * nb - gamma(u) + b + 1);
                        end
                    end
                end
                M3 = M3 * sqrt(2) * sqrt((2 * a + 1) * (2 * b + 1)) / pi;
                
                row1 = e^2 + e + epsilon(p) + 1;
                null_mat(row1,ind_start:ind_stop) = M3(:)';
            end
        end
    end
end
S = spnull(null_mat);

mat_name = [data_fname0, ...
    'Sarr_gbnull_nmax_',num2str(Nmax),'.mat'];
save(mat_name,'S');

rmpath(genpath(util_dir));