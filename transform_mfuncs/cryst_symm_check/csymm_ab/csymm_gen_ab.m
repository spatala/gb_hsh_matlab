function [] = csymm_gen_ab()
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

s1 = set_vars(); Nmax = s1.Nmax; pt_grp = s1.pt_grp;

data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/cryst_symm/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ga_s, gb_s, num_gen] = get_symmgen_mats(pt_grp);
symm_orders = zeros(Nmax^2,2);
ct3 = 1;

for ct1=0:Nmax
    for ct2=0:Nmax
        a_val = ct1; b_val = ct2;
        %         [ct1, ct2]
        for ct4 = 1:2*num_gen
            [v0, d0] = compute_eigen(ga_s,gb_s,ct4,a_val,b_val);
            col0 = (abs(imag(diag(d0)))<1e-5 & abs(real(diag(d0))-1)<1e-5);
            if any(col0)
                if (ct4 == 1)
                    X0 = eye(size(v0,1), size(v0,1));
                else
                    X0 = orth(v1(:,col1));
                end
                [v1, d1] = combine_XY_symms(X0, v0, col0);
                col1 = (abs(imag(diag(d1)))<1e-5 & abs(real(diag(d1))-1)<1e-5);
                if ~any(col1)
                    break;
                end
            end
        end
        if any(col1)
%             save_symm_arr(a_val, b_val, v1, col1, data_fname0);
            symm_orders(ct3,:) = [ct1, ct2]; ct3 = ct3 + 1;
        end
    end
end

symm_orders(ct3:end,:) = [];
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
% save(mat_name,'symm_orders');
rmpath(genpath(util_dir));
end

function [v, d] = compute_eigen(ga_s,gb_s,n_gen, a_val,b_val)
gs1 = ga_s{n_gen}; gs2 = gb_s{n_gen};
R1 = full(generate_csymm(gs1, gs2, [a_val,b_val]));
[v,d] = eig(R1);
end

function [v, d] = combine_XY_symms(X0, v, col)
P0 = X0*X0';
Y1 = orth(v(:,col)); Q1 = Y1*Y1';
[v, d] = eig(P0*Q1*P0);
end

function save_symm_arr(a_val, b_val, v, col, fname)
S = orth(v(:,col));
% size(S,2)
mat_name = [fname,'Sarr_',num2str(a_val),'_',num2str(b_val),'.mat'];
save(mat_name,'S');
end