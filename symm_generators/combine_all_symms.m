function [] = combine_all_symms(top_dir, pt_grp, Nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 
%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, num_gen, Laue] = get_symmgen_mats(pt_grp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Laue
    num_symm = 2*num_gen+2;
else
    num_symm = 2*num_gen+1;
end

nmax = max(symm_orders(:));

for ct4 = 1:num_symm
    ct4
    [v0, d0] = compute_eigen(ct4, nmax, data_fname0);
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
    save_symm_arr(Nmax, v1, col1, data_fname0);
end

end

function [v, d] = combine_XY_symms(X0, v, col)
P0 = X0*X0';
Y1 = orth(v(:,col)); Q1 = Y1*Y1';
[v, d] = eig(P0*Q1*P0);
end

function [v, d] = compute_eigen(ngen, nmax, fname)
mat_name = [fname,'symm_mat_full_ngen_',num2str(ngen)...
    ,'_nmax_',num2str(nmax),'.mat'];
s1 = load(mat_name); R1 = full(s1.symm_mat); [v,d] = eig(R1);
end


function save_symm_arr(Nmax, v, col, fname)
S = orth(v(:,col));
mat_name = [fname,'Sarr_combined_nmax_',num2str(Nmax),'.mat'];
save(mat_name,'S');
end




