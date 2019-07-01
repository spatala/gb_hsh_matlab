function [] = gen_symm_orders(top_dir, pt_grp, Nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 
%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ga_s, gb_s, num_gen] = get_symmgen_mats(pt_grp);
symm_orders = zeros(Nmax^2,2);
ct3 = 1;
for ct1=0:Nmax
    for ct2=0:Nmax
        a_val = ct1; b_val = ct2;
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
            symm_orders(ct3,:) = [ct1, ct2]; ct3 = ct3 + 1;
        end
    end
end

symm_orders(ct3:end,:) = [];
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
save(mat_name,'symm_orders');

end

function [v, d] = compute_eigen(ga_s,gb_s,n_gen, a_val,b_val)
gs1 = ga_s{n_gen}; gs2 = gb_s{n_gen};
Na = 2*a_val; Nb = 2*b_val; R1 = full(so4_irrep(gs1,gs2,Na,Nb));
[v,d] = eig(R1);
end

function [v, d] = combine_XY_symms(X0, v, col)
P0 = X0*X0';
Y1 = orth(v(:,col)); Q1 = Y1*Y1';
[v, d] = eig(P0*Q1*P0);
end