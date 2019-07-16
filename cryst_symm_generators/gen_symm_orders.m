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
    ct1
    for ct2=0:Nmax
        a_val = ct1; b_val = ct2;
        for ct4 = 1:2*num_gen
            S0 = compute_eigen(ga_s,gb_s,ct4,a_val,b_val);
            if (size(S0,2) > 0)
                if (ct4 == 1)
                    X0 = eye(size(S0,1), size(S0,1));
                else
                    X0 = S1;
                end
                S1 = combine_XY_symms(X0, S0);
                if (size(S1,2) == 0)
                    break;
                end
            end
        end
        if (size(S1,2) > 0)
            symm_orders(ct3,:) = [ct1, ct2]; ct3 = ct3 + 1;
        end
    end
end

symm_orders(ct3:end,:) = [];
data_fname1 = [data_fname0,'nmax_',num2str(Nmax),'/'];
mat_name = [data_fname1,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
save(mat_name,'symm_orders');

end

function S = compute_eigen(ga_s,gb_s,n_gen, a_val,b_val)
gs1 = ga_s{n_gen}; gs2 = gb_s{n_gen};
Na = 2*a_val; Nb = 2*b_val; R1 = so4_irrep(gs1,gs2,Na,Nb);
nsz = size(R1,1); R2 = R1 - speye(nsz,nsz); S = spnull(R2);
end

function S = combine_XY_symms(X0, Y1)
P0 = X0*X0'; Q1 = Y1*Y1'; R1 = P0*Q1*P0;
nsz = size(R1); R2 = R1 - speye(nsz); S = spnull(R2);
end