function [] = generate_ab_symms(top_dir, pt_grp, Nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function to compute symmetrized eigen-vectors for crystal symmetry
%%%% operations. The eigen-vectors are for a given (a,b) order.
%%%% 
%%%% Input:
%%%% top_dir:
%%%%    String. Directory location for saving the symmetrized the
%%%%    eigen-vectors.
%%%% pt_grp:
%%%%    The underlying point-group of the crystal.
%%%% Nmax:
%%%%    The maximum order for (a,b) for generating symmetrized
%%%%    eigen-vectors.
%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname1 = [data_fname0,'Sarr_ab/'];
data_fname2 = [data_fname0,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ga_s, gb_s, num_gen, Laue] = get_symmgen_angs(pt_grp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname2,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncryst_symm = 2*num_gen;
nosymm_inds = zeros(nsymm,1);
nosymm_ct = 1;

for ns_ord = 1:nsymm
    ns_ord
    
    a_val = symm_orders(ns_ord,1); b_val = symm_orders(ns_ord,2);
    Na = 2*a_val; Nb = 2*b_val;
    
    for ct1=1:ncryst_symm
        mat_name = [data_fname1,...
            'Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
        if ~(isfile(mat_name))
            symm_mat = generate_c1symm(ct1, ga_s, gb_s, Na, Nb);
            save_mat = save_Sarr(symm_mat, mat_name);
            if ~(save_mat)
                nosymm_inds(nosymm_ct) = ns_ord;
                nosymm_ct = nosymm_ct + 1;
            end
        end
    end
    
    if Laue
        ct1 = ct1 + 1;
        mat_name = [data_fname1,...
            'Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
        if ~(isfile(mat_name))
            symm_mat = generate_ypi_left_ab(a_val,b_val);
            save_mat = save_Sarr(symm_mat, mat_name);
            if ~(save_mat)
                nosymm_inds(nosymm_ct) = ns_ord;
                nosymm_ct = nosymm_ct + 1;
            end
        end
    end
end

if nosymm_ct > 1
    nosymm_inds(nosymm_ct:end) = [];
    symm_orders(nosymm_inds,:) = [];
end

mat_name = [data_fname2,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
save(mat_name,'symm_orders');
end


function symm_mat = generate_c1symm(ct1, ga_s, gb_s, a_val, b_val)
gs1 = ga_s{ct1}; gs2 = gb_s{ct1};

U_a = rotation_mat(a_val, gs1);
U_b = rotation_mat(b_val, gs2);

symm_mat = sparse((kron(U_a, U_b)).');

% Rr_ab_12 = so4_irrep(gs1,gs2,Na,Nb);
% symm_mat = transpose(Rr_ab_12);
end

function save_mat = save_Sarr(R1, mat_name)
nsz = size(R1,1);
R2 = R1 - speye(nsz,nsz);
S = spnull(R2);
if size(S,2) > 0
    save(mat_name, 'S');
    save_mat = 1;
else
    save_mat = 0;
end
end

