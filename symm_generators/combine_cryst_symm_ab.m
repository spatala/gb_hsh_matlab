function [] = combine_cryst_symm_ab(top_dir, pt_grp, Nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
data_fname1 = [data_fname,'Sarr_ab/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, num_gen, Laue] = get_symmgen_mats(pt_grp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Laue
    num_symm = 2*num_gen + 1;
else
    num_symm = 2*num_gen;
end

Sab_arr = cell(nsymm,1);

num_rows = 0; num_cols = 0;

for ns_ord = 1:nsymm
    ns_ord
    X1 = cell(num_symm,1);
    a_val = symm_orders(ns_ord,1); b_val = symm_orders(ns_ord,2);
    c_val = min(a_val, b_val); Nc = 2*c_val;
    
    ct1 = 1;
    mat_name = [data_fname1,...
        'Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
    s1 = load(mat_name); X1{ct1} = s1.S;
    if Laue
        for ct1 = 2:num_symm-1
            ct0 = ct1-1;
            mat_name = [data_fname1,...
                'Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
            X1{ct1} = get_Xarr_proj(X1{ct0}, mat_name);
        end
        X1{ct1} = kron(eye(Nc+1), X1{ct1});
        ct1 = ct1 + 1; ct0 = ct1-1;
        mat_name = [data_fname1,...
            'Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
        X1{ct1} = get_Xarr_proj(X1{ct0}, mat_name);
    else
        for ct1 = 2:num_symm
            ct0 = ct1-1;
            mat_name = [data_fname1,...
                'Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
            X1{ct1} = get_Xarr_proj(X1{ct0}, mat_name);
        end
        X1{ct1} = kron(eye(Nc+1), X1{ct1});
    end
    
    Sab = X1{ct1};
    
    Na = 2*a_val; Nb = 2*b_val;
    c_val = min(a_val, b_val); Nc = 2*c_val;
    num_rows = num_rows + (Na+1)*(Nb+1)*(Nc+1);
    num_cols = num_cols + size(Sab,2);
    
    Sab_arr{ns_ord} = Sab;
end

S1_arr = sparse(zeros(num_rows, num_cols));
row_ind_start = 1;
col_ind_start = 1;
for ns_ord = 1:nsymm    
    Sab = Sab_arr{ns_ord};
    row_ind_stop = row_ind_start + size(Sab,1) - 1;
    col_ind_stop = col_ind_start + size(Sab,2) - 1;    
    S1_arr(row_ind_start:row_ind_stop, col_ind_start:col_ind_stop) = Sab;
    row_ind_start = row_ind_stop+1;
    col_ind_start = col_ind_stop+1;
end
S = clean(S1_arr);
mat_name = [data_fname0,'Sarr_abc_combined_csymm_nmax_',num2str(Nmax),'.mat'];
save(mat_name, 'S');
end


function Xarr = get_Xarr_proj(X1, mat_name)
P1 = X1*X1';
s1 = load(mat_name); Y1 = s1.S;
Q1 = Y1*Y1';
R1 = P1*Q1*P1;

nsz = size(R1,1);
R2 = R1 - speye(nsz,nsz);
Xarr = spnull(R2);
end
