clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Set the directory for Utilifity functions
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

%%%% Sets variables (Nmax and point-group)
s1 = set_vars();
Nmax = s1.Nmax; pt_grp = s1.pt_grp;
%%%% File path for storing data-files (.mat files)
fname = [top_dir,'data_files', '/ptgrp_',pt_grp];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Open the file symm_ab....m, which contains the (a,b) indices that are
%%%% satisfied by the crystal symmetries of the pt-grp
mat_name = [fname,'/cryst_symm/symm_ab_',...
pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1, b1);
num_rows = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The variable tot_inds contains the mapping between the column-index
%%%% for the M-function and the (a,b,gamma,alpha, beta) values.
tot_inds = mbp_inds_ab_array(symm_orders, num_rows);

%%%% ges_mat is the $A$ matrix that specifies the grain-exchange-symmetry
%%%% for the coefficients.
num_cols = num_rows;
ges_mat = sparse(num_rows, num_cols);
for ct1=1:num_cols
%%%% Get the (a,b,gamma,alpha, beta) corresponding to the column.
    a1 = tot_inds(ct1,3); b1 = tot_inds(ct1,4);
    gamma1 = tot_inds(ct1,5); 
    alpha1 = tot_inds(ct1,6); beta1 = tot_inds(ct1,7);
    
%%%% Find the index for (b,a,gamma,beta, alpha) row.
    a2 = b1; b2 = a1;
    gamma2 = gamma1;
    alpha2 = beta1; beta2 = alpha1;
    ind1 = find((tot_inds(:,3) == a2) & ...
        (tot_inds(:,4) == b2) & ...
        (tot_inds(:,5) == gamma2) & ...
        (tot_inds(:,6) == alpha2) & ...
        (tot_inds(:,7) == beta2));
    
% %%%% The (a,b,?,?,?)th column contains a 1 in the (b,a,?,?,?)th row
%     ges_mat(ind1,ct1) = 1;
% %%%% The (a,b,?,?,?)th column contains a (?1)^{a+b} in the 
% %%%% (a,b,?,?,?)th row
%     ges_mat(ct1,ct1) = (-1)^(a1+b1);

% % Set 1:
%     ges_mat(indct1) = 1;
%     ges_mat(ct1,ct1) = (-1)^(a1+b1);
% % Set 2:
%     ges_mat(ind1,ct1) = (-1)^(a1+b1);
%     ges_mat(ct1,ct1) = 1;
% Set 3:
    ges_mat(ct1,ind1) = 1;
    ges_mat(ct1,ct1) = (-1)^(a1+b1);
% % Set 4:
%     ges_mat(ct1,ind1) = (-1)^(a1+b1);
%     ges_mat(ct1,ct1) = 1;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Finding the column-space for a sparse matrix
%%%% Using QR decomposition and the stepform of R matrix. 
%%%% https://math.stackexchange.com/questions/748500/how-to-find-linearly-independent-columns-in-a-matrix/748538#748538
[Q,R] = qr(ges_mat);
j_inds = zeros(num_rows, 1); ct2 = 1; st1 = 0;
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
j_inds(ct2:end) =[]; col1 = ges_mat; col1(:,j_inds) = [];
%%%% Save the column-space in .mat file.
mat_name = [fname,'/ge_symm/Y_ges_Nmax_',...
    num2str(Nmax),'.mat']; save(mat_name, 'col1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check that the column-space calculation is indeed correct.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y1=col1;
%%%% Get an orthogonal basis for the column-space
Y1=orth(full(col1));
mat_name = [fname,'/ge_symm/orthY_ges_Nmax_',...
    num2str(Nmax),'.mat']; save(mat_name, 'Y1');
%%%% Get the eigenvectors of eigen-value one for the projection matrix
Q_mat = Y1*Y1'; [v, d] = eig(Q_mat);
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
S = orth(v(:,col));
if any(col) S = orth(v(:,col)); end
mat_name = [fname,'/ge_symm/Sarr_ges_Nmax_',...
    num2str(Nmax),'.mat'];
save(mat_name,'S');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmpath(genpath(util_dir));