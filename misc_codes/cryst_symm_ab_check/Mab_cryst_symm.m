%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Code to check the rotation of SO(4) Functions R^(a,b)
%%%%%%%

clear all; clc;

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

s1 = set_vars();
Nmax = s1.Nmax; 
% Nmax = 1;
pt_grp = s1.pt_grp;

data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/cryst_symm/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symm_orders = zeros(Nmax^2,2);
ct3 = 1;
for ct1=0:Nmax
    for ct2=0:Nmax
        symm_orders(ct3,:) = [ct1,ct2];
        ct3 = ct3 + 1;
    end
end
% symm_orders = [4,3];
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_cols = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
nsymm = size(symm_orders,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Symmetry generator
Zpi = vrrotvec2mat([0,0,1,pi]); gs1 = Zpi; gs2 = eye(3);
csymm_mat = generate_csymm(gs1, gs2, symm_orders);

% mat_name = ['CS_Nmax_',num2str(Nmax),'.mat']; 
% save(mat_name,'csymm_mat');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath(genpath(util_dir));