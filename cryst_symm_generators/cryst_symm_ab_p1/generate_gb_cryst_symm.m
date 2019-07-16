function [] = generate_gb_cryst_symm()
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

s1 = set_vars(); Nmax = s1.Nmax; pt_grp = s1.pt_grp;

% %%%%
% gen_symm_orders(top_dir, pt_grp, Nmax);
% %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ga_s, gb_s, ~, ~] = get_symmgen_mats(pt_grp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ns_ord = 3;
a_val = symm_orders(ns_ord,1); b_val = symm_orders(ns_ord,2);
c_val = min(a_val, b_val);
Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val;

ncryst_symm = 4;

for ct1=1:ncryst_symm
    S = generate_c1symm(ct1, ga_s, gb_s, Na, Nb, 'ab');
    mat_name = ['Sarr_ab_',num2str(a_val),'_',...
        num2str(b_val),'_',num2str(ct1),'.mat'];
    save(mat_name, 'S');

    S = generate_c1symm(ct1, ga_s, gb_s, Na, Nb, 'abc');
    mat_name = ['Sarr_abc_',num2str(a_val),'_',...
        num2str(b_val),'_',num2str(ct1),'.mat'];
    save(mat_name, 'S');
%     check_equi_basis(kron(eye(Nc+1),S1), S2)
end

X1 = cell(4,1); X2 = cell(4,1);

ct1 = 1;
mat_name = ['Sarr_ab_',num2str(a_val),'_',...
    num2str(b_val),'_',num2str(ct1),'.mat'];
s1 = load(mat_name); X1{ct1} = s1.S;
mat_name = ['Sarr_abc_',num2str(a_val),'_',...
    num2str(b_val),'_',num2str(ct1),'.mat'];
s1 = load(mat_name); X2{ct1} = s1.S;


for ct1 = 2:4
    % ct1 = 2; ct0 = ct1-1;
    ct0 = ct1-1;
    mat_name = ['Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
    X1{ct1} = get_Xarr_proj(X1{ct0}, mat_name);
    mat_name = ['Sarr_abc_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
    X2{ct1} = get_Xarr_proj(X2{ct0}, mat_name);
    tX1 = kron(eye(Nc+1), X1{ct1}); check_equi_basis(tX1,X2{ct1})
end


rmpath(genpath(util_dir));

end


function Sarr = generate_c1symm(ct1, ga_s, gb_s, Na, Nb, t1_typ)
gs1 = ga_s{ct1}; gs2 = gb_s{ct1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rr_ab_12 = so4_irrep(gs1,gs2,Na,Nb);
trRr_ab_12 = transpose(Rr_ab_12);

if strcmp(t1_typ, 'ab')
    [v,d] = eig(trRr_ab_12);
    col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
    if any(col)
        Sarr = orth(v(:,col));
    end
end

if strcmp(t1_typ, 'abc')
    Nc = min(Na,Nb);
    csymm_mat = kron(eye(Nc+1), trRr_ab_12);
    [v,d] = eig(csymm_mat);
    col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
    if any(col)
        Sarr = orth(v(:,col));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function Xarr = get_Xarr_proj(X1, mat_name)
P1 = X1*X1';
s1 = load(mat_name); Y1 = s1.S;
Q1 = Y1*Y1';
R1 = P1*Q1*P1;
[v,d] = eig(R1);
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
if any(col)
    Xarr = orth(v(:,col));
end
end
