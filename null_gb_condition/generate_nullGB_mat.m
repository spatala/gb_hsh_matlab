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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pt_grp = 'C1'; Nmax = 1;

data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];

mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];

s1 = load(mat_name);
symm_orders = s1.symm_orders;
tot_inds = mbp_inds_ab_array(symm_orders);

nsz = size(tot_inds,1);

e_max = 2*Nmax;
e_eps_vals = zeros((e_max+1)^2,2);
ct3 = 1;
for ct1=0:e_max
    for ct2 = -ct1:ct1
        e_eps_vals(ct3,1) = ct1;
        e_eps_vals(ct3,2) = ct2;
        ct3 = ct3 + 1;
    end
end
gbnull_mat = sparse((e_max+1)^2,nsz);

for ct2 = 1:(e_max+1)^2
    e1 = e_eps_vals(ct2,1);
    eps1 = e_eps_vals(ct2,2);
    for ct1=1:nsz
        a1 = tot_inds(ct1,2);
        b1 = tot_inds(ct1,3);
        
        g1 = tot_inds(ct1,4);
        al_1 = tot_inds(ct1,5);
        be_1 = tot_inds(ct1,6);
        try
            cleb1 = clebschgordan(a1, al_1, b1, be_1, e1, -eps1);
            cleb2 = clebschgordan(a1, g1, b1, -g1, e1, 0);
            
            if a1 == b1
                PI_ab = sqrt(2*a1+1);
            else
                if (a1 > b1)
                    PI_ab = sqrt(prod(2*b1+1:2*a1+1));
                else
                    PI_ab = sqrt(prod(2*a1+1:2*b1+1));
                end
            end
            
            PI_e = sqrt(2*e1+1);
            
            gbnull_mat(ct2,ct1) = (sqrt(2)*PI_ab)/(pi*PI_e)*cleb1*cleb2;
        catch
            gbnull_mat(ct2,ct1) = 0;
        end
    end
end

S = spnull(gbnull_mat);

mat_name = [data_fname0, 'Sarr_gbnull_nmax_',num2str(Nmax),'.mat'];
save(mat_name,'S','gbnull_mat');

nsz1 = size(S,2);
rand_vec = 100*rand(nsz1,1);

C_vec = zeros(nsz,1);
for ct1 = 1:nsz1
    C_vec = C_vec + rand_vec(ct1)*S(:,ct1);
end

norm(gbnull_mat*C_vec)

rmpath(genpath(util_dir));
% tot_inds = 