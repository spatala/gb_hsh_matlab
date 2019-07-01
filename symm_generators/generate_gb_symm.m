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

%%%% 
gen_symm_orders(top_dir, pt_grp, Nmax);
%%%% 
generate_csymm_nmax(top_dir, pt_grp, Nmax);
%%%%
generate_ges(top_dir, pt_grp, Nmax)
%%%% 
generate_symm_combos(top_dir, pt_grp, Nmax);
%%%% 
rmpath(genpath(util_dir));

