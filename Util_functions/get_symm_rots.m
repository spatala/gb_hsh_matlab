function symm_rots = get_symm_rots(g1,g2, pt_grp, data_fname)
% clear all; clc;
% 
% curr_pwd = split(pwd,'/');
% top_dir = '';
% for ct1=1:length(curr_pwd)
%     top_dir = strcat(top_dir,curr_pwd{ct1},'/');
%     if (strcmp(curr_pwd{ct1},'gb_hsh_matlab'))
%         break;
%     end
% end
% util_dir = strcat(top_dir,'Util_functions','/');
% addpath(genpath(util_dir));
% 
% s1 = set_vars();
% Nmax = s1.Nmax; pt_grp = s1.pt_grp;
% 
% data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
% data_fname0 = [data_fname,'ge_symm/'];
% data_fname1 = [data_fname,'cryst_symm/'];


mat_name = [data_fname, 'SymmMat_', pt_grp, '.mat'];
s1 = load(mat_name);
SymmMat = s1.SymmMat;
nsz = size(SymmMat,1);

symm_rots = zeros(3,6,2*nsz*nsz);
ct3 = 1;
for ct1=1:nsz
    gs1 = g1*SymmMat{ct1};
    for ct2=1:nsz
        gs2 = g2*SymmMat{ct1};
        symm_rots(:,:,ct3) = [gs1,gs2];
        ct3 = ct3 + 1;
    end
end

ypi = vrrotvec2mat([0,1,0,pi]);
g1p = ypi*g2; g2p = ypi*g1;

for ct1=1:nsz
    gs1 = g1p*SymmMat{ct1};
    for ct2=1:nsz
        gs2 = g2p*SymmMat{ct1};
        symm_rots(:,:,ct3) = [gs1,gs2];
        ct3 = ct3 + 1;
    end
end

end