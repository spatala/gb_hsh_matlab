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
addpath(genpath(top_dir));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;

rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand())); r1 = rots1(:,1:3);


a_val = 5; Na = 2*a_val;
ax_ang = vrrotmat2vec(r1);
U_a = rotation(ax_ang(1:3), ax_ang(4), Na);

trots = rot_mats(:,:,floor(size(rot_mats,3)*rand())); 
tl = trots(:,1:3); tr = trots(:,4:6);
ax_ang_tl = vrrotmat2vec(tl);
Ul = rotation(ax_ang_tl(1:3), ax_ang_tl(4), Na);
ax_ang_tr = vrrotmat2vec(tr);
Ur = rotation(ax_ang_tr(1:3), ax_ang_tr(4), Na);

tmat = tl*r1*tr;
ax_ang_t = vrrotmat2vec(tmat);
tU = rotation(ax_ang_t(1:3), ax_ang_t(4), Na);
norm(Ul*U_a*Ur - tU)
