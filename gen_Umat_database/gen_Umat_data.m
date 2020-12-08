%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate SO(3) irreducible representatives for a large number of 
clear all; clc;

top_dir = get_top_dir();
util_dir = [top_dir, 'Util_functions/'];
addpath(genpath(util_dir));

num = 100;
wm1 = linspace(0,pi,num);
th1 = linspace(0,pi,num);
ph1 = linspace(0,2*pi,2*num);

[wm2, th2, ph2] = meshgrid(wm1, th1, ph1);

wm = wm2(:); th = th2(:); ph = ph2(:);

du1 = [top_dir,'data_files/Umat_grid/Umat_disc_',num2str(num),'/'];

a = 0:16;
num_a = size(a,2);
n_angs = size(wm,1);

for ct1 =1:num_a
    % av = a(num_a+1-ct1);
    av = a(ct1);
    av
    U_mat = zeros(2*av+1, 2*av+1, n_angs);
    for ct2=1:n_angs
        w_1 = wm(ct2);
        th_1 = th(ct2);
        ph_1 = ph(ct2);
        U_mat(:,:,ct2) = rotation_mat(av, [-w_1, th_1, ph_1]);
    end
    mat_name = [du1,'Umat_',num2str(av),'_disc_',...
        num2str(num),'.mat'];
    save(mat_name,'U_mat','-v7.3');
end

rmpath(genpath(util_dir));