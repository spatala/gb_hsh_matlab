%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate SO(3) irreducible representatives for a large number of
clear all; clc;

top_dir = get_top_dir();
util_dir = [top_dir, 'Util_functions/'];
du1 = [top_dir,'data_files/Umat_grid/'];
addpath(genpath(util_dir));

% num = 10;
% num = 25;
% num = 150;


wm1 = linspace(0,pi,num);
th1 = linspace(0,pi,num);
ph1 = linspace(0,2*pi,2*num);

[wm2, th2, ph2] = meshgrid(wm1, th1, ph1);
wm_sz = size(wm2);



a = 1:2;
num_a = size(a,2);
n_angs = size(wm2,1);

for ct1 =1:num_a
    % av = a(num_a+1-ct1);
    av = a(ct1);
    av
    U_mat = zeros(2*av+1, 2*av+1, wm_sz(1),wm_sz(2),wm_sz(3));
    for i=1:wm_sz(1)
        disp(i)
        for j=1:wm_sz(2)
            for k=1:wm_sz(3)
                w_1 = wm2(i,j,k);
                th_1 = th2(i,j,k);
                ph_1 = ph2(i,j,k);
                U_mat(:,:,i,j,k) = rotation_mat(av, [2*w_1, th_1, ph_1]);
            end
        end
    end
        mat_name = [du1,'Umat_',num2str(av),'_disc_',...
            num2str(num),'.mat'];
    save(mat_name,'U_mat','-v7.3');
end

rmpath(genpath(util_dir));