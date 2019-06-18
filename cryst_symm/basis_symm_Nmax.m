clear all; clc;

addpath(genpath('../Util_functions/'));

N = [0,4];

num_cols = 0;
for ct1=N
    for ct2=N
        a_val = ct1; b_val = ct2;
        c_val = min(a_val, b_val);
        
        mat_name = ['Sarr_',num2str(a_val),'_',num2str(b_val),'.mat'];        
        s1 = load(mat_name);
        Svec = s1.S;
        num_cols = num_cols + size(Svec,2);
    end
end
tot_inds = mbp_inds(N);
num_rows = size(tot_inds,1);
Sarr_nmax = zeros(num_rows, num_cols);

start_row_ind = 1;
col_ind = 1;
for ct1=N
    for ct2=N
        a_val = ct1; b_val = ct2;
        c_val = min(a_val, b_val);
        
        stop_row_ind = start_row_ind + (2*a_val+1)*(2*b_val+1)*(2*c_val+1)-1;

        mat_name = ['Sarr_',num2str(a_val),'_',num2str(b_val),'.mat'];
        s1 = load(mat_name);
        Svec = s1.S;

        Sarr_nmax(start_row_ind:stop_row_ind, col_ind) = Svec;
        start_row_ind = stop_row_ind + 1;
        col_ind = col_ind + 1;
        
    end
end
save('Sarr_nmax_4.mat','Sarr_nmax');

rmpath(genpath('../Util_functions/'));