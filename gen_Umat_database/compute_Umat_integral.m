function [] = compute_Umat_integral()
top_dir = get_top_dir();

num = 150;
diff_val = s3_diff_elem(num);

N = 2;

du0 = '/data_files/Umat_grid/';
du1 = ['Umat_',num2str(N),'_disc_',num2str(num),'.mat'];

mat_name = [top_dir, du0, du1];

s1 = load(mat_name);
U_mat = s1.U_mat;
U_mat1 = zeros(num, num, 2*num);

diff_Atot = zeros((2*N+1)^2,1);

ct1 = 1;

for a1 = 1:2*N+1
    for a2 = 1:2*N+1
        disp(ct1)
        U_mat1(:,:,:) = U_mat(a1, a2, :, :, :);
        U_avg = compute_grid_avg(U_mat1);
        UU_star = U_avg.*conj(U_avg);
        diff_val1 = UU_star.*diff_val;
        Atot = sum(sum(sum(diff_val1)));
        diff_Atot(ct1) = (abs(Atot - 2*pi*pi/(2*N+1)));
        ct1 = ct1 + 1;
    end
end
disp(max(diff_Atot));

end