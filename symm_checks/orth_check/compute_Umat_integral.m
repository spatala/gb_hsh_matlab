function [] = compute_Umat_integral()
top_dir = get_top_dir();

num = 25;
i1 = 1; i3 = 2;
i2 = num-1; i4 = num;
i5 = 2*num-1; i6 = 2*num;

ind0 = i1:i2;
ind1 = i3:i4;
ind2 = i1:i5;
ind3 = i3:i6;

diff_val = s3_diff_elem(num);

N = 2;

mat_name1 = [top_dir,'/data_files/',...
    'Umat_grid/Umat_disc_',num2str(num),...
    '/Umat_',num2str(N),...
    '_disc_',num2str(num),'.mat'];
s1 = load(mat_name1);
U_mat = full(s1.U_mat);

nshape = [num,num,2*num];


for a1 = 1:2*N+1
    for a2 = 1:2*N+1
        U_avg = compute_Uavg(U_vec0(ct1,:), nshape, ind0, ind1, ind2, ind3);
        diff_val1 = abs(U_avg).*diff_val;
        Atot = sum(sum(sum(diff_val1)));
        disp(Atot);
    end
end

end

% function U_vec0 = reshape_Uvec(U_vec1, N1)
% nsz0 = (2*N1+1); nsz1 = size(U_vec1,3);
% 
% U_vec0 = zeros(nsz0^2, nsz1);
% ct3 = 1;
% for ct1 = 1:nsz0
%     for ct2 = 1:nsz0
%         U_vec0(ct3,:) = U_vec1(ct1,ct2,:);
%         ct3 = ct3 + 1;
%     end
% end
% end



function U_avg = compute_Uavg(U_vec0, nshape, ind0, ind1, ind2, ind3)

% U1 = reshape(transpose(U_vec0), nshape);
U_avg = (U1(ind0,ind0,ind2) + ...
    U1(ind1,ind0,ind2) + ...
    U1(ind0,ind1,ind2) + ...
    U1(ind0,ind0,ind3) + ...
    U1(ind1,ind1,ind3) + ...
    U1(ind0,ind1,ind3) + ...
    U1(ind1,ind0,ind3) + ...
    U1(ind1,ind1,ind2))/8;

end



