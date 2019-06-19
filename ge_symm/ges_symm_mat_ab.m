clear all; clc;

addpath(genpath('../Util_functions/'));
addpath(genpath('../GB_Parameters/'));


a_val = 4; b_val = 4;

tot_inds = mbp_ab(4,4);
num_inds = size(tot_inds,1);

% ges_mat = sparse(num_inds, num_inds);
ges_mat = zeros(num_inds, num_inds);

for ct1=1:num_inds
%     ct1
    a1 = tot_inds(ct1,2);
    b1 = tot_inds(ct1,3);
    gamma1 = tot_inds(ct1,4);
    alpha1 = tot_inds(ct1,5);
    beta1 = tot_inds(ct1,6);
    
    a2 = b1;
    b2 = a1;
    gamma2 = gamma1;
    alpha2 = beta1;
    beta2 = alpha1;
    
    ind1 = find((tot_inds(:,2) == a2) & ...
        (tot_inds(:,3) == b2) & ...
        (tot_inds(:,4) == gamma2) & ...
        (tot_inds(:,5) == alpha2) & ...
        (tot_inds(:,6) == beta2));
    
    ges_mat(ct1,ind1) = 1;
    ges_mat(ct1,ct1) = (-1)^(a1+b1);
    
end

rmpath(genpath('../Util_functions/'));
rmpath(genpath('../GB_Parameters/'));


% col1 = double(colspace(sym(ges_mat)));

rank1 = rank(ges_mat);
[Q,R] = qr(ges_mat);

st1 = 0;
j_inds = zeros(num_inds, 1);
ct2 = 1;
for ct1=1:num_inds
    ct1
    ind1 = max(find(abs(R(:,ct1))));
    if ct1 == 1
        st1 = ind1;
    else
        if (ind1 > st1)
            st1 = ind1;
        else
            j_inds(ct2) = ct1;
            ct2 = ct2 + 1;
        end
    end
end
j_inds(ct2:end) =[];
col1 = ges_mat;
col1(:,j_inds) = [];

save('Y_ges.mat', 'col1');

Q_mat = col1*col1';
[v, d] = eig(Q_mat);
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
S = orth(v(:,col));
if any(col)
    S = orth(v(:,col));
end
save('Sarr.mat','S');