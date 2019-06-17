clear all; clc;

N = 1;
tot_inds = mbp_inds(N);
num_inds = size(tot_inds,1);

ges_mat = zeros(num_inds, num_inds);

for ct1=1:num_inds
    ct1
    a1 = tot_inds(ct1,3);
    b1 = tot_inds(ct1,4);
    gamma1 = tot_inds(ct1,5);
    alpha1 = tot_inds(ct1,6);
    beta1 = tot_inds(ct1,7);
    
    a2 = b1;
    b2 = a1;
    gamma2 = gamma1;
    alpha2 = beta1;
    beta2 = alpha1;
    
    ind1 = find((tot_inds(:,3) == a2) & ...
        (tot_inds(:,4) == b2) & ...
        (tot_inds(:,5) == gamma2) & ...
        (tot_inds(:,6) == alpha2) & ...
        (tot_inds(:,7) == beta2));
    
    ges_mat(ct1,ind1) = (-1)^(a1+b1);
    ges_mat(ct1,ct1) = 1;
    
    % ind1
    % tot_inds([ct1,ind1],3:end)
    
end

col1 = double(colspace(sym(ges_mat)));


Q_mat = col1*col1';
[v, d] = eig(Q_mat);
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
S = orth(v(:,col));
if any(col)
    S = orth(v(:,col));
end
save('Sarr.mat','S');