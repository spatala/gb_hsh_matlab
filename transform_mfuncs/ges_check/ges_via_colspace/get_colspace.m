function col1 = get_colspace(ges_mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Finding the column-space for a sparse matrix
%%%% Using QR decomposition and the stepform of R matrix. 
%%%% https://math.stackexchange.com/questions/748500/how-to-find-linearly-independent-columns-in-a-matrix/748538#748538
[Q,R] = qr(ges_mat);
num_rows = size(ges_mat,1);
j_inds = zeros(num_rows, 1); ct2 = 1; st1 = 0;
for ct1=1:num_rows
    ct1
    ind1 = find(abs(R(:,ct1)), 1, 'last');
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
j_inds(ct2:end) =[]; col1 = ges_mat; col1(:,j_inds) = [];
end