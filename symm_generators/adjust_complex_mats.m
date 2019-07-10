function M1 = adjust_complex_mats(M1, tol1)
if nargin<2
    tol1 = 1e-14;
end
M1(abs(M1(:)) < tol1) =  0;

M1_re = adjust_vals(real(M1), tol1);
M1_im = adjust_vals(imag(M1), tol1);

M1 = M1_re + 1i*M1_im;
end

function M1 = adjust_vals(M1, tol1)
if nargin<2
    tol1 = 1e-14;
end
ct1= 0; M1(abs(M1(:)-ct1) < tol1) = ct1;
ct1= 1; M1(abs(M1(:)-ct1) < tol1) = ct1;
ct1=-1; M1(abs(M1(:)-ct1) < tol1) = ct1;
end
