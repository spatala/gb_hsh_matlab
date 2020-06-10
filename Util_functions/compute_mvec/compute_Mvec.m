function Mvec = compute_Mvec(s, nrots, num_cols, vec_inds)


a1_rots  = s.a1;
a2_rots  = s.a2;
lb1_rots = s.lb1;
lb2_rots = s.lb2;
q1_rots  = s.q1;
q2_rots  = s.q2;
Q1_rots  = s.Q1;
Q2_rots  = s.Q2;


a_val   = vec_inds(:,3);
alp_val = vec_inds(:,6);
al_val  = vec_inds(:,5);
b_val   =  vec_inds(:,4);
bep_val =  vec_inds(:,7);
be_val  = -vec_inds(:,5);


Mvec = sparse(nrots,num_cols);

for ct2 = 1:size(vec_inds,1)
    a = a_val(ct2); b = b_val(ct2);
    U1 = rotation_wo_svd(a1_rots,lb1_rots, ...
        q1_rots, Q1_rots, a, alp_val(ct2), al_val(ct2));
    U2 = rotation_wo_svd(a2_rots,lb2_rots, ...
        q2_rots, Q2_rots, b, bep_val(ct2), be_val(ct2));
    PI_ab = realsqrt((2*a+1)*(2*b+1));
    % M = PI_ab*M/(2*pi*pi);
    Mvec(:,vec_inds(ct2,1)) = PI_ab*(U1.*U2)/(sqrt(2*(pi^3)));
end
end