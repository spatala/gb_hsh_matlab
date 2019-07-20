function Mvec = compute_Mvec(s, nrots, num_rows, symm_num, ...
                                        tot_Uprops, ind_ranges)


a1_rots  = s.a1;  a2_rots  = s.a2;
lb1_rots = s.lb1; lb2_rots = s.lb2;
q1_rots  = s.q1;  q2_rots  = s.q2;
Q1_rots  = s.Q1;  Q2_rots  = s.Q2;


st1 = ind_ranges(symm_num,1); st2 = ind_ranges(symm_num,2);
U_props = tot_Uprops(st1:st2,:); vec_inds = U_props(:,1);
a_val   = U_props(:,2); alp_val = U_props(:,5); al_val  = U_props(:,4);
b_val   =  U_props(:,3); bep_val =  U_props(:,6); be_val  = -U_props(:,4);


Mvec = sparse(nrots,num_rows);

for ct2 = 1:size(vec_inds,1)
    a = a_val(ct2); b = b_val(ct2);
    U1 = rotation_wo_svd(a1_rots,lb1_rots, q1_rots, Q1_rots, a, alp_val(ct2), al_val(ct2));
    U2 = rotation_wo_svd(a2_rots,lb2_rots, q2_rots, Q2_rots, b, bep_val(ct2), be_val(ct2));
    if a == b
        PI_ab = sqrt(2*a+1);
    else
        if (a > b)
            PI_ab = sqrt(prod(2*b+1:2*a+1));
        else
            PI_ab = sqrt(prod(2*a+1:2*b+1));
        end
    end

    %     M = PI_ab*M/(2*pi*pi);
    Mvec(:,vec_inds(ct2)) = PI_ab*(U1.*U2)/(2*(pi^2));
end
end
