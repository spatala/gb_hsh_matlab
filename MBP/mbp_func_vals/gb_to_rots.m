function [g1,g2] = gb_to_rots(M1, bpn1)
zvec = [0,0,1];
x = cross(bpn1, zvec);
if (norm(x) < 1e-14)
    if (dot(bpn1, [0,0,1]) > 0)
        g1 = eye(3); bpn1 = zvec;
    else
        g1 = vrrotvec2mat([1,0,0,pi]); bpn1 = -zvec;
    end
else
    v_mat=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
    c = dot(bpn1, zvec);
    R = eye(3) + v_mat + (v_mat^2)*(1/(1+c));
    g1 = R;
end
g2 = g1*M1;

end