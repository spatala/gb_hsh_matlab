function rots = mbp_to_rots(mbp)
%
% Convert GB parameters from Misorientation-Boundaryplane (MBP) to
% SO3xSO3 parametersation (rots).
%
% Input
% mbp: 3 X 4 matrix
%         First three columns are the misorientation matrix and the
%         last column is the boundary-plane.
%
% Output
% rots: 3 X 6 matrix
%         First three columns are the first orientation matrix (g1)
%         and the second three are the second orientaiton matrix (g2)
%


M1 = mbp(:,1:3); bpn1 = mbp(:,4);
zvec = [0,0,1];
x = cross(bpn1, zvec);
if (norm(x) < 1e-14)
    if (dot(bpn1, [0,0,1]) > 0)
        g1 = eye(3);
    else
        g1 = vrrotvec2mat([1,0,0,pi]);
    end
else
    v_mat=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
    c = dot(bpn1, zvec);
    R = eye(3) + v_mat + (v_mat^2)*(1/(1+c));
    g1 = R;
end
g2 = g1*M1;
rots = [g1, g2];

end