%%% Convert random grain boundary to SO(3)XSO(3) parameterization

clear all; clc;
s1 = load('rand_gb.mat');
gbs = s1.gbs;
bpn = gbs(:,5:7);
rot_quats = gbs(:,1:4);

num_gbs = size(gbs,1);
check_err = zeros(num_gbs, 2);
for ct1=1:num_gbs
    bpn1 = bpn(ct1,:);
    x = cross(bpn1, [0,0,1]);
    if (norm(x) < 1e-14)
        if (dot(bpn1, [0,0,1]) > 0)
            g1 = eye(3);
            bpn1 = [0,0,1];
        else
            g1 = vrrotvec2mat([1,0,0,pi]);
            bpn1 = [0,0,-1];
        end
    else
        v_mat=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
        c = dot(bpn1, [0,0,1]);
        R = eye(3) + v_mat + (v_mat^2)*(1/(1+c));
        g1 = R;
        M1 = quat2mat(rot_quats(ct1,:));
        g2 = g1*M1;
    end
    
    %%% Check
    M1 = (g1^(-1))*g2;
    quat1 = Mat2Quat(M1);
    check_err(ct1, :) = [min(norm(quat1 - rot_quats(ct1,:)), norm(quat1 + rot_quats(ct1,:))), norm(bpn1 - ((g1^(-1))*([0,0,1]'))')];
end

max(max(check_err))



