%%% Convert random grain boundary to SO(3)XSO(3) parameterization

clear all; clc;
s1 = load('rand_gb.mat');
gbs = s1.gbs;
bpn = gbs(:,5:7);
rot_quats = gbs(:,1:4);

num_gbs = size(gbs,1);
check_err = zeros(num_gbs, 2);

rot_mats = zeros(3,6,4);
for ct1=1:num_gbs
% for ct1=1:4
    bpn1 = bpn(ct1,:);
    M1 = quat2mat(rot_quats(ct1,:));
    
    [g1,g2] = gb_to_rots(M1, bpn1);
    rot_mats(:,:,ct1) = [g1, g2];
%     %%% Check
%     M1 = (g1^(-1))*g2;
%     quat1 = Mat2Quat(M1);
%     check_err(ct1, :) = [min(norm(quat1 - rot_quats(ct1,:)), ...
%         norm(quat1 + rot_quats(ct1,:))), ...
%         norm(bpn1 - ((g1^(-1))*([0,0,1]'))')];
end

% max(max(check_err))

save('rand_gb_rots.mat','rot_mats');



