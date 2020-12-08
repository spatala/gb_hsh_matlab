function s = convert_gbrots(rot_mats)
%
% Convert (g1,g2) rotation matrices (of GBs) into 
% quaternions and Cayley-Klein parameters
% 
% - Input:
%   + rot_mats: Array of size 3X6XN
%   		N GB parameters in SO3XSO3 parameterization
% 
% - Output:
%   + s: struct of N-arrays (8 elements)
%   	- ((s.a1, s.lb1, s.q1, s.Q1), (s.a2, s.lb2, s.q2, s.Q2))
%   	- s.a  :   q_0 - i*q_3
%   	- s.lb : -q_2 - i*q_1
%   	- s.q  : q_0
%   	- s.Q  : (q_1, q_2, q_3)
% 
% - Notes:
% 	+ This struct is used to compute MBP basis functions
%     
nrots = size(rot_mats,3);

a1_rots = zeros(nrots,1);
lb1_rots = zeros(nrots,1);
q1_rots = zeros(nrots,1);
Q1_rots  = zeros(nrots,3);
a2_rots = zeros(nrots,1);
lb2_rots = zeros(nrots,1);
q2_rots = zeros(nrots,1);
Q2_rots  = zeros(nrots,3);

for ct1=1:nrots
    rots1 = rot_mats(:,:,ct1);
    g1 = rots1(:,1:3);
    ax_ang1 = vrrotmat2vec(g1');
    ax1 = ax_ang1(1:3);
    an1 = ax_ang1(4);
    quat1 = [cos(an1/2), ax1*sin(an1/2)];
    q1 = quat1(1);
    Q1=quat1(2:4);
    a1 = q1-1i*Q1(3);
    lb1 = -Q1(2)-1i*Q1(1);
    
    a1_rots(ct1,:) = a1;
    lb1_rots(ct1,:) = lb1;
    q1_rots(ct1,:) = q1;
    Q1_rots(ct1,:) = Q1;
    
    g2 = rots1(:,4:6);
    ax_ang2 = vrrotmat2vec(g2');
    ax2 = ax_ang2(1:3);
    an2 = ax_ang2(4);
    quat2 = [cos(an2/2), ax2*sin(an2/2)];
    q2 = quat2(1);
    Q2=quat2(2:4);
    a2 = q2-1i*Q2(3);
    lb2 = -Q2(2)-1i*Q2(1);
    
    a2_rots(ct1,:) = a2;
    lb2_rots(ct1,:) = lb2;
    q2_rots(ct1,:) = q2;
    Q2_rots(ct1,:) = Q2;
end

s = struct();
s.a1 = a1_rots;
s.lb1 = lb1_rots;
s.q1 = q1_rots;
s.Q1 = Q1_rots;

s.a2 = a2_rots;
s.lb2 = lb2_rots;
s.q2 = q2_rots;
s.Q2 = Q2_rots;

end
