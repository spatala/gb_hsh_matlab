clear all; clc;

curr_pwd = split(pwd,'/');
top_dir = '';
for ct1=1:length(curr_pwd)
    top_dir = strcat(top_dir,curr_pwd{ct1},'/');
    if (strcmp(curr_pwd{ct1},'gb_hsh_matlab'))
        break;
    end
end
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(util_dir));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation along the Z-axis.
ind_start = 1; 
num = 10; sc = 3;

a_tot = zeros(sc*num,1); lb_tot = zeros(sc*num,1);
q_tot = zeros(sc*num,1); Q_tot = zeros(sc*num,3);

ax = [0,0,1]; ax = ax./sum(ax.^2,2);
an = rand(num,1)*2*pi;

q = cos(an/2); Q = sin(an/2).*ax;
a = q - 1i*Q(:,3); lb = -Q(:,2)-1i*Q(:,1);

ind_stop = ind_start - 1 + num;

a_tot(ind_start:ind_stop,:) = a; lb_tot(ind_start:ind_stop,:) = lb;
q_tot(ind_start:ind_stop,:) = q; Q_tot(ind_start:ind_stop,:) = Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation of pi and axis in the X-Y plane.
ind_start = ind_stop + 1;
phi = rand(num,1)*2*pi; 
ax = [cos(phi), sin(phi), 0*phi]; ax = ax./sum(ax.^2,2);
an = pi;

q = cos(an/2); Q = sin(an/2).*ax;
a = q - 1i*Q(:,3); lb = -Q(:,2)-1i*Q(:,1);

ind_stop = ind_start - 1 + num;

a_tot(ind_start:ind_stop,:) = a; lb_tot(ind_start:ind_stop,:) = lb;
q_tot(ind_start:ind_stop,:) = q; Q_tot(ind_start:ind_stop,:) = Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation of pi and axis in the X-Y plane.
ind_start = ind_stop + 1;
u1 = rand(num,1); u2 = rand(num,1); u3 = rand(num,1);
q = sqrt(1-u1).*sin(2*pi*u2);
Q = [sqrt(1-u1).*cos(2*pi*u2),...
    sqrt(u1).*sin(2*pi*u3),...
    sqrt(u1).*cos(2*pi*u3),];

a = q - 1i*Q(:,3); lb = -Q(:,2)-1i*Q(:,1);

ind_stop = ind_start - 1 + num;

a_tot(ind_start:ind_stop,:) = a; lb_tot(ind_start:ind_stop,:) = lb;
q_tot(ind_start:ind_stop,:) = q; Q_tot(ind_start:ind_stop,:) = Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_val = 4;
Na = 2*a_val;
U1 = zeros(Na+1,Na+1,sc*num);
for alp_val = -a_val:a_val
    ri = alp_val + a_val + 1;
    for al_val = -a_val:a_val
        cj = al_val + a_val + 1;
        U1(ri,cj,:) = rotation_wo_svd(a_tot,lb_tot,q_tot, ...
            Q_tot, a_val, alp_val, al_val);
    end
end

for ct1 = 1:sc*num
    tU1 = U1(:,:,ct1);
    
    quat = [q_tot(ct1), Q_tot(ct1,:)];
    ax_ang = vrrotmat2vec(quat2mat(quat));
    ax = ax_ang(1:3); an = ax_ang(4);

    U2 = rotation(ax,an,Na);
    norm(tU1 - U2)

end
