function [] = check_ypi_left()

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

s1 = set_vars(); Nmax = s1.Nmax; pt_grp = s1.pt_grp;

a_max = 5;
diff1 = zeros(a_max^2,1);
diff2 = zeros(a_max^2,1);
diff3 = zeros(a_max^2,1);
tct1 = 1;
for a_val = 0:a_max
    for b_val = 0:a_max
        % a_val = 4; b_val = 3;
        Na = 2*a_val; Nb = 2*b_val;
        Mrot = generate_ypi_left_ab(a_val,b_val);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Only works (so far) for
        %%%%% 1) X,Y,Z-pi rotations!; 2) Z-theta rotations!
        Ypi = [0,1,0,pi]; gl1 = vrrotvec2mat(Ypi); gl2 = vrrotvec2mat(Ypi);
        % Zth = [0,0,1,rand()]; gl1 = vrrotvec2mat(Zth); gl2 = vrrotvec2mat(Zth);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gr1 = eye(3); gr2 = eye(3);
        Rr_ab_12 = so4_irrep(gr1,gr2,Na,Nb);
        Rl_ab_12 = so4_irrep(gl1,gl2,Na,Nb);
        trRr_ab_12 = transpose(Rr_ab_12);
        
        tr_Rrot = kron(Rl_ab_12,trRr_ab_12);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mrot1 = calc_Mrot_mat_Rab(tr_Rrot, a_val, b_val);
        Mrot1(abs(Mrot1)<1e-10) = 0; Mrot1 = sparse(Mrot1);
        
        Mdiff = abs(Mrot-Mrot1);
        diff1(tct1) = norm(full(Mrot-Mrot1));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
        rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
        g1 = rots1(:,1:3); g2 = rots1(:,4:6);
        Mvec = calc_Mvec(g1, g2, [a_val, b_val]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ypi = vrrotvec2mat([0,1,0,pi]);
        g1p = ypi*g1; g2p = ypi*g2;
        SMvec = calc_Mvec(g1p, g2p, [a_val, b_val]);
        
        diff2(tct1) = norm(Mvec*Mrot  - SMvec);
        diff3(tct1) = norm(Mvec*Mrot1 - SMvec);
        tct1 = tct1 + 1;
    end
end
end


% function Ypi_symm_mat = check_ypi_left(a_val, b_val)
% % function Ypi_symm_mat = generate_ypi_left(symm_orders)
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
% % num_cols = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
% % nsymm_ab = size(symm_orders,1);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Ypi_symm_mat = sparse(num_cols, num_cols);
%
% % ind_start = 1;
% % for ct1 = 1:nsymm_ab
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     a_val = symm_orders(ct1,1); b_val = symm_orders(ct1,2);
%     Mrot = generate_ypi_left_ab(a_val,b_val);
%
% %     c_val = min(a_val, b_val);
% %     Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val; nsz = (Na+1)*(Nb+1);
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%% Only works (so far) for
% %     %%%%% 1) X,Y,Z-pi rotations!; 2) Z-theta rotations!
% %     Ypi = [0,1,0,pi]; gl1 = vrrotvec2mat(Ypi); gl2 = vrrotvec2mat(Ypi);
% %     % Zth = [0,0,1,rand()]; gl1 = vrrotvec2mat(Zth); gl2 = vrrotvec2mat(Zth);
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     gr1 = eye(3); gr2 = eye(3);
% %     Rr_ab_12 = so4_irrep(gr1,gr2,Na,Nb);
% %     Rl_ab_12 = so4_irrep(gl1,gl2,Na,Nb);
% %     trRr_ab_12 = transpose(Rr_ab_12);
% %     tr_Rrot = kron(Rl_ab_12,trRr_ab_12);
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     Mrot = calc_Mrot_mat_Rab(tr_Rrot, a_val, b_val);
% %     Mrot(abs(Mrot)<1e-10) = 0; Mrot = sparse(Mrot);
%
%
% %     norm(full(Mrot)-full(Mrot1))
%
% %     nsz1 = nsz*(Nc+1);
% %     ind_stop = ind_start + nsz1 - 1;
% %     Ypi_symm_mat(ind_start:ind_stop, ind_start:ind_stop) = Mrot;
% %     ind_start = ind_stop + 1;
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % end
% end
