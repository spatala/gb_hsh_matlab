function [a_ax, a_ang, b_ax, b_ang, num_gen] = get_symm_gens(pt_grp)
if (strcmp(pt_grp,'O'))
    num_gen = 2;
    
    a_ax1 = [0 0 1]; a_ang1 = pi/2;
    b_ax1 = [0 0 1]; b_ang1 = 0;
    
    a_ax2 = [1 0 0]; a_ang2 = pi/2;
    b_ax2 = [1 0 0]; b_ang2 = 0;
    
    a_ax3 = [0 0 1]; a_ang3 = 0;
    b_ax3 = [0 0 1]; b_ang3 = pi/2;
    
    a_ax4 = [1 0 0]; a_ang4 = 0;
    b_ax4 = [1 0 0]; b_ang4 = pi/2;
    
    a_ax = [a_ax1; a_ax2; a_ax3; a_ax4];
    a_ang = [a_ang1; a_ang2; a_ang3; a_ang4];
    b_ax = [b_ax1; b_ax2; b_ax3; b_ax4];
    b_ang = [b_ang1; b_ang2; b_ang3; b_ang4];
end

if (strcmp(pt_grp,'C2'))
    num_gen = 1;
    
    a_ax1 = [0 0 1]; a_ang1 = pi;
    b_ax1 = [0 0 1]; b_ang1 = 0;
    
    a_ax2 = [0 0 1]; a_ang2 = 0;
    b_ax2 = [0 0 1]; b_ang2 = pi;
    
    a_ax = [a_ax1; a_ax2];
    a_ang = [a_ang1; a_ang2];
    b_ax = [b_ax1; b_ax2];
    b_ang = [b_ang1; b_ang2];
end
end