function [ga_s, gb_s, num_gen] = get_symmgen_mats(pt_grp)
if (strcmp(pt_grp,'O'))
    num_gen = 2;
    
    a_ax1 = [0 0 1]; a_ang1 = pi/2;
    b_ax1 = [0 0 1]; b_ang1 = 0;
    ga_s1 = vrrotvec2mat([a_ax1, a_ang1]);
    gb_s1 = vrrotvec2mat([b_ax1, b_ang1]);
    
    a_ax2 = [1 0 0]; a_ang2 = pi/2;
    b_ax2 = [1 0 0]; b_ang2 = 0;
    ga_s2 = vrrotvec2mat([a_ax2, a_ang2]);
    gb_s2 = vrrotvec2mat([b_ax2, b_ang2]);
    
    a_ax3 = [0 0 1]; a_ang3 = 0;
    b_ax3 = [0 0 1]; b_ang3 = pi/2;
    ga_s3 = vrrotvec2mat([a_ax3, a_ang3]);
    gb_s3 = vrrotvec2mat([b_ax3, b_ang3]);
    
    a_ax4 = [1 0 0]; a_ang4 = 0;
    b_ax4 = [1 0 0]; b_ang4 = pi/2;
    ga_s4 = vrrotvec2mat([a_ax4, a_ang4]);
    gb_s4 = vrrotvec2mat([b_ax4, b_ang4]);
    
    ga_s = cell(num_gen^2,1);
    ga_s{1} = ga_s1; ga_s{2} = ga_s2;
    ga_s{3} = ga_s3; ga_s{4} = ga_s4;
    gb_s = cell(num_gen^2,1);
    gb_s{1} = gb_s1; gb_s{2} = gb_s2;
    gb_s{3} = gb_s3; gb_s{4} = gb_s4;
end

if (strcmp(pt_grp,'D2'))
    num_gen = 2;
    
    a_ax1 = [0 0 1]; a_ang1 = pi;
    b_ax1 = [0 0 1]; b_ang1 = 0;
    ga_s1 = vrrotvec2mat([a_ax1, a_ang1]);
    gb_s1 = vrrotvec2mat([b_ax1, b_ang1]);
    
    a_ax2 = [1 0 0]; a_ang2 = pi;
    b_ax2 = [1 0 0]; b_ang2 = 0;
    ga_s2 = vrrotvec2mat([a_ax2, a_ang2]);
    gb_s2 = vrrotvec2mat([b_ax2, b_ang2]);
    
    a_ax3 = [0 0 1]; a_ang3 = 0;
    b_ax3 = [0 0 1]; b_ang3 = pi;
    ga_s3 = vrrotvec2mat([a_ax3, a_ang3]);
    gb_s3 = vrrotvec2mat([b_ax3, b_ang3]);
    
    a_ax4 = [1 0 0]; a_ang4 = 0;
    b_ax4 = [1 0 0]; b_ang4 = pi;
    ga_s4 = vrrotvec2mat([a_ax4, a_ang4]);
    gb_s4 = vrrotvec2mat([b_ax4, b_ang4]);
    
    ga_s = cell(num_gen^2,1);
    ga_s{1} = ga_s1; ga_s{2} = ga_s2;
    ga_s{3} = ga_s3; ga_s{4} = ga_s4;
    gb_s = cell(num_gen^2,1);
    gb_s{1} = gb_s1; gb_s{2} = gb_s2;
    gb_s{3} = gb_s3; gb_s{4} = gb_s4;
end

if (strcmp(pt_grp,'C2'))
    num_gen = 1;
    
    a_ax1 = [0 0 1]; a_ang1 = pi;
    b_ax1 = [0 0 1]; b_ang1 = 0;
    ga_s1 = vrrotvec2mat([a_ax1, a_ang1]);
    gb_s1 = vrrotvec2mat([b_ax1, b_ang1]);
    
    a_ax2 = [0 0 1]; a_ang2 = 0;
    b_ax2 = [0 0 1]; b_ang2 = pi;
    ga_s2 = vrrotvec2mat([a_ax2, a_ang2]);
    gb_s2 = vrrotvec2mat([b_ax2, b_ang2]);
    
    ga_s = cell(num_gen^2,1);
    ga_s{1} = ga_s1; ga_s{2} = ga_s2;
    gb_s = cell(num_gen^2,1);
    gb_s{1} = gb_s1; gb_s{2} = gb_s2;
end
end