function [ga_s, gb_s, num_gen, Laue] = get_symmgen_angs(pt_grp)
%
% Function to set the symmetry generators for each point group.
% The symmetries are for grain boundaries in the $SO(3) \times SO(3)$
% parameterization, i.e. $(O_a, O_b) ~ (O_b S_i, O_a S_j)$
% 
% - Input:
%   + pt_grp:	(string) 
%   			Point-group of the underlying crystal
% 
% - Output:
%   + ga_s, gb_s: (Cell arrays)
%    		grain boundary symmetry generators.
%   + num_gen: (integer)
%    		Number of generators for the point-group.
%   + Laue: (string)
%    		0 for non-Laue and 1 for Laue group.
% 
% - Notes:
% + Only proper rotation elements are set.
%   + The symmetry operation is provided using three angles
%   	- w : rotation angle
%   	- th: polar angle of the axis
%   	- ph: azimuthal angle of the axis
%


if (strcmp(pt_grp,'O') || strcmp(pt_grp,'Oh'))
    num_gen = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ga_s1 = [pi/2, 0, 0]; gb_s1 = [0,0,0];
    ga_s2 = [pi/2, pi/2, 0]; gb_s2 = [0,0,0];
    ga_s3 = [0,0,0]; gb_s3 = [pi/2, 0, 0];
    ga_s4 = [0,0,0]; gb_s4 = [pi/2, pi/2, 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ga_s = cell(num_gen^2,1);
    ga_s{1} = ga_s1; ga_s{2} = ga_s2;
    ga_s{3} = ga_s3; ga_s{4} = ga_s4;
    gb_s = cell(num_gen^2,1);
    gb_s{1} = gb_s1; gb_s{2} = gb_s2;
    gb_s{3} = gb_s3; gb_s{4} = gb_s4;
    
    if strcmp(pt_grp,'O')
        Laue = 0;
    else
        Laue = 1;
    end
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
    
    Laue = 0;
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
    
    Laue = 0;
end

if (strcmp(pt_grp,'C1'))
    num_gen = 1;
    
    a_ax1 = [0 0 1]; a_ang1 = 0;
    b_ax1 = [0 0 1]; b_ang1 = 0;
    ga_s1 = vrrotvec2mat([a_ax1, a_ang1]);
    gb_s1 = vrrotvec2mat([b_ax1, b_ang1]);
    
    a_ax2 = [0 0 1]; a_ang2 = 0;
    b_ax2 = [0 0 1]; b_ang2 = 0;
    ga_s2 = vrrotvec2mat([a_ax2, a_ang2]);
    gb_s2 = vrrotvec2mat([b_ax2, b_ang2]);
    
    ga_s = cell(num_gen^2,1);
    ga_s{1} = ga_s1; ga_s{2} = ga_s2;
    gb_s = cell(num_gen^2,1);
    gb_s{1} = gb_s1; gb_s{2} = gb_s2;
    
    Laue = 0;
end

end