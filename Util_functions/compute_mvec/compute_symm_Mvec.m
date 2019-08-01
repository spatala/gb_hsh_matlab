function symm_Mvec = compute_symm_Mvec(rots, S, symm_num, data_fname0, Nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input:
%%%% 1) rots:
%%%%    All the N GB rotation matrices (3 X 6 X N)
%%%% 2) pt_grp:
%%%%    Point-group symmetry of the underlying crystal    
%%%% 3) Nmax:
%%%%    Maximum of the allowed (a,b) values.
%%%% 

num_rows = size(S,1);
nrots = size(rots,3); s = convert_gbrots(rots);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_MabInds_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
tot_Uprops=s1.tot_Uprops; ind_ranges=s1.ind_ranges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Mvec = compute_Mvec(s, nrots, num_rows, symm_num, ...
                                        tot_Uprops, ind_ranges);
symm_Mvec = Mvec*S(:,symm_num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end