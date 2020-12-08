function SMvec=calc_Mvec_symm(top_dir, pt_grp, Nmax, coeffs_typ, rots)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,coeffs_typ,'_',num2str(Nmax),'/'];
mat_name = [data_fname0,'symm_ab_',pt_grp,'_',...
    coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);

mat_name = [top_dir, 'data_files/','ptgrp_',pt_grp,'/', ...
                    coeffs_typ,'_',num2str(Nmax),'/',...
            'Sarr_MabInds_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); tot_Uprops = s1.tot_Uprops;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_inds = mbp_inds_ab_array(symm_orders);
vec_inds = tot_inds(unique(tot_Uprops(:,1)),:);
nrots = size(rots,3);
s = convert_gbrots(rots);
num_rows = size(S,1);
Mvec = compute_Mvec(s, nrots, num_rows, vec_inds);
SMvec = Mvec*S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end