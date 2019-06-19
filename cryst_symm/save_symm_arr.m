function save_symm_arr(a_val, b_val, v, col, pt_grp,fname)
S = orth(v(:,col));
mat_name = [fname,'/ptgrp_',pt_grp,'/cryst_symm/Sarr_',...
    num2str(a_val),'_',num2str(b_val),'.mat'];
save(mat_name,'S');
end