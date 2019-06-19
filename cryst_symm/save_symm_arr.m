function save_symm_arr(a_val, b_val, v, col, pt_grp,fname)
S = orth(v(:,col));
% if size(S,2)>1
%     [ct1, ct2, size(S,2)]
% end
c_val = min(a_val, b_val);
S = repmat(S, 2*c_val+1, 1);

mat_name = [fname,'/ptgrp_',pt_grp,'/cryst_symm/Sarr_',...
    num2str(a_val),'_',num2str(b_val),'.mat'];
save(mat_name,'S');
end