function null_mat = generate_gb_null_ret(top_dir,pt_grp, Nmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_inds = mbp_inds_ab_array(symm_orders);

a_val = symm_orders(:,1)'; b_val = symm_orders(:,2)';
c_val = min(a_val, b_val);
num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));

null_mat = sparse((2*Nmax+1)^2,num_cols);

for ct1 = 1:size(symm_orders,1)
    a = a_val(ct1); b = b_val(ct1);
    disp("++++++++++++++++++++++++++++++++++++++++++++++++++")
    disp([a,b]);
    
    c = min(a,b);
    na = 2*a+1; nb = 2*b+1; ng = 2*c+1;
    alpha = -a:a; beta = -b:b; gamma = -c:c;
    
    cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
    ind_start = find(cond1,1);
    ind_stop  = find(cond1,1,'last');
    
    for e = abs(a - b):(a + b)
        epsilon = -e:e;
        ne = 2 * e + 1;
        
        [C, m1, m2] = clebsch_gordan(a, b, e, 0);
        ind = (m1 + a) * nb + m2 + b + 1;
        Cz = sparse(ind, 1, C, na * nb, 1);
        for p = 1:ne % \epsilon
            [C, m1, m2] = clebsch_gordan(a, b, e, -epsilon(p));
            ind = (m1 + a) * nb + m2 + b + 1;
            Ce = sparse(ind, 1, C, na * nb, 1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            M3 = zeros(na * nb, ng);
            for s = 1:na % \alpha
                for t = 1:nb % \beta
                    tbeta = -(alpha(s) + epsilon(p));
                    row1 = (alpha(s) + a) * nb + tbeta + b + 1;
                    row = (alpha(s) + a) * nb + beta(t) + b + 1;
                    
                    if (abs(Ce(row)) > 0)
                        disp([tbeta-beta(t),s,t, row-row1]);
%                         disp(alpha(s) + beta(t) + epsilon(p))
                    end
                    
                    for u = 1:ng % \gamma
                    

                        M3(row, u) = M3(row, u) + ...
                            (Ce(row)/ sqrt(2 * e + 1)) * ...
                            Cz((gamma(u) + a) * nb - gamma(u) + b + 1);
                    end
                end
            end
            M3 = M3 * sqrt(2) * sqrt((2 * a + 1) * (2 * b + 1)) / pi;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            row1 = e^2 + e + epsilon(p) + 1;
            null_mat(row1,ind_start:ind_stop) = transpose(M3(:));
        end
    end
end
% S = spnull(null_mat);
% 
% mat_name = [data_fname0, ...
%     'Sarr_gbnull_nmax_',num2str(Nmax),'.mat'];
% save(mat_name,'S');

end