function [] = check_mbp_ab(Nmax)
% 
% For given values of (a,b), the set of basis functions have the
% folloing indices
%   + mbp: N X 6 matrix. The columns are as follows:
%     - Index number for fixed a, b.
%     - Value of 'a'
%     - Value of 'b'
%     - Value of 'gamma' - range [-c, c], where $c = min(a,b)$
%     - Value of 'alpha' - range [-a, a]
%     - Value of 'beta'  - range [-b, b]
%   
% In the following code, we wrote two functions:
%   + ind_to_gab: Converts the Index number to (gamma, alpha, beta) for
%   given mbp_inds.
%   + gab_to_inds: Convers (a, b, gamma, alpha, beta) to index.
% 
% 


a = randi([0,Nmax]);
b = randi([0,Nmax]);

inds = mbp_ab(a,b);


disp([Nmax, a, b])
[diff_g, diff_a, diff_b] = ind_to_gab(inds);
disp(max(diff_g))
disp(max(diff_a))
disp(max(diff_b))
diff_i1 = gab_to_inds(gamma_val, alpha_val, beta_val, a, b);
disp(max(diff_i1))
end


function [diff_g, diff_a, diff_b] = ind_to_gab(inds)
n1 = size(inds,1);

beta_val = zeros(n1,1);
alpha_val = zeros(n1,1);
gamma_val = zeros(n1,1);


for i1 = 1:n1
    %%%% Convert i1 to (gamma, alpha, beta)
    bi = rem(i1-1,2*b+1);
    beta_val(i1) = bi - b;
    tbi = beta_val(i1) + b + 1;
    
    t1_val = (i1 - tbi)/(2*b+1);
    j1 = rem(t1_val, (2*a+1));
    alpha_val(i1) = j1+1 - (a+1);
    tai = alpha_val(i1) + a + 1;
    
    l1 = t1_val - (tai-1);
    tgi = l1/(2*a+1)+1;
    gamma_val(i1) = tgi - c - 1;
end


diff_g = (abs(inds(:,4) - gamma_val));
diff_a = (abs(inds(:,5) - alpha_val));
diff_b = (abs(inds(:,6) - beta_val));

end

function diff_i1 = gab_to_inds(gamma_val, alpha_val, beta_val, a, b)
c = min(a,b);
n1 = size(gamma_val,1);
diff_i1 = zeros(n1,1);
for i1 = 1:n1
    %%%% Convert (gamma, alpha, beta) to i1
    t_g = c + gamma_val(i1) + 1;
    t_a = a + alpha_val(i1) + 1;
    t_b = b + beta_val(i1) + 1;
    
    t_i1 = (t_g - 1)*(2*a+1)*(2*b+1) + (t_a - 1)*(2*b+1) + t_b;
    diff_i1(i1) = abs(t_i1 - i1);
end
end