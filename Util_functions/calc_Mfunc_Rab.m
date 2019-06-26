% function M = calc_Mfunc_Rab(Rab, a,b)
function [M, inds, Rab_inds] = calc_Mfunc_Rab(Rab, a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Get the value of MBP functions for a given $(a,b)$ and R^(ab) 
%%%%%%% matrix
%%%%%%%
%%%%%%% Input
%%%%%%% Rab: (2*a+1)*(2*b+1) square matrix
%%%%%%%         Is given by kron(Ua(g1), Ub(g2)
%%%%%%% 
%%%%%%% a, b: Integers
%%%%%%%         order of MBP function (gamma_MBP^{a,b}_{alpha, beta})
%%%%%%%
%%%%%%%
%%%%%%% Output
%%%%%%% M: 1 X N row vector.
%%%%%%%         The value of the function MBP function 
%%%%%%%         (gamma_MBP^{a,b}_{alpha, beta})
%%%%%%%         for different values of a, b, gamma, alpha, beta
%%%%%%%         are calculated for the GB parameter (given by rot)
%%%%%%%

c = min(a, b);

Na = 2*a; Nb = 2*b; Nc = 2*min(a,b);
tf1 = (Na+1)*(Nb+1); num_inds = tf1*(Nc+1);

% M = zeros(1, num_inds);
% for ct1=1:Nc+1
%     g1 = ct1-(c+1); ai = g1 + a + 1; bi = -g1 + b + 1;
%     ind_row = (ai-1)*(2*b+1)+bi;
%     M((ct1-1)*tf1+1:ct1*tf1) = Rab(ind_row, :);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_inds = (Na+1)*(Nb+1)*(Nc+1);
ct1 = 1;
M = zeros(1, num_inds);
inds = zeros(num_inds,6);
Rab_inds = zeros(num_inds,3);
for g1=-c:c
    al_j = g1; be_j = -g1;
    for al_i=-a:a
        for be_i=-b:b
            ai = al_i + a + 1; aj = al_j + a + 1;
            bi = be_i + b + 1; bj = be_j + b + 1;
            ind_row = (ai-1)*(2*b+1)+bi;
            ind_col = (aj-1)*(2*b+1)+bj;
            inds(ct1,:) = [ct1, a, b, g1, al_i, be_i];
            Rab_inds(ct1,:) = [ct1, ind_row, ind_col];
            M(ct1) = Rab(ind_row, ind_col);
            ct1 = ct1 + 1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if a == b
%     PI_ab = sqrt(2*a+1);
% else
%     if (a > b)
%         PI_ab = sqrt(prod(2*b+1:2*a+1));
%     else
%         PI_ab = sqrt(prod(2*a+1:2*b+1));
%     end
% end
% 
% M = PI_ab*M/(2*pi*pi);
end