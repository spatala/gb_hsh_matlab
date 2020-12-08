function [] = check_mbp_ab()
clear all; clc;

Nmax = 10;
% a = randi([0,Nmax]);
% b = randi([0,Nmax]);
a = 3; b = 4;
c = min(a,b);

inds = mbp_ab(a,b);
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

%%%% Convert (gamma, alpha, beta) to i1
end
end



