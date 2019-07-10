clear all; clc;


t_ax = [1/sqrt(2),1/sqrt(2),0];
t_ang = pi;

syms ta real;

q = 0;
Q = [ta,ta,0];

a = q-1i*Q(3);
lb = -Q(2)-1i*Q(1);

for a_val = 1:5
    N = 2*a_val;
    R = sym(zeros(N+1, N+1));
    
    
    sd = (Q(2)-1i*Q(1))/(-Q(2)-1i*Q(1));
    R(1,N+1) = lb^N;
    for d = 2:N+1
        R(d,N-d+2) = R(d-1,N-d+3)*sd;
    end
    
    R1 = double(subs(R, ta, 1/sqrt(2)));
    R2 = rotation(t_ax,t_ang,N);
    norm(R1(:) - R2(:))
end

%     % Force matrix to be orthogonal
%     R = R^r;
%     [U,~,V] = svd(R);
%     R = U*V';
% return;