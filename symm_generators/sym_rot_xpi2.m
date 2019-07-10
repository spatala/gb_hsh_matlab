clear all; clc;


t_ax = [1,0,0];
t_ang = pi/2;

for a_val = 1:20
    N = 2*a_val;
    R2 = rotation(t_ax,t_ang,N);
    
    display(R2)
end


% syms ta real;
% 
% q = ta;
% Q = [ta,0,0];
% 
% a = q-1i*Q(3);
% lb = -Q(2)-1i*Q(1);
% 
% for a_val = 1:5
%     N = 2*a_val;
%     R = sym(zeros(N+1, N+1));
%     
%  
%     % General rotation
%     aa = q^2+Q(3)^2;
%     bb = Q(2)^2+Q(1)^2;
%     
%     R(1,1) = a^N;
%     R(N+1,N+1) = lb^N;
%     j = N/2;
%     for mup = j-(1:2*j)
%         %         display(j-mup+1);
%         R(1,j-mup+1) = realsqrt(prod(1:2*j)/(prod(1:(j+mup))*prod(1:(j-mup))))*a^(j+mup)*lb^(j-mup);
%     end
%     for mpr = j-(1:2*j)
%         for mup = j-(0:2*j-1)
%             rmin = max([0,mup-mpr]);
%             rmax = min([j-mpr,j+mup]);
%             c = j+mup-rmin;
%             d = j-mpr-rmin;
%             e = mpr-mup+rmin;
%             f = rmin;
%             ltau = 1;
%             for n = rmax-(rmin:(rmax-1))
%                 ltau = 1-((c-n+1)*(d-n+1)*bb)/((e+n)*(f+n)*aa)*ltau;
%             end
%             R(j-mpr+1,j-mup+1) = (realsqrt(prod(1:(j+mpr))*prod(1:(j-mpr))*...
%                 prod(1:(j+mup))*prod(1:(j-mup)))*(-1)^rmin*(a^c*a'^d*lb^e*lb'^f)...
%                 /(prod(1:c)*prod(1:d)*prod(1:e)*prod(1:f))*ltau);
%         end
%         R(j-mpr+1,N+1) = realsqrt(prod(1:2*j)/(prod(1:(j+mpr))*prod(1:(j-mpr))))*a'^(j-mpr)*lb^(j+mpr);
%     end
%     display(R)
%     R1 = double(subs(R, ta, 1/sqrt(2)));
%     R2 = rotation(t_ax,t_ang,N);
%     norm(R1(:) - R2(:))
% end
