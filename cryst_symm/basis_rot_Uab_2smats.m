clear all; clc;

addpath(genpath('../Util_functions/'));
addpath(genpath('../GB_Parameters/'));

Nmax = 4;
for ct1=0:Nmax
    for ct2=0:Nmax
        a_val = ct1; Na = 2*a_val;
        b_val = ct2; Nb = 2*b_val;
        
        % Symmetry generators
        a_ax1 = [0 0 1]; a_ang1 = pi;
        b_ax1 = [0 0 1]; b_ang1 = 0;
        
        a_ax2 = [0 0 1]; a_ang2 = 0;
        b_ax2 = [0 0 1]; b_ang2 = pi;
        
        
        R1 = kron(rotation(a_ax1,a_ang1,Na),rotation(b_ax1,b_ang1,Nb));
        [v,d] = eig(R1);
        col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
        if any(col)
            X0 = eye(size(R1)); P0 = X0*X0';
            Y1 = orth(v(:,col)); Q1 = Y1*Y1';
            [v, d] = eig(P0*Q1*P0);
            col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
            if any(col)
                X1 = orth(v(:,col)); P1 = X1*X1';
                R2 = kron(rotation(a_ax2,a_ang2,Na),rotation(b_ax2,b_ang2,Nb));
                [v,d] = eig(R2);
                col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
                if any(col)
                    Y2 = orth(v(:,col)); Q2 = Y2*Y2';
                    [v, d] = eig(P1*Q2*P1);
                    col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
                    
                    if any(col)
                        
                        S = orth(v(:,col));
                        if size(S,2)>1
                            [ct1, ct2, size(S,2)]
                        end
                        c_val = min(a_val, b_val);
                        S = repmat(S, 2*c_val+1, 1);
                        mat_name = ['Sarr_',num2str(a_val),'_',num2str(b_val),'.mat'];
                        save(mat_name,'S');
                        
                    end
                end
            end
        end
    end
end
%


rmpath(genpath('../Util_functions/'));
rmpath(genpath('../GB_Parameters/'));