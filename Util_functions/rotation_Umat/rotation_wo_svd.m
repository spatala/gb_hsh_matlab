function U_val = rotation_wo_svd(a,lb,q, Q, a_val, alp_val, al_val)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function to compute elements of SO(3) Irreducible representative 
%%%% U^{a_val}_{alp_val, al_val} for a rotation given by (q,Q) quaterion.
%%%%
%%%% U^{a_val}_{alp_val, al_val}:
%%%%    a_val is the order of the rotation matrix. 
%%%%
%%%% Optimized for array of rotations as input (suppose N rotations).
%%%%
%%%% Input:
%%%% q: 
%%%%    q0 component of the rotation quaternion. Array of size N X 1. N is
%%%%    the number rotations.
%%%% Q: 
%%%%    [q1,q2,q3] components array of the rotation quaternion. Array of
%%%%    size N X 3.
%%%%
%%%% alp_val: Corresponds to row-index.    Range -alp_val:alp_val
%%%% al_val:  Corresponds to column-index. Range -al_val:al_val
%%%%
%%%%
%%%% a = q-1i*Q(3); lb = -Q(2)-1i*Q(1);
%%%% 
%%%% 
%%%% Output:
%%%% U_val:
%%%%    Value of U in the alp_val row and al_val column. Size N X 1, where
%%%%    N is the number of rotations.
%%%%    
%%%% 


%%%%%%
N = 2*a_val;
%%%%%% alp_val represents rows
ri = alp_val + a_val + 1;
%%%%%% al_val represents columns
cj = al_val + a_val + 1;
%%%%%%

ind1 = find(abs(lb) < 1e-10);
ind2 = find(abs(a) < 1e-10);
ind3 = find((abs(a) >= 1e-10) & (abs(lb) >= 1e-10));

ngbs = size(q,1);

U_val = zeros(ngbs,1);

%%%%%%% Rotation about z axis
if size(ind1,1) > 0
    q1 = q(ind1,:); Q1 = Q(ind1,:); a1 = a(ind1,:);
    U_val(ind1) = zaxis_rot(ri, cj, q1, Q1, a1, N);
end

%%%%%%% Rotation about axis in x-y plane
if size(ind2,1) > 0
    Q2 = Q(ind2,:); lb2 = lb(ind2,:);
    U_val(ind2) = xyplane_pi_rot(ri, cj, Q2, lb2, N);
end

if size(ind3,1) > 0
    q3 = q(ind3,:); Q3 = Q(ind3,:); 
    a3 = a(ind3,:); lb3 = lb(ind3,:);
    U_val(ind3) = gen_rot(ri, cj, q3, Q3, a3, lb3, N);
end

end

function U_val = zaxis_rot(ri, cj, q, Q, a, N)
if (ri == cj)
    du = (q+1i*Q(:,3))./(q-1i*Q(:,3));
    if (ri == 1)
        U_val = a.^N; return;
    else
        tR1 = a.^N;
        for d = 2:ri
            tR1 = tR1.*du;
        end
        U_val = tR1; return;
    end
else
    U_val = zeros(size(q)); return;
end
end

function U_val = xyplane_pi_rot(ri, cj, Q, lb, N)
sd = (Q(:,2)-1i*Q(:,1))./(-Q(:,2)-1i*Q(:,1));
if (ri + cj == N + 2)
    if (ri == 1 && cj == N+1)
        U_val = lb.^N; return;
    else
        tR1 = lb.^N;
        for d = 2:ri
            tR1 = tR1.*sd;
        end
        U_val = tR1; return;
    end
else
    U_val = zeros(size(lb)); return;
end
end

function U_val = gen_rot(ri, cj, q, Q, a, lb, N)
%%%%%%% General rotation
j = N/2;
aa = q.^2+Q(:,3).^2; bb = Q(:,2).^2+Q(:,1).^2;

if (ri == 1)
    if (cj == 1)
        U_val = a.^N; return;
    else
        mup = j - cj + 1;
        U_val = realsqrt(prod(1:2*j)/(prod(1:(j+mup))*prod(1:(j-mup))))...
            *a.^(j+mup).*lb.^(j-mup);
        return;
    end
end
if (cj == N+1)
    mpr = j - ri + 1;
    U_val = realsqrt(prod(1:2*j)/(prod(1:(j+mpr))*prod(1:(j-mpr))))*...
        ((transpose(a')).^(j-mpr)).*lb.^(j+mpr);
    return;
end
if ((ri ~= 1) || (cj ~= N+1))
    mpr = j - ri + 1;
    mup = j - cj + 1;
    
    rmin = max([0,mup-mpr]);
    rmax = min([j-mpr,j+mup]);
    
    c = j+mup-rmin;
    d = j-mpr-rmin;
    e = mpr-mup+rmin;
    f = rmin;
    ltau = ones(size(aa));
    
    for n = rmax-(rmin:(rmax-1))
%         n
        ltau = 1-((((c-n+1)*(d-n+1)*bb)./((e+n)*(f+n)*aa)).*ltau);
    end
    
%     [ri, cj]
    
    U_val = realsqrt(prod(1:(j+mpr))*prod(1:(j-mpr))*prod(1:(j+mup))*prod(1:(j-mup)))*(-1)^rmin*...
        (a.^c.*transpose(a').^d.*lb.^e.*transpose(lb').^f)...
        /(prod(1:c)*prod(1:d)*prod(1:e)*prod(1:f)).*ltau;
    return;
end
end
