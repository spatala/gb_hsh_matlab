function diff_val = s3_diff_elem_v1(pts)

n1 = 2;
c = zeros(n1*pts,pts,pts);
b = zeros(n1*pts,pts,pts);
a = zeros(n1*pts,pts,pts);

A = linspace(0,pi,pts);
B = linspace(0,pi,pts);
C = linspace(0,2*pi,n1*pts);

for j = 1:pts
    a(:,:,j) = A(j);
end

for j = 1:pts
    b(:,j,:) = B(j);
end

for j = 1:n1*pts
    c(j,:,:) = C(j);
end

wm3 = compute_Uavg(a);
th3 = compute_Uavg(b);
% ph3 = compute_Uavg(c);

dw1 = A(2)-A(1);
dth1 = B(2)-B(1);
dph1 = C(2)-C(1);
diff_val = (dw1.*(sin(wm3).^2)).*(sin(th3).*dth1).*(dph1);

end