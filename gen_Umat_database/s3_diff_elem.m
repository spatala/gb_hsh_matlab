function diff_val = s3_diff_elem(num)

wm1 = linspace(0,pi,num);
th1 = linspace(0,pi,num);
ph1 = linspace(0,2*pi,2*num);

[wm2, th2, ph2] = meshgrid(wm1, th1, ph1);

i1 = 1; i3 = 2;
i2 = num-1; i4 = num;
i5 = 2*num-1; i6 = 2*num;

ind0 = i1:i2;
ind1 = i3:i4;
ind2 = i1:i5;
ind3 = i3:i6;

%%% wm is changing in 2nd dimension (:,i,:)
dw1 = wm2(1,2,1) - wm2(1,1,1);
wm3 = compute_grid_avg(wm2);

%%% th2 is changing in 1st dimension (i,:,:)
dth1 = th2(2,1,1) - th2(1,1,1);
th3 = compute_grid_avg(th2);

%%% ph2 is changing in 3rd dimension (:,:,i)
dph1 = ph2(1,1,2) - ph2(1,1,1);

diff_val = (dw1.*(sin(wm3).^2)).*(sin(th3).*dth1).*(dph1);

end
