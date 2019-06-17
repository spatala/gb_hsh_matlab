function mbp = rots_to_mbp(rots)
g1 = rots(:,1:3); g2 = rots(:,4:6);

g1_inv = g1^(-1);
zvec = [0,0,1]';
m1 = g1_inv*g2;
n1 = g1_inv*zvec;
mbp = [m1, n1];
end