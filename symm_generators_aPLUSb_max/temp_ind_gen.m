clear all; clc;

Nmax = 3;

ab_inds = zeros((Nmax*(Nmax+1)/2),2);
tct3 = 1;
for tct1=0:Nmax
    for tct2=0:tct1
        ab_inds(tct3,:) = [tct2,tct1-tct2];
        tct3 = tct3 + 1;
    end
end

tct1 = 0; tct2 = 0; tct3 = 0;