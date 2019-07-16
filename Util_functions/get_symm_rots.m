function [symm_rots, Laue] = get_symm_rots(g1,g2, pt_grp, data_fname)
Laue = 0;
switch pt_grp
    case 'Oh'
        pt_grp = 'O';
        Laue = 1;
    case 'D2h'
        pt_grp = 'D2';
        Laue = 1;
end
% if (strcmp(pt_grp, 'Oh'))
%     Laue = 1; pt_grp = 'Oh';
% else
%     Laue = 0;
% end


mat_name = [data_fname, 'SymmMat_', pt_grp, '.mat'];
s1 = load(mat_name);
SymmMat = s1.SymmMat;
nsz = size(SymmMat,1);
if Laue
    symm_rots = zeros(3,6,4*nsz*nsz);
else
    symm_rots = zeros(3,6,2*nsz*nsz);
end
ct3 = 1;
for ct1=1:nsz
    gs1 = g1*SymmMat{ct1};
    for ct2=1:nsz
        gs2 = g2*SymmMat{ct1};
        symm_rots(:,:,ct3) = [gs1,gs2];
        ct3 = ct3 + 1;
    end
end

ypi = vrrotvec2mat([0,1,0,pi]);

g1p = ypi*g2; g2p = ypi*g1;
for ct1=1:nsz
    gs1 = g1p*SymmMat{ct1};
    for ct2=1:nsz
        gs2 = g2p*SymmMat{ct1};
        symm_rots(:,:,ct3) = [gs1,gs2];
        ct3 = ct3 + 1;
    end
end

if Laue
    g1p = ypi*g1; g2p = ypi*g2;
    for ct1=1:nsz
        gs1 = g1p*SymmMat{ct1};
        for ct2=1:nsz
            gs2 = g2p*SymmMat{ct1};
            symm_rots(:,:,ct3) = [gs1,gs2];
            ct3 = ct3 + 1;
        end
    end
    g1p = g2; g2p = g1;
    for ct1=1:nsz
        gs1 = g1p*SymmMat{ct1};
        for ct2=1:nsz
            gs2 = g2p*SymmMat{ct1};
            symm_rots(:,:,ct3) = [gs1,gs2];
            ct3 = ct3 + 1;
        end
    end
end
end