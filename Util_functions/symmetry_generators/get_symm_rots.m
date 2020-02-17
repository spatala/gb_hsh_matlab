function [symm_rots, Laue] = get_symm_rots(g1,g2, pt_grp, data_fname, opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 
%%%% Inputs:
%%%%    1) g1, g2: $3 \times 3$ rotationa matrices of the 
%%%%        $SO(3) \times SO(3)$ parameterization
%%%%    2) pt_grp: string
%%%%        Crystal point group to apply symmetries
%%%%    3) data_fname: string
%%%%        Location of the symmetry elements (matrices)
%%%%    4) opt:
%%%%        1: crystal symmetries only (if Laue (M,n) = (M,-n))
%%%%        2: crystal symmetries + Grain Exchange symmetry
%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Laue = 0;
switch pt_grp
    case 'Oh'
        pt_grp = 'O'; Laue = 1;
    case 'D2h'
        pt_grp = 'D2'; Laue = 1;
end

mat_name = [data_fname, 'SymmMat_', pt_grp, '.mat'];
s1 = load(mat_name);
SymmMat = s1.SymmMat;
nsz = size(SymmMat,1);
nrot_csymm = nsz*nsz;

if opt == 1
    if Laue
        symm_rots = zeros(3,6,2*nrot_csymm);
    else
        symm_rots = zeros(3,6,1*nrot_csymm);
    end
elseif  opt == 2
    if Laue
        symm_rots = zeros(3,6,4*nrot_csymm);
    else
        symm_rots = zeros(3,6,2*nrot_csymm);
    end
end

ct3 = 1;
for ct1=1:nsz
    gs1 = g1*SymmMat{ct1};
    for ct2=1:nsz
        gs2 = g2*SymmMat{ct2};
        symm_rots(:,:,ct3) = [gs1,gs2];
        ct3 = ct3 + 1;
    end
end

ypi = vrrotvec2mat([0,1,0,pi]);
if Laue
    for ct1 = 1:nrot_csymm
        gs1 = ypi*symm_rots(:,1:3,ct1);
        gs2 = ypi*symm_rots(:,4:6,ct1);
        symm_rots(:,:,ct3) = [gs1,gs2];
        ct3 = ct3 + 1;
    end
end

if opt == 2
    if Laue
        nsymm = 2*nrot_csymm;
    else
        nsymm = nrot_csymm;
    end
    for ct1 = 1:nsymm
        gs1 = ypi*symm_rots(:,4:6,ct1);
        gs2 = ypi*symm_rots(:,1:3,ct1);
        symm_rots(:,:,ct3) = [gs1,gs2];
        ct3 = ct3 + 1;
    end
end
end