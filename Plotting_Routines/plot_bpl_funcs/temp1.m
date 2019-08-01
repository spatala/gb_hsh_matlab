symm_num = 14;
symm_Mvec = Mvecs*S(:,symm_num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[1200,50,1200,1200]); hold on;
fval = abs(symm_Mvec);
fval = reshape(fval, [num+1, num+1]);
surf(X,Y,Z, fval);
shading interp;
axis equal; axis off;
% view([1,0,0])
view([1,1,1]);
colorbar;

