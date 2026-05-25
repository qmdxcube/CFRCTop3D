%%%%%%%%%% Preparation for filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H,Hs]=PreFILTER3D(nelx,nely,nelz,rmin)
[dy,dx,dz] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
H = max(0,rmin-sqrt(dx.^2+dy.^2+dz.^2));
Hs = convn(ones(nelz,nely,nelx),H,'same');