function verifyOptimality3D(hx,hy,hz,U,D,S0,xPhys,yPhys,pxPhys,pyPhys,pzPhys)
d = S0(1,1,:) - 2*S0(1,2,:)+ S0(2,2,:) - S0(6,6,:);
y = yPhys(:);
d(y<1e-3)=0;
if any(d>0)
    warning('This program is exclusively designed for the analysis of weak shear materials.');
    return;
end
%% Compute collinearity metric between the dominant principal direction and fiber direction
[nelz,nely,nelx] = size(xPhys);
nElem = nelx*nely*nelz; % Number of elements
nNode = (nelx+1)*(nely+1)*(nelz+1); % number of nodes
nodenrs = int32(reshape(1:nNode,1+nelz,1+nely,1+nelx));
edofVec = reshape(nodenrs(1:end-1,1:end-1,1:end-1),nElem,1);
nids = edofVec+int32([0,(nely+1)*(nelz+1),(nely+2)*(nelz+1),nelz+1,1+[0,(nely+1)*(nelz+1),(nely+2)*(nelz+1),nelz+1]]);
edofMat = 3*repelem(nids,1,3)+repmat(int32([-2,-1,0]),1,8);
B=zeros(6,24);
a = hx/2; b = hy/2; c = hz/2;
dNx = [-1 1 1 -1 -1 1 1 -1]/8/a;
dNy = [-1 -1 1 1 -1 -1 1 1]/8/b;
dNz = [-1 -1 -1 -1 1 1 1 1]/8/c;
B(1,1:3:24)=dNx;B(2,2:3:24)=dNy;B(3,3:3:24)=dNz;%delta1, delta2, delta3
B(4,2:3:24)=dNz;B(4,3:3:24)=dNy;%taoyz
B(5,1:3:24)=dNz;B(5,3:3:24)=dNx;%taoxz
B(6,1:3:24)=dNy;B(6,2:3:24)=dNx;%taoxy
Ud = U(edofMat');
strain = B*Ud;
strain_3d = reshape(strain, 6, 1, []);
stress_3d = pagemtimes(D, strain_3d);
stress = squeeze(stress_3d);
collinearity = zeros(nElem,1);
l1 = pxPhys(:); m1 = pyPhys(:); n1 = pzPhys(:);
for e = 1:nElem
    sigma = [stress(1,e), stress(6,e), stress(5,e); % 3x3 stress tensor
        stress(6,e), stress(2,e), stress(4,e);
        stress(5,e), stress(4,e), stress(3,e)];
    [pVectors, pStress] = eig(sigma); % Compute principal stresses and directions
    pStress = diag(pStress);
    [~, idx] = sort(abs(pStress), 'descend');
    pVectors_sorted = pVectors(:, idx);
    dPDirection = pVectors_sorted(:, 1);% Find the dominant principal direction
    fiberDirection = [l1(e); m1(e); n1(e)]; % Fiber direction
    collinearity(e, 1) = abs(dot(dPDirection, fiberDirection));
end
collinearity = reshape(collinearity,nelz,nely,nelx);

xthreshold = 0.5; ythreshold = 0.01;
activeElems = (xPhys >= xthreshold) & (yPhys > ythreshold);
faceMask = repelem(activeElems(:), 6);
[Z, Y, X] = ndgrid(0:nelz, 0:nely, 0:nelx);
coords = [X(:)*hx, Y(:)*hy, Z(:)*hz];
facePattern = [1 4 8 5; 4 8 7 3; 3 7 6 2; 2 6 5 1; 1 4 3 2; 8 5 6 7];
elementFace = reshape(nids(:, facePattern')', 4, [])';

figure('Color', 'w', 'Name', 'Colinearity');
patch('Faces', elementFace(faceMask, :), 'Vertices', coords, ...
    'FaceVertexCData', repelem(collinearity(activeElems), 6), ...
    'FaceColor', 'flat', 'EdgeAlpha', 0.1, 'EdgeColor', 'k', 'LineWidth', 0.1);

setup_3d_style('Colinearity', [0, 1], jet);
hold off;