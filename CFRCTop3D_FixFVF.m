function CFRCTop3D_FixFVF(nelx,nely,nelz,he,volfrac,fvf,rxmin,rpmin)
% Fixed fiber volume fraction version.
%% Material properties
Em = 1.4e9; vm = 0.394;                                      % Matrix material properties
Ef1 = 230e9; Ef2 = 15e9; vf12 = 0.2; vf23 = 0.25; Gf12 = 24e9; % Fiber material properties
[D0,~,S0] = computeConstitutiveMatrix3D(Em,vm,Ef1,Ef2,vf12,vf23,Gf12,fvf); % Local stiffness for fixed fvf
%% Define support and load
nElem = nelx*nely*nelz;             % Number of elements
nNode = (nelx+1)*(nely+1)*(nelz+1); % Number of nodes
nDof  = 3*nNode;                    % Number of degrees of freedom
% USER-DEFINED SUPPORT FIXED DOFs
ifs = 0;
[jfs,kfs] = meshgrid(0:nely,0:nelz);
fixednid = ifs*(nely+1)*(nelz+1) + jfs*(nelz+1) + kfs + 1;
fixeddofs = [3*fixednid(:)-2; 3*fixednid(:)-1; 3*fixednid(:)];
alldofs = 1:nDof;
freedofs = setdiff(alldofs,fixeddofs);
% USER-DEFINED LOAD DOFs
il = nelx; jl = 0:nely; kl = 0;
loadnid = il*(nely+1)*(nelz+1) + jl*(nelz+1) + kl + 1;
loaddofs = 3*loadnid(:);
F = sparse(loaddofs,1,-20000/(nely+1),nDof,1);
U = zeros(nDof,1);
%% Finite element analysis preparation
nodenrs = int32(reshape(1:nNode,1+nelz,1+nely,1+nelx));
edofVec = reshape(nodenrs(1:end-1,1:end-1,1:end-1),nElem,1);
nids = edofVec + int32([0,(nely+1)*(nelz+1),(nely+2)*(nelz+1),nelz+1, ...
                        1+[0,(nely+1)*(nelz+1),(nely+2)*(nelz+1),nelz+1]]);
edofMat = 3*repelem(nids,1,3) + repmat(int32([-2,-1,0]),1,8);
[sI,sII] = find(tril(ones(24)));
[iK,jK] = deal(edofMat(:,sI)',edofMat(:,sII)');
numD = 21; % Number of independent stiffness components in a 6x6 symmetric matrix
KET  = zeros(24,24,numD); % Template stiffness matrices (TSMs) in full matrix form
KETc = zeros(300,numD);   % Lower-triangular entries of the TSMs
trilIndex = tril(true(24));
[rowD,colD] = find(tril(ones(6)));
for k = 1:numD
    DT = zeros(6,6);
    i = rowD(k); j = colD(k);
    DT(i,j) = 1;
    DT(j,i) = 1;
    KEt = H8BasicKe(he,he,he,DT);
    KET(:,:,k) = KEt;
    KETc(:,k) = KEt(trilIndex);
end
%% Preparation for filtering and projection
[Hx,Hxs] = PreFILTER3D(nelx,nely,nelz,rxmin);
[Hp,Hps] = PreFILTER3D(nelx,nely,nelz,rpmin);
filt   = @(v,H,Hs) convn(v,H,'same')./Hs; % Filter for design variables, input/output are 3D arrays
ifiltV = @(v,H,Hs) reshape(convn(reshape(v,nelz,nely,nelx)./Hs,H,'same'),[],1); % Sensitivity filter, input/output are vectors
prj  = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v-eta))) ./ ...
                      (tanh(beta*eta)+tanh(beta*(1-eta)));
dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2) ./ ...
                      (tanh(beta*eta)+tanh(beta*(1-eta)));
%% Initial designs for topology and fiber orientation
x  = volfrac*ones(nelz,nely,nelx);
px = ones(nelz,nely,nelx);
py = zeros(nelz,nely,nelx);
pz = zeros(nelz,nely,nelx);
%% Parameters
% Optimization parameters
iter = 0;
change = 1;
maxiter = 200;
tol = 0.01;
% Heaviside projection parameters
loopbeta = 0;
beta = 1;
eta = 0.5;
% MMA parameters for density variables
Nc = 1;
Nd = nElem;
aa0 = 1;
aa = zeros(Nc,1);
cc = 1e4*ones(Nc,1);
dd = ones(Nc,1);
xval = x(:);
xold1 = xval;
xold2 = xval;
xmin = zeros(nElem,1);
xmax = ones(nElem,1);
xL = xmin;
xU = xmax;
% Treat orientation variables as one coupled orientation vector
pval = [px(:); py(:); pz(:)];
compliance = zeros(maxiter,1);
volumefrac = zeros(maxiter,1);
contentfracHis = zeros(maxiter,1);
%% Start iteration
while (change > tol && iter < maxiter)
    iter = iter + 1;
    loopbeta = loopbeta + 1;
    %% Filtering, projection and orientation normalization
    xTilde = filt(x,Hx,Hxs);
    xPhys = prj(xTilde,eta,beta);
    pxTilde = filt(px,Hp,Hps);
    pyTilde = filt(py,Hp,Hps);
    pzTilde = filt(pz,Hp,Hps);
    Lp = sqrt(pxTilde.^2 + pyTilde.^2 + pzTilde.^2);
    Lp = max(Lp,1e-12); % Avoid division by zero during orientation normalization
    pxPhys = pxTilde ./ Lp;
    pyPhys = pyTilde ./ Lp;
    pzPhys = pzTilde ./ Lp;
    xPhysV = xPhys(:);
    yPhys = fvf*ones(size(xPhys));
    yPhysV = yPhys(:);
    %% FE analysis
    l1 = pxPhys(:); m1 = pyPhys(:); n1 = pzPhys(:);
    [xmp,dxmp] = MaterialInterpolation(xPhys);
    xmpV  = xmp(:); dxmpV = dxmp(:);
    % Compute 21 lower-triangular stiffness components directly.
    [Dcomp,dDlcomp,dDmcomp,dDncomp] = ...
        ComputeElasticityMatrices3DLowerTri21(D0,l1,m1,n1);
    % Assemble element stiffness entries using 21 lower-triangular components.
    sK = KETc * (Dcomp .* xmpV.');
    K = sparse(iK,jK,sK,nDof,nDof);
    clear sK
    K = K + K' - spdiags(diag(K),0,nDof,nDof);
    U(freedofs,1) = K(freedofs,freedofs) \ F(freedofs,1);
    compliance(iter,1) = F' * U;
    volumefrac(iter,1) = mean(xPhysV);
    contentfracHis(iter,1) = mean(xPhysV .* yPhysV);
    %% Sensitivity analysis
    if mod(iter,20) == 1
        scale = 10/compliance(iter);
    end
    cex  = zeros(nElem,1);
    cepx = zeros(nElem,1);
    cepy = zeros(nElem,1);
    cepz = zeros(nElem,1);
    Ud = U(edofMat);
    for k = 1:numD
        KE = KET(:,:,k);
        uku = sum((Ud*KE).*Ud,2);
        cex  = cex  + Dcomp(k,:).'   .* uku;
        cepx = cepx + dDlcomp(k,:).' .* uku;
        cepy = cepy + dDmcomp(k,:).' .* uku;
        cepz = cepz + dDncomp(k,:).' .* uku;
    end
    dcdxPhysV  = -dxmpV .* cex;
    dcdpxPhysV = -xmpV  .* cepx;
    dcdpyPhysV = -xmpV  .* cepy;
    dcdpzPhysV = -xmpV  .* cepz;
    %% Filtering of sensitivities
    dx = dprj(xTilde(:),eta,beta);
    % Objective function sensitivities
    dcdx = ifiltV(dcdxPhysV .* dx, Hx, Hxs);
    G = [dcdpxPhysV, dcdpyPhysV, dcdpzPhysV];
    P = [l1, m1, n1];
    dcdpTiles = (G - P.*sum(P.*G,2)) ./ Lp(:);
    dcdpx = ifiltV(dcdpTiles(:,1), Hp, Hps);
    dcdpy = ifiltV(dcdpTiles(:,2), Hp, Hps);
    dcdpz = ifiltV(dcdpTiles(:,3), Hp, Hps);
    f0val = scale*compliance(iter);
    df0dx = scale*dcdx;
    % Volume constraint
    dgvdx = ifiltV(dx/nElem, Hx, Hxs);
    fval(1,1) = volumefrac(iter,1)/volfrac - 1;
    dfvaldx(1,:) = dgvdx.'/volfrac;
    %% Design variable update
    % Update density variables by MMA.
    [xmma,~,~,~,~,~,~,~,~,xL,xU] = mmasub(Nc,Nd,iter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfvaldx,xL,xU,aa0,aa,cc,dd);
    xold2 = xold1;
    xold1 = xval;
    xval = xmma;
    x = reshape(xval,nelz,nely,nelx);
    change = max(abs(xval - xold1));
    % Update fiber orientation variables.
    df0dp = scale*[dcdpx; dcdpy; dcdpz];
    pval = updateOrientationVariables(df0dp,pval,iter);
    px = reshape(pval(1:nElem),nelz,nely,nelx);
    py = reshape(pval(nElem+(1:nElem)),nelz,nely,nelx);
    pz = reshape(pval(2*nElem+(1:nElem)),nelz,nely,nelx);
    %% Print results
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f FVF.:%7.3f ch.:%7.3f\n', ...
        iter,compliance(iter,1),volumefrac(iter,1),fvf,change);
    %% Update Heaviside regularization parameter
    if beta < 128 && (loopbeta >= 20 || change < tol)
        beta = 2*beta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end
end
compliance = compliance(1:iter,1);
volumefrac = volumefrac(1:iter,1);
contentfracHis = contentfracHis(1:iter,1);
%% Visualization of the optimized design
yPhys = fvf*ones(size(xPhys));
plotOptimizedDesign(he,he,he,xPhys,yPhys,pxPhys,pyPhys,pzPhys,U)
plotHistory(compliance,volumefrac,contentfracHis)
%% Validation of the optimality of fiber orientation
% Reconstruct the full 6-by-6-by-nElem stiffness matrix only for validation and saving.
D = zeros(6,6,nElem);
for k = 1:numD
    Dk = reshape(Dcomp(k,:),1,1,nElem);
    D(rowD(k),colD(k),:) = Dk;
    D(colD(k),rowD(k),:) = Dk;
end
verifyOptimality3D(he,he,he,U,D,S0,xPhys,yPhys,pxPhys,pyPhys,pzPhys)
filename = sprintf('CFRCTop3D_FixFVF_%dx%dx%d_vf%.1ffvf%.2f.mat',nelx,nely,nelz,volfrac,fvf);
save(filename,'compliance','volumefrac','contentfracHis','xPhys','yPhys','pxPhys','pyPhys','pzPhys','U','D','S0','he');
end
