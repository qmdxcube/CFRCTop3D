function CFRCTop3D(nelx,nely,nelz,he,volfrac,contentfrac,rxmin,rymin,rpmin)
%% Material properties
Em = 1.4e9; vm = 0.394; % Material properties of matrix
Ef1 = 230e9; Ef2 = 15e9; vf12 = 0.2; vf23 = 0.25; Gf12 = 24e9; % Material properties of fiber
ymax = 0.5; % Upper bound for fiber volume fraction
%% Define support and load
nElem = nelx*nely*nelz; % Number of elements
nNode = (nelx+1)*(nely+1)*(nelz+1); % number of nodes
nDof = 3*nNode; % Number of dofs
% USER-DEFINED SUPPORT FIXED DOFs
ifs = 0;[jfs,kfs] = meshgrid(0:nely,0:nelz);% Coordinates
fixednid = ifs*(nely+1)*(nelz+1)+jfs*(nelz+1)+kfs+1; % Fixed Node IDs
fixeddofs = [3*fixednid(:)-2;3*fixednid(:)-1;3*fixednid(:)]; % Fixed DOFs
alldofs = 1:nDof; % All dofs
freedofs = setdiff(alldofs,fixeddofs); % Free DOFs
% USER-DEFINED LOAD DOFs
il=nelx;jl=0:nely;kl = 0;% Coordinates
loadnid = il*(nely+1)*(nelz+1)+jl*(nelz+1)+kl+1; % Loading Node IDs
loaddofs = 3*loadnid(:);% Loading DOFs
F = sparse(loaddofs,1,-20000/(nely+1),nDof,1);
U = zeros(nDof,1);
%% Finite element analysis preparation
nodenrs = int32(reshape(1:nNode,1+nelz,1+nely,1+nelx));
edofVec = reshape(nodenrs(1:end-1,1:end-1,1:end-1),nElem,1);
nids = edofVec+int32([0,(nely+1)*(nelz+1),(nely+2)*(nelz+1),nelz+1,1+[0,(nely+1)*(nelz+1),(nely+2)*(nelz+1),nelz+1]]);
edofMat = 3*repelem(nids,1,3)+repmat(int32([-2,-1,0]),1,8);
[sI, sII] = find(tril(ones(24)));
[iK,jK] = deal(edofMat(:,sI)',edofMat(:,sII)');
numD = 21; % Number of independent stiffness components in a 6x6 symmetric matrix
KET  = zeros(24,24,numD); % Template stiffness matrices (TSMs) in full matrix form
KETc = zeros(300,numD);   % Lower-triangular entries of the TSMs
trilIndex = tril(true(24));
[rowD, colD] = find(tril(ones(6)));
for n = 1:numD
    DT = zeros(6,6);
    i = rowD(n); j = colD(n);
    DT(i,j) = 1; DT(j,i) = 1;
    KEt = H8BasicKe(he,he,he,DT);
    KET(:,:,n) = KEt;
    KETc(:,n) = KEt(trilIndex);
end
%% Preparation for filtering and projection
[Hx,Hxs] = PreFILTER3D(nelx,nely,nelz,rxmin);
[Hy,Hys] = PreFILTER3D(nelx,nely,nelz,rymin);
[Hp,Hps] = PreFILTER3D(nelx,nely,nelz,rpmin);
filt = @(v,H,Hs)convn(v,H,'same')./Hs;% filtering for design variables
ifiltV = @(v,H,Hs) reshape(convn(reshape(v,nelz,nely,nelx)./Hs,H,'same'),[],1); % adjoint operation for sensitivities
prj = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v-eta)))./(tanh(beta*eta)+tanh(beta*(1-eta))); % projection
dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));% proj. x-derivative
%% Initialization of topology, fiber volume fraction and fiber orientation
yinit = contentfrac/volfrac;
x = volfrac*ones(nelz,nely,nelx);
y = yinit*ones(nelz,nely,nelx);
px = ones(nelz,nely,nelx);
py = zeros(nelz,nely,nelx);
pz = zeros(nelz,nely,nelx);
%% Parameters
% Optimization parameters
iter = 0;
change = 1;
maxiter = 2;
tol = 0.01;
% Heaviside projection parameters
loopbeta = 0;
beta = 1;
eta = 0.5;
% MMA parameters for density variables
Nc = 200;
Nd = 2*nElem;
aa0 = 1;
aa = zeros(Nc,1);
cc = 1e4*ones(Nc,1);
dd = ones(Nc,1);
xyval = [x(:);y(:)];
xyold1 = xyval;
xyold2 = xyval;
xymin = [zeros(nElem,1);zeros(nElem,1)];
xymax = [ones(nElem,1);ymax*ones(nElem,1)];
xyL = xymin;
xyU = xymax;
% Treat orientation variables as one coupled vector
pval = [px(:); py(:); pz(:)];
compliance = zeros(maxiter,1);
volumefrac = zeros(maxiter,1);
contentfracHis = zeros(maxiter,1);
%% START ITERATION
while (change>tol && iter<maxiter)
    iter = iter+1;
    loopbeta = loopbeta+1;
    %% FILTERING AND PROJECTION/NORMALIZATION
    xTilde = filt(x,Hx,Hxs); yPhys = filt(y,Hy,Hys);
    xPhys = prj(xTilde,eta,beta);
    pxTilde = filt(px,Hp,Hps); pyTilde = filt(py,Hp,Hps); pzTilde = filt(pz,Hp,Hps);
    Lp = sqrt(pxTilde.^2 + pyTilde.^2 + pzTilde.^2);
    Lp = max(Lp,1e-12); % Avoid division by zero
    pxPhys = pxTilde./Lp; pyPhys = pyTilde./Lp; pzPhys = pzTilde./Lp;
    xPhysV = xPhys(:); yPhysV = yPhys(:);
    %% FE-ANALYSIS
    l1 = pxPhys(:); m1 = pyPhys(:); n1 = pzPhys(:);
    [xmp,dxmp] = MaterialInterpolation(xPhys); % Material interpolation
    xmpV  = xmp(:); dxmpV = dxmp(:);
    [D0,dD0y,~] = computeConstitutiveMatrix3D(Em,vm,Ef1,Ef2,vf12,vf23,Gf12,yPhys); % Compute stiffness matrices in material coordinates
    [Dcomp,dDlcomp,dDmcomp,dDncomp,dDycomp] = ...
        ComputeElasticityMatrices3DLowerTri21(D0,l1,m1,n1,dD0y); % 21 lower-triangular global stiffness components 
    sK = KETc * (Dcomp .* xmpV.'); % Assemble element stiffness entries using the 21 lower-triangular stiffness components directly. 
    K = sparse(iK,jK,sK,nDof,nDof);
    clear sK
    K = K+K' - spdiags(diag(K),0,nDof,nDof);
    U(freedofs,1) = K(freedofs,freedofs)\F(freedofs,1);
    compliance(iter,1) = F'*U;
    volumefrac(iter,1) = mean(xPhysV);
    contentfracHis(iter,1) = mean(xPhysV.*yPhysV);
    %% SENSITIVITY ANALYSIS
    if mod(iter, 20) == 1
        scale = 10/compliance(iter); 
    end
    cex = zeros(nElem,1);
    cey = zeros(nElem,1);
    cepx = zeros(nElem,1);
    cepy = zeros(nElem,1);
    cepz = zeros(nElem,1);
    Ud = U(edofMat);
    for k = 1:numD
        KE = KET(:,:,k);
        uku = sum((Ud*KE).*Ud,2);
        cex  = cex  + Dcomp(k,:).'   .* uku;
        cey  = cey  + dDycomp(k,:).' .* uku;
        cepx = cepx + dDlcomp(k,:).' .* uku;
        cepy = cepy + dDmcomp(k,:).' .* uku;
        cepz = cepz + dDncomp(k,:).' .* uku;
    end
    dcdxPhysV  = -dxmpV .* cex;
    dcdyPhysV  = -xmpV  .* cey;
    dcdpxPhysV = -xmpV  .* cepx;
    dcdpyPhysV = -xmpV  .* cepy;
    dcdpzPhysV = -xmpV  .* cepz;
    %% FILTERING OF SENSITIVITIES
    % Objective function
    dx = dprj(xTilde(:),eta,beta);
    dcdx = ifiltV(dcdxPhysV .* dx, Hx, Hxs);
    dcdy = ifiltV(dcdyPhysV, Hy, Hys);
    G = [dcdpxPhysV, dcdpyPhysV, dcdpzPhysV]; 
    P = [l1,m1,n1]; 
    dcdpTiles = (G-P.*sum(P.*G,2))./Lp(:);
    dcdpx = ifiltV(dcdpTiles(:,1), Hp, Hps);
    dcdpy = ifiltV(dcdpTiles(:,2), Hp, Hps);
    dcdpz = ifiltV(dcdpTiles(:,3), Hp, Hps);
    f0val = scale*compliance(iter);
    df0dxy = scale*[dcdx; dcdy];
    % Constraint functions
    dgvdx = ifiltV(dx/nElem, Hx, Hxs);
    dgvdy = zeros(nElem,1);
    dgfdx = ifiltV(yPhysV.*dx/nElem, Hx, Hxs);
    dgfdy = ifiltV(xPhysV/nElem, Hy, Hys);
    fval = zeros(Nc,1);
    fval(1,1) = volumefrac(iter,1)/volfrac-1;
    fval(2,1) = contentfracHis(iter,1)/contentfrac-1;
    dfvaldxy = zeros(Nc,Nd);
    dfvaldxy(1,:) = [dgvdx.', dgvdy.'] / volfrac;
    dfvaldxy(2,:) = [dgfdx.', dgfdy.'] / contentfrac;
    %% DESIGN VARIABLE UPDATE
    % Update density variables and fiber volume fractions
    [xymma,~,~,~,~,~,~,~,~,xyL,xyU] = mmasub(Nc,Nd,iter,xyval,xymin,xymax,xyold1,xyold2,f0val,df0dxy,fval,dfvaldxy,xyL,xyU,aa0,aa,cc,dd);
    xyold2 = xyold1; xyold1 = xyval; xyval = xymma;
    x = reshape(xymma(1:nElem,1),nelz,nely,nelx); 
    y = reshape(xymma(nElem+(1:nElem),1),nelz,nely,nelx);
    change = max(abs(x(:)-xyold1(1:nElem)));
    % Update fiber orientation variables.
    df0dp = scale*[dcdpx; dcdpy; dcdpz];
    pval = updateOrientationVariables(df0dp,pval,iter);
    px = reshape(pval(1:nElem),nelz,nely,nelx);
    py = reshape(pval(nElem+(1:nElem)),nelz,nely,nelx);
    pz = reshape(pval(2*nElem+(1:nElem)),nelz,nely,nelx);
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f fVf.:%7.3f ch.:%7.3f\n',iter,compliance(iter,1),volumefrac(iter,1),contentfracHis(iter,1),change);
    %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if  beta < 128 && (loopbeta >= 20 || change < tol)
        beta = 2*beta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end
end
compliance = compliance(1:iter,1);
volumefrac = volumefrac(1:iter,1);
contentfracHis = contentfracHis(1:iter,1);
%% VISUALIZATION OF THE OPTIMIZED DESIGN
plotOptimizedDesign(he,he,he,xPhys,yPhys,pxPhys,pyPhys,pzPhys,U)
plotHistory(compliance,volumefrac,contentfracHis)
%% VALIDATION OF THE OPTIMALITY OF FIBER ORIENTATION
[~,~,S0] = computeConstitutiveMatrix3D(Em,vm,Ef1,Ef2,vf12,vf23,Gf12,yPhys);
%% Reconstruct the full 6-by-6-by-nElem stiffness matrix only for validation and saving.
D = zeros(6,6,nElem);
for k = 1:numD
    Dk = reshape(Dcomp(k,:),1,1,nElem);
    D(rowD(k),colD(k),:) = Dk;
    D(colD(k),rowD(k),:) = Dk;
end
verifyOptimality3D(he,he,he,U,D,S0,xPhys,yPhys,pxPhys,pyPhys,pzPhys)
filename = sprintf('CFRCTop3D_%dx%dx%d_vf%.1fcf%.2f.mat',nelx,nely,nelz,volfrac,contentfrac);
save(filename,'compliance','volumefrac','contentfracHis','xPhys','yPhys','pxPhys','pyPhys','pzPhys','U','D','S0','he','timeHis');
