function [Dcomp, dDlcomp, dDmcomp, dDncomp, dDycomp] = ...
    ComputeElasticityMatrices3DLowerTri21(D0, l, m, n, dD0y)
% ComputeElasticityMatrices3DLowerTri21
% Compute the 21 lower-triangular components of the transformed 3D
% elasticity matrix and, when requested, their derivatives.
%
% Lower-triangular Voigt order:
% [11,21,31,41,51,61,22,32,42,52,62,33,43,53,63,44,54,64,55,65,66]
%
% Usage:
%   [Dcomp,dDlcomp,dDmcomp,dDncomp] = ...
%       ComputeElasticityMatrices3DLowerTri21(D0,l,m,n)
%
%   [Dcomp,dDlcomp,dDmcomp,dDncomp,dDycomp] = ...
%       ComputeElasticityMatrices3DLowerTri21(D0,l,m,n,dD0y)
%
% Inputs:
%   D0    - Local material stiffness matrices, size 6 x 6 x nElem
%   l,m,n - Direction cosines of the fiber direction, each size nElem x 1
%   dD0y  - Optional. Derivative of D0 with respect to fiber content y,
%           size 6 x 6 x nElem. Required only when dDycomp is requested.
%
% Outputs:
%   Dcomp   - Lower-triangular stiffness components, size 21 x nElem
%   dDlcomp - Derivative of Dcomp with respect to l, size 21 x nElem
%   dDmcomp - Derivative of Dcomp with respect to m, size 21 x nElem
%   dDncomp - Derivative of Dcomp with respect to n, size 21 x nElem
%   dDycomp - Optional. Derivative of Dcomp with respect to y,
%             size 21 x nElem

% Ensure direction variables use column-vector layout.
l = l(:);
m = m(:);
n = n(:);

% -------------------------------------------------------------------------
% 1. Extract independent local stiffness components.
%    For the transversely isotropic local material, these five terms are
%    sufficient to build the transformed 3D stiffness matrix.
% -------------------------------------------------------------------------
C11 = squeeze(D0(1,1,:)); C11 = C11(:);
C22 = squeeze(D0(2,2,:)); C22 = C22(:);
C12 = squeeze(D0(1,2,:)); C12 = C12(:);
C23 = squeeze(D0(2,3,:)); C23 = C23(:);
C44 = squeeze(D0(5,5,:)); C44 = C44(:);

% -------------------------------------------------------------------------
% 2. Convert local stiffness components to invariant coefficients.
% -------------------------------------------------------------------------
a1 = C23;
a2 = (C22 - C23) / 2;
a3 = C12 - C23;
a4 = C44 - (C22 - C23) / 2;
a5 = C11 + C22 - 2*C12 - 4*C44;

% -------------------------------------------------------------------------
% 3. Build the 21 lower-triangular stiffness components.
% -------------------------------------------------------------------------
Dcomp = fillStiffnessLowerTri21(a1, a2, a3, a4, a5, l, m, n);

% -------------------------------------------------------------------------
% 4. Directional derivatives with respect to l, m, and n.
%    a1 and a2 are direction-independent, so only a3, a4, and a5 are needed.
% -------------------------------------------------------------------------
[dDlcomp, dDmcomp, dDncomp] = ...
    computeDirectionalDerivsLowerTri21(a3, a4, a5, l, m, n);

% -------------------------------------------------------------------------
% 5. Derivative with respect to fiber content y.
%    This is computed only when the caller requests dDycomp.
% -------------------------------------------------------------------------
if nargout >= 5
    if nargin < 5 || isempty(dD0y)
        error('dD0y is required when dDycomp is requested.');
    end

    dC11y = squeeze(dD0y(1,1,:)); dC11y = dC11y(:);
    dC22y = squeeze(dD0y(2,2,:)); dC22y = dC22y(:);
    dC12y = squeeze(dD0y(1,2,:)); dC12y = dC12y(:);
    dC23y = squeeze(dD0y(2,3,:)); dC23y = dC23y(:);
    dC44y = squeeze(dD0y(5,5,:)); dC44y = dC44y(:);

    da1y = dC23y;
    da2y = (dC22y - dC23y) / 2;
    da3y = dC12y - dC23y;
    da4y = dC44y - (dC22y - dC23y) / 2;
    da5y = dC11y + dC22y - 2*dC12y - 4*dC44y;

    dDycomp = fillStiffnessLowerTri21(da1y, da2y, da3y, da4y, da5y, l, m, n);
end

end
