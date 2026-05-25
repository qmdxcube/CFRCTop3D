function [D0,dD0df,S0] = computeConstitutiveMatrix3D(Em,vm,Ef1,Ef2,vf12,vf23,Gf12,f)
% computeConstitutiveMatrix3D
% -------------------------------------------------------------------------
% Computes the local 3D constitutive stiffness matrix D0 of a transversely
% isotropic composite material and its derivative with respect to the fiber
% volume fraction f.
%
% Inputs:
%   Em    - Matrix Young's modulus
%   vm    - Matrix Poisson's ratio
%   Ef1   - Fiber Young's modulus in the longitudinal direction
%   Ef2   - Fiber Young's modulus in the transverse direction
%   vf12  - Fiber major Poisson's ratio
%   vf23  - Fiber transverse Poisson's ratio
%   Gf12  - Fiber longitudinal shear modulus
%   f     - Fiber volume fraction field, size: nelz x nely x nelx
%
% Outputs:
%   D0     - Local stiffness matrix, size: 6 x 6 x nElem
%   dD0df  - Derivative of D0 with respect to f, size: 6 x 6 x nElem
%   S0     - Local compliance matrix, size: 6 x 6 x nElem
%
% Notes:
%   The material is described in its local coordinate system. The first
%   material direction is the fiber direction. The compliance matrix S0 is
%   first assembled from effective engineering constants and then inverted
%   page-wise to obtain the stiffness matrix D0.
% -------------------------------------------------------------------------

% Number of elements in each direction and total number of elements
[nelz,nely,nelx] = size(f);
nElem = nelx*nely*nelz;
% Reshape f into page format so that each element corresponds to one page
% in the 6 x 6 x nElem constitutive matrices.
f = reshape(f,1,1,nElem);
% Allocate compliance matrix and its derivative with respect to f
S0 = zeros(6,6,nElem);
dSdf = zeros(6,6,nElem);
% Matrix shear modulus and fiber transverse shear modulus
Gm = Em/2/(1+vm);
Gf23 = Ef2/2/(1+vf23);
% Effective engineering constants based on the selected mixture rules
% E1:  longitudinal Young's modulus
% E2:  transverse Young's modulus
% G12: longitudinal shear modulus
% G23: transverse shear modulus
E1 = f*Ef1 + (1-f)*Em;
E2 = Em/(1-sqrt(f)*(1-Em/Ef2));
G12 = Gm/(1-sqrt(f)*(1-Gm/Gf12));
G23 = Gm/(1-sqrt(f)*(1-Gm/Gf23));
% Effective Poisson's ratios
v12 = f*vf12 + (1-f)*vm;
v23 = E2/2./G23 - 1;
% Assemble the local compliance matrix S0.
% Voigt notation order: [11, 22, 33, 23, 13, 12].
% The normal strain block is 3 x 3, and the shear terms are diagonal.
S0(1:3,1:3,:) = [ 1./E1,     -v12./E1,  -v12./E1; ...
                  -v12./E1,   1./E2,     -v23./E2; ...
                  -v12./E1,  -v23./E2,    1./E2 ];
S0(4,4,:) = 1./G23;
S0(5,5,:) = 1./G12;
S0(6,6,:) = 1./G12;
% Convert compliance matrix to stiffness matrix for every element/page
D0 = pageinv(S0);
% Derivatives of the effective engineering constants with respect to f
dE1 = Ef1 - Em;
dE2 = (0.5*Em*(1-Em/Ef2).*f.^(-0.5))./(1-sqrt(f)*(1-Em/Ef2)).^2;
dG12 = (0.5*Gm*(1-Gm/Gf12).*f.^(-0.5))./(1-sqrt(f)*(1-Gm/Gf12)).^2;
dG23 = (0.5*Gm*(1-Gm/Gf23).*f.^(-0.5))./(1-sqrt(f)*(1-Gm/Gf23)).^2;
% Derivatives of the effective Poisson's ratios with respect to f
dv12 = vf12 - vm;
dv23 = 0.5*(dE2.*G23 - E2.*dG23)./G23.^2;
% Assemble derivative of the compliance matrix with respect to f.
dSdf(1:3,1:3,:) = [ -1./(E1.^2).*dE1, ...
                    -1./E1.*dv12 + v12./E1.^2.*dE1, ...
                    -1./E1.*dv12 + v12./E1.^2.*dE1; ...
                    -1./E1.*dv12 + v12./E1.^2.*dE1, ...
                    -1./E2.^2.*dE2, ...
                    -1./E2.*dv23 + v23./E2.^2.*dE2; ...
                    -1./E1.*dv12 + v12./E1.^2.*dE1, ...
                    -1./E2.*dv23 + v23./E2.^2.*dE2, ...
                    -1./E2.^2.*dE2 ];
dSdf(4,4,:) = -1./G23.^2.*dG23;
dSdf(5,5,:) = -1./G12.^2.*dG12;
dSdf(6,6,:) = -1./G12.^2.*dG12;
% Derivative of stiffness matrix using D = inv(S):
%   dD/df = -D * (dS/df) * D
dD0df = pagemtimes(-pagemtimes(D0,dSdf),D0);
end
