function C = fillStiffnessLowerTri21(a1, a2, a3, a4, a5, l, m, n)
% fillStiffnessLowerTri21
% Construct the 21 lower-triangular components of a 6x6 transformed
% stiffness matrix in Voigt notation.
%
% Lower-triangular Voigt order:
% [11,21,31,41,51,61,22,32,42,52,62,33,43,53,63,44,54,64,55,65,66]
%
% Inputs:
%   a1-a5 - Invariant stiffness coefficients, each size nElem x 1
%   l,m,n - Direction cosines, each size nElem x 1
%
% Output:
%   C     - Lower-triangular stiffness components, size 21 x nElem

% Ensure all input arrays use row-vector layout for 21 x nElem output assignment.
a1 = a1(:).'; a2 = a2(:).'; a3 = a3(:).'; a4 = a4(:).'; a5 = a5(:).';
l  = l(:).';  m  = m(:).';  n  = n(:).';

nElem = numel(l);
C = zeros(21,nElem);

l2 = l.^2;
m2 = m.^2;
n2 = n.^2;

lm = l.*m;
ln = l.*n;
mn = m.*n;

% Column 1 of the lower triangle: (1,1), (2,1), ..., (6,1)
C(1,:) = a1 + 2*a2 + 2*a3.*l2 + 4*a4.*l2 + a5.*l2.^2;
C(2,:) = a1 + a3.*(l2 + m2) + a5.*l2.*m2;
C(3,:) = a1 + a3.*(l2 + n2) + a5.*l2.*n2;
C(4,:) = a3.*mn + a5.*l2.*mn;
C(5,:) = a3.*ln + 2*a4.*ln + a5.*l2.*ln;
C(6,:) = a3.*lm + 2*a4.*lm + a5.*l2.*lm;

% Column 2 of the lower triangle: (2,2), ..., (6,2)
C(7,:)  = a1 + 2*a2 + 2*a3.*m2 + 4*a4.*m2 + a5.*m2.^2;
C(8,:)  = a1 + a3.*(m2 + n2) + a5.*m2.*n2;
C(9,:)  = a3.*mn + 2*a4.*mn + a5.*m2.*mn;
C(10,:) = a3.*ln + a5.*m2.*ln;
C(11,:) = a3.*lm + 2*a4.*lm + a5.*m2.*lm;

% Column 3 of the lower triangle: (3,3), ..., (6,3)
C(12,:) = a1 + 2*a2 + 2*a3.*n2 + 4*a4.*n2 + a5.*n2.^2;
C(13,:) = a3.*mn + 2*a4.*mn + a5.*n2.*mn;
C(14,:) = a3.*ln + 2*a4.*ln + a5.*n2.*ln;
C(15,:) = a3.*lm + a5.*n2.*lm;

% Column 4 of the lower triangle: (4,4), (5,4), (6,4)
C(16,:) = a2 + a4.*(m2 + n2) + a5.*m2.*n2;
C(17,:) = a4.*lm + a5.*lm.*n2;
C(18,:) = a4.*ln + a5.*ln.*m2;

% Column 5 of the lower triangle: (5,5), (6,5)
C(19,:) = a2 + a4.*(l2 + n2) + a5.*l2.*n2;
C(20,:) = a4.*mn + a5.*mn.*l2;

% Column 6 of the lower triangle: (6,6)
C(21,:) = a2 + a4.*(l2 + m2) + a5.*l2.*m2;

end
