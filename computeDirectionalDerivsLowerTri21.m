function [dDlcomp, dDmcomp, dDncomp] = computeDirectionalDerivsLowerTri21(a3, a4, a5, l, m, n)
% computeDirectionalDerivsLowerTri21
% Compute directional derivatives of the 21 lower-triangular stiffness
% components with respect to the direction cosines l, m, and n.
%
% Lower-triangular Voigt order:
% [11,21,31,41,51,61,22,32,42,52,62,33,43,53,63,44,54,64,55,65,66]
%
% Inputs:
%   a3-a5 - Direction-dependent invariant stiffness coefficients,
%           each size nElem x 1
%   l,m,n - Direction cosines, each size nElem x 1
%
% Outputs:
%   dDlcomp - dDcomp/dl, size 21 x nElem
%   dDmcomp - dDcomp/dm, size 21 x nElem
%   dDncomp - dDcomp/dn, size 21 x nElem
%
% Note:
%   a1 and a2 do not appear in this function because they are independent
%   of the direction cosines l, m, and n.

% Use row-vector layout for direct assignment to 21 x nElem arrays.
a3 = a3(:).';
a4 = a4(:).';
a5 = a5(:).';
l  = l(:).';
m  = m(:).';
n  = n(:).';

nElem = numel(l);
dDlcomp = zeros(21,nElem);
dDmcomp = zeros(21,nElem);
dDncomp = zeros(21,nElem);

l2 = l.^2;
m2 = m.^2;
n2 = n.^2;

lm = l.*m;
ln = l.*n;
mn = m.*n;

l3 = l.^3;
m3 = m.^3;
n3 = n.^3;

% Frequently appearing coefficient combination.
term_a34 = a3 + 2*a4;

% -------------------------------------------------------------------------
% Column 1 of the lower triangle: (1,1), (2,1), ..., (6,1)
% -------------------------------------------------------------------------
% D11 = a1 + 2a2 + (2a3+4a4)l^2 + a5*l^4
dDlcomp(1,:) = 4*term_a34.*l + 4*a5.*l3;

% D21 = D12 = a1 + a3(l^2+m^2) + a5*l^2*m^2
dDlcomp(2,:) = 2*a3.*l + 2*a5.*l.*m2;
dDmcomp(2,:) = 2*a3.*m + 2*a5.*m.*l2;

% D31 = D13 = a1 + a3(l^2+n^2) + a5*l^2*n^2
dDlcomp(3,:) = 2*a3.*l + 2*a5.*l.*n2;
dDncomp(3,:) = 2*a3.*n + 2*a5.*n.*l2;

% D41 = D14 = a3*m*n + a5*l^2*m*n
dDlcomp(4,:) = 2*a5.*l.*mn;
dDmcomp(4,:) = a3.*n + a5.*l2.*n;
dDncomp(4,:) = a3.*m + a5.*l2.*m;

% D51 = D15 = (a3+2a4)*l*n + a5*l^3*n
dDlcomp(5,:) = term_a34.*n + 3*a5.*l2.*n;
dDncomp(5,:) = term_a34.*l + a5.*l3;

% D61 = D16 = (a3+2a4)*l*m + a5*l^3*m
dDlcomp(6,:) = term_a34.*m + 3*a5.*l2.*m;
dDmcomp(6,:) = term_a34.*l + a5.*l3;

% -------------------------------------------------------------------------
% Column 2 of the lower triangle: (2,2), ..., (6,2)
% -------------------------------------------------------------------------
% D22 = a1 + 2a2 + (2a3+4a4)m^2 + a5*m^4
dDmcomp(7,:) = 4*term_a34.*m + 4*a5.*m3;

% D32 = D23 = a1 + a3(m^2+n^2) + a5*m^2*n^2
dDmcomp(8,:) = 2*a3.*m + 2*a5.*m.*n2;
dDncomp(8,:) = 2*a3.*n + 2*a5.*n.*m2;

% D42 = D24 = (a3+2a4)*m*n + a5*m^3*n
dDmcomp(9,:) = term_a34.*n + 3*a5.*m2.*n;
dDncomp(9,:) = term_a34.*m + a5.*m3;

% D52 = D25 = a3*l*n + a5*m^2*l*n
dDlcomp(10,:) = a3.*n + a5.*m2.*n;
dDmcomp(10,:) = 2*a5.*m.*ln;
dDncomp(10,:) = a3.*l + a5.*m2.*l;

% D62 = D26 = (a3+2a4)*l*m + a5*l*m^3
dDlcomp(11,:) = term_a34.*m + a5.*m3;
dDmcomp(11,:) = term_a34.*l + 3*a5.*m2.*l;

% -------------------------------------------------------------------------
% Column 3 of the lower triangle: (3,3), ..., (6,3)
% -------------------------------------------------------------------------
% D33 = a1 + 2a2 + (2a3+4a4)n^2 + a5*n^4
dDncomp(12,:) = 4*term_a34.*n + 4*a5.*n3;

% D43 = D34 = (a3+2a4)*m*n + a5*m*n^3
dDmcomp(13,:) = term_a34.*n + a5.*n3;
dDncomp(13,:) = term_a34.*m + 3*a5.*n2.*m;

% D53 = D35 = (a3+2a4)*l*n + a5*l*n^3
dDlcomp(14,:) = term_a34.*n + a5.*n3;
dDncomp(14,:) = term_a34.*l + 3*a5.*n2.*l;

% D63 = D36 = a3*l*m + a5*n^2*l*m
dDlcomp(15,:) = a3.*m + a5.*n2.*m;
dDmcomp(15,:) = a3.*l + a5.*n2.*l;
dDncomp(15,:) = 2*a5.*n.*lm;

% -------------------------------------------------------------------------
% Column 4 of the lower triangle: (4,4), (5,4), (6,4)
% -------------------------------------------------------------------------
% D44 = a2 + a4(m^2+n^2) + a5*m^2*n^2
dDmcomp(16,:) = 2*a4.*m + 2*a5.*m.*n2;
dDncomp(16,:) = 2*a4.*n + 2*a5.*n.*m2;

% D54 = D45 = a4*l*m + a5*l*m*n^2
dDlcomp(17,:) = a4.*m + a5.*m.*n2;
dDmcomp(17,:) = a4.*l + a5.*l.*n2;
dDncomp(17,:) = 2*a5.*lm.*n;

% D64 = D46 = a4*l*n + a5*l*n*m^2
dDlcomp(18,:) = a4.*n + a5.*n.*m2;
dDmcomp(18,:) = 2*a5.*ln.*m;
dDncomp(18,:) = a4.*l + a5.*l.*m2;

% -------------------------------------------------------------------------
% Column 5 of the lower triangle: (5,5), (6,5)
% -------------------------------------------------------------------------
% D55 = a2 + a4(l^2+n^2) + a5*l^2*n^2
dDlcomp(19,:) = 2*a4.*l + 2*a5.*l.*n2;
dDncomp(19,:) = 2*a4.*n + 2*a5.*n.*l2;

% D65 = D56 = a4*m*n + a5*m*n*l^2
dDlcomp(20,:) = 2*a5.*mn.*l;
dDmcomp(20,:) = a4.*n + a5.*n.*l2;
dDncomp(20,:) = a4.*m + a5.*m.*l2;

% -------------------------------------------------------------------------
% Column 6 of the lower triangle: (6,6)
% -------------------------------------------------------------------------
% D66 = a2 + a4(l^2+m^2) + a5*l^2*m^2
dDlcomp(21,:) = 2*a4.*l + 2*a5.*l.*m2;
dDmcomp(21,:) = 2*a4.*m + 2*a5.*m.*l2;

end
