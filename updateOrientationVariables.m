%% -------------------------------------- MMA UPDATE SCHEME (UNCONSTRAINED)
function znew = updateOrientationVariables(dfdz,z,Iter) 
zMin = -1.0; zMax = 1.0; move = 1/8*0.95^(Iter-1);
AsymInit = 0.5;
albefa = 0.1;
% AsymInc = 1.2; AsymDecr = 0.7;
% Compute asymptotes L and U:
L = z - AsymInit*(zMax-zMin);  U = z + AsymInit*(zMax-zMin);    
% Compute bounds alpha and beta
zzz1 = L + albefa*(z-L);   zzz2 = z-move*(zMax-zMin);  zzz = max(zzz1,zzz2);  alpha = max(zMin,zzz);
zzz1 = U - albefa*(U-z);  zzz2 = z+move*(zMax-zMin);  zzz = min(zzz1,zzz2);    beta = min(zMax,zzz);
% Solve unconstrained subproblem
% feps = 0.00001; 
p = (U-z).^2.*(max(dfdz,0));
q = (z-L).^2.*(-min(dfdz,0));
p = sqrt(p);
q = sqrt(q);
zCnd = (L.*p + U.*q)./(p + q);
znew = max(alpha,min(beta,zCnd));