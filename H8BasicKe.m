function [KE] = H8BasicKe(lx,ly,lz,D)
% half-length of each edge
a=lx/2;b=ly/2;c=lz/2;
% Two Gauss points in all directions
xx=[-1/sqrt(3), 1/sqrt(3)]; yy = xx;zz=xx;
ww=[1,1];
% Initialize
KE = zeros(24,24);
B=zeros(6,24);
detJ=a*b*c;
xic=  [-1 1 1 -1 -1 1 1 -1];
etac= [-1 -1 1 1 -1 -1 1 1];
zetac=[-1 -1 -1 -1 1 1 1 1];
for ii=1:length(xx)
    for jj=1:length(yy)
        for kk=1:length(zz)
            xi = xx(ii);eta = yy(jj);zeta = zz(kk);
            dNx=xic.*(1+eta*etac).*(1+zeta*zetac)/8/a;
            dNy=(1+xi*xic).*etac.*(1+zeta*zetac)/8/b;
            dNz=(1+xi*xic).*(1+eta*etac).*zetac/8/c;
            B(1,1:3:24)=dNx;B(2,2:3:24)=dNy;B(3,3:3:24)=dNz;
            B(4,2:3:24)=dNz;B(4,3:3:24)=dNy;
            B(5,1:3:24)=dNz;B(5,3:3:24)=dNx;
            B(6,1:3:24)=dNy;B(6,2:3:24)=dNx;
            weight = ww(ii)*ww(jj)*ww(kk)*detJ;
            KE = KE + weight*(B' * D * B);
        end
    end
end