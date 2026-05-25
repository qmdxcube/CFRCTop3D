function [xp,dxp]=MaterialInterpolation(x)
epsilon=1e-6; penal=3;
xp=epsilon+(1-epsilon)*x.^penal;
dxp=penal*(1-epsilon)*x.^(penal-1);
end