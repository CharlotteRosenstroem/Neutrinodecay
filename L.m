function LZ=L(OL,OM,LH)
syms z 

f=(2.*sqrt(OM+OL).*atanh(sqrt(1+OM/OL)))./(3.*sqrt(1+OM/OL).*OL)-(2.*sqrt(OL+OM.*(1+z).^3).*atanh(sqrt(1+(OM.*(1+z).^3)./OL)))./(3.*OL.*sqrt(1+(OM.*(1+z).^3)./OL));
LZ=f.*LH;

end

