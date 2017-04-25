function D=Di(LH,OM,OL,k)
syms z E0
Z=Z2(LH,OM,OL); 

D=Z.^(-k*LH/E0);
end