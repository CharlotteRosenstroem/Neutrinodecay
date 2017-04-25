function fz=f(a,b,c,B,C,eta,OM,OL,LH)
syms z
LZ=L(OL,OM,LH)

fz=((1+z).^(a*eta)+((1+z)./B).^(b*eta)+((1+z)./C).^(c*eta)).^(1/eta)/(4*pi*LZ)

end 