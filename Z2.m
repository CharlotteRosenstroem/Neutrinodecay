function Z=Z2(LH,OM,OL)
syms x z 

f=1./((1+x)^2*sqrt(OM*(1+x)^3+OL));

g(x)=int(f,x); 
Z=exp(g(z)-g(0))
end