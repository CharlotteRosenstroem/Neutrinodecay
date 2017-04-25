LH=3.86, OM=0.27, OL=0.73, k=10^(-2), a=3.4, b=-0.3, c=-3.5, eta=-10, B=5000, C=9

syms z E0 

D=Di(LH,OM,OL,k)
fz=f(a,b,c,B,C,eta,OM,OL,LH)


g= @(z,E0) fz*D

t(z,E0)=int(g,E0);

n(z,E0)=int(t,z)

n_zE0=n(6,10000)-n(0,1)


