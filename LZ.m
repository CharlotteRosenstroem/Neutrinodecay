function L=LZ(LH,OM,OL,z)
L=[];

for i=1:length(z);
    z=i*0.001;
f= @(x) 1./((1+x).*sqrt(OM*(1+x).^3+OL));

q=integral(f,0,z);

L(i)=LH*q;
end
