function D=Di(LH,OM,OL,z,E0,k)
Z2=[]
for i=1:length(z);
z=i*0.001;
f= @(x) 1./((1+x).^2.*sqrt(OM.*(1+x).^3+OL));

q=integral(f,0,z);

Z2(i)=exp(q);

end
D=Z2.^(-k*LH/E0);
end