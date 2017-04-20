function fz=f(rho0,z,a,b,c,B,C,eta)

fz=rho0.*((1+z).^(a*eta)+((1+z)./B).^(b*eta)+((1+z)./C).^(c*eta)).^(1/eta)

end 