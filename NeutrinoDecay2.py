import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math
import cmath


#returns LZ in MPc, default values for OM and OL are taken from wikipedia page for Lambda CDM, and LH from 1208.4600. Make everything consistent later
def LZ(z, LH=3.89e3, OM=0.3089, OL=0.6911):
        def integrand(x):
            return 1./((1+x)*np.sqrt(OM*(1.+x)**3. + OL))
        def integrator(z):
            #The following line is how you evaluate an integral, numerically
            return integrate.quad(integrand, 0., z, limit=50)[0]
        if np.isscalar(z):
            #if z is a single redshift
            return integrator(z)*LH
        else:
            # if z is a large array of redshifts, quickly evaluate for all in the array. 'map' is a functional programming option of python that does that much quicker than a manual loop.
            return np.asarray(map(integrator, z))*LH
        



def Z2(z, OM=0.3089, OL=0.6911):
        def integrand(x):
            return 1./((1+x)**2.*np.sqrt(OM*(1.+x)**3. + OL))
        def integrator(z):
            #The following line is how you evaluate an integral, numerically
            return integrate.quad(integrand, 0., z, limit=50)[0]
        if np.isscalar(z):
            #if z is a single redshift
            return np.exp(integrator(z))
        else:
            # if z is a large array of redshifts, quickly evaluate for all in the array. 'map' is a functional programming option of python that does that much quicker than a manual loop.
            return np.exp(np.asarray(map(integrator, z)))   



def Rc(z, LH=3.89e3, OM=0.3089, OL=0.6911):
        def integrand(x):
            return 1./(np.sqrt(OM*(1.+x)**3. + OL))
        def integrator(z):
            return integrate.quad(integrand, 0., z, limit=50)[0]
        if np.isscalar(z):
            return integrator(z)*LH
        else: 
            return np.asarray(map(integrator,z))*LH   


def RL(z, LH=3.89e3, OM=0.3089, OL=0.6911):
        def integrand(x):
            return 1/(np.sqrt(OM*(1.+x)**3. + OL))
        def integrator(z):
            return integrate.quad(integrand, 0., z, limit=50)[0]
        if np.isscalar(z):
            return integrator(z)*LH*(1+z)
        else: 
            return np.asarray(map(integrator, z))*LH*(1+z)

def dRc(z,LH=3.89e3, OM=0.3089, OL=0.6911):
	return LH/np.sqrt(OM*(1.+z)**3. + OL)


def Di(z, k, E, LH=3.89e3, OM=0.3089, OL=0.6911):
	return Z2(z)**(-k*LH/E)



def f(z, rho=5e-5, a=3.4 , b=-0.3 , c=-3.5 , B=5000 , C=9 , eta=-10 , LH=3.89e3, OM=0.3089, OL=0.6911):
	return rho*((1+z)**(a*eta)+((1+z)/B)**(b*eta)+((1+z)/C)**(c*eta))**(1/eta)



def nv(E,k,L,zmax):
	def integrand(z,k,E):
		return f(z)*Di(z,k,E)*Rc(z)**2*dRc(z)*(L*((1+z)*E)**(-2))/RL(z)**2
	def integrator(E,k=k,zmax=zmax):
		return integrate.quad(integrand,0.,zmax, args=(k,E), limit=50)[0]
	if np.isscalar(E):
		return integrator(E,k,zmax)
	else:
		return np.asarray(map(integrator, E))


def c(theta):
	return math.cos(theta)


def s(theta):
	return math.sin(theta)


theta12=0.5764
theta13=0.1468
theta23=0.7222
delta = math.pi*1.35

def v1(ve,vm,vt):
	return ve*c(theta13)*c(theta12)-vm*(c(theta23)*s(theta12)+s(theta23)*c(theta12)*s(theta13)*cmath.exp(-1j*delta))+vt*(s(theta23)*s(theta12)-c(theta23)*c(theta12)*s(theta13)*cmath.exp(-1j*delta)) 

def v2(ve,vm,vt):
	return ve*c(theta13)*s(theta12)+vm*(c(theta23)*c(theta12)-s(theta23)*s(theta12)*s(theta13)*cmath.exp(-1j*delta))-vt*(s(theta23)*s(theta12)+c(theta23)*c(theta12)*s(theta13)*cmath.exp(-1j*delta))#here there was a bug

def v3(ve,vm,vt):
	return ve*s(theta13)*cmath.exp(1j*delta)+vm*c(theta13)*s(theta23)+vt*c(theta13)*c(theta23)


m1=abs(v1(1.,1.,1.))
m2=abs(v2(1.,1.,1.))
m3=abs(v3(1.,1.,1.))

print m1, m2, m3

Ls=3e41


def v1E(E):
	return m1*nv(E,1e-17,Ls,3.)
	

def v2E(E,k):
	return m2*nv(E,k,Ls,3.)
	
def v3E(E,k):
	return m3*nv(E,k,Ls,3.)





def veE(E,k):	
	return abs(v1E(E)*c(theta13)*c(theta12)+v2E(E,k)*c(theta13)*s(theta12)+v3E(E,k)*s(theta13)*cmath.exp(-1j*delta))


def vmE(E,k):
	return abs(-v1E(E)*(c(theta23)*s(theta12)+s(theta23)*c(theta12)*s(theta13)*cmath.exp(1j*delta))+v2E(E,k)*(c(theta23)*c(theta12)-s(theta23)*s(theta12)*s(theta13)*cmath.exp(1j*delta))+v3E(E,k)*c(theta13)*s(theta23))


def vtE(E,k):
	return abs(v1E(E)*(s(theta23)*s(theta12)-c(theta23)*c(theta12)*s(theta13)*cmath.exp(1j*delta))-v2E(E,k)*(s(theta23)*c(theta12)+c(theta23)*s(theta12)*s(theta13)*cmath.exp(1j*delta))+v3E(E,k)*c(theta13)*c(theta23))


def veE1(k,Emin,Emax):
	def integrand(E,k):
		return veE(E,k)
	def integrator(k,Emin=Emin,Emax=Emax):
		return integrate.quad(integrand,Emin,Emax, args=(k), limit=50)[0]
	if np.isscalar(k):
		return integrator(k,Emin,Emax)
	else:
		return np.asarray(map(integrator, k))

def vmE1(k,Emin,Emax):
	def integrand(E,k):
		return vmE(E,k)
	def integrator(k,Emin=Emin,Emax=Emax):
		return integrate.quad(integrand,Emin,Emax, args=(k), limit=50)[0]
	if np.isscalar(k):
		return integrator(k,Emin,Emax)
	else:
		return np.asarray(map(integrator, k))

def vtE1(k,Emin,Emax):
	def integrand(E,k):
		return vtE(E,k)
	def integrator(k,Emin=Emin,Emax=Emax):
		return integrate.quad(integrand,Emin,Emax, args=(k), limit=50)[0]
	if np.isscalar(k):
		return integrator(k,Emin,Emax)
	else:
		return np.asarray(map(integrator, k))



ks=np.power(10.,np.linspace(-5.,5.,100))

vme=vmE1(ks,10.,10000.)/veE1(ks,10.,10000.)
vte=vtE1(ks,10.,10000.)/veE1(ks,10.,10000.)


ks1=np.power(ks,-1.)

plt.figure(1)
R_1,=plt.plot(ks1,vme,color='red',label='vm/ve')
R_2,=plt.plot(ks1,vte,color='green', label='vt/ve')
plt.xscale('log')
plt.xlabel('$k^{-1}[s/eV]$', fontsize=15)
plt.ylabel('Flavor ratio', fontsize=15)
plt.legend((R_1, R_2,),('vm/ve','vt/ve'))
plt.title('v1 stable, v2 and v3 lifetime: $k^{-1}$')
plt.savefig('flavorratio.png')
plt.show()

#
"""
ks1=np.power(ks,-1.)

veEE=veE(ks)
vmEE=vmE(ks)
vtEE=vtE(ks)

plt.figure(1)
plt.plot(ks1,veEE,color='red')
plt.plot(ks1,vmEE,color='green')
plt.plot(ks1,vtEE,color='blue')
plt.xscale('log')
plt.yscale('log')
plt.show()
"""
