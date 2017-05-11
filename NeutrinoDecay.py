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
            return integrate.quad(integrand, 0., z)[0]
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
            return integrate.quad(integrand, 0., z)[0]
        if np.isscalar(z):
            #if z is a single redshift
            return np.exp(integrator(z))
        else:
            # if z is a large array of redshifts, quickly evaluate for all in the array. 'map' is a functional programming option of python that does that much quicker than a manual loop.
            return np.exp(np.asarray(map(integrator, z)))   



def Di(z, k, E, LH=3.89e3, OM=0.3089, OL=0.6911):
	return Z2(z)**(-k*LH/E)



def f(z, a=3.4 , b=-0.3 , c=-3.5 , B=5000 , C=9 , eta=-10 , LH=3.89e3, OM=0.3089, OL=0.6911):
	return ((1+z)**(a*eta)+((1+z)/B)**(b*eta)+((1+z)/C)**(c*eta))**(1/eta)/(4*math.pi*LZ(z))



#
"""
def nv(E,k,zmax):
	def integrand(z,k,E):
		return f(z)*Di(z,k,E)*(1+z)**(-2)*E**(-2)
	def integrator(E,k=k,zmax=zmax):
		return integrate.quad(integrand,0,zmax, args=(k,E))[0]
	if np.isscalar(E):
		return integrator(E,k,zmax)
	else:
		return np.asarray(map(integrator, E))
		
Es=np.power(10., np.linspace(0.,5.,100))
"""


def nv(k,Emin,Emax,zmax):
	def integrand(z,k,E):
		return f(z)*Di(z,k,E)*(1+z)**(-2)*E**(-2)
	def integrator(k,E,zmax=zmax):
		return integrate.quad(integrand,0,zmax, args=(k,E))[0]
	def integratorE(k,Emin=Emin,Emax=Emax,zmax=zmax):
		return integrate.quad(integrator,Emin, Emax, args=(k,zmax))[0]	
	if np.isscalar(k):
		return integratorE(k,Emin,Emax,zmax)
	else:
		return np.asarray(map(integratorE, k))

print nv(1000,10.,10000.,3.)

def c(theta):
	return math.cos(theta)


def s(theta):
	return math.sin(theta)


theta12=0.5764
theta13=0.1463
theta23=0.7223

def v1(ve,vm,vt):
	return ve*c(theta13)*c(theta12)-vm*abs(c(theta23)*s(theta12)+s(theta23)*c(theta12)*s(theta13)*cmath.exp(1j*1.35*math.pi))+vt*abs(s(theta23)*s(theta12)-c(theta23)*c(theta12)*s(theta13)*cmath.exp(1j*1.35*math.pi)) 

def v2(ve,vm,vt):
	return ve*c(theta13)*s(theta12)+vm*abs(c(theta23)*c(theta12)-s(theta23)*s(theta12)*s(theta13)*cmath.exp(1j*1.35*math.pi))-vt*abs(s(theta23)*c(theta12)+c(theta23)*c(theta12)*s(theta13)*cmath.exp(1j*1.35*math.pi))

def v3(ve,vm,vt):
	return -ve*abs(s(theta13)*cmath.exp(-1j*1.35*math.pi))+vm*c(theta13)*s(theta23)+vt*c(theta13)*c(theta23)


m1=v1(1.,2.,0.)
m2=v2(1.,2.,0.)
m3=v3(1.,2.,0.)

#print m1, m2, m3


nvs=nv(1e-17,10.,10000.,3.)

v1E=m1*nvs

def v2E(k):
	return m2*nv(k,10.,10000.,3.)

def v3E(k):
	return m3*nv(k,10.,10000.,3.)

def veE(k):
	return v1E*c(theta13)*c(theta12)+v2E(k)*c(theta13)*s(theta12)-v3E(k)*abs(s(theta13)*cmath.exp(-1j*1.35*math.pi))


def vmE(k):
	return -v1E*abs(c(theta23)*s(theta12)+s(theta23)*c(theta12)*s(theta13)*cmath.exp(1j*1.35*math.pi))+v2E(k)*abs(c(theta23)*c(theta12)-s(theta23)*s(theta12)*s(theta13)*cmath.exp(1j*1.35*math.pi))+v3E(k)*c(theta13)*s(theta23)


def vtE(k):
	return v1E*abs(s(theta23)*s(theta12)-c(theta23)*c(theta12)*s(theta13)*cmath.exp(1j*1.35*math.pi))-v2E(k)*abs(s(theta23)*c(theta12)+c(theta23)*s(theta12)*s(theta13)*cmath.exp(1j*1.35*math.pi))+v3E(k)*c(theta13)*c(theta23)

ks=np.power(10.,np.linspace(-4.,4.,10))

veE1=veE(ks)
vmE1=vmE(ks)
vtE1=vtE(ks)

ks1=np.power(ks,-1.)

plt.figure(1)
R_1,=plt.plot(ks1,veE1,color='red',label='ve')
R_2,=plt.plot(ks1,vmE1,color='green', label='vm')
R_3,=plt.plot(ks1,vtE1,color='blue', label='vt')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$k^{-1}[s/eV]$', fontsize=15)
plt.ylabel('$Flux[s^{-1}Mpc^{-2}TeV^{-1}]$', fontsize=15)
plt.legend((R_1, R_2, R_3),('ve','vm','vt'))
plt.title('v1 stable, v2 and v3 lifetime: k^{-1}')
plt.show()



#
"""
veE1=veE(1)
vmE1=vmE(1)
vtE1=vtE(1)

plt.figure(1)
v_e,=plt.plot(Es,veE1,color='red',label='ve')
v_m,=plt.plot(Es,vmE1,color='green', label='vm')
v_t,=plt.plot(Es,vtE1,color='blue', label='vt')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$E[TeV]$', fontsize=15)
plt.ylabel('$Flux[s^{-1}Mpc^{-2}TeV^{-1}]$', fontsize=15)
plt.legend((v_e, v_m, v_t),('ve','vm','vt'))
plt.title('v1 stable, v2 and v3 lifetime:1 s/ev')
plt.savefig('k1.png')

veE2=veE(1e-2)
vmE2=vmE(1e-2)
vtE2=vtE(1e-2)

plt.figure(2)
v_e2,=plt.plot(Es,veE2,color='red',label='ve')
v_m2,=plt.plot(Es,vmE2,color='green', label='vm')
v_t2,=plt.plot(Es,vtE2,color='blue', label='vt')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$E[TeV]$', fontsize=15)
plt.ylabel('$Flux[s^{-1}Mpc^{-2}TeV^{-1}]$', fontsize=15)
plt.legend((v_e2, v_m2, v_t2),('ve','vm','vt'))
plt.title('v1 stable, v2 and v3 lifetime: 100 s/eV')
plt.savefig('k100.png')

veE3=veE(1e-4)
vmE3=vmE(1e-4)
vtE3=vtE(1e-4)

plt.figure(3)
v_e3,=plt.plot(Es,veE3,color='red',label='ve')
v_m3,=plt.plot(Es,vmE3,color='green', label='vm')
v_t3,=plt.plot(Es,vtE3,color='blue', label='vt')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$E[TeV]$', fontsize=15)
plt.ylabel('$Flux[s^{-1}Mpc^{-2}TeV^{-1}]$', fontsize=15)
plt.legend((v_e3, v_m3, v_t3),('ve','vm','vt'))
plt.title('v1 stable, v2 and v3 lifetime: 10000 s/eV')
plt.savefig('k10000.png')


veE4=veE(1e4)
vmE4=vmE(1e4)
vtE4=vtE(1e4)


plt.figure(4)
v_e4,=plt.plot(Es,veE4,color='red',label='ve')
v_m4,=plt.plot(Es,vmE4,color='green', label='vm')
v_t4,=plt.plot(Es,vtE4,color='blue', label='vt')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$E[TeV]$', fontsize=15)
plt.ylabel('$Flux[s^{-1}Mpc^{-2}TeV^{-1}]$', fontsize=15)
plt.legend((v_e4, v_m4, v_t4),('ve','vm','vt'))
plt.title('v1 stable, v2 and v3 lifetime: 0.0001 s/eV')
plt.savefig('k00001.png')
plt.show()

"""

#for testing, generate 100 redshifts at equal logarithmic intervals between 1.-3 and 1.e+1
#zs = np.power(10., np.linspace(-3, 1, 100))
#calculate corresponding LZ
#LZs = LZ(zs)
#calculate corresponding Z2s
#Z2s = Z2(zs)
#calculate corresponding Di
#Diz=Di(zs, 10e-2, 1000)
#calculate corresponding f
#fz=f(zs)
#print fz


#Uncomment the following code fragment to generate Figure 1.
"""
plt.plot(zs, LZs)
plt.xscale('log')
plt.xlabel('z', fontsize=15)
plt.yscale('log')
plt.ylabel('L(z)[MPc]', fontsize=15)
plt.show()
"""
    
#Uncomment the following code fragment to generate Figure 2.
"""
plt.plot(zs, Z2s)
plt.xscale('log')
plt.xlabel('z', fontsize=15)
plt.ylabel('Z2', fontsize=15)
plt.show()
"""

#Uncomment the following code fragment to generate Figure 3.
"""
plt.plot(zs,Diz)
plt.xlabel('z',fontsize=15)
plt.ylabel('Di(z,E0)', fontsize=15)
plt.show()
"""
