import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math


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


def nv(E,k, zmax):
	def integrand(z,k,E):
		return f(z)*Di(z,k,E)*(1+z)**(-2)*E**(-2)
	def integrator(E,k=k,zmax=zmax):
		return integrate.quad(integrand,0,zmax, args=(k,E))[0]
	if np.isscalar(E):
		return integrator(E,k,zmax)
	else:
		return np.asarray(map(integrator, E))
	

	
Es=np.power(10., np.linspace(0.,4.,100))

nvs=nv(Es,1e-2,3.)
nvs2=nv(Es,1.,3.)
nvs3=nv(Es,1e4,3)
nvs4=nv(Es,1e-17,3)
#plt.plot(Es,nvs,color='red')
plt.plot(Es,nvs2, color='blue')
#plt.plot(Es,nvs3, color='green')
plt.plot(Es, nvs4, color='purple')
plt.xscale('log')
plt.yscale('log')
plt.show()




def c(theta):
	return math.cos(theta)


def s(theta):
	return math.sin(theta)


theta12=0.5764
theta13=0.1463
theta23=0.7223

def v1(ve,vm,vt):
	return ve*c(theta13)*c(theta12)+vm*(-c(theta23)*s(theta12)-s(theta23)*c(theta12)*s(theta13))+vt*(s(theta23)*s(theta12)-c(theta23)*c(theta12)*s(theta13)) 

def v2(ve,vm,vt):
	return ve*c(theta13)*s(theta12)+vm*(c(theta23)*c(theta12)-s(theta23)*s(theta12)*s(theta13))+vt*(-s(theta23)*s(theta12)-c(theta23)*c(theta12)*s(theta13))

def v3(ve,vm,vt):
	return ve*s(theta13)+vm*c(theta13)*s(theta23)+vt*c(theta13)*c(theta23)


m1=abs(v1(1.,2.,0.))
m2=abs(v2(1.,2.,0.))
m3=abs(v3(1.,2.,0.))

print m1, m2, m3

v1E=m1*nvs4
v2E=m2*nvs2
v3E=m3*nvs2

veE=abs(v1E*c(theta13)*c(theta12)+v2E*c(theta13)*s(theta12)+v3E*s(theta13))
vmE=abs(v1E*(-c(theta23)*s(theta12)-s(theta23)*c(theta12)*s(theta13))+v2E*(c(theta23)*c(theta12)-s(theta23)*s(theta12)*s(theta13))+v3E*c(theta13)*s(theta23))
vtE=abs(v1E*(s(theta23)*s(theta12)-c(theta23)*c(theta12)*s(theta13))+v2E*(-s(theta23)*c(theta12)-c(theta23)*s(theta12)*s(theta13))+v3E*c(theta13)*c(theta23))



v_e,=plt.plot(Es,veE,color='red',label='ve')
v_m,=plt.plot(Es,vmE,color='green', label='vm')
v_t,=plt.plot(Es,vtE,color='blue', label='vt')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$E[TeV]$', fontsize=15)
plt.ylabel('$Flux[s^{-1}Mpc^{-2}TeV^{-1}]$', fontsize=15)
plt.legend((v_e, v_m, v_t),('ve','vm','vt'))
plt.show()

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
