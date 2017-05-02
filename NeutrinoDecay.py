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
	

	
Es=np.power(10., np.linspace(2,10,100))

nvs=nv(3,10e-2,Es)

plt.plot(Es,nvs)
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
