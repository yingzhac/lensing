import numpy as np
import time


h = 0.7
M_200 = 1e10
dndV = 2.37e-3  #if 1e9.9-1e10.1, dndV=2.37e-3
z_S = 1.2
G = 4.79*1e-20
sigma_sh = 0.25
c = 5.0
H_0 = 100 * h
#rho_cr = 2.77*1e11 * h**2
n_S_LSST = 50.0
n_S_DES = 5.0
A_DES = 5e3*(np.pi/180.0)**2
A_LSST = 2e4*(np.pi/180.0)**2

sample_fg = 10000
sample_bg = 2000
sigma = 0.25


def hubble(z):
    return np.sqrt(0.3*(1+z)**3+0.7)

def comoving_distance(z):
    z_prime = np.linspace(0,z,1000)
    integrand = 1.0/hubble(z_prime)
    integral = np.trapz(integrand,z_prime)
    return integral*3e3/h

def gamma(x):
    if x<1:
        g = 8*np.arctanh(np.sqrt((1.0-x)/(1.0+x)))/(x**2*np.sqrt(1.0-x**2))  +  4*np.log(x/2.0)/x**2  -  2.0/(x**2-1.0)  +   4*np.arctanh(np.sqrt((1.0-x)/(1.0+x)))/((x**2-1.0)*np.sqrt(1.0-x**2))######
    elif x==1:
        g = 10.0/3.0 +4*np.log(1.0/2.0)
    else:
        g = 8*np.arctan(np.sqrt((x-1.0)/(1.0+x)))/(x**2*np.sqrt(x**2-1.0))  +  4*np.log(x/2.0)/x**2  -  2.0/(x**2-1.0)  +   4*np.arctan(np.sqrt((x-1.0)/(1.0+x)))/pow(x**2-1.0,1.5)######
    return g

def delta_c(c):
    return 200.0*c**3/( 3.0*(np.log(1+c)-c/(1.0+c)) )
    
def rho_cr(z):
    return 2.77*1e11 * h**2 * hubble(z)**2

def R_S(M_200,z):
    return pow(M_200*4.302*1e-9/(100.0 * c**2 * (100.0*h*hubble(z))**2 ),1.0/3.0)

def N_BG(theta_in_arcmin,n_S):
    return n_S * 2*np.pi * theta_in_arcmin


z_L = np.random.rand(sample_fg)*0.1
D_L = np.array([comoving_distance(z)/(1+z) for z in z_L])
theta_x = np.load('theta_x.npy')
theta_y = np.load('theta_y.npy')


t1 = time.time()
M_200s = np.logspace(8,12,101)
shear_catalog = np.zeros((101,sample_fg,sample_bg))

for mm in range(101):
    np.save('sim_doing.npy',mm)
    M_200 = M_200s[mm]
    r_x = np.outer(D_L,theta_x*np.pi/(180.0*60))
    r_y = np.outer(D_L,theta_y*np.pi/(180.0*60))
    r_p = np.sqrt(r_x**2+r_y**2)
    D_S = comoving_distance(z_S)/(1+z_S)
    D_SL = D_S - D_L
    Sigma_cr = D_S / (4 * np.pi * G * D_SL * D_L)  # in units: M_solar Mpc^{-2}
    r_s = np.array([R_S(M_200,z) for z in z_L])
    rho_crit = np.array([rho_cr(z) for z in z_L])
    gamma_unit = r_s*delta_c(c)*rho_crit/Sigma_cr
    #theta_in_rad = np.array([r_p[i]/D_L[i] for i in range(sample_fg)])
    x = np.array([r_p[i]/r_s[i] for i in range(sample_fg)])


    gamma_in_unit = np.zeros((sample_fg,sample_bg))
    for i in range(sample_fg):
        for j in range(sample_bg):
            gamma_in_unit[i,j] = gamma(x[i,j])

    gamma_unitless = np.array([gamma_in_unit[i] * gamma_unit[i] for i in range(sample_fg)])
    shear_catalog[mm] = gamma_unitless


    
np.save('shear_ctlg.npy',shear_catalog)


t2 = time.time()

print t2-t1
