from matplotlib import pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import matplotlib as mpl
from functions import ln_funct
import libspheroids

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"})

def display_solution(c, p, xopt):
    """docstring for display_solution"""
    grid = plt.GridSpec(3,2)
    plt.subplot(grid[0,0])
    plt.semilogy(c.interp_wavelengths, p.calc_vector[:6], 'b')
    plt.semilogy(c.interp_wavelengths, c.meas_vector_interp[:6],'bo')
    plt.ylabel(r'Back. coeff. 1/km $\times$ sr', color='b')
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    
    #print(c.wavelengths[:2], p.calc_vector[-2:], c.meas_vector[-2:])
    plt.subplot(grid[1,0])
    plt.semilogy(c.interp_wavelengths, p.calc_vector[6:],'g')
    plt.semilogy(c.interp_wavelengths, c.meas_vector_interp[6:], 'go')
    plt.ylabel(r'Ext. coeff. 1/km', color='g')
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    
    
    plt.subplot(grid[2,0])
    plt.plot(c.interp_wavelengths, p.calc_vector[6:]/p.calc_vector[:6], 'r')
    plt.plot(c.interp_wavelengths, c.meas_vector_interp[6:]/c.meas_vector_interp[:6], 'ro')
    plt.ylabel('LR ratio sr', color='r')
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    plt.xlabel('Wavelength, nm')
    
    plt.subplot(grid[0,1])
    plt.plot(c.interp_wavelengths, p.lidar_depol_ratio, 'kd-')
    plt.ylabel(r'lid. dep. rat., \%')
    plt.xlabel('Wavelength, nm')
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    
    rr = np.logspace(np.log10(c.r_min), np.log10(c.r_max), 50)
    plt.subplot(grid[1:,1])
    plt.semilogx(rr, ln_funct(rr, xopt.params['lnN1'].value, xopt.params['sigma1'].value, xopt.params['rm1'].value), 'b-')
    plt.semilogx(rr, ln_funct(rr, xopt.params['lnN2'].value, xopt.params['sigma2'].value, xopt.params['rm2'].value), 'g-')
    plt.semilogx(libspheroids.mo_dls.rrr[:c.knots_count], libspheroids.mo_dls.sd[:c.knots_count], 'ro', mfc=None)
    plt.ylabel(r'$\frac{dV}{d\ln r}$, $\frac{\mu m^3}{\mu m^3}$')
    plt.xlabel(r'Radius, $\mu m$')
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    
    plt.suptitle('Lidar measurements processing results', fontsize="x-large")
    plt.tight_layout()
    print(c.interp_wavelengths)
    print(p.lidar_depol_ratio)
    print(p.back_lidar_ratio)
    plt.show()
