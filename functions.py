import numpy as np

def ln_funct(x, lnN, sigma, rm):
    return np.exp(lnN)/(np.sqrt(2*np.pi)*sigma)*np.exp(-0.5*((np.log(x)-np.log(rm))/sigma)**2)
