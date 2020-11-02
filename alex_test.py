import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.special import erf, erfc

###############################################################################

def linear(x, m, c):
    return np.add(np.multiply(m, x), c)

def gauss(x, a, u, s):
    return np.multiply(a,
                       np.exp(np.multiply(np.divide(np.power(np.subtract(x, u), 2),
                                                    2 * s ** 2),
                                          -1)))
    return y
    
def gauss_with_linear(x, a, u, s, m, c):
    return np.add(gauss(x, a, u, s), linear(x, m, c))
    
def asymmetric_gauss(x, a, u, s, alpha):
    amp = (a / (s * np.sqrt(2 * np.pi)))
    spread = np.exp((-(x - u) ** 2.0) / (2 * s ** 2.0))
    skew = (1 + erf((alpha * (x - u)) / (s * np.sqrt(2))))
    return np.multiply(amp * spread, skew)
    
def asymmetric_gauss_with_linear(x, a, u, s, m, c, alpha):
    return np.add(linear(x, m, c), asymmetric_gauss(x, a, u, s, alpha))
    
def double_gauss_with_linear(x, a, a2, u, u2, s, s2, m, c):
    return np.add(linear(x, m, c), np.add(gauss(x, a, u, s), gauss(x, a2, u2, s2)))
    
def RMSE(x, y, p_opt, fitting_func):
    return np.sqrt(np.divide(np.sum([(y0 - y) ** 2 for y, y0 in zip(y, fitting_func(x, *p_opt))]), len(y)))

###############################################################################

x, y = np.loadtxt("spring_dips_clipped_2.txt", dtype=float, skiprows=1, unpack=True)
data_x = x[2050:2300]
data_y = y[2050:2300]

fit_x = x[2050:2300]    # this does not need to be different but its good for testing ...
fit_y = y[2050:2300]    # ... fitting of data sets to more constrained boundaries

fn = double_gauss_with_linear
p0 = (-3.6, -1.6, 4.54, 4.56, 0.01, -0.02, -6, 64)

###############################################################################

popt, pcov = curve_fit(fn, fit_x, fit_y, p0=p0, maxfev=10000)
plt.plot(data_x, data_y, "k.", markersize=2)
plt.plot(fit_x, fn(fit_x, *popt))
plt.fill_between(
    fit_x, fn(fit_x, *np.subtract(popt, np.sqrt(np.diag(pcov)))), fn(fit_x, *np.add(popt, np.sqrt(np.diag(pcov)))), alpha=0.2)
plt.show()
for opt, cov in zip(popt, np.sqrt(np.diag(pcov))):
    print("{} pm {}".format(opt, cov))
print(RMSE(fit_x, fit_y, popt, fn))
#plt.plot(fit_x, fn(fit_x, *p0))
plt.show()