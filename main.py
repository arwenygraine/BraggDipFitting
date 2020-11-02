import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.special import erf, erfc

def gauss(x, a, u, s, m, c):
    y = np.add(np.multiply(a, np.exp(np.multiply(np.divide(np.power(np.subtract(x, u), 2),
                                               2 * s ** 2), -1))),
               np.add(np.multiply(m, x), c))
    return y

def double_gauss(x, a, a2, u, u2, s, s2, m, m2, c, c2):
    y = np.add(gauss(x, a, u, s, m, c), gauss(x, a2, u2, s2, m2, c2))
    return y

def super_gauss(x, a, u, s, m, c):
    y = np.add(np.multiply(a, np.exp(np.multiply(np.power(np.divide(np.power(np.subtract(x, u), 2),
                                                           2 * s ** 2), 2), -1))),
               np.add(np.multiply(m, x), c))
    return y

def asymmetric_gauss(x, a, u, s, m, c, alpha):
    density_fn = np.multiply(np.divide(1, np.sqrt(np.multiply(2, np.pi))), np.exp(np.multiply(np.divide(np.power(x, 2), 2), -1)))
    distribution_fn = np.multiply(0.5, np.add(1, erf(np.divide(np.multiply(alpha, x), np.sqrt(2)))))
    y = np.multiply(np.multiply(2, gauss(x, a, u, s, m, c)), distribution_fn)
    print(distribution_fn, y)
    return y

def exp_mod_gauss(x, amp, mu, sig, tau, m, c):
    const_term = np.multiply(np.divide(np.multiply(-1*amp, sig), tau), np.sqrt(np.divide(np.pi, 2)))
    exp_term = np.exp(np.subtract(np.multiply(0.5, np.power(np.divide(sig, tau), 2)), np.divide(np.subtract(x, mu), tau)))
    erf_term = erf(np.multiply(np.divide(1, np.sqrt(2)), np.subtract(np.divide(sig, tau), np.divide(np.subtract(x, mu), sig))))
    y = np.add(np.multiply(np.multiply(const_term, exp_term), erf_term), np.add(np.multiply(m, x), c))
    print("const: " + str(const_term) + "  exp: " + str(exp_term) + "  erf: " + str(erf_term))
    return y

def RMSE(x, y, p_opt, fitting_func):
    return np.sqrt(np.divide(np.sum([(y0 - y) ** 2 for y, y0 in zip(y, fitting_func(x, *p_opt))]), len(y)))


with open('spring_dips_clipped_2.txt', 'r') as f:
    lines = f.readlines()
    x = []
    y = []
    for index, line in enumerate(lines):
        if index != 0:
            x.append(float(line.split()[0]))
            y.append(float(line.split()[1]))

        else:
            x_title = line.split()[0]
            y_title = line.split()[1]

#set data range
x = x[2050:2300]
y = y[2050:2300]

# x = np.linspace(0, 10, 10000)

troughs, _ = find_peaks(np.negative(x))

#average transmission for e
transmission_average = np.average(y)
print(transmission_average)

#find min y value, find index, map onto x values, find position of minimum
minimum = np.amin(y)
y_index = 0
for index, value in enumerate(y):
    if value == minimum:
        y_index = index
    else:
        pass
x_min_pos = x[y_index]

#average slope across data range
slope = (y[len(y)-1] - y[0]) / (x[len(x)-1] - x[0])

#height of peak (difference between straight line value and peak value at xmin)
#line_x_min = slope * x_min_pos + transmission_average
#print(line_x_min)
peak_height = minimum - transmission_average
print(peak_height)

fn = double_gauss

p0 = [-10, x_min_pos, 0.1, slope, transmission_average]
p0 = [-10, x_min_pos, 0.1, slope, transmission_average, -11, 4.65, 0.03, -7, 70]
# p0 = [-4, 4.55, 0.02, -6, 66, 20]       #'a, u, s, m, c, alpha'
# p0 = [-20, 3.4, 0.6, 0.8, -0.1, 39]     #exp_mod_gauss
print(p0)

#test, fitting
# popt, pcov = curve_fit(fn, x, y, p0=p0)
# print(popt, pcov)
#
# plt.xlabel(x_title)
# plt.ylabel(y_title)
# plt.plot(x, y)
# #plt.plot(troughs, x[troughs], "x")
# plt.plot(x, fn(x, *popt))
# plt.fill_between(
#     x, fn(x, *np.subtract(popt, np.sqrt(np.diag(pcov)))), fn(x, *np.add(popt, np.sqrt(np.diag(pcov)))), alpha=0.2)
# plt.show()
#
# for opt, cov in zip(popt, np.sqrt(np.diag(pcov))):
#     print("{} pm {}".format(opt, cov))
# print(RMSE(x, y, popt, fn))


#test, no fitting
plt.plot(x, fn(x, *p0))
#plt.plot(x, y)
plt.show()