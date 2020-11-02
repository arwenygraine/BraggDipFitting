import numpy as np
import matplotlib.pyplot as plt

#f(x) = a * exp(-((x-b)**2)/2c**2)

def gauss(x, a, b, c, d, e):
    y = a * np.exp(-((x-b)**2)/2*c**2) + (d*x)+e
    return y

x = np.linspace(-10, 10, num=1000)

y = gauss(x, -3, 0, 1, -0.5, 0)

plt.plot(x, y)
plt.show()