import numpy as np
import matplotlib.pyplot as plt
import orix
from orix.quaternion.orientation import Orientation, Misorientation
from orix.quaternion.rotation import Rotation
from orix.quaternion.symmetry import D6
from orix.quaternion.orientation_region import OrientationRegion
from orix.vector.neo_euler import AxAngle
from orix.vector import Vector3d
from orix import plot
from itertools import product, combinations

import math

from scipy.spatial.transform import Rotation as R


def eqn4(a, h, k, l, a11, a12, a13):
    top = np.abs(np.add(np.add(np.multiply(h, a11), np.multiply(k, a12)), np.multiply(l, a13)))
    denom = np.add(np.add(np.power(h, 2), np.power(k, 2)), np.power(l, 2))
    lambda_calc = 2*a*top/denom
    return lambda_calc


h, k, l, m, d, f = np.loadtxt("Ni.hkl", dtype=float, skiprows=1, unpack=True)

# ori = Orientation.from_euler(np.radians([10, 20, 30]), convention="Krakow_Hielscher")
# print(ori)
# print(ori.size)

# r = R.from_quat([0, 0, np.sin(np.pi/4), np.cos(np.pi/4)])
b_r = R.from_quat([1, 0, 0, 0])
b_matrix = b_r.as_matrix()
r = R.from_quat([1, 0, 0, 0])
matrix = r.as_matrix()
print(matrix)

lambda_max = np.multiply(2, d)
lambda_result = eqn4(3.5195, h=h, k=k, l=l, a11=matrix[0, 0], a12=matrix[0, 1], a13=matrix[0, 2])
cos_alpha = np.divide(lambda_result, lambda_max)
alpha_rad = np.arccos(cos_alpha)
alpha_deg = np.multiply(alpha_rad, (180 / np.pi))
for index, value in enumerate(alpha_deg):
    if math.isnan(value):
        print("aligned with {} {} {} direction".format(int(h[index]), int(k[index]), int(l[index])))
    else:
        pass


print(lambda_result)


x, y = np.loadtxt("spring_dips_clipped_2.txt", dtype=float, skiprows=1, unpack=True)

fig = plt.figure(figsize=plt.figaspect(0.5))
ax11 = fig.add_subplot(121, projection="3d")

#####################################################################
# quiver for incident beam
len_beam = 5
b_start = [0, 0, 0]
try:
    x_to_x = 1/int(b_matrix[0,0])
except ZeroDivisionError:
    x_to_x = 0
try:
    y_to_y = 1/int(b_matrix[1,1])
except ZeroDivisionError:
    y_to_y = 0
try:
    z_to_x = 1/int(b_matrix[2,0])
except ZeroDivisionError:
    z_to_x = 0
# b_inc = [(1 / matrix[0, 0]), (1 / matrix[1, 1]), (1 / matrix[2, 0])]
ax11.quiver(b_start[0], b_start[1], b_start[2], x_to_x, y_to_y, z_to_x, length=len_beam, color="red")
#####################################################################
r_cube = [-1, 1]
for s, e in combinations(np.array(list(product(r_cube, r_cube, r_cube))), 2):
    print(s, e)
    if np.sum(np.abs(s-e)) == r_cube[1]-r_cube[0]:

        ax11.plot3D(*zip(s, e), color="b")

ax11.plot([-4,4],[0,0],zs=[0,0], color="black", markersize=5)
ax11.plot([0,0],[-4,4],zs=[0,0], color="black", markersize=5)
ax11.plot([0,0],[0,0],zs=[-4,4], color="black", markersize=5)

ax12 = fig.add_subplot(122)
ax12.plot(x, y, "k.", markersize=2)
for index, value in enumerate(lambda_result):
    ax12.axvline(x=value, color="r")
for index, value in enumerate(lambda_max):
    ax12.axvline(x=value, color="g")
plt.show()
