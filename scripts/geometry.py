
# Standard library
import math

# Third party
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

degrees = 10


def rectangle_slopes(radians):
    '''Return the 4 slopes of a rectangle.'''
    a1, a3 = np.tan(radians), np.tan(radians - math.radians(90))
    a2, a4 = a1, a3
    return a1, a2, a3, a4


def intersect(a, x, y):
    '''Intersection of straight line with y-axis.'''
    return y - a * x


def find_point_on_line(a, b, d, m):
  '''Find new point on straight line with slope m, with a known point (a,b),
  at a distance to that point, d.'''
  k = d / np.sqrt(1 + m ** 2)
  x = a + k
  y = b + k * m
  return x, y


def transform_rectangle(i, ec, height, coords=[0, 0], width=30):
    ts = eval("ax{i}.transData")
    tr = mpl.transforms.Affine2D().rotate_deg_around(coords[0], coords[1],
                                                     degrees)
    t = tr + ts
    rect = mpl.patches.Rectangle((coords[0], coords[1]), width, height,
                                 linewidth=1, edgecolor=ec, facecolor='none',
                                 transform=t)
    exec("ax{i}.add_patch(rect)")


def add_circle(j, arr, x, xC, y, yC, c, z):
    '''Add a circle patch to figure.'''
    for i in range(len(arr)):
        circ = plt.Circle((x[i] - xC, y[i] - yC), radius=arr[i], color=c,
                          fill=False, zorder=z)
        exec("ax{j}.add_patch(circ)")


if __name__ == '__main__':
    print(rectangle_slopes(math.radians(70)))
    print(intersect(.5, 3, 12))
    print(find_point_on_line(1, 1, 5, .5))
