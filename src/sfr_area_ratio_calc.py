
import h5py
import math
import os
from typing import Final
# from matplotlib.colors import LogNorm
# import matplotlib.pyplot as plt
import numpy as np
import TrackGalaxy
from src.rectangle import (rectangle_slopes,
                      intersect,
                      find_point_on_line,
                      transform_rectangle,
                      add_circle
                      )

path = os.getcwd() + '/'
T = TrackGalaxy.TrackGalaxy(np.array([67]), '1330-3', Dir=path)

# Read in position, SFR, stellar mass, gas mass, dark matter mass of all the
# galaxies (and subhalos) from the simulation
Attrs = T.GetGroups(67, Attrs=['/Subhalo/SubhaloPos', '/Subhalo/SubhaloSFR',
                    '/Subhalo/SubhaloMassType', '/Subhalo/SubhaloHalfmassRad'])
Pos = Attrs['/Subhalo/SubhaloPos']  # in comoving kpc/h
SFR = Attrs['/Subhalo/SubhaloSFR']  # in Msun/yr

factor: float = 1e10 / 0.7

Mstar = Attrs['/Subhalo/SubhaloMassType'][:, 4] * factor  # in Msun
Mgas = Attrs['/Subhalo/SubhaloMassType'][:, 0] * factor  # in Msun
Mdm = Attrs['/Subhalo/SubhaloMassType'][:, 1] * factor  # in Msun
Rhalf = Attrs['/Subhalo/SubhaloHalfmassRad']  # in comoving kpc/h

# there is a lot of subhalos (thousands), but most of them don't have any stars.
# So let us just pick out the galaxies with at least (Mstar > 1e7 Msun):
GoodIDs = np.where(Mstar > 1e7)
Pos = Pos[GoodIDs]
SFR = SFR[GoodIDs]
Mstar = Mstar[GoodIDs]
Mgas = Mgas[GoodIDs]
Mdm = Mdm[GoodIDs]

x_galaxy = Pos[:, 0]
y_galaxy = Pos[:, 1]
z_galaxy = Pos[:, 2]
r_galaxy = Rhalf[GoodIDs]

galaxy_path = os.path.abspath(os.path.join(os.getcwd(), '..')) + '/'
figure_path = galaxy_path + 'figures/'
text_files_path = '.'
Filename = 'snapshot_067.hdf5'
SnapshotFile = h5py.File(Filename,'r')

# Center coordinates of galaxy
xC: Final = 39937.98
yC: Final = 34857.863
zC: Final = 37441.234

# Gas
Pos0 = SnapshotFile['PartType0/Coordinates'].value
Vel0 = SnapshotFile['PartType0/Velocities'].value
SFR_gas = SnapshotFile['PartType0/StarFormationRate'].value
V0 = SnapshotFile['PartType0/Potential'].value
Masses_gas = SnapshotFile['PartType0/Masses'].value
x_gas = Pos0[:, 0]
y_gas = Pos0[:, 1]
z_gas = Pos0[:, 2]
vx_gas = Vel0[:, 0]
vy_gas = Vel0[:, 1]
vz_gas = Vel0[:, 2]

# Dark Matter
Pos1 = SnapshotFile['PartType1/Coordinates'].value
Vel1 = SnapshotFile['PartType1/Velocities'].value
V1 = SnapshotFile['PartType1/Potential'].value
x_dm = Pos1[:, 0]
y_dm = Pos1[:, 1]
z_dm = Pos1[:, 2]
vx_dm = Vel1[:, 0]
vy_dm = Vel1[:, 1]
vz_dm = Vel1[:, 2]

# Stars
Pos4 = SnapshotFile['PartType4/Coordinates'].value
Vel4 = SnapshotFile['PartType4/Velocities'].value
V4 = SnapshotFile['PartType4/Potential'].value
Masses_star = SnapshotFile['PartType4/Masses'].value
x_star = Pos4[:, 0]
y_star = Pos4[:, 1]
z_star = Pos4[:, 2]
vx_star = Vel4[:, 0]
vy_star = Vel4[:, 1]
vz_star = Vel4[:, 2]

# AGN
Pos5 = SnapshotFile['PartType5/Coordinates'].value
Vel5 = SnapshotFile['PartType5/Velocities'].value
V5 = SnapshotFile['PartType5/Potential'].value
Masses_agn = SnapshotFile['PartType5/Masses'].value
x_agn = Pos5[:, 0]
y_agn = Pos5[:, 1]
z_agn = Pos5[:, 2]
vx_agn = Vel5[:, 0]
vy_agn = Vel5[:, 1]
vz_agn = Vel5[:, 2]

x: float = x_star - xC
y: float = y_star - yC
x_g: float = x_gas - xC
y_g: float = y_gas - yC

degrees: float = 10
radians: float = math.radians(degrees)

# Bridge rectangle ----------------------------------------------
Ba1, Ba2, Ba3, Ba4 = rectangle_slopes(radians)

# BL, BR, TL, TR: Bottom Left, Bottom Right, Top Left, Top Right
BxBL, ByBL, BxBR, ByBR = -10, -25, 19.544, -19.791
BxTL, ByTL, BxTR, ByTR = -12, -13, 17.5, -8

# intersections
Bb1 = intersect(Ba1, BxBL, ByBL)
Bb2 = intersect(Ba2, BxTR, ByTR)
Bb3 = intersect(Ba3, BxTL, ByTL)
Bb4 = intersect(Ba4, BxTR, ByTR)

Gas_RectangleIDs = np.where((y_g > Ba1 * x_g + Bb1)
                            * (y_g < Ba2 * x_g + Bb2)
                            * (y_g > Ba3 * x_g + Bb3)
                            * (y_g < Ba4 * x_g + Bb4)
                            )
RectangleIDs = np.where((y > Ba1 * x + Bb1)
                        * (y < Ba2 * x + Bb2)
                        * (y > Ba3 * x + Bb3)
                        * (y < Ba4 * x + Bb4)
                        )
hB, wB = 12, 30  # height and width
AB: float = hB * wB  # area = 360
ratioB = sum(SFR_gas[Gas_RectangleIDs]) / AB  # 0.000293586789828
# print(f'{ratioB =}')

# areas of galaxy1 and galaxy2 rectangles
hG, wG = 24, 30  # heights and widths
AG = hG * wG  # areas (each) = 720

# Galaxy 1 rectangle (lower)
G1xBL, G1yBL = 0, -47
G1xBR, G1yBR = find_point_on_line(G1xBL, G1yBL, wG, Ba1)
G1xTL, G1yTL = find_point_on_line(G1xBL, G1yBL, -hG, Ba3)
G1xTR, G1yTR = find_point_on_line(G1xTL, G1yTL, wG, Ba2)

# intersections
G1b1 = intersect(Ba1, G1xBL, G1yBL)
G1b2 = intersect(Ba2, G1xTR, G1yTR)
G1b3 = intersect(Ba3, G1xTL, G1yTL)
G1b4 = intersect(Ba4, G1xTR, G1yTR)

G1IDs = np.where((y > Ba1 * x + G1b1)
                 * (y < Ba2 * x + G1b2)
                 * (y > Ba3 * x + G1b3)
                 * (y < Ba4 * x + G1b4)
                 )
Gas_G1IDs = np.where((y_g > Ba1 * x_g + G1b1)
                     * (y_g < Ba2 * x_g + G1b2)
                     * (y_g > Ba3 * x_g + G1b3)
                     * (y_g < Ba4 * x_g + G1b4)
                     )
ratioG1 = sum(SFR_gas[Gas_G1IDs]) / AG  # 0.000649032024383
# print(f'{ratioG1 =}')

# Galaxy 2 rectangle (higher)
G2xBL, G2yBL = -12,-13
G2xBR, G2yBR = find_point_on_line(G2xBL, G2yBL, wG, Ba1)
G2xTL, G2yTL = find_point_on_line(G2xBL, G2yBL, -hG, Ba3)
G2xTR, G2yTR = find_point_on_line(G2xTL, G2yTL, wG, Ba2)

# intersections
G2b1 = intersect(Ba1, G2xBL, G2yBL)
G2b2 = intersect(Ba2, G2xTR, G2yTR)
G2b3 = intersect(Ba3, G2xTL, G2yTL)
G2b4 = intersect(Ba4, G2xTR, G2yTR)

G2IDs = np.where((y > Ba1 * x + G2b1)
                 * (y < Ba2 * x + G2b2)
                 * (y > Ba3 * x + G2b3)
                 * (y < Ba4 * x + G2b4)
                 )
Gas_G2IDs = np.where((y_g > Ba1 * x_g + G2b1)
                     * (y_g < Ba2 * x_g + G2b2)
                     * (y_g > Ba3 * x_g + G2b3)
                     * (y_g < Ba4 * x_g + G2b4)
                     )
ratioG2 = sum(SFR_gas[Gas_G2IDs]) / AG  # 0.0190599170674
# print(f'{ratioG2 =}')

Bins: tuple = (100, 100)
