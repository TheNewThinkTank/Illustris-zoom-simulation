
# Standard library
import h5py
import math
import os
from typing import Final, Tuple

# Third party
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

# Local application
import TrackGalaxy
from geometry import rectangle_slopes, intersect, find_point_on_line, transform_rectangle, add_circle


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

# Switches for figures
panel_1 = 0
panel = 0
rectangles_and_stars_2dhist = 0
centered_gas_2dhist = 0
centered_galaxies_stars = 0
centered_galaxies_gas = 0
rectangular_region_stars_2dhist = 0

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

Bins: Tuple = (100, 100)

if panel_1:
    plt.figure(5, figsize=(19, 8))
    plt.subplots_adjust(left=0.06, right=0.97, top=0.96, wspace=0.25,
                        hspace=0.1, bottom=0.08)

    factor_1 = 1e10 / 0.7 / ((70.0 / 100.0) ** 2)
    factor_2 = (70.0 / 100.0) ** 2
    Range_1 = np.array([(-35, 35), (-50, 20)])

    plt.subplot(2, 4, 1)
    w_241 = Masses_star * factor_1
    plt.hist2d(x_star - xC, y_star - yC, bins=Bins, range=Range_1,
               weights=w_241, norm=LogNorm(), vmax=1e10, vmin=1e3)
    # plt.colorbar()
    # plt.ylabel('$y$ [kpc]', fontsize=24)
    plt.gca().set_xticks([])

    plt.subplot(2, 4, 2)
    w_242 = Masses_star[RectangleIDs] * factor_1
    plt.hist2d(x_star[RectangleIDs] - xC, y_star[RectangleIDs] - yC,
               bins=Bins, range=Range_1, weights=w_242, norm=LogNorm(),
               vmax=1e10, vmin=1e3)
    plt.gca().set_yticks([])
    plt.gca().set_xticks([])

    plt.subplot(2, 4, 3)
    w_243 = Masses_star[G1IDs] * factor_1
    plt.hist2d(x_star[G1IDs] - xC, y_star[G1IDs] - yC, bins=Bins,
               range=Range_1, weights=w_243, norm=LogNorm(), vmax=1e10,
               vmin=1e3)
    plt.gca().set_yticks([])
    plt.gca().set_xticks([])

    plt.subplot(2, 4, 4)
    w_244 = Masses_star[G2IDs] * factor_1
    plt.hist2d(x_star[G2IDs] - xC, y_star[G2IDs] - yC, bins=Bins,
               range=Range_1, weights=w_244, norm=LogNorm(), vmax=1e10,
               vmin=1e3)
    cbar = plt.colorbar()
    cbar.set_label('Surface stellar Density [M$_\odot$ kpc$^{-2}$]')
    plt.gca().set_yticks([])
    plt.gca().set_xticks([])

    plt.subplot(2, 4, 5)
    w_245 = SFR_gas / factor_2
    plt.hist2d(x_g, y_g, bins=Bins, range=Range_1, weights=w_245,
               norm=LogNorm(), vmax=1, vmin=1e-5)
    plt.ylabel('\t\t\t\t $y$ [kpc]', fontsize=20)

    plt.subplot(2, 4, 6)
    w_246 = SFR_gas[Gas_RectangleIDs] / factor_2
    plt.hist2d(x_g[Gas_RectangleIDs], y_g[Gas_RectangleIDs], bins=Bins,
               range=Range_1, weights=w_246, norm=LogNorm(), vmax=1, vmin=1e-5)
    plt.xlabel('\t\t\t\t $x$ [kpc]', fontsize=20)
    plt.gca().set_yticks([])

    plt.subplot(2, 4, 7)
    w_247 = SFR_gas[Gas_G1IDs] / factor_2
    plt.hist2d(x_g[Gas_G1IDs], y_g[Gas_G1IDs], bins=Bins, range=Range_1,
               weights=w_247, norm=LogNorm(), vmax=1, vmin=1e-5)
    # cbar=plt.colorbar()
    # cbar.set_label('Surface SFR [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]')
    plt.gca().set_yticks([])

    plt.subplot(2, 4, 8)
    w_248 = SFR_gas[Gas_G2IDs] / factor_2
    plt.hist2d(x_g[Gas_G2IDs], y_g[Gas_G2IDs], bins=Bins, range=Range_1,
               weights=w_248, norm=LogNorm(), vmax=1, vmin=1e-5)
    cbar = plt.colorbar()
    cbar.set_label('Surface SFR [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]')
    plt.gca().set_yticks([])
    # plt.axes().set_aspect('equal')

    plt.savefig(figure_path + 'panel_1.png')
    # plt.show()
    # sys.exit()

Range = np.array([(-80, 80), (-80, 80)])
denom = (160.0 / 100.0) ** 2

if panel:
    f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(6, 6))
    f.subplots_adjust(hspace=0, wspace=0)

    Z = x_star - xC, y_star - yC
    im = ax1.imshow(Z, extent = (-80, 80, -80, 80), interpolation="nearest")
    cb = plt.colorbar(im)
    plt.setp(cb.ax.get_yticklabels(), visible=False)

    w1 = Masses_star / denom
    ax1.hist2d(x_star - xC, y_star - yC, bins=Bins, range=Range, weights=w1,
               normed=LogNorm())
    ax1.set_ylabel('y_star-yC')

    w2 = Masses_star[RectangleIDs] / denom
    ax2.set_title('Bridge')
    ax2.hist2d(x_star[RectangleIDs] - xC, y_star[RectangleIDs] - yC, bins=Bins,
               range=Range, weights=w2, normed=LogNorm())
    ax2.tick_params(axis='y', which='both', left='off', right='off',
                    labelleft='off')

    w3 = Masses_star[G1IDs] / denom
    ax3.set_title('Galaxy 1')
    ax3.hist2d(x_star[G1IDs] - xC, y_star[G1IDs] - yC, bins=Bins, range=Range,
               weights=w3, normed=LogNorm())
    ax3.tick_params(axis='y', which='both', left='off', right='off',
                    labelleft='off')
    ax3.set_xlabel('x_star-xC')

    w4 = Masses_star[G2IDs] / denom
    ax4.set_title('Galaxy 2')
    ax4.hist2d(x_star[G2IDs] - xC, y_star[G2IDs] - yC, bins=Bins, range=Range,
               weights=w4, normed=LogNorm())
    ax4.tick_params(axis='y', which='both', left='off', right='off',
                    labelleft='off')

    transform_rectangle(1, 'b', 12, [-10, -25])
    transform_rectangle(1, 'r', 24, [-12, -13])
    transform_rectangle(1, 'g', 24, [0, -47])
    transform_rectangle(2, 'b', 12, [-10, -25])
    transform_rectangle(3, 'g', 24, [0, -47])
    transform_rectangle(4, 'r', 24, [-12, -13])

    f.savefig(figure_path + 'panel.png')

if rectangles_and_stars_2dhist:
    f, (ax1) = plt.subplots(1, 1, figsize=(6, 6))
    w = Masses_star[RectangleIDs] / denom
    plt.hist2d(x_star[RectangleIDs] - xC, y_star[RectangleIDs] - yC,
               bins=Bins, range=Range, weights=w, normed=LogNorm())
    # plt.hist2d(x_star[G1IDs] - xC, y_star[G1IDs] - yC, bins=Bins,
    #            range=Range, weights=Masses_star[G1IDs] / denom,
    #            normed=LogNorm())
    # plt.hist2d(x_star[G2IDs] - xC, y_star[G2IDs] - yC, bins=Bins,
    #            range=Range, weights=Masses_star[G2IDs] / denom,
    #            normed=LogNorm())
    plt.colorbar()

    transform_rectangle(1, 'b', 12, [-10, -25])
    transform_rectangle(1, 'r', 24, [-12, -13])
    transform_rectangle(1, 'g', 24, [0, -47])

    ax1.set_xlabel(r'$x-x_c$', fontsize=10)
    ax1.set_ylabel(r'$y-y_c$', fontsize=10)
    ax1.set_title('Centralized stars', fontsize=15)
    ax1.set_xlim(-80, 80)
    ax1.set_ylim(-80, 80)
    f.savefig(figure_path + 'rectangles_and_stars_2dhist.png')

if centered_gas_2dhist:
    f, (ax1) = plt.subplots(1, 1, figsize=(6, 6))
    plt.hist2d(x_gas - xC, y_gas - yC, bins=Bins, range=Range,
               weights=SFR_gas / denom, norm=LogNorm())
    plt.colorbar()
    add_circle(1, r_galaxy, x_galaxy, xC, y_galaxy, yC, 'r', 2)
    ax1.set_xlabel(r'$x-x_c$', fontsize=10)
    ax1.set_ylabel(r'$y-y_c$', fontsize=10)
    ax1.set_title('Centralized gas', fontsize=15)
    ax1.set_xlim(-80, 80)
    ax1.set_ylim(-80, 80)
    f.savefig(figure_path + 'centered_gas_2dhist.png')

if centered_galaxies_stars:
    f, (ax1) = plt.subplots(1, 1, figsize=(6, 6))
    ax1.plot(x_star - xC, y_star - yC, 'ob', ms=1, zorder=1)
    ax1.plot(x_galaxy - xC, y_galaxy - yC, 'or', ms=1, zorder=2)
    add_circle(1, r_galaxy, x_galaxy, xC, y_galaxy, yC, 'k', 3)
    ax1.set_xlabel(r'$x-x_c$', fontsize=10)
    ax1.set_ylabel(r'$y-y_c$', fontsize=10)
    ax1.set_title('Centralized galaxies, Rhalf and stars', fontsize=15)
    ax1.set_xlim(-300, 200)
    ax1.set_ylim(-300, 200)
    f.savefig(figure_path + 'centered_galaxies_stars.png')

if centered_galaxies_gas:
    f, (ax1) = plt.subplots(1, 1, figsize=(6, 6))
    ax1.plot(x_gas - xC, y_gas - yC, 'ob', ms=1, zorder=1)
    ax1.plot(x_galaxy - xC, y_galaxy - yC, 'or', ms=1, zorder=2)
    add_circle(1, r_galaxy, x_galaxy, xC, y_galaxy, yC, 'k', 3)
    ax1.set_xlabel(r'$x-x_c$', fontsize=10)
    ax1.set_ylabel(r'$y-y_c$', fontsize=10)
    ax1.set_title('Centralized galaxies, Rhalf and gas', fontsize=15)
    ax1.set_xlim(-300, 200)
    ax1.set_ylim(-300, 200)
    f.savefig(figure_path + 'centered_galaxies_gas.png')

# RectangleIDs = np.where((y_star - yC > a1 * (x_star - xC) + b1)
#                         * (y_star - yC < a2 * (x_star - xC) + b2)
#                         * (y_star - yC > a3 * (x_star - xC) + b3)
#                         * (y_star - yC < a4 * (x_star - xC) + b4))

if rectangular_region_stars_2dhist:
    f, (ax) = plt.subplots(1, 1, figsize=(6, 6))
    plt.hist2d(x_star[RectangleIDs] - xC, y_star[RectangleIDs] - yC,
               bins=Bins, range=Range,
               weights=Masses_star[RectangleIDs] / denom, normed=True)  # norm=LogNorm()
    plt.colorbar()
    # transform_rectangle('', 'b', 12, [-10, -25])
    ax.set_xlabel(r'$x-x_c$', fontsize=10)
    ax.set_ylabel(r'$y-y_c$', fontsize=10)
    ax.set_title('Stars between merging galaxies', fontsize=15)
    ax.set_xlim(-80, 80)
    ax.set_ylim(-80, 80)
    f.savefig(figure_path + 'rectangular_region_stars_2dhist.png')

plt.show()
