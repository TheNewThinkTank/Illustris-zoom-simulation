
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

# Switches for figures
panel_1 = 0
panel = 0
rectangles_and_stars_2dhist = 0
centered_gas_2dhist = 0
centered_galaxies_stars = 0
centered_galaxies_gas = 0
rectangular_region_stars_2dhist = 0

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
