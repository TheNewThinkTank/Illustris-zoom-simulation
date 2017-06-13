import h5py, IPython, astropy, pylab, os, PIL, math, sys, copy, TrackGalaxy
import matplotlib.pyplot as plt
import modules.Image as Image
import modules.Stars as Stars
import modules.Gas as Gas
import modules.BH as BH
import matplotlib.patches as patches
import matplotlib as mpl
import numpy as np
from astropy.cosmology import WMAP9
from matplotlib.colors import LogNorm
from modules.Gas import CalcTandDens
from matplotlib import rc

path = os.getcwd() + '/'
T    = TrackGalaxy.TrackGalaxy(np.array([67] ),'1330-3',Dir = path)

#Read in position, SFR, stellar mass, gas mass, dark matter mass of all the galaxies (and subhalos) from the simulation
Attrs = T.GetGroups(67,Attrs=['/Subhalo/SubhaloPos','/Subhalo/SubhaloSFR','/Subhalo/SubhaloMassType','/Subhalo/SubhaloHalfmassRad'])
Pos   = Attrs['/Subhalo/SubhaloPos']#in comoving kpc/h
SFR   = Attrs['/Subhalo/SubhaloSFR']#in Msun/yr
Mstar = Attrs['/Subhalo/SubhaloMassType'][:,4]*1e10/0.7#in Msun
Mgas  = Attrs['/Subhalo/SubhaloMassType'][:,0]*1e10/0.7#in Msun
Mdm   = Attrs['/Subhalo/SubhaloMassType'][:,1]*1e10/0.7#in Msun
Rhalf = Attrs['/Subhalo/SubhaloHalfmassRad']#in comoving kpc/h

#there is a lot of subhalos (thousands), but most of them don't have any stars. So let us just pick out the galaxies with at least Mstar>1e7 Msun:
GoodIDs = np.where(Mstar > 1e7)
Pos     = Pos[GoodIDs]
SFR     = SFR[GoodIDs]
Mstar   = Mstar[GoodIDs]
Mgas    = Mgas[GoodIDs]
Mdm     = Mdm[GoodIDs]

x_galaxy = Pos[:,0]
y_galaxy = Pos[:,1]
z_galaxy = Pos[:,2]
r_galaxy = Rhalf[GoodIDs]

galaxy_path     = os.path.abspath(os.path.join(os.getcwd(), '..')) + '/'
figure_path     = galaxy_path + 'figures/'
text_files_path = galaxy_path + 'text_files/'
Filename        = 'snapshot_067.hdf5'
SnapshotFile    = h5py.File(Filename,'r')

# Center coordinates of galaxy
xC         = 39937.98       
yC         = 34857.863
zC         = 37441.234

# Gas
Pos0       = SnapshotFile['PartType0/Coordinates'].value 
Vel0       = SnapshotFile['PartType0/Velocities'].value  
SFR_gas    = SnapshotFile['PartType0/StarFormationRate'].value
V0         = SnapshotFile['PartType0/Potential'].value  
Masses_gas = SnapshotFile['PartType0/Masses'].value    
x_gas      = Pos0[:,0]
y_gas      = Pos0[:,1]
z_gas      = Pos0[:,2]
vx_gas     = Vel0[:,0]
vy_gas     = Vel0[:,1]
vz_gas     = Vel0[:,2] 

# Dark Matter
Pos1       = SnapshotFile['PartType1/Coordinates'].value 
Vel1       = SnapshotFile['PartType1/Velocities'].value  
V1         = SnapshotFile['PartType1/Potential'].value  
x_dm       = Pos1[:,0]
y_dm       = Pos1[:,1]
z_dm       = Pos1[:,2]
vx_dm      = Vel1[:,0]
vy_dm      = Vel1[:,1]
vz_dm      = Vel1[:,2] 

# Stars
Pos4        = SnapshotFile['PartType4/Coordinates'].value 
Vel4        = SnapshotFile['PartType4/Velocities'].value  
V4          = SnapshotFile['PartType4/Potential'].value  
Masses_star = SnapshotFile['PartType4/Masses'].value   
x_star      = Pos4[:,0]
y_star      = Pos4[:,1]
z_star      = Pos4[:,2]
vx_star     = Vel4[:,0]
vy_star     = Vel4[:,1]
vz_star     = Vel4[:,2] 

# AGN
Pos5        = SnapshotFile['PartType5/Coordinates'].value 
Vel5        = SnapshotFile['PartType5/Velocities'].value  
V5          = SnapshotFile['PartType5/Potential'].value  
Masses_agn  = SnapshotFile['PartType5/Masses'].value   
x_agn       = Pos5[:,0]
y_agn       = Pos5[:,1]
z_agn       = Pos5[:,2]
vx_agn      = Vel5[:,0]
vy_agn      = Vel5[:,1]
vz_agn      = Vel5[:,2] 

# Switches for figures
Fig_rectangle_IDs 						 = 0
centered_stars_2dhist     			     = 0
panel 									 = 1
rectangles_and_stars_2dhist				 = 0
centered_gas_2dhist          	 		 = 0
centered_galaxies_stars		  			 = 0
centered_galaxies_gas        			 = 0
rectangular_region_stars_2dhist 		 = 0

x = x_star-xC
y = y_star-yC

degrees = 10
radians = math.radians(degrees)

# bridge rectangle
# slopes
Ba1, Ba3 = np.tan(radians), np.tan(radians - math.radians(90))
Ba2, Ba4 = Ba1, Ba3

# BL, BR, TL, TR: Bottom Left, Bottom Right, Top Left, Top Right
BxBL, ByBL, BxBR, ByBR = -10, -25, 19.544, -19.791
BxTL, ByTL, BxTR, ByTR = -12, -13, 17.5, -8

# intersections
Bb1 = ByBL - Ba1*BxBL
Bb2 = ByTR - Ba2*BxTR
Bb3 = ByTL - Ba3*BxTL
Bb4 = ByTR - Ba4*BxTR

RectangleIDs = np.where((y > Ba1*x + Bb1)*(y < Ba2*x + Bb2)*(y > Ba3*x + Bb3)*(y < Ba4*x + Bb4))
hB, wB       = 12, 30 # height and width
AB           = hB*wB  # area = 360
ratioB       = sum(SFR_gas[RectangleIDs]) / AB # 0.000293586789828


# areas of galaxy1 and galaxy2 rectangles
hG, wG = 24, 30 # heights and widths
AG     = hG*wG  # areas (each) = 720

def find_point_on_line(a,b,d,m):
	'''
	find new point on straight line with slope m, with a known point (a,b), at a distance to that point, d.
	'''
	k = d/np.sqrt(1 + m**2)
	x = a + k
	y = b + k*m
	return x, y

# Galaxy 1 rectangle (lower)
G1xBL, G1yBL = 0, -47
G1xBR, G1yBR = find_point_on_line(G1xBL, G1yBL, wG, Ba1)
G1xTL, G1yTL = find_point_on_line(G1xBL, G1yBL, -hG, Ba3)
G1xTR, G1yTR = find_point_on_line(G1xTL, G1yTL, wG, Ba2)

# intersections
G1b1 = G1yBL - Ba1*G1xBL
G1b2 = G1yTR - Ba2*G1xTR
G1b3 = G1yTL - Ba3*G1xTL
G1b4 = G1yTR - Ba4*G1xTR

G1IDs = np.where((y > Ba1*x + G1b1)*(y < Ba2*x + G1b2)*(y > Ba3*x + G1b3)*(y < Ba4*x + G1b4))
ratioG1 = sum(SFR_gas[G1IDs]) / AG # 0.000649032024383

# Galaxy 2 rectangle (higher)
G2xBL, G2yBL = -12,-13
G2xBR, G2yBR = find_point_on_line(G2xBL, G2yBL, wG, Ba1)
G2xTL, G2yTL = find_point_on_line(G2xBL, G2yBL, -hG, Ba3)
G2xTR, G2yTR = find_point_on_line(G2xTL, G2yTL, wG, Ba2)

# intersections
G2b1 = G2yBL - Ba1*G2xBL
G2b2 = G2yTR - Ba2*G2xTR
G2b3 = G2yTL - Ba3*G2xTL
G2b4 = G2yTR - Ba4*G2xTR

G2IDs = np.where((y > Ba1*x+G2b1)*(y < Ba2*x+G2b2)*(y > Ba3*x+G2b3)*(y < Ba4*x+G2b4))
ratioG2 = sum(SFR_gas[G2IDs]) / AG # 0.0190599170674

# Find the mass of stars in the bridge and in the galaxies:
'''
print 'Masses_star[RectangleIDs] = ', Masses_star[RectangleIDs] 
print 'Masses_star[G1IDs] = ', Masses_star[G1IDs] 
print 'Masses_star[G2IDs] = ', Masses_star[G2IDs] 

print 'len(Masses_star[RectangleIDs]) = ', len(Masses_star[RectangleIDs]) 
print 'len(Masses_star[G1IDs]) = ', len(Masses_star[G1IDs]) 
print 'len(Masses_star[G2IDs]) = ', len(Masses_star[G2IDs]) 
'''



if Fig_rectangle_IDs:
	f,(ax) = plt.subplots(1,1,figsize=(6,6))
	plt.plot([-10, 19.544], [-25, -19.791],color='k',marker='o')
	plt.plot([-12, 17.5], [-13, -8], color='b',marker='o')
	plt.plot([-12, -10], [-13, -25],color='r',marker='o')
	plt.plot([17.5, 19.544], [-8, -19.791], color='g',marker='o')
	a1 = 0.176313295424
	a2 = 0.169491525424
	a3 = -6
	a4 = -5.76859099804
	b1 = -23.2368670458
	b2 = -10.9661016949
	b3 = -85
	b4 = 92.9503424657
	x = np.random.uniform(-13,20,1000)
	y = np.random.uniform(-26,-7,1000)
	RectangleIDs = np.where((y > a1*x+b1)*(y < a2*x+b2)*(y > a3*x+b3)*(y < a4*x+b4))
	plt.plot(x[RectangleIDs],y[RectangleIDs],marker='.')
	plt.grid(True)
	ax.set_xlabel('x',fontsize=10)
	ax.set_ylabel('y',fontsize=10)
	ax.set_title(r'Tilted Rectangle and selected points',fontsize=15)
	f.savefig(figure_path + 'Fig_rectangle_IDs.png')

if centered_stars_2dhist:
	f,(ax1) = plt.subplots(1,1,figsize=(6,6))
	f.subplots_adjust(hspace=0,wspace=0)
	plt.hist2d(x_star-xC,y_star-yC, bins=(100,100), range=np.array([(-80,80), (-80,80)]), weights = Masses_star/((160.0/100.0)**2), norm=LogNorm())
	plt.colorbar()
	for i in range(len(r_galaxy)):
		circ = plt.Circle((x_galaxy[i]-xC,y_galaxy[i]-yC), radius=r_galaxy[i], color='black', fill=False,zorder=3)
		ax1.add_patch(circ)

	ts = ax1.transData
	coords = [-10,-25]
	tr = mpl.transforms.Affine2D().rotate_deg_around(coords[0],coords[1], 10)
	t = tr + ts
	rect1 = mpl.patches.Rectangle((-10,-25),30,12,linewidth=1,edgecolor='b',facecolor='none',transform=t)
	ax1.add_patch(rect1)

	tsg1 = ax1.transData
	coordsg1 = [-12,-13]
	trg1 = mpl.transforms.Affine2D().rotate_deg_around(coordsg1[0],coordsg1[1], 10)
	tg1 = trg1 + tsg1
	rectg1 = mpl.patches.Rectangle((-12,-13),30,24,linewidth=1,edgecolor='r',facecolor='none',transform=tg1)
	ax1.add_patch(rectg1)

	tsg2 = ax1.transData
	coordsg2 = [0,-47]
	trg2 = mpl.transforms.Affine2D().rotate_deg_around(coordsg2[0],coordsg2[1], 10)
	tg2 = trg2 + tsg2
	rectg2 = mpl.patches.Rectangle((0,-47),30,24,linewidth=1,edgecolor='g',facecolor='none',transform=tg2)
	ax1.add_patch(rectg2)

	ax1.set_xlabel(r'$x-x_c$',fontsize=10)
	ax1.set_ylabel(r'$y-y_c$',fontsize=10)
	ax1.set_title(r'Centralized stars',fontsize=15)
	ax1.set_xlim(-80,80)
	ax1.set_ylim(-80,80)
	f.savefig(figure_path + 'centered_stars_2dhist.png')


if panel:
	f,(ax1, ax2, ax3, ax4) = plt.subplots(1,4,figsize=(6,6))
	f.subplots_adjust(hspace=0,wspace=0)

	ax1.hist2d(x_star-xC,y_star-yC, bins=(100,100), range=np.array([(-80,80), (-80,80)]), weights = Masses_star/((160.0/100.0)**2), normed = LogNorm())
	ax1.set_ylabel('y_star-yC')

	ax2.set_title('Bridge')
	ax2.hist2d(x_star[RectangleIDs]-xC,y_star[RectangleIDs]-yC, bins=(100,100), range=np.array([(-80,80), (-80,80)]), weights = Masses_star[RectangleIDs]/((160.0/100.0)**2), normed = LogNorm())
	ax2.tick_params(
    	axis='y',       
    	which='both',    
    	left='off',      
    	right='off',         
    	labelleft='off')

	ax3.set_title('Galaxy 1')
	ax3.hist2d(x_star[G1IDs]-xC,y_star[G1IDs]-yC, bins=(100,100), range=np.array([(-80,80), (-80,80)]), weights = Masses_star[G1IDs]/((160.0/100.0)**2), normed = LogNorm())
	ax3.tick_params(
    	axis='y',       
    	which='both',    
    	left='off',      
    	right='off',         
    	labelleft='off')
	ax3.set_xlabel('x_star-xC')

	ax4.set_title('Galaxy 2')
	ax4.hist2d(x_star[G2IDs]-xC,y_star[G2IDs]-yC, bins=(100,100), range=np.array([(-80,80), (-80,80)]), weights = Masses_star[G2IDs]/((160.0/100.0)**2), normed = LogNorm())
	ax4.tick_params(
    	axis='y',       
    	which='both',    
    	left='off',      
    	right='off',         
    	labelleft='off')

	ts = ax1.transData
	coords = [-10,-25]
	tr = mpl.transforms.Affine2D().rotate_deg_around(coords[0],coords[1], 10)
	t = tr + ts
	rect1 = mpl.patches.Rectangle((-10,-25),30,12,linewidth=1,edgecolor='b',facecolor='none',transform=t)
	ax1.add_patch(rect1)

	tsg1 = ax1.transData
	coordsg1 = [-12,-13]
	trg1 = mpl.transforms.Affine2D().rotate_deg_around(coordsg1[0],coordsg1[1], 10)
	tg1 = trg1 + tsg1
	rectg1 = mpl.patches.Rectangle((-12,-13),30,24,linewidth=1,edgecolor='r',facecolor='none',transform=tg1)
	ax1.add_patch(rectg1)

	tsg2 = ax1.transData
	coordsg2 = [0,-47]
	trg2 = mpl.transforms.Affine2D().rotate_deg_around(coordsg2[0],coordsg2[1], 10)
	tg2 = trg2 + tsg2
	rectg2 = mpl.patches.Rectangle((0,-47),30,24,linewidth=1,edgecolor='g',facecolor='none',transform=tg2)
	ax1.add_patch(rectg2)


	ts = ax2.transData
	coords = [-10,-25]
	tr = mpl.transforms.Affine2D().rotate_deg_around(coords[0],coords[1], 10)
	t = tr + ts
	rect1 = mpl.patches.Rectangle((-10,-25),30,12,linewidth=1,edgecolor='b',facecolor='none',transform=t)
	ax2.add_patch(rect1)

	tsg2 = ax3.transData
	coordsg2 = [0,-47]
	trg2 = mpl.transforms.Affine2D().rotate_deg_around(coordsg2[0],coordsg2[1], 10)
	tg2 = trg2 + tsg2
	rectg2 = mpl.patches.Rectangle((0,-47),30,24,linewidth=1,edgecolor='g',facecolor='none',transform=tg2)
	ax3.add_patch(rectg2)

	tsg1 = ax4.transData
	coordsg1 = [-12,-13]
	trg1 = mpl.transforms.Affine2D().rotate_deg_around(coordsg1[0],coordsg1[1], 10)
	tg1 = trg1 + tsg1
	rectg1 = mpl.patches.Rectangle((-12,-13),30,24,linewidth=1,edgecolor='r',facecolor='none',transform=tg1)
	ax4.add_patch(rectg1)
	
	f.savefig(figure_path + 'panel.png')




if rectangles_and_stars_2dhist:
	f,(ax1) = plt.subplots(1,1,figsize=(6,6))
	f.subplots_adjust(hspace=0,wspace=0)
	plt.hist2d(x_star[RectangleIDs]-xC,y_star[RectangleIDs]-yC, bins=(100,100), range=np.array([(-80,80), (-80,80)]), weights = Masses_star[RectangleIDs]/((160.0/100.0)**2), normed = LogNorm())
	#plt.hist2d(x_star[G1IDs]-xC,y_star[G1IDs]-yC, bins=(100,100), range=np.array([(-80,80), (-80,80)]), weights = Masses_star[G1IDs]/((160.0/100.0)**2), normed = LogNorm())
	#plt.hist2d(x_star[G2IDs]-xC,y_star[G2IDs]-yC, bins=(100,100), range=np.array([(-80,80), (-80,80)]), weights = Masses_star[G2IDs]/((160.0/100.0)**2), normed = LogNorm())
	plt.colorbar()

	ts = ax1.transData
	coords = [-10,-25]
	tr = mpl.transforms.Affine2D().rotate_deg_around(coords[0],coords[1], 10)
	t = tr + ts
	rect1 = mpl.patches.Rectangle((-10,-25),30,12,linewidth=1,edgecolor='b',facecolor='none',transform=t)
	ax1.add_patch(rect1)

	tsg1 = ax1.transData
	coordsg1 = [-12,-13]
	trg1 = mpl.transforms.Affine2D().rotate_deg_around(coordsg1[0],coordsg1[1], 10)
	tg1 = trg1 + tsg1
	rectg1 = mpl.patches.Rectangle((-12,-13),30,24,linewidth=1,edgecolor='r',facecolor='none',transform=tg1)
	ax1.add_patch(rectg1)

	tsg2 = ax1.transData
	coordsg2 = [0,-47]
	trg2 = mpl.transforms.Affine2D().rotate_deg_around(coordsg2[0],coordsg2[1], 10)
	tg2 = trg2 + tsg2
	rectg2 = mpl.patches.Rectangle((0,-47),30,24,linewidth=1,edgecolor='g',facecolor='none',transform=tg2)
	ax1.add_patch(rectg2)

	ax1.set_xlabel(r'$x-x_c$',fontsize=10)
	ax1.set_ylabel(r'$y-y_c$',fontsize=10)
	ax1.set_title(r'Centralized stars',fontsize=15)
	ax1.set_xlim(-80,80)
	ax1.set_ylim(-80,80)
	f.savefig(figure_path + 'rectangles_and_stars_2dhist.png')

if centered_gas_2dhist:
	f,(ax1) = plt.subplots(1,1,figsize=(6,6))
	f.subplots_adjust(hspace=0,wspace=0)
	plt.hist2d(x_gas-xC,y_gas-yC, bins=(100,100), range=np.array([(-80,80), (-80,80)]), weights = SFR_gas/((160.0/100.0)**2), norm=LogNorm())
	plt.colorbar()
	for i in range(len(r_galaxy)):
		circ = plt.Circle((x_galaxy[i]-xC,y_galaxy[i]-yC), radius=r_galaxy[i], color='red', fill=False,zorder=2)
		ax1.add_patch(circ)
	ax1.set_xlabel(r'$x-x_c$',fontsize=10)
	ax1.set_ylabel(r'$y-y_c$',fontsize=10)
	ax1.set_title(r'Centralized gas',fontsize=15)
	ax1.set_xlim(-80,80)
	ax1.set_ylim(-80,80)
	f.savefig(figure_path + 'centered_gas_2dhist.png')

if centered_galaxies_stars:
	f,(ax1) = plt.subplots(1,1,figsize=(6,6))
	f.subplots_adjust(hspace=0,wspace=0)
	ax1.plot(x_galaxy-xC,y_galaxy-yC,'or',ms=1,zorder=2)
	ax1.plot(x_star-xC,y_star-yC,'ob',ms=1,zorder=1)
	for i in range(len(r_galaxy)):
		circ = plt.Circle((x_galaxy[i]-xC,y_galaxy[i]-yC), radius=r_galaxy[i], color='black', fill=False,zorder=3)
		ax1.add_patch(circ)	
	ax1.set_xlabel(r'$x-x_c$',fontsize=10)
	ax1.set_ylabel(r'$y-y_c$',fontsize=10)
	ax1.set_title(r'Centralized galaxies, Rhalf and stars',fontsize=15)
	ax1.set_xlim(-300,200)
	ax1.set_ylim(-300,200)
	f.savefig(figure_path + 'centered_galaxies_stars.png')

if centered_galaxies_gas:
	f,(ax1) = plt.subplots(1,1,figsize=(6,6))
	f.subplots_adjust(hspace=0,wspace=0)
	ax1.plot(x_galaxy-xC,y_galaxy-yC,'or',ms=1,zorder=2)
	ax1.plot(x_gas-xC,y_gas-yC,'ob',ms=1,zorder=1)
	for i in range(len(r_galaxy)):
		circ = plt.Circle((x_galaxy[i]-xC,y_galaxy[i]-yC), radius=r_galaxy[i], color='black', fill=False,zorder=3)
		ax1.add_patch(circ)	
	ax1.set_xlabel(r'$x-x_c$',fontsize=10)
	ax1.set_ylabel(r'$y-y_c$',fontsize=10)
	ax1.set_title(r'Centralized galaxies, Rhalf and gas',fontsize=15)
	ax1.set_xlim(-300,200)
	ax1.set_ylim(-300,200)
	f.savefig(figure_path + 'centered_galaxies_gas.png')

#RectangleIDs = np.where((y_star-yC > a1*(x_star-xC)+b1)*(y_star - yC < a2*(x_star-xC)+b2)*(y_star-yC > a3*(x_star-xC)+b3)*(y_star-yC < a4*(x_star-xC)+b4))

if rectangular_region_stars_2dhist:
	f,(ax) = plt.subplots(1,1,figsize=(6,6))
	#plt.hist2d(x_star[RectangleIDs]-xC,y_star[RectangleIDs]-yC, bins=(100,100), range=np.array([(-80,80), (-80,80)]), weights = Masses_star/((160.0/100.0)**2), norm=LogNorm())
	plt.hist2d(x_star[RectangleIDs]-xC,y_star[RectangleIDs]-yC, bins=(100,100), range=np.array([(-80,80), (-80,80)]), weights = Masses_star[RectangleIDs]/((160.0/100.0)**2), normed = True)	
	plt.colorbar()
	#ts = ax.transData
	#coords = [-10,-25]
	#tr = mpl.transforms.Affine2D().rotate_deg_around(coords[0],coords[1], 10)
	#t = tr + ts
	#rect1 = mpl.patches.Rectangle((-10,-25),30,12,linewidth=1,edgecolor='b',facecolor='none',transform=t)
	#ax.add_patch(rect1)
	ax.set_xlabel(r'$x-x_c$',fontsize=10)
	ax.set_ylabel(r'$y-y_c$',fontsize=10)
	ax.set_title(r'Stars between merging galaxies',fontsize=15)
	ax.set_xlim(-80,80)
	ax.set_ylim(-80,80)
	f.savefig(figure_path + 'rectangular_region_stars_2dhist.png')

plt.show()


