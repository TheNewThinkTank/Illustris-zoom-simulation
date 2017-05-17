#import matplotlib
import numpy as np
import h5py, IPython, astropy, pylab
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9
import TrackGalaxy
import sys,copy
from matplotlib.colors import LogNorm
import PIL
from modules.Gas import CalcTandDens
import modules.Image as Image
import modules.Stars as Stars
import modules.Gas as Gas
import modules.BH as BH
from matplotlib import rc
import matplotlib.patches as patches
import matplotlib as mpl
import math

path = '/Users/gustavcollinrasmussen/Desktop/galaxy/mergerzoomanalysis/'
T = TrackGalaxy.TrackGalaxy(np.array([67] ),'1330-3',Dir = path)

#Read in position, SFR, stellar mass, gas mass, dark matter mass of all the galaxies (and subhalos) from the simulation
Attrs = T.GetGroups(67,Attrs=['/Subhalo/SubhaloPos','/Subhalo/SubhaloSFR','/Subhalo/SubhaloMassType','/Subhalo/SubhaloHalfmassRad'])
Pos = Attrs['/Subhalo/SubhaloPos']#in comoving kpc/h
SFR = Attrs['/Subhalo/SubhaloSFR']#in Msun/yr
Mstar = Attrs['/Subhalo/SubhaloMassType'][:,4]*1e10/0.7#in Msun
Mgas = Attrs['/Subhalo/SubhaloMassType'][:,0]*1e10/0.7#in Msun
Mdm = Attrs['/Subhalo/SubhaloMassType'][:,1]*1e10/0.7#in Msun
Rhalf = Attrs['/Subhalo/SubhaloHalfmassRad']#in comoving kpc/h

#there is a lot of subhalos (thousands), but most of them don't have any stars. So let us just pick out the galaxies with at least Mstar>1e7 Msun:
GoodIDs = np.where(Mstar > 1e7)
Pos = Pos[GoodIDs]
SFR = SFR[GoodIDs]
Mstar = Mstar[GoodIDs]
Mgas = Mgas[GoodIDs]
Mdm = Mdm[GoodIDs]

x_galaxy = Pos[:,0]
y_galaxy = Pos[:,1]
z_galaxy = Pos[:,2]
r_galaxy = Rhalf[GoodIDs]

galaxy_path     = '/Users/gustavcollinrasmussen/Desktop/galaxy/'
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
rotate_rectangle_example 				 = 0
rotate_rectangle 						 = 0
points_and_rectangle 					 = 0
points_in_rectangle 					 = 0
Fig_rectangle_IDs 						 = 0

centered_stars_2dhist     			     = 0

centered_gas_2dhist          	 		 = 0
centered_galaxies_stars		  			 = 0
centered_galaxies_gas        			 = 0
rectangular_region_stars_2dhist 		 = 0
rectangular_region_stars_plot 			 = 0


x = x_star-xC
y = y_star-yC

degrees = 10
radians = degrees*math.pi/180

# bridge rectangle
# slopes
Ba1, Ba3 = np.tan(radians), np.tan(radians - 90 * math.pi / 180)
Ba2, Ba4 = Ba1, Ba3

# BL, BR, TL, TR: Bottom Left, Bottom Right, Top Left, Top Right
BxBL, ByBL, BxBR, ByBR = -10, -25, 19.544, -19.791
BxTL, ByTL, BxTR, ByTR = -12, -13, 17.5, -8

# intersections
Bb1 = ByBL - Ba1*BxBL
Bb2 = ByBR - Ba2*BxBR
Bb3 = ByTL - Ba3*BxTL
Bb4 = ByTR - Ba4*BxTR

RectangleIDs = np.where((y > Ba1*x + Bb1)*(y < Ba2*x + Bb2)*(y > Ba3*x + Bb3)*(y < Ba4*x + Bb4))
hB, wB       = 12, 30 # height and width
AB           = hB*wB # area = 360
ratioB       = sum(SFR_gas[RectangleIDs]) / AB # SFR/area = 0.000293586789828

# areas of galaxy1 and galaxy2 rectangles
hG, wG = 24, 30 # heights and widths
AG = hG*wG # areas (each) = 720

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
#print 'G1xBR, G1yBR = ', G1xBR, G1yBR
G1xTL, G1yTL = find_point_on_line(G1xBL, G1yBL, -hG, Ba3)
#print 'G1xTL, G1yTL = ', G1xTL, G1yTL
G1xTR, G1yTR = find_point_on_line(G1xTL, G1yTL, wG, Ba2)
#print 'G1xTR, G1yTR = ', G1xTR, G1yTR

# intersections
G1b1 = G1yBL - Ba1*G1xBL
G1b2 = G1yBR - Ba2*G1xBR
G1b3 = G1yTL - Ba3*G1xTL
G1b4 = G1yTR - Ba4*G1xTR

G1IDs = np.where((y > Ba1*x+G1b1)*(y < Ba2*x+G1b2)*(y > Ba3*x+G1b3)*(y < Ba4*x+G1b4))
ratioG1 = sum(SFR_gas[G1IDs]) / AG # SFR/area = 
print 'ratioG1 = ', ratioG1 # 0

#print type(RectangleIDs) # tuple
#print type(G1IDs) # tuple 

# Galaxy 2 rectangle (higher)
G2xBL, G2yBL = -12,-13
G2xBR, G2yBR = find_point_on_line(G2xBL, G2yBL, wG, Ba1)
G2xTL, G2yTL = find_point_on_line(G2xBL, G2yBL, -hG, Ba3)
G2xTR, G2yTR = find_point_on_line(G1xTL, G1yTL, wG, Ba2)

# intersections
G2b1 = G2yBL - Ba1*G2xBL
G2b2 = G2yBR - Ba2*G2xBR
G2b3 = G2yTL - Ba3*G2xTL
G2b4 = G2yTR - Ba4*G2xTR

G2IDs = np.where((y > Ba1*x+G2b1)*(y < Ba2*x+G2b2)*(y > Ba3*x+G2b3)*(y < Ba4*x+G2b4))
ratioG2 = sum(SFR_gas[G2IDs]) / AG # SFR/area = 
print 'ratioG2 = ', ratioG2



if rotate_rectangle_example:
	f,(ax1) = plt.subplots(1,1,figsize=(6,6))
	f.subplots_adjust(hspace=0,wspace=0)
	ts = ax1.transData
	coords = [1,1]
	tr = mpl.transforms.Affine2D().rotate_deg_around(coords[0],coords[1], 10)
	t = tr + ts
	rec0 = mpl.patches.Rectangle((1,1),3,2,linewidth=1,edgecolor='r',facecolor='none')
	ax1.add_patch(rec0)
	#Rotated rectangle patch
	rect1 = mpl.patches.Rectangle((1,1),3,2,linewidth=1,edgecolor='b',facecolor='none',transform=t)
	ax1.add_patch(rect1)
	# Rectangles lower left corner
	plt.plot([1], [1], marker='o', markersize=3, color="green")
	ax1.set_xlim(0,6)
	ax1.set_ylim(-1,4)
	# Rectangles lower left corner
	plt.plot([1], [1], marker='o', markersize=3, color="green")
	plt.grid(True)
	f.savefig(figure_path + 'rotate_rectangle_example.png')

if rotate_rectangle:
	f,(ax1) = plt.subplots(1,1,figsize=(6,6))
	ts = ax1.transData
	coords = [-10,-25]
	tr = mpl.transforms.Affine2D().rotate_deg_around(coords[0],coords[1], 10)
	t = tr + ts
	rec0 = mpl.patches.Rectangle((-10,-25),30,12,linewidth=1,edgecolor='r',facecolor='none')
	ax1.add_patch(rec0)
	#Rotated rectangle patch
	rect1 = mpl.patches.Rectangle((-10,-25),30,12,linewidth=1,edgecolor='b',facecolor='none',transform=t)
	ax1.add_patch(rect1)
	#plt.plot([150,360],[260,50],'k-')
	plt.plot([19.544],[-19.791], marker='o', markersize=7, color="green")
	plt.plot([-10],[-25], marker='o', markersize=7, color="green")
	# Guesses
	plt.plot([17.5],[-8], marker='o', markersize=7, color="green")
	plt.plot([-12],[-13], marker='o', markersize=7, color="green")
	plt.grid(True)
	ax1.set_xlim(-20,30)
	ax1.set_ylim(-30,-5)
	ax1.set_xlabel('x',fontsize=10)
	ax1.set_ylabel('y',fontsize=10)
	ax1.set_title(r'Rectangle transformation',fontsize=15)
	f.savefig(figure_path + 'rotate_rectangle.png')

if points_and_rectangle:
	f,(ax) = plt.subplots(1,1,figsize=(6,6))
	plt.plot([-10, 19.544], [-25, -19.791],color='k',marker='o',zorder=1)
	plt.plot([-12, 17.5], [-13, -8], color='b',marker='o',zorder=1)
	plt.plot([-12, -10], [-13, -25],color='r',marker='o',zorder=1)
	plt.plot([17.5, 19.544], [-8, -19.791], color='g',marker='o',zorder=1)
	x = np.random.uniform(-13,20,1000)
	y = np.random.uniform(-26,-7,1000)
	plt.plot(x,y,marker='.',zorder=0)
	plt.grid(True)
	ax.set_xlabel('x',fontsize=10)
	ax.set_ylabel('y',fontsize=10)
	ax.set_title(r'Tilted Rectangle and random points',fontsize=15)
	f.savefig(figure_path + 'points_and_rectangle.png')

if points_in_rectangle:
	#corners = [[-10, -25], [19.544, -19.791], [-12, -13], [17.5, -8]]
	f,(ax) = plt.subplots(1,1,figsize=(6,6))
	plt.plot([-10, 19.544], [-25, -19.791],color='k',marker='o')
	plt.plot([-12, 17.5], [-13, -8], color='b',marker='o')
	plt.plot([-12, -10], [-13, -25],color='r',marker='o')
	plt.plot([17.5, 19.544], [-8, -19.791], color='g',marker='o')
	# Find slopes
	# print 'slope of 1st (black) line : ', (-19.791 -(-25))/(19.544 -(-10)) # 0.176313295424
	a1 = 0.176313295424
	#print 'slope of 2nd (blue) line : ', (-8-(-13))/(17.5 -(-12))
	a2 = 0.169491525424
	#print 'slope of 3rd (red) line : ', (-25-(-13))/(-10 -(-12))
	a3 = -6
	#print 'slope of 4th (green) line : ', (-19.791-(-8))/(19.544 -17.5)
	a4 = -5.76859099804
	# Find intersections
	#print 'intersection of 1st (black) line : ', -19.791 - (0.176313295424)*(19.544) # -23.2368670458
	b1 = -23.2368670458
	#print 'intersection of 2nd (blue) line : ', -8 - (0.169491525424)*(17.5)
	b2 = -10.9661016949
	#print 'intersection of 3rd (red) line : ', -25 - (-6)*(-10)
	b3 = -85
	#print 'intersection of 4th (green) line : ', -19.791 - (-5.76859099804)*(19.544)
	b4 = 92.9503424657
	x = np.random.uniform(-13,20,1000)
	y = np.random.uniform(-26,-7,1000)
	#x1 = np.arange(-12., 20., .1)
	for i in range(len(x)):
		if y[i] > a1*x[i]+b1 and y[i] < a2*x[i]+b2 and y[i] > a3*x[i]+b3 and y[i] < a4*x[i]+b4: 
			plt.plot(x[i],y[i],marker='.')
	plt.grid(True)
	ax.set_xlabel('x',fontsize=10)
	ax.set_ylabel('y',fontsize=10)
	ax.set_title(r'Tilted Rectangle and selected points',fontsize=15)
	f.savefig(figure_path + 'points_in_rectangle.png')

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

if rectangular_region_stars_plot:
	f,(ax) = plt.subplots(1,1,figsize=(6,6))
	plt.plot(x_star[RectangleIDs]-xC,y_star[RectangleIDs]-yC)	
	ax.set_xlabel(r'$x-x_c$',fontsize=10)
	ax.set_ylabel(r'$y-y_c$',fontsize=10)
	ax.set_title(r'Stars between merging galaxies',fontsize=15)
	ax.set_xlim(-80,80)
	ax.set_ylim(-80,80)
	f.savefig(figure_path + 'rectangular_region_stars_plot.png')

plt.show()












