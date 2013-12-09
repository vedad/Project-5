import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def readInfo():
	"""
	Reads information file to get parameters to set as plot title.
	"""
	infoFile = open("../data/info.dat", 'r')
	global N, R0, dt, Tmax, n
	info = infoFile.readline().split()
	N, R0, dt, Tmax = int(info[0]), float(info[1]), float(info[2]), float(info[3])
	n = int(round(Tmax / dt)) 
	infoFile.close()
	global properties; properties = '$N = %g$, $R_0 = %g$, $dt = %g$' 
	return None
readInfo()

def readTime():
	"""
	Reads the time and makes an array.
	"""
	global time; time = np.zeros(n)
	inFile = open("../data/energy/cluster/clusterEnergy.dat", 'r')
	j = 0
	for line in inFile:
		info = line.split()
		time[j] = float(info[0])
		j += 1
	inFile.close()
	return None
readTime()
	
def readPositions():
	"""
	Reads the positions of every object.
	"""
	global positions
	positions = np.zeros((N,n,3))
	
	for i in range(N):
		inFile = open("../data/objects/obj%s.dat" % i, 'r')
		inFile.readline()
		j = 0
		for line in inFile:
			info = line.split()
			positions[i,j,0] = float(info[1])
			positions[i,j,1] = float(info[2])
			positions[i,j,2] = float(info[3])
			j += 1
		inFile.close()	
	return None
readPositions()

def readEnergy():
	"""
	Reads the energy of every object.
	"""
	global averageKineticEnergy, averagePotentialEnergy
	kineticEnergy = np.zeros(n)
	potentialEnergy = np.zeros(n)
	boundKinEn = np.zeros(n)
	boundPotEn = np.zeros(n)
	totalEnergy = np.zeros(n)
	boundCounter = np.zeros(n)

	for i in range(N):
		inFile = open("../data/objects/obj%s.dat" % i, 'r')
		inFile.readline()
		j = 0
		for line in inFile:
			info = line.split()
			objectKinEn = float(info[4])
			objectPotEn = float(info[5])
			objectTotalEn = float(info[6])
			bound[i,j] = int(info[7])	
			if (objectTotalEn < 0):
				boundKinEn[j] += objectKinEn
				boundPotEn[j] += objectPotEn
				boundCounter[j] += 1
			if (objectTotalEn >= 0):
				boundPotEn[j] -= objectPotEn
			kineticEnergy[j] += objectKinEn
			potentialEnergy[j] += objectPotEn
			objectKinEn = 0
			objectPotEn = 0
			j += 1

	boundPotEn /= 2.
	averageKineticEnergy = boundKinEn / boundCounter
	averagePotentialEnergy = boundPotEn / boundCounter
	
	return None
readEnergy()
	
def readClusterEnergy():

	inFile = open("../data/energy/cluster/clusterEnergy.dat", 'r')
	global kineticEnergy, potentialEnergy, clusterEnergy, boundClusterEnergy
	kineticEnergy = np.zeros(n)
	potentialEnergy = np.zeros(n)
	clusterEnergy = np.zeros(n)
	boundClusterEnergy = np.zeros(n)
		
	i = 0;
	for line in inFile:
		info = line.split()
		kineticEnergy[i] = info[1]
		potentialEnergy[i] = info[2]
		clusterEnergy[i] = info[3]
		boundClusterEnergy[i] = info[4]
		i += 1
	
	inFile.close()

	return None
readClusterEnergy()

def readBoundObjects():

	bound = np.zeros((N,n))
	global boundObjects; boundObjects = np.zeros(n)

	for i in range(N):
		inFile = open("../data/objects/obj%s.dat" % i, 'r')
		inFile.readline()
		j = 0
		for line in inFile:
			info = line.split()
			bound[i,j] = int(info[7])
			j += 1
	
	inFile.close()
	for i in range(n):
		boundObjects[i] = np.sum(bound[:,i])
	
	return None
readBoundObjects()
		
def animateCluster():

	fig = plt.figure()
	fig.suptitle('$\mathrm{Evolution \ of \ star \ cluster}$. $N = %g$, $R_0 = %g$, $dt = %g$' % (N,R0,dt))
	
	ax3 = fig.add_subplot(111, projection='3d')
	dots = []
	
	for i in range(N):
		x,y,z = positions[i,0]
		dots.append(ax3.plot([x],[y],[z],'o')[0])

	ax3.set_xlim([-1.5 * R0, 1.5 * R0])
	ax3.set_ylim([-R0,R0])
	ax3.set_zlim([-R0,R0])

	plt.draw()

	for i in range(1,n):
		ax3.set_title('$t = %g$' % time[i])
		for j in range(N):
			x,y,z = positions[j,i]
			dots[j].set_data([x],[y])
			dots[j].set_3d_properties([z])
		plt.draw()
	return None

def plotEnergy():

	fig = plt.figure()
	fig.suptitle('$\mathrm{Total \ energy \ for \ the \ cluster}$', fontsize='14')

	ax = fig.add_subplot(111)
	ax.set_title('$N = %g$, $R_0 = %g$, $dt = %g$' % (N,R0,dt))
	ax.set_xlabel('$t \\ [\\tau_{\mathrm{crunch}}]$', fontsize='14')
	ax.set_ylabel('$\mathrm{Cluster \ energy} \ [\mathrm{unknown}]$', fontsize='14')
	ax.plot(time,clusterEnergy, label='$E$')
	ax.hold('on')
	ax.plot(time,kineticEnergy, '--', label='$E_k$')
	ax.plot(time,potentialEnergy, '--', label='$E_p$')
	ax.legend(loc='best')
	ax.hold('off')
	ax.grid('on')

	fig2 = plt.figure()
	fig2.suptitle('$\mathrm{Cluster \ energy \ with \ all \ and \ bound \ objects}$', fontsize='14')
	ax2 = fig2.add_subplot(111)
	ax2.set_title('$N = %g$, $R_0 = %g$, $dt = %g$' % (N,R0,dt))
	ax2.set_xlabel('$t \\ [\\tau_{\mathrm{crunch}}]$', fontsize='14')
	ax2.set_ylabel('$\mathrm{Cluster \ energy} \ [\mathrm{unknown}]$', fontsize='14')
	ax2.plot(time,clusterEnergy, label='$\mathrm{All \ objects}$')
	ax2.hold('on')
	ax2.plot(time,boundClusterEnergy, label='$\mathrm{Bound \ objects}$')
	ax2.hold('off')
	ax2.legend(loc='best')
	ax2.grid('on')

	fig3 = plt.figure()
#	fig3.suptitle('$\mathrm{Testing \ the \ virial \ theorem}$', fontsize='14')
	ax3 = fig3.add_subplot(111)
	ax3.set_title(properties % (N,R0,dt), fontsize='24')
	ax3.set_xlabel('$t \\ [{\\tau_{\mathrm{crunch}}}]$', fontsize='24')
	ax3.set_ylabel('$\langle K \\rangle/\langle V \\rangle$', fontsize='24')
	ax3.plot(time,averageKineticEnergy/averagePotentialEnergy)
	ax3.grid('on')

	return None


def plotBoundObjects():
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	fig.suptitle('$\mathrm{Number \ of \ bound \ objects \ over \ time}$', fontsize='14')
	ax.set_title('$N = %g$, $R_0 = %g$, $dt = %g$' % (N,R0,dt))
	ax.set_xlabel('$t \ [\\tau_{\mathrm{crunch}}]$', fontsize='14')
	ax.set_ylabel('$\mathrm{Number \ of \ bound \ objects}$', fontsize='14')
	ax.plot(time,boundObjects)
	ax.grid('on')
	
	return None
	
plt.ion()
plotEnergy()
#animateCluster()
plotBoundObjects()
plt.ioff()

plt.show()
