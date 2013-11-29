import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def readInfo():
	
	infoFile = open("../data/info.dat", 'r')
	global N, R0, dt, Tmax, n
	info = infoFile.readline().split()
	N, R0, dt, Tmax = int(info[0]), float(info[1]), float(info[2]), float(info[3])
	n = int(round(Tmax / dt))

	infoFile.close()
	return None
readInfo()

def readPositions():
	
	global positions, time
	positions = np.zeros((N,n,3))
	time = np.zeros(n)
	
	for i in range(N):
		inFile = open("../data/objects/obj%s.dat" % i, 'r')
		inFile.readline()
		j = 0
		for line in inFile:
			info = line.split()
			time[j] = float(info[0])
			positions[i,j,0] = float(info[1])
			positions[i,j,1] = float(info[2])
			positions[i,j,2] = float(info[3])
			j += 1
		inFile.close()	
	print time
	return None
readPositions()

def readEnergy():

	global totalEnergy
#	totalEnergy = np.zeros(n)
	kineticEnergy = np.zeros(n)
	potentialEnergy = np.zeros(n)
	
	for i in range(N):
		inFile = open("../data/conservations/energy/obj%s.dat" % i, 'r')
		inFile.readline()
		j = 0
		for line in inFile:
			info = line.split()
			objectKinEn = float(info[1])
			objectPotEn = float(info[2])
			kineticEnergy[j] += objectKinEn
			potentialEnergy[j] += objectPotEn
			objectKinEn = 0
			objectPotEn = 0
			j += 1

	potentialEnergy /= 2.
	totalEnergy = kineticEnergy + potentialEnergy
	return None
readEnergy()
	
		
def animateCluster():

	plt.ion()
	fig = plt.figure()
	fig.suptitle('$\mathrm{Evolution \ of \ star \ cluster}$. $N = %g$, $R_0 = %g$, $dt = %g$' % (N,R0,dt))
	
	ax3 = fig.add_subplot(111, projection='3d')
	dots = []
	
	for i in range(N):
		x,y,z = positions[i,0]
		dots.append(ax3.plot([x],[y],[z],'o')[0])

	ax3.set_xlim([-R0,R0])
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
	plt.ioff()
	plt.show()
	return None

def plotEnergy():

	fig = plt.figure()
	fig.suptitle('$\mathrm{Total \ energy \ for \ the \ cluster}$', fontsize='14')

	ax = fig.add_subplot(111)
	ax.set_title('$N = %g$, $R_0 = %g$, $dt = %g$' % (N,R0,dt))
	ax.set_xlabel('$t \\ [\\tau_{\mathrm{crunch}}]$', fontsize='14')
	ax.set_ylabel('$\mathrm{Total \ energy}$', fontsize='14')
	ax.plot(time,totalEnergy)
	ax.grid('on')
	plt.show()
	return None
plotEnergy()
animateCluster()
