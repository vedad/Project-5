import matplotlib.pyplot as plt
import numpy as np

def readInfo():

	inFile = open("../data/info.dat", 'r')
	global dt, tMax
	for line in inFile:
		info = line.split()
		dt = float(info[2])
		tMax = float(info[3])

	return None
readInfo()

def read():

	inFile = open("../data/objects/Star#2.dat",'r')

	global _t; _t = []
	global _xB; _xB = []
	global _yB; _yB = []
	global _zB; _zB = []

	inFile.readline()
	for line in inFile:
		columns = line.split()
		_t.append(float(columns[0]))
		_xB.append(float(columns[4]))
		_yB.append(float(columns[5]))
		_zB.append(float(columns[6]))

	inFile.close()

	_xB, _yB, _zB = np.asarray(_xB), np.asarray(_yB), np.asarray(_zB)

	_yB = _yB / 2.
	_zB = _xB + _yB

def read2():

	inFile = open("../data/objects/Star#1.dat", 'r')

	global _xA; _xA = []
	global _yA; _yA = []
	global _zA; _zA = []

	inFile.readline()
	for line in inFile:
		columns = line.split()
		_xA.append(float(columns[4]))
		_yA.append(float(columns[5]))
		_zA.append(float(columns[6]))
	inFile.close()

	_xA, _yA, _zA = np.asarray(_xA), np.asarray(_yA), np.asarray(_zA)
	_yA = _yA / 2.
	_zA = _xA + _yA
	return None

def read3():
	
	global totalEnergy; totalEnergy = []
	global energyError
	inFile = open("../data/energy/cluster/clusterEnergy.dat", 'r')
	counter = 0
	for line in inFile:
		columns = line.split()
		if (counter == 0):
			realEnergy = float(columns[4])
		totalEnergy.append(float(columns[4]))
		counter += 1
	inFile.close()
	totalEnergy = np.asarray(totalEnergy)
	energyError = realEnergy - totalEnergy
	return None
	

read()
read2()
read3()
_z = _zA + _zB
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('$T = %g, \\ dt = %g$' % (tMax,dt), fontsize='26')
ax.plot(_t,energyError)
ax.set_xlabel('$t \\ \mathrm{[year]}$', fontsize='26')
ax.set_ylabel('$\mathrm{Energy \ error}$', fontsize='26')
ax.grid('on')

plt.show()
