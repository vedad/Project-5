#!/usr/bin/env python
"""
Created on Tir 22 Oct 2013

Script for reading and plotting data from celestial objects outputted from C++
simulations of the solar system.

"""
from matplotlib import pyplot as plt
import numpy as np
import os,sys
from mpl_toolkits.mplot3d import Axes3D

"""Constants"""

OBJECTS_PATH = '../data/objects/'

"""Classes"""

class ObjectData(object):
    """
    Contains the data for one single object.
    """
    def __init__(self, datafile):
        """
        Name for object is found in datafile.

        @param datafile The file containing positional data.
        """
        self.name = self._findName(datafile)
        self._read(datafile)

    def _read(self, datafile):
        """
        Reads datafile and puts into arrays.
        """
        inData = open(datafile, 'r')
        x = []
        y = []
        z = []
		
        # First read till end of header
        inData.readline()
        for line in inData:
            columns = line.split()
            x.append(float(columns[1]))
            y.append(float(columns[2]))
            z.append(float(columns[3]))

        self.x,self.y,self.z = np.asarray(x), np.asarray(y), np.asarray(z)
        inData.close()

    def _findName(self, datafile):
        """
        Finds the name of the object from the datafile.

        @param datafile The file containing data.
        """
        inData = open(datafile, 'r')

        for line in inData:
            if line.startswith('Positions'):
                inData.close()
                return line.strip().split(':')[1].strip()
        
        # Rapid err msg, no time dude
        print 'If you see this the datafile has wrong syntax.'
        sys.exit(1)

"""Methods"""

if __name__ == '__main__':
    """
    Running script will loop through all files in objects folder and add them
    to the plot along with name as label.
    """
    files = os.listdir(OBJECTS_PATH)
    labels = ["Sirius A", "Sirius B"]

    fig = plt.figure()
    fig.suptitle('$\mathrm{Newtonian \ two-body \ problem \ with \ Leapfrog}$', fontsize='26')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title('$T = 200$, $dt=1.0$', fontsize='26')
    ax.set_xlabel('$x \\ \mathrm{[AU]}$', fontsize='26')
    ax.set_ylabel('$y \\ \mathrm{[AU]}$', fontsize='26')
    ax.set_zlabel('$z \\ \mathrm{[AU]}$', fontsize='26')
    ax.grid('on')

    i=0	
    for datafile in files:
        data = ObjectData(os.path.join(OBJECTS_PATH,datafile))
        ax.plot(data.x,data.y,data.z,label=labels[i])
        i += 1

    ax.legend(loc='best')

    plt.show()
