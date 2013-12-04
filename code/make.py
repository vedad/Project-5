"How does one use this crap?"

import os,sys
os.system('c++ -c -O3 *.cpp')
os.system('c++ -o %s -larmadillo -lblas -llapack *.o' % sys.argv[1])
