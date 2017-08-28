import numpy as np
import pylab as pl
import math
data=np.loadtxt('WHM_Results.txt')#Load output from C++ model
b=0
for t in range(0, 399):
    if data[t,1] > 0.00001:
    	T = t;
    if data[t,1] > b:
    	b = data[t,1]
    if data[t,1] < 0.00001:
    	break

#print(b)
#print(T)
pl.plot(data[:,0],data[:,1])
pl.plot([0,(2 * T) + 6],[10,10])#Adds the microscopy threshold.
pl.xlabel('Time (Days)')
pl.ylabel('Parasitaemia')
pl.yscale('log')
pl.xlim(0.0, 6 + (2 * T))
pl.ylim(0.00001, b * 1.2) #Scaling of 1.2 just provides some space above 1st peak
pl.savefig("WHM_one_run.pdf", format='pdf')
pl.show()