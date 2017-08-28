import numpy as np
import pylab as pl
import math
data=np.loadtxt('WHM_Results_with_drugs.txt')
b=0
#Loop to check how long episode lasts, and get max parasitaemia 
for t in range(0, 399):
    if data[t,1] > 0.00001:
    	T = t;
    if data[t,2] > b:
    	b = data[t,2]
    if data[t,1] < 0.00001:
    	break

#now load PK data 
dataPK=np.loadtxt('WHM_Results_with_PKPD.txt')
dimm=dataPK.shape[0] - 1 #I've placed some useful variables in the last row
#print(dimm)

fev = dataPK[dimm,3] #Timing of first fever, in increments of dt
dt = dataPK[dimm,2] #dt (in hours)
delta = dt * (1/24.0) * dataPK[dimm,4] #Start of treatment, converted into days
#print(fev)
#print(delta)
pl.plot(data[:,0],data[:,1], marker ='o')
pl.plot(data[:,0],data[:,2], marker ='o')

#Microscopy threshold, higher (more realistic?) than used in Molineaux et al., Parasitology 122 379 (2001)
pl.plot([0,(2 * T) + 6],[40,40])
pl.plot([0,(2 * T) + 6],[fev,fev],linestyle='-.')#Fever Threshold
pl.plot([delta,delta],[0.00001 , b * 1.4])#Treatment starts
pl.plot([delta+1,delta+1],[0.00001 , b * 1.4])#24 hrs after treatment starts
pl.plot([delta+28,delta+28],[0.00001 , b * 1.4])#28 day check, if we get that far!
pl.plot(dataPK[0:dimm,0],dataPK[0:dimm,1])#Artemether
pl.plot(dataPK[0:dimm,0],dataPK[0:dimm,3])#Lumefantrine

pl.xlabel('Time (Days)')
pl.ylabel('Parasitaemia [PRBCs / $\mu L$]')
pl.title('Falciparum infection treated with Artemether & Lumefantrine')
pl.yscale('log')

pl.xlim(0.0, 6 + (2 * T))
pl.ylim(0.00001, b * 1.5) #Scaling of 1.5 just provides some space above 1st peak

pl.savefig("WHM_one_run_withPK.pdf", format='pdf')
pl.show()