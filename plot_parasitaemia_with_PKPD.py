import numpy as np
import pylab as pl
import math
data=np.loadtxt('WHM_Results_with_drugs.txt')
b=0
for t in range(0, 399):
    if data[t,1] > 0.00001:
    	T = t;
    if data[t,2] > b:
    	b = data[t,2]
    if data[t,1] < 0.00001:
    	break

#now load PK data 
dataPK=np.loadtxt('WHM_Results_with_PKPD.txt')
dimm=dataPK.shape[0] - 1
print(dimm)
#T=29 # Not sure what you want for 'T' here. Too long, one can't see drugs well
#print(b)
#print(T)
#test=np.log10(40)
#print(test)
fev = dataPK[dimm,3]
dt = dataPK[dimm,2]
delta = dt * (1/24.0) * dataPK[dimm,4] #Start of treatment, converted into days
print(fev)
print(delta)
pl.plot(data[:,0],data[:,1], marker ='o')
pl.plot(data[:,0],data[:,2], marker ='o')
#pl.plot([0,(2 * T) + 6],[10,10])#Adds the microscopy threshold.
pl.plot([0,(2 * T) + 6],[40,40])#Adds the (more realistic?) microscopy threshold.
pl.plot([0,(2 * T) + 6],[fev,fev],linestyle='-.')#Fever Threshold
pl.plot([delta,delta],[0.00001 , b * 1.4])#Treatment starts
pl.plot([delta+1,delta+1],[0.00001 , b * 1.4])#24 hrs after treatment starts
pl.plot([delta+28,delta+28],[0.00001 , b * 1.4])#28 day check, if we get that far!
pl.plot(dataPK[0:dimm,0],dataPK[0:dimm,1])#Artemether
pl.plot(dataPK[0:dimm,0],dataPK[0:dimm,3])#Lumefantrine
pl.xlabel('Time (Days)')
pl.ylabel('Parasitaemia [PRBCs / $\mu L$]')
#pl.title('Falciparum infection treated with Artemether & slow-release Lumefantrine')
pl.title('Falciparum infection treated with Artemether & Lumefantrine')
pl.yscale('log')
pl.xlim(0.0, 6 + (2 * T))
pl.ylim(0.00001, b * 1.5) #Scaling of 1.5 just provides some space above 1st peak
#pl.Circle((fev,b), radius = 51, color ='g', fill = True)
pl.savefig("WHM_one_run_withPK.pdf", format='pdf')
pl.show()