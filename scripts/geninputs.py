import os
import numpy as np




x, y, z, beta, theta, phi = map(float,input("Enter x y z beta theta phi: ").split())




v = beta * 2.997925E8
vx = v * np.sin(theta)*np.cos(phi)
vy = v * np.sin(theta)*np.sin(phi)
vz = v * np.cos(theta)

isproton = "1"
dt = "1.0e-11"
tmax = "2.0e-7"
stepout = "100"
radius = "14.1421"


line1 = str(x) + "\t" + str(y) + "\t" +  str(z) + "\n"
line2 = str(vx) + "\t" + str(vy) + "\t" +  str(vz) + "\n"
line3 = isproton + "\t" + dt + tmax + "\n"
line4 = stepout + "\n"
line5 = radius + "\n"

file = open("testinput.dat","w") 
 
for line in line1,line2,line3, line4, line5:
	file.write(line)
file.close() 



# xlist = np.arange(0.0,.1,0.1)
# ylist = np.arange(0.0,.1,0.1)
# zlist = np.arange(0.0,.1,0.1)
# vxlist = np.arange(0.0,.1,0.1)
# vylist = np.arange(0.0,.1,0.1)
# vzlist = np.arange(0.0,1,.1)

# filecounter=0
# isproton=1
# dt = 1.0e-11
# tmax = 2.0e-4
# stepout = 1




# for x in xlist:
# 	for y in ylist:
# 		for z in zlist:
# 			for vx in vxlist:
# 				for vy in vylist:
# 					for vz in vzlist:
# 						print(x,y,z)
# 						print(vx,vy,vz)
# 						print("")
# 						currfile = open("input"+str(filecounter),"w")
# 						currfile.write(str(x) + "\t" + str(y) + "\t" + str(z))
# 						currfile.write("\n")
# 						currfile.write(str(vx) + "\t" + str(vy) + "\t" + str(vz))
# 						currfile.write("\n")
# 						currfile.write(str(isproton) + "\t" + str(dt) + "\t" + str(tmax))
# 						currfile.write("\n")
# 						currfile.write(str(stepout))
# 						filecounter+=1









# while(x <= xmax):
# 	while(y <= ymax):
# 		while(z <= zmax):
# 			while(vx <= vxmax):
# 				while(vy <= vymax):
# 					while(vz <= vzmax):
# 						print(x,y,z)
# 						print(vx,vy,vz)

# 						x += xstep
# 						y += ystep
# 						z += zstep
# 						vx += vxstep
# 						vy += vystep
# 						vz += vzstep




