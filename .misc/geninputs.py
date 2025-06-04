import os
import numpy as np



xlist = np.arange(0.0,.1,0.1)
ylist = np.arange(0.0,.1,0.1)
zlist = np.arange(0.0,.1,0.1)
vxlist = np.arange(0.0,.1,0.1)
vylist = np.arange(0.0,.1,0.1)
vzlist = np.arange(0.0,1,.1)

filecounter=0
isproton=1
dt = 1.0e-11
tmax = 2.0e-4
stepout = 1




for x in xlist:
	for y in ylist:
		for z in zlist:
			for vx in vxlist:
				for vy in vylist:
					for vz in vzlist:
						print(x,y,z)
						print(vx,vy,vz)
						print("")
						currfile = open("input"+str(filecounter),"w")
						currfile.write(str(x) + "\t" + str(y) + "\t" + str(z))
						currfile.write("\n")
						currfile.write(str(vx) + "\t" + str(vy) + "\t" + str(vz))
						currfile.write("\n")
						currfile.write(str(isproton) + "\t" + str(dt) + "\t" + str(tmax))
						currfile.write("\n")
						currfile.write(str(stepout))
						filecounter+=1









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




