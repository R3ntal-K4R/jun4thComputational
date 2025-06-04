import sys

if len(sys.argv)-1 == 0:
    print("Give SRIM file as argument during terminal call")
    exit()
elif len(sys.argv)-1 == 1:
    SRIM = open(sys.argv[1], 'r')
    outFile = open("elossfile.dat", "w")
else:
    SRIM = open(sys.argv[1], 'r')
    outFile = open(sys.argv[2], 'w')

energy = 0
elec = 0
nuclear = 0
ran = 0
longitudinal = 0
lateral = 0

for i in range(25):
    SRIM.readline()

# outFile.write("energy[eV]  elec[eV/m]    nuclear[eV/m]  range[A]   longitudinal[A]    lateral[A]\n")
outFile.write("energy[eV]  elec[eV/m]    nuclear[eV/m]  range[m]   longitudinal[m]    lateral[m]\n")
line = SRIM.readline()
while(line[0]!="-"):
    line = line.split(" ")
    while '' in line:
        line.remove('')

    if line[1] == 'keV':
        energy = ((float(line[0])*1000))
    elif line[1] == 'MeV':
        energy = ((float(line[0])*1000000))
    elif line[1] == "GeV":
        energy = ((float(line[0])*1000000000))
    else:
        energy = (line[0])
    elec = (float(line[2])*10e9)    #for keV/micron stopping units
    nuclear = (float(line[3])*10e9)

    #Commented out section puts everything in Angstroms, boo.
    """
    if line[5] == 'um':
        ran = ((float(line[4])*10000))
    elif line[5] == 'mm':
        ran = ((float(line[4])*10000000))
    elif line[5] == 'm':
        ran = ((float(line[4])*10000000000))
    elif line[5] == 'km':
        ran = ((float(line[4])*10000000000000))
    else:
        ran = (line[4])

    if line[7] == 'um':
        longitudinal = ((float(line[6])*10000))
        # print(longitudinal)
    elif line[7] == 'mm':
        longitudinal = ((float(line[6])*10000000))
    elif line[7] == 'm':
        longitudinal = ((float(line[6])*10000000000))
    elif line[7] == 'km':
        longitudinal = ((float(line[6])*10000000000000))
    else:
        longitudinal = (line[6])

    if line[9] == 'um':
        lateral = ((float(line[8])*10000))
    elif line[9] == 'mm':
        lateral = ((float(line[8])*10000000))
    elif line[9] == 'm':
        lateral = ((float(line[8])*10000000000))
    else:
        lateral = (line[8])   
    """
    #SI Unit conversion (Chad)
    if line[5] == 'A':
        ran = ((float(line[4])*10e-10))
    elif line[5] == 'um':
        ran = ((float(line[4])*10e-6))
    elif line[5] == 'mm':
        ran = ((float(line[4])*10e-3))
    elif line[5] == 'm':
        ran = ((float(line[4])))
    elif line[5] == 'km':
        ran = ((float(line[4])*1000))
    else:
        ran = (line[4])

    if line[7] == 'A':
        longitudinal = ((float(line[6])*10e-10))
    elif line[7] == 'um':
        longitudinal = ((float(line[6])*10e-6))
        # print(longitudinal)
    elif line[7] == 'mm':
        longitudinal = ((float(line[6])*10e-3))
    elif line[7] == 'm':
        longitudinal = ((float(line[6])))
    elif line[7] == 'km':
        longitudinal = ((float(line[6])*1000))
    else:
        longitudinal = (line[6])

    if line[9] == 'A':
        lateral = ((float(line[8])*10e-10))
    elif line[9] == 'um':
        lateral = ((float(line[8])*10e-6))
    elif line[9] == 'mm':
        lateral = ((float(line[8])*10e-3))
    elif line[9] == 'm':
        lateral = ((float(line[8])))
    else:
        lateral = (line[8]) 
    # print(line) 
    line = SRIM.readline()
    # print(longitudinal)
    # print(energy + "  " + elec + "  " + nuclear + "  " + ran + "  " + longitudinal + "  " + lateral)
    outFile.write(str(energy) + "  " + str(elec) + "  " + str(nuclear) + "  " + str(ran) + "  " + str(longitudinal) + "  " + str(lateral) + "\n")

SRIM.close()
outFile.close()