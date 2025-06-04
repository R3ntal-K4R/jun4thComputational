import sys
units = ""

if len(sys.argv)-1 == 0:
    print("Give SRIM file as argument during terminal call")
    exit()
elif len(sys.argv)-1 == 1:
    SRIM = open(sys.argv[1], 'r')
    outFile = open("elossfile.dat", "w")
elif len(sys.argv)-1 == 2:
    SRIM = open(sys.argv[1], 'r')
    outFile = open(sys.argv[2], 'w')
else:
    SRIM = open(sys.argv[1], 'r')
    outFile = open(sys.argv[2], 'w')
    units = sys.argv[3]

if units == "SI" or units == "si":
    factorE = 1.60218e-19
    factorD = 1e-10
    factorS = factorE/factorD
else:
    factorE = 1
    factorD = 1
    factorS = 1

energy = 0
elec = 0
nuclear = 0
ran = 0
longitudinal = 0
lateral = 0

for i in range(25):
    SRIM.readline()

if units == "SI" or units == "si":
    outFile.write("energy(J)  elec(J/m)    nuclear(J/m)  range(m)   longitudinal(m)    lateral(m)\n")
else:
    outFile.write("energy(eV)  elec(eV/A)    nuclear(eV/A)  range(A)   longitudinal(A)    lateral(A)\n")

line = SRIM.readline()
while(line[0]!="-"):
    line = line.split(" ")
    while '' in line:
        line.remove('')

    if line[1] == 'keV':
        energy = ((float(line[0])*1e3))
    elif line[1] == 'MeV':
        energy = ((float(line[0])*1e6))
    else:
        energy = (float(line[0]))

    elec = float(line[2])
    nuclear = float(line[3])

    if line[5] == 'um':
        ran = ((float(line[4])*1e4))
    elif line[5] == 'mm':
        ran = ((float(line[4])*1e7))
    elif line[5] == 'm':
        ran = ((float(line[4])*1e10))
    else:
        ran = (float(line[4]))

    if line[7] == 'um':
        longitudinal = ((float(line[6])*1e4))
    elif line[7] == 'mm':
        longitudinal = ((float(line[6])*1e7))
    elif line[7] == 'm':
        longitudinal = ((float(line[6])*1e10))
    else:
        longitudinal = float(line[6])

    if line[9] == 'um':
        lateral = ((float(line[8])*1e4))
    elif line[9] == 'mm':
        lateral = ((float(line[8])*1e7))
    elif line[9] == 'm':
        lateral = ((float(line[8])*1e10))
    else:
        lateral = float(line[8])   
    line = SRIM.readline()
    outFile.write(str(factorE*energy) + "  " + str(factorS*elec) + "  " + str(factorS*nuclear) + "  " + str(factorD*ran) + "  " + str(factorD*longitudinal) + "  " + str(factorD*lateral) + "\n")

SRIM.close()
outFile.close()