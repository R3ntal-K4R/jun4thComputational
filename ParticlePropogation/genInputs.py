import os, random
import numpy as np

def createFile():
    random.seed()

    while True:
        x = random.randrange(1,200) * 0.1
        y = random.randrange(1,200) * 0.1
        if (np.sqrt(x**2 + y**2) >= 15 or np.sqrt(x**2 + y**2) <= 5):
            break
    z = random.randrange(1,200) * 0.1

    beta = 2.997925e8 * random.random()
    theta = random.random() * 2
    phi = random.random() * 2

    vx = beta * np.sin(theta)*np.cos(phi)
    vy = beta * np.sin(theta)*np.sin(phi)
    vz = beta * np.cos(theta)

    charge = 1
    mass = 1
    dt = 1.0e-13
    tmax = 1.0e-7

    line1 = str(x) + "\t" + str(y) + "\t" +  str(z) + "\n"
    line2 = str(vx) + "\t" + str(vy) + "\t" +  str(vz) + "\n"
    line3 = str(charge) + "\t" + str(mass) + "\n"
    line4 = str(dt) + "\t" + str(tmax) + "\n"
    line5 = str(500) + "\n"
    line6 = str(5) + "\n"

    lines = [line1, line2, line3, line4, line5, line6]
    with open("input.dat", "w") as f:
        for line in lines:
            f.write(line)
    return

def main():
    for run in range(21,51):
        try:
            os.mkdir("runs/"+str(run)+"/")
        except:
            pass
        os.chdir("runs/"+str(run)+"/")
        createFile()
        os.system("cp ../../particlePropagation .")
        os.system("./particlePropagation ../../Density__1.2510E_min_03_g_cm3__Target_N_7_100.00_100.00__Particle_Hydrogen.txt ../../Bfield.dat input.dat")
        os.chdir("../../")

if(__name__ == "__main__"):
    main()
