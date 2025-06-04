import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

#note requires ffmpeg and imagemagick to save gif
input_x=10.0
input_y=10.0
input_z=10.0
height_cylin=5.0
height_sphere= height_cylin
radius_cylin=5.0




def update_lines(num, dataLines, lines):
    for line, data in zip(lines, dataLines):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
        
    return lines

coords= []
array_x=[]
array_y=[]
array_z=[]
i=1
unique= open("output1.dat","r")
next(unique)
for line in unique:
    i=i+1
    if (i%1000==0):
        t,x,y,z,vx,vy,vz,ax,ay,az,rmag,Bx,By,Bz,eloss,blah =line.split()
        coords.append([float(x),float(y),float(z)])
        array_x.append(float(x))
        array_y.append(float(y))
        array_z.append(float(z))
        
    
unique.close()

data=np.array([[array_x,array_y,array_z]])
print(data)

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)





# NOTE: Can't pass empty arrays into 3d version of plot()
lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

# Setting the axes properties
ax.set_xlim3d([0.0, 20.0])
ax.set_xlabel('X')

ax.set_ylim3d([0.0, 20.0])
ax.set_ylabel('Y')

ax.set_zlim3d([0.0, 20.0])
ax.set_zlabel('Z')

ax.set_title('3D Test')

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, update_lines, 50, fargs=(data, lines),
                                   interval=100, repeat=True)




#working cylinder
# plot_3D_cylinder(radius, height, elevation=elevation, resolution=resolution, color=color, x_center=x_center, y_center=y_center)
def data_for_cylinder_along_z(center_x,center_y,radius,height_z):
    z = np.linspace(0, height_z, 50)+input_z
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

from mpl_toolkits.mplot3d import Axes3D



Xc,Yc,Zc = data_for_cylinder_along_z(input_x,input_y,radius_cylin,height_cylin)
ax.plot_surface(Xc, Yc, Zc, alpha=0.5)

from itertools import product, combinations



scale_sp=radius_cylin


# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = scale_sp*np.cos(u)*np.sin(v)+input_x
y = scale_sp*np.sin(u)*np.sin(v)+input_y
z = scale_sp*np.cos(v)+input_z
ax.plot_wireframe(x, y, z, color="r")





# from matplotlib import rcParams
# rcParams['animation.convert_path'] = r'C:\Program Files\ImageMagick\convert'
# rcParams['animation.ffmpeg_path'] = r'C:\Users\wscot\Downloads\ffmpeg-20200617-0b3bd00-win64-static\ffmpeg-20200617-0b3bd00-win64-static\bin'


plt.show()
# line_ani.save('myAnimation.gif', writer="imagemagick", fps=30)

