import pynbody
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import pandas as pd

# Extract and showing Data
file='/mnt/data0/jillian/samuel/test5stampede/blackholes.pkl'
# file='/mnt/data0/jillian/samuel/test7/blackholes.pkl'
with open(file,'rb') as f: #loads pickle files 
    data = pickle.load(f)
# print(data)
mtx=data['data']
# BlackHole 1 data
x1 = mtx[0]['x']
y1 = mtx[0]['y']
z1 = mtx[0]['z']
vx1 = mtx[0]['vx']
vy1 = mtx[0]['vy']
vz1 = mtx[0]['vz']
t1 = mtx[0]['Time']
m1 = mtx[0]['mass']
#BlackHole 2 data
x2 = mtx[1]['x']
y2 = mtx[1]['y']
z2 = mtx[1]['z']
vx2 = mtx[1]['vx']
vy2 = mtx[1]['vy']
vz2 = mtx[1]['vz']
t2 = mtx[1]['Time']
m2 = mtx[1]['mass']

# Finding the Time of Merging of BH1 and BH2
# Special Data
R_t1 = (t1*38844.49)-7160 #Converts time "t1" of BH1 to Megayears (Myr)
R_t2 = (t2*38844.49)-7160 #Converts time "t2" of BH2 to Megayears (Myr)
t1_bcol = R_t1[:len(R_t2)] #Time before colission
t1_acol = R_t1[len(R_t2):] #Time after colission
x1_bcol = x1[:len(R_t2)] 
y1_bcol = y1[:len(R_t2)]
z1_bcol = z1[:len(R_t2)]
x1_acol = x1[len(R_t2):] 
y1_acol = y1[len(R_t2):]
z1_acol = z1[len(R_t2):]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x1_bcol[100:], y1_bcol[100:], z1_bcol[100:], label='BH1 Trajectory',color = 'red', linewidth = 1.0)
ax.scatter([x1_bcol[len(R_t2)-1]], [y1_bcol[len(R_t2)-1]], [z1_bcol[len(R_t2)-1]], color='black', s=50)
# ax.scatter([x1_bcol[0]], [y1_bcol[0]], [z1_bcol[0]], color='green', s=25)
ax.plot(x2[100:], y2[100:], z2[100:], label='BH2 Trajectory',color = 'blue', linewidth = 1.0)
ax.scatter([x2[len(R_t2)-1]], [y2[len(R_t2)-1]], [z2[len(R_t2)-1]], color='black', s=50)
# ax.scatter([x2[0]], [y2[0]], [z2[0]], color='green', s=25)
ax.plot(x1_acol, y1_acol, z1_acol, label='SBH Trajectory',color = 'purple', linewidth = 1.75,)
ax.scatter([x1_acol[0]], [y1_acol[0]], [z1_acol[0]], color='black', s=250)
ax.scatter([x1_acol[len(t1_acol)-1]], [y1_acol[len(t1_acol)-1]], [z1_acol[len(t1_acol)-1]], color='black', marker = 'x', s=50)
ax.set_xlabel('X [kpc]')
ax.set_ylabel('Y [kpc]')
ax.set_zlabel('Z [kpc]')
ax.grid(True)
ax.legend()

#Animation for Black Holes merging

all_x = np.concatenate((x1_bcol, x2, x1_acol))
all_y = np.concatenate((y1_bcol, y2, y1_acol))
all_z = np.concatenate((z1_bcol, z2, z1_acol))

x_limits = (min(all_x), max(all_x))
y_limits = (min(all_y), max(all_y))
z_limits = (min(all_z), max(all_z))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Initialize the line objects (Describes and graph the trajectories)
line1, = ax.plot([], [], [], 'r-', label='BH1 Trajectory')
line2, = ax.plot([], [], [], 'b-', label='BH2 Trajectory')
line_merge, = ax.plot([], [], [], 'g-', label='Merged BH Trajectory')

# Initialize the scatter point (Black Holes points as particles)
point1, = ax.plot([], [], [], 'go', label='Black Hole 1')
point2, = ax.plot([], [], [], 'mo', label='Black Hole 2')
point_merge, = ax.plot([], [], [], 'co', label='Merged Black Hole')

ax.set_xlim(x_limits)
ax.set_ylim(y_limits)
ax.set_zlim(z_limits)

ax.legend()

# Initialization function: plot the background of each frame
def init():
    line1.set_data([], [])
    line1.set_3d_properties([])
    line2.set_data([], [])
    line2.set_3d_properties([])
    line_merge.set_data([], [])
    line_merge.set_3d_properties([])
    point1.set_data([], [])
    point1.set_3d_properties([])
    point2.set_data([], [])
    point2.set_3d_properties([])
    point_merge.set_data([], [])
    point_merge.set_3d_properties([])
    return line1, line2, line_merge, point1, point2, point_merge

# Animation function: update the plot for each frame
def update(frame):
    if frame < len(x1_bcol):
        line1.set_data(x1_bcol[:frame], y1_bcol[:frame])
        line1.set_3d_properties(z1_bcol[:frame])
        point1.set_data(x1_bcol[frame], y1_bcol[frame])
        point1.set_3d_properties(z1_bcol[frame])
    else:
        line1.set_data(x1_bcol, y1_bcol)
        line1.set_3d_properties(z1_bcol)
        point1.set_data(x1_bcol[-1], y1_bcol[-1])
        point1.set_3d_properties(z1_bcol[-1])
    
    if frame < len(x2):
        line2.set_data(x2[:frame], y2[:frame])
        line2.set_3d_properties(z2[:frame])
        point2.set_data(x2[frame], y2[frame])
        point2.set_3d_properties(z2[frame])
    else:
        line2.set_data(x2, y2)
        line2.set_3d_properties(z2)
        point2.set_data(x2[-1], y2[-1])
        point2.set_3d_properties(z2[-1])
    if frame >= len(x1_bcol) + len(x2):
        # Post-merger trajectory
        merge_index = frame - (len(x1_bcol) + len(x2))
        if merge_index < len(x1_acol):
            line_merge.set_data(x1_acol[:merge_index], y1_acol[:merge_index])
            line_merge.set_3d_properties(z_merge[:merge_index])
            point_merge.set_data(x1_acol[merge_index], y1_acol[merge_index])
            point_merge.set_3d_properties(z1_acol[merge_index])
        else:
            line_merge.set_data(x1_acol, y1_acol)
            line_merge.set_3d_properties(z1_acol)
            point_merge.set_data(x1_acol[-1], y1_acol[-1])
            point_merge.set_3d_properties(z1_acol[-1])

    return line1, line2, line_merge, point1, point2, point_merge


ani = FuncAnimation(fig, update, frames=len(x1_bcol) + len(x2) + len(x1_acol), init_func=init, blit=True, interval=2.5)

plt.show()