import pynbody
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from scipy.constants import gravitational_constant as G

# Extract and showing Data
file='/mnt/data0/jillian/samuel/test7/blackholes.pkl'
with open(file,'rb') as f: #loads pickle files 
    data = pickle.load(f)
# print(data)
mtx = data['data']
# BlackHole 1 data
x1 = mtx[0]['x']
y1 = mtx[0]['y']
z1 = mtx[0]['z']
vx1 = mtx[0]['vx']
vy1 = mtx[0]['vy']
vz1 = mtx[0]['vz']
t1 = mtx[0]['Time']
m1 = mtx[0]['mass']
# BlackHole 2 data
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
t1_bcol = R_t1[:len(R_t2)] #Time before collision
t1_acol = R_t1[len(R_t2):] #Time after collision
'''
# Showing the new time conversion also the time before and after colission
print("Time for BH1", t1, len(t1))
print("Time for BH2", t2, len(t2))
print("t1 in Myr is:", R_t1, len(R_t1))
print("t2 in Myr is:", R_t2, len(R_t2))
print("Time before colission:",t1_bcol ,len(t1_bcol))
print("Time after colission:", t1_acol ,len(t1_acol))
'''
print("The merging of BH1 and BH2 is aproximately at:", "{:.2f}".format(R_t2[len(R_t2)-1]),"Megayears")

# Finding the last distance for merging
# Divide the BH1 position before and after merging
x1_bcol = x1[:len(R_t2)] 
y1_bcol = y1[:len(R_t2)]
z1_bcol = z1[:len(R_t2)]
x1_acol = x1[len(R_t2):] 
y1_acol = y1[len(R_t2):]
z1_acol = z1[len(R_t2):]
dist_BH = np.sqrt((x2-x1_bcol)**2+(y2-y1_bcol)**2+(z2-z1_bcol)**2)
print("The merging of BH1 and BH2 is aproximately at:", "{:.2f}".format(dist_BH[len(R_t2)-1]),"kpc")

# Making a 3D plot of BHs position before and after
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x1_bcol, y1_bcol, z1_bcol, label='BH1 Trajectory',color = 'red', linewidth = 1.5)
ax.scatter([x1_bcol[len(R_t2)-1]], [y1_bcol[len(R_t2)-1]], [z1_bcol[len(R_t2)-1]], color='black', s=175)
#ax.scatter([x1_bcol[100]], [y1_bcol[100]], [z1_bcol[100]], color='green', s=25)
ax.plot(x2, y2, z2, label='BH2 Trajectory',color = 'blue', linewidth = 1.5)
ax.scatter([x2[len(R_t2)-1]], [y2[len(R_t2)-1]], [z2[len(R_t2)-1]], color='black', s=175)
#ax.scatter([x2[100]], [y2[100]], [z2[100]], color='green', s=25)
ax.plot(x1_acol, y1_acol, z1_acol, label='SBH Trajectory',color = 'green', linewidth = 1.75)
ax.scatter([x1_acol[len(t1_acol)-1]], [y1_acol[len(t1_acol)-1]], [z1_acol[len(t1_acol)-1]], color='black', marker = 'o', s=350)
# ax.scatter([x1_acol[0]], [y1_acol[0]], [z1_acol[0]], color='black', s=250)
ax.set_xlabel('X [kpc]')
ax.set_ylabel('Y [kpc]')
ax.set_zlabel('Z [kpc]')
ax.grid(True)
ax.legend()

# Making a 3D plot of BHs position before
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x1_bcol, y1_bcol, z1_bcol, label='BH1 Trajectory',color = 'red', linewidth = 2.0, marker = 'o')
ax.scatter([x1_bcol[len(R_t2)-1]], [y1_bcol[len(R_t2)-1]], [z1_bcol[len(R_t2)-1]], color='black', s=150)
# ax.scatter([x1_bcol[0]], [y1_bcol[0]], [z1_bcol[0]], color='green', s=25)
ax.plot(x2, y2, z2, label='BH2 Trajectory',color = 'blue', linewidth = 2.0, marker = 'o')
ax.scatter([x2[len(R_t2)-1]], [y2[len(R_t2)-1]], [z2[len(R_t2)-1]], color='black', s=150)
# ax.scatter([x2[0]], [y2[0]], [z2[0]], color='green', s=25)
ax.set_xlabel('X [kpc]')
ax.set_ylabel('Y [kpc]')
ax.set_zlabel('Z [kpc]')
ax.grid(True)
ax.legend()

'''
# Distance of BH1 and BH2 before merging
fig = plt.figure(figsize = (10,8))
ax1 = fig.add_subplot(2, 2, 1)  
ax1.plot(R_t2, x1_bcol, label = 'BH1 x-position', color = 'red')
ax1.plot(R_t2, x2, label = 'BH2 x-position', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax1.set_xlim(left = 0)
plt.xlabel("Time [Megayears]")
plt.ylabel("BH distance [kpc]")
ax1.set_title('BH1 & BH2 x-position')
ax1.legend()

ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(R_t2, y1_bcol, label = 'BH1 y-position', color = 'red')
ax2.plot(R_t2, y2, label = 'BH2 y-position', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax2.set_xlim(left = 0)
plt.xlabel("Time [Megayears]")
plt.ylabel("BH distance [kpc]")
ax2.set_title('BH1 & BH2 y-position')
ax2.legend()

ax3 = fig.add_subplot(2, 2, 3) 
ax3.plot(R_t2, z1_bcol, label = 'BH1 z-position', color = 'red')
ax3.plot(R_t2, z2, label = 'BH2 z-position', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax3.set_xlim(left = 0)
plt.xlabel("Time [Megayears]")
plt.ylabel("BH distance [kpc]")
ax3.set_title('BH1 & BH2 z-position')
ax3.legend()

ax4 = fig.add_subplot(2, 2, 4) 
ax4.plot(R_t2 , dist_BH , label = 'BH1 & BH2 Distance', linestyle = '-', color = 'purple')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax4.set_xlim(left = 0)
plt.xlabel("Time [Megayears]")
plt.ylabel("BH distance [kpc]")
ax4.set_title('Distance between BH1 & BH2 around the time')
ax4.legend()
plt.tight_layout()

# Making a plot for the change of mass of BBHS
fig = plt.figure()
plt.plot(R_t1 , m1 , label = 'BH1 mass', linestyle = '--', linewidth = 2.0 , color = 'orange')
plt.plot(R_t2 , m2 , label = 'BH2 mass', linestyle = '-', color = 'purple')
plt.xlabel("Time [Megayears]")
plt.ylabel("BH mass [Solar Mass]")
plt.title("BH1 and BH2 mass around the time")
plt.legend()
'''
# Plot for BH1 & BH2 xyz-velocities before merging
fig = plt.figure(figsize = (10,8))
ax1 = fig.add_subplot(2, 2, 1) 
ax1.plot(R_t2, vx1[:len(R_t2)], label = 'BH1 x-velocity', color = 'red')
ax1.plot(R_t2, vx2, label = 'BH2 x-velocity', color = 'blue')
ax1.set_xlim([0 , R_t2[len(R_t2)-1]])
plt.xlabel("Time [Megayears]")
plt.ylabel("Velocity [km/s]")
ax1.set_title('BH1 & BH2 x-velocity')
ax1.legend()

ax2 = fig.add_subplot(2, 2, 2) 
ax2.plot(R_t2, vy1[:len(R_t2)], label = 'BH1 y-velocity', color = 'red')
ax2.plot(R_t2, vy2, label = 'BH2 y-velocity', color = 'blue')
ax2.set_xlim([0 , R_t2[len(R_t2)-1]])
plt.xlabel("Time [Megayears]")
plt.ylabel("Velocity [km/s]")
ax2.set_title('BH1 & BH2 y-velocity')
ax2.legend()

ax3 = fig.add_subplot(2, 2, 3) 
ax3.plot(R_t2, vz1[:len(R_t2)], label = 'BH1 z-velocity', color = 'red')
ax3.plot(R_t2, vz2, label = 'BH2 z-velocity', color = 'blue')
ax3.set_xlim([0 , R_t2[len(R_t2)-1]])
plt.xlabel("Time [Megayears]")
plt.ylabel("Velocity [km/s]")
ax3.set_title('BH1 & BH2 z-velocity')
ax3.legend()

# Plot of BH1 & BH2 velocity before merging
v1_mag = np.sqrt((vx1**2)+(vy1**2)+(vz1**2))
v2_mag = np.sqrt((vx2**2)+(vy2**2)+(vz2**2))

ax4 = fig.add_subplot(2, 2, 4) 
ax4.plot(R_t2, v1_mag[:len(R_t2)], label = 'BH1 velocity', color = 'red')
ax4.plot(R_t2, v2_mag, label = 'BH2 velocity', color = 'blue')
ax4.set_xlim([0 , R_t2[len(R_t2)-1]])
plt.xlabel("Time [Megayears]")
plt.ylabel("Velocity [km/s]")
ax4.set_title('BH1 & BH2 velocity')
ax4.legend()
plt.tight_layout()

fig = plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax1.plot(R_t1, vx1, label = 'SBH X-velocity', color = 'orange', linestyle = '--')
ax1.plot(R_t1, vy1, label = 'SBH Y-velocity', color = 'blue', linestyle = '--')
ax1.plot(R_t1, vz1, label = 'SBH Z-velocity', color = 'green', linestyle = '--')
ax1.plot(t1_bcol, vx1[:len(R_t2)], label = 'BH1 X-velocity', color = 'orange')
ax1.plot(t1_bcol, vy1[:len(R_t2)], label = 'BH1 Y-velocity', color = 'blue')
ax1.plot(t1_bcol, vz1[:len(R_t2)], label = 'BH1 Z-velocity', color = 'green')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax1.set_xlim([0,80])
plt.xlabel("Time [Megayears]")
plt.ylabel("Velocity [km/s]")
ax1.set_title('XYZ velocity of BH1 and SBH')
ax1.legend()

ax2 = fig.add_subplot(2,1,2)
#plt.plot(t1_bcol, vx1[:len(R_t2)], label = 'BH1 X-velocity', color = 'red')
#plt.plot(t1_bcol, vy1[:len(R_t2)], label = 'BH1 Y-velocity', color = 'blue')
#plt.plot(t1_bcol, vz1[:len(R_t2)], label = 'BH1 Z-velocity', color = 'green')
ax2.plot(R_t2, vx2, label = 'BH2 X-velocity', color = 'orange')
ax2.plot(R_t2, vy2, label = 'BH2 Y-velocity', color = 'blue')
ax2.plot(R_t2, vz2, label = 'BH2 Z-velocity', color = 'green')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax2.set_xlim([0,80])
plt.xlabel("Time [Megayears]")
plt.ylabel("Velocity [km/s]")
ax2.set_title('XYZ velocity of BH2')
ax2.legend()
plt.tight_layout()

'''
# Determining the xyz-momentum of BH1 & BH2
lmx1 = m1*vx1
lmy1 = m1*vy1
lmz1 = m1*vz1
lmx2 = m2*vx2
lmy2 = m2*vy2
lmz2 = m2*vz2

# BH1 & BH2 xyz-momentum
fig = plt.figure(figsize = (10,8))
plt.plot(R_t1, lmx1, label = 'x-momentum BH1' , color = 'orange')
plt.plot(R_t1, lmy1, label = 'y-momentum BH1' , color = 'purple')
plt.plot(R_t1, lmz1, label = 'z-momentum BH1' , color = 'cyan')
plt.plot(R_t2, lmx2, label = 'x-momentum BH2' , linestyle = '--', color = 'red' )
plt.plot(R_t2,lmy2, label = 'y-momentum BH2' , linestyle = '--', color = 'blue' )
plt.plot(R_t2,lmz2, label = 'z-momentum BH2' , linestyle = '--', color = 'green' )
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
plt.xlabel("Time [Megayears]")
plt.ylabel("BH momentum [Msun*km/s]")
plt.title("BH1 & BH2 xyz-momentum around the time")
plt.legend()

# Total xyz-momentum before and after merging
tot_lmx_bcol = m1[:len(R_t2)] *vx1[:len(R_t2)] + m2*vx2
tot_lmy_bcol = m1[:len(R_t2)] *vy1[:len(R_t2)] + m2*vy2
tot_lmz_bcol = m1[:len(R_t2)] *vz1[:len(R_t2)] + m2*vz2

tot_lmx_acol = m1[len(R_t2):]* vx1[len(R_t2):]
tot_lmy_acol = m1[len(R_t2):]* vy1[len(R_t2):]
tot_lmz_acol = m1[len(R_t2):]* vz1[len(R_t2):]

# Plot for BH1 & BH2 xyz-TOTAL momentum before and after merging
fig = plt.figure( figsize = (10,8))
ax1 = fig.add_subplot(2, 2, 1)
ax1.plot(t1_bcol , tot_lmx_bcol , label = 'Total x-momentum before merging', color = 'red')
ax1.plot(t1_acol , tot_lmx_acol , label = 'Total x-momentum after merging', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
plt.xlabel("Time [Megayears]")
plt.ylabel("BH momentum [Msun*km/s]")
ax1.set_title("Total x-momentum around time")
ax1.legend()

ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(t1_bcol , tot_lmy_bcol , label = 'Total y-momentum before merging', color = 'red')
ax2.plot(t1_acol , tot_lmy_acol , label = 'Total y-momentum after merging', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
plt.xlabel("Time [Megayears]")
plt.ylabel("BH momentum [Msun*km/s]")
ax2.set_title("Total y-momentum around time")
ax2.legend()

ax3 = fig.add_subplot(2, 2, 3)
ax3.plot(t1_bcol , tot_lmz_bcol , label = 'Total z-momentum before merging', color = 'red')
ax3.plot(t1_acol , tot_lmz_acol , label = 'Total z-momentum after merging', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
plt.xlabel("Time [Megayears]")
plt.ylabel("BH momentum [Msun*km/s]")
ax3.set_title("Total z-momentum around time")
ax3.legend()

# Determine the Linear Momentum
lm1 = m1*v1_mag
lm2 = m2*v2_mag
tot_lm_bcol = lm2 + lm1[:len(R_t2)]
tot_lm_acol = lm1[len(R_t2):]
ax4 = fig.add_subplot(2, 2, 4)
ax4.plot(t1_bcol , tot_lm_bcol , label = 'Total momentum before merging', color = 'red')
ax4.plot(t1_acol , tot_lm_acol , label = 'Total momentum after merging', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
plt.xlabel("Time [Megayears]")
plt.ylabel("BH momentum [Msun*km/s]")
ax4.set_title("Total z-momentum around time")
ax4.legend()
plt.tight_layout()

print(len(R_t1))
print(len(R_t2))

print("BH1 initial xyz position")
print("{:.2e}".format(x1_bcol[0]))
print("{:.2e}".format(y1_bcol[0]))
print("{:.2e}".format(z1_bcol[0]))

print("BH1 final xyz position")
print("{:.2e}".format(x1_bcol[len(R_t2)-1]))
print("{:.2e}".format(y1_bcol[len(R_t2)-1]))
print("{:.2e}".format(z1_bcol[len(R_t2)-1]))

print("BH2 initial xyz position")
print("{:.2e}".format(x2[0]))
print("{:.2e}".format(y2[0]))
print("{:.2e}".format(z2[0]))

print("BH2 final xyz position")
print("{:.2e}".format(x2[len(R_t2)-1]))
print("{:.2e}".format(y2[len(R_t2)-1]))
print("{:.2e}".format(z2[len(R_t2)-1]))

print("SBH initial xyz position")
print("{:.2e}".format(x1_acol[0]))
print("{:.2e}".format(y1_acol[0]))
print("{:.2e}".format(z1_acol[0]))

print("SBH final xyz position")
print("{:.2e}".format(x1[len(R_t1)-1]))
print("{:.2e}".format(y1[len(R_t1)-1]))
print("{:.2e}".format(z1[len(R_t1)-1]))

print("BH1 initial xyz velocity")
print("{:.2e}".format(vx1[0]))
print("{:.2e}".format(vy1[0]))
print("{:.2e}".format(vz1[0]))

print("BH2 initial xyz velocity")
print("{:.2e}".format(vx2[0]))
print("{:.2e}".format(vy2[0]))
print("{:.2e}".format(vz2[0]))

print("SBH initial xyz velocity")
print("{:.2e}".format(vx1[(len(R_t2)-1):][0]))
print("{:.2e}".format(vy1[(len(R_t2)-1):][0]))
print("{:.2e}".format(vz1[(len(R_t2)-1):][0]))

print("BH1 & BH2 xyz Momentum")
print("xBefore:","{:.2e}".format(lmx1[len(R_t2)-1]))
print("xBefore:","{:.2e}".format(lmx2[len(R_t2)-1]))
print("xAfter:","{:.2e}".format(lmx1[len(R_t2)]))

print("yBefore:","{:.2e}".format(lmy1[len(R_t2)-1]))
print("yBefore:","{:.2e}".format(lmy2[len(R_t2)-1]))
print("yAfter:","{:.2e}".format(lmy1[len(R_t2)]))

print("zBefore:","{:.2e}".format(lmz1[len(R_t2)-1]))
print("zBefore:","{:.2e}".format(lmz2[len(R_t2)-1]))
print("zAfter:","{:.2e}".format(lmz1[len(R_t2)]))

print("The average velocity of BH1:","{:.2f}".format(np.mean(v1_mag[:len(R_t2)])),"km/s")
print("The average velocity of BH2:","{:.2f}".format(np.mean(v2_mag)),"km/s")
print("The average velocity of SBH:","{:.2f}".format(np.mean(v1_mag[len(R_t2):])),"km/s")
'''
'''
# Determine the Kinetic Energy and Gravitational Potential Energy
KE1 = 0.5*m1*(v1_mag**2)
KE2 = 0.5*m2*(v2_mag**2)
U = G * (m2 * m1[:len(R_t2)])/dist_BH

fig = plt.figure()
plt.plot(R_t1, KE1, label = 'BH1 kinetic energy' , color = 'red')
plt.plot(R_t2, KE2, label = 'BH2 kinetic energy' , color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
plt.xlabel("Time [Megayears]")
plt.ylabel("BH kinetic energy [Msun*(km/s)^2]")
plt.title("BH1 & BH2 kinetic energy around the time")
plt.legend()

fig = plt.figure()
plt.plot(R_t2, U, label = 'BH2 Potential energy' , color = 'purple')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
plt.xlabel("Time [Megayears]")
plt.ylabel("BH kinetic energy [Msun*(km/s)^2]")
plt.title("BH1 & BH2 kinetic energy around the time")
plt.legend()
'''

plt.show()
