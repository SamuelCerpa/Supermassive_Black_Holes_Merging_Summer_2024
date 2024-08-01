import pynbody
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Oldest, does not work well
# file = "/mnt/data0/jillian/samuel/testj/testj.BlackHoles"
# s = pynbody.load("/mnt/data0/jillian/samuel/testj/testj.000012")

# Strange results, strange time
# file = "/mnt/data0/jillian/samuel/test5stampede/dennis.BlackHoles"
# s = pynbody.load("/mnt/data0/jillian/samuel/test5stampede/dennis.001104")

# Used on presentation
file = "/mnt/data0/jillian/samuel/newmerge/testj.BlackHoles"
s = pynbody.load("/mnt/data0/jillian/samuel/newmerge/testj.000012")

munits = float(s["mass"].units)
posunits = float(s["x"].units)
velunits = float(s["vx"].units)
potunits = float(s["phi"].units)
tunits = 1./7.587954066206311e-19  # seconds?
#print("t unit = ",tunits)  #  super unhelpful
ttwo = 7.587954066206311e-19  # actual tunit into seconds
Eunits = munits*potunits
# this is in terrible units
efactor = 4.36e+65 / 2.18e22   #  used units command
  #  "2.18e+22 km**2 Msol a**-1 s**-2"  is the written unit
def read_orbit_file(file):
    # read in the file
    print("reading file")
    df = pd.read_table(file, sep="\s+")
    df.columns = ["iord", "time", "step", "mass[su]", "x[su]", "y[su]", "z[su]", "vx[su]", "vy[su]", "vz[su]", "pot", "Mdot", "dM", "E", "dt_eff", "a","zero","thing"]
    print("file read")
    # converting units
    df["time"] = df["time"] * tunits / 3.1556926e+16#  gigayear
    df["mass[su]"] *= munits
    df["x[su]"] *= posunits
    df["y[su]"] *= posunits
    df["z[su]"] *= posunits
    df["vx[su]"] *= velunits
    df["vy[su]"] *= velunits
    df["vz[su]"] *= velunits
    df["pot"] *= potunits
    df["Mdot"] *= munits/tunits  # msun / sec
    df["dM"] *= munits
    df["E"] *= Eunits*efactor # now in erg
    df["dt_eff"] = df["dt_eff"] * tunits #  seconds?    tunits * 3.1556926e+16 # gigayear to second
    df = df.rename(columns={"mass[su]": "mass[Msun]", "x[su]": "x[kpc]", "y[su]": "y[kpc]","z[su]": "z[kpc]", "vx[su]": "vx[km/s]","vy[su]": "vy[km/s]","vz[su]": "vz[km/s]", "E": "E [erg]", "time": "time[Gyr]"})
#    print(df.columns)
    return df

df = read_orbit_file(file)
print(df.columns)

BH1 = np.where(df["iord"] == 8388608)[0]
BH2 = np.where(df["iord"] == 8388609)[0]
# 194 elements
x1 = np.array(df["x[kpc]"][BH1])
y1 = np.array(df["y[kpc]"][BH1])
z1 = np.array(df["z[kpc]"][BH1])
vx1 = np.array(df["vx[km/s]"][BH1])
vy1 = np.array(df["vy[km/s]"][BH1]) 
vz1 = np.array(df["vz[km/s]"][BH1]) 
t1 = np.array(df["time[Gyr]"][BH1])
m1 = np.array(df["mass[Msun]"][BH1])
# 195 elements #To work in the same number array [1:]
x2 = np.array(df["x[kpc]"][BH2])[1:]
y2 = np.array(df["y[kpc]"][BH2])[1:]
z2 = np.array(df["z[kpc]"][BH2])[1:]
vx2 = np.array(df["vx[km/s]"][BH2])[1:]
vy2 = np.array(df["vy[km/s]"][BH2])[1:]
vz2 = np.array(df["vz[km/s]"][BH2])[1:]
t2 = np.array(df["time[Gyr]"][BH2])[1:]
m2 = np.array(df["mass[Msun]"][BH2])[1:]

# Finding the Time of Merging of BH1 and BH2
# Changing from Gigayears to Megayears ?
R_t1 = (t1-7.69864571)*(10**3) # 4628 elements
R_t2 = (t2-7.69864571)*(10**3)
t1_bcol = R_t1[:len(R_t2)] # 194 elements
t1_acol = R_t1[len(R_t2):] # 4434 elements
print("The merging of BH1 and BH2 is aproximately at:", "{:.2f}".format(R_t2[len(R_t2)-1]),"Megayears")

# Finding the last distance for merging
# Divide the BH1 position before and after merging
x1_bcol = x1[:len(R_t2)] 
y1_bcol = y1[:len(R_t2)]
z1_bcol = z1[:len(R_t2)]
x1_acol = x1[len(R_t2):]
y1_acol = y1[len(R_t2):]
z1_acol = z1[len(R_t2):]
dist_BH = np.sqrt((x2 - x1_bcol)**2 + (y2 - y1_bcol)**2 + (z2 - z1_bcol)**2)
print("The merging of BH1 and BH2 is aproximately at:", "{:.2f}".format(dist_BH[len(R_t2)-1]),"kpc")

# Making a 3D plot of BHs position before
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x1_bcol, y1_bcol, z1_bcol, label='BH1 Trajectory',color = 'red', linewidth = 0.5)
ax.scatter([x1_bcol[len(R_t2)-1]], [y1_bcol[len(R_t2)-1]], [z1_bcol[len(R_t2)-1]], color='black', s=50)
ax.scatter([x1_bcol[0]], [y1_bcol[0]], [z1_bcol[0]], color='green', s=25)
ax.plot(x2, y2, z2, label='BH2 Trajectory',color = 'blue', linewidth = 0.5)
ax.scatter([x2[len(R_t2)-1]], [y2[len(R_t2)-1]], [z2[len(R_t2)-1]], color='black', s=50)
ax.scatter([x2[0]], [y2[0]], [z2[0]], color='green', s=25)
ax.plot(x1_acol, y1_acol, z1_acol, label='SBH Trajectory',color = 'purple', linewidth = 0.5)
ax.set_xlabel('X [kpc]')
ax.set_ylabel('Y [kpc]')
ax.set_zlabel('Z [kpc]')
ax.grid(True)
ax.legend()

# Distance of BH1 and BH2 before merging
fig = plt.figure(figsize = (10,8))
ax1 = fig.add_subplot(2, 2, 1)  
ax1.plot(t1_bcol, x1_bcol, label = 'BH1 x-position', color = 'red')
ax1.plot(R_t2, x2, label = 'BH2 x-position', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1], label = 'Time of merging', linestyle = '--', color ='black')
ax1.set_xlim([0 , R_t2[len(R_t2)-1]])
plt.xlabel("Time [Megayears]")
plt.ylabel("BH distance [kpc]")
ax1.set_title('BH1 & BH2 x-position')
ax1.legend()

ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(t1_bcol, y1_bcol, label = 'BH1 y-position', color = 'red')
ax2.plot(R_t2, y2, label = 'BH2 y-position', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax2.set_xlim([0 , R_t2[len(R_t2)-1]])
plt.xlabel("Time [Megayears]")
plt.ylabel("BH distance [kpc]")
ax2.set_title('BH1 & BH2 y-position')
ax2.legend()

ax3 = fig.add_subplot(2, 2, 3) 
ax3.plot(t1_bcol, z1_bcol, label = 'BH1 z-position', color = 'red')
ax3.plot(R_t2, z2, label = 'BH2 z-position', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax3.set_xlim([0 , R_t2[len(R_t2)-1]])
plt.xlabel("Time [Megayears]")
plt.ylabel("BH distance [kpc]")
ax3.set_title('BH1 & BH2 z-position')
ax3.legend()

ax4 = fig.add_subplot(2, 2, 4) 
ax4.plot(t1_bcol , dist_BH , label = 'BH1 & BH2 Distance', linestyle = '-', color = 'purple')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax4.set_xlim([0 , R_t2[len(R_t2)-1]])
plt.xlabel("Time [Megayears]")
plt.ylabel("BH distance [kpc]")
plt.title("Distance between BH1 & BH2 around the time")
ax4.set_title('BH1 & BH2 position')
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

# Finding xyz-velocity before merging
fig = plt.figure(figsize = (10,8))
ax1 = fig.add_subplot(2, 2, 1)  
ax1.plot(t1_bcol, vx1[:len(R_t2)], label = 'BH1 x-velocity', color = 'red')
ax1.plot(R_t2, vx2, label = 'BH2 x-velocity', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1], label = 'Time of merging', linestyle = '--', color ='black')
ax1.set_xlim(left = 0)
plt.xlabel("Time [Megayears]")
plt.ylabel("BH velocity [km/s]")
ax1.set_title('BH1 & BH2 x-velocity')
ax1.legend()

ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(t1_bcol, vy1[:len(R_t2)], label = 'BH1 y-velocity', color = 'red')
ax2.plot(R_t2, vy2, label = 'BH2 y-velocity', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax2.set_xlim(left = 0)
plt.xlabel("Time [Megayears]")
plt.ylabel("BH velocity [km/s]")
ax2.set_title('BH1 & BH2 y-velocity')
ax2.legend()

ax3 = fig.add_subplot(2, 2, 3) 
ax3.plot(t1_bcol, vz1[:len(R_t2)], label = 'BH1 z-velocity', color = 'red')
ax3.plot(R_t2, vz2, label = 'BH2 z-velocity', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax3.set_xlim(left = 0)
plt.xlabel("Time [Megayears]")
plt.ylabel("BH velocity [km/s]")
ax3.set_title('BH1 & BH2 z-velocity')
ax3.legend()

# Plot of BH1 & BH2 velocity
v1_mag = np.sqrt((vx1**2)+(vy1**2)+(vz1**2))
v2_mag = np.sqrt((vx2**2)+(vy2**2)+(vz2**2))

ax4 = fig.add_subplot(2, 2, 4)
ax4.plot(t1_bcol , v1_mag[:len(R_t2)] , label = 'BH1 velocity', color = 'red')
ax4.plot(R_t2 , v2_mag , label = 'BH2 velocity', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
plt.xlabel("Time [Megayears]")
plt.ylabel("BH velocity [km/s]")
ax4.set_title("BH1 & BH2 velocity around time")
ax4.legend()
plt.tight_layout()

print("The average velocity of BH1:","{:.2f}".format(np.mean(v1_mag[:len(R_t2)])),"km/s")
print("The average velocity of BH2:","{:.2f}".format(np.mean(v2_mag)),"km/s")
print("The average velocity of SBH:","{:.2f}".format(np.mean(v1_mag[len(R_t2):])),"km/s")

# Finding xyz-velocity 
fig = plt.figure(figsize = (10,8))
ax1 = fig.add_subplot(2, 2, 1)  
ax1.plot(R_t1, vx1, label = 'BH1 x-velocity', color = 'red')
ax1.plot(R_t2, vx2, label = 'BH2 x-velocity', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1], label = 'Time of merging', linestyle = '--', color ='black')
ax1.set_xlim(left = 0)
plt.xlabel("Time [Megayears]")
plt.ylabel("BH velocity [km/s]")
ax1.set_title('BH1 & BH2 x-velocity')
ax1.legend()

ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(R_t1, vy1, label = 'BH1 y-velocity', color = 'red')
ax2.plot(R_t2, vy2, label = 'BH2 y-velocity', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax2.set_xlim(left = 0)
plt.xlabel("Time [Megayears]")
plt.ylabel("BH velocity [km/s]")
ax2.set_title('BH1 & BH2 y-velocity')
ax2.legend()

ax3 = fig.add_subplot(2, 2, 3) 
ax3.plot(R_t1, vz1, label = 'BH1 z-velocity', color = 'red')
ax3.plot(R_t2, vz2, label = 'BH2 z-velocity', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax3.set_xlim(left = 0)
plt.xlabel("Time [Megayears]")
plt.ylabel("BH velocity [km/s]")
ax3.set_title('BH1 & BH2 z-velocity')
ax3.legend()

# Plot of BH1 & BH2 velocity
ax4 = fig.add_subplot(2, 2, 4)
ax4.plot(R_t1 , v1_mag, label = 'BH1 velocity', color = 'red')
ax4.plot(R_t2 , v2_mag , label = 'BH2 velocity', color = 'blue')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
plt.xlabel("Time [Megayears]")
plt.ylabel("BH velocity [km/s]")
ax4.set_title("BH1 & BH2 velocity around time")
ax4.legend()
plt.tight_layout()

# Determining the xyz-momentum of BH1 & BH2
lmx1 = m1*vx1
lmy1 = m1*vy1
lmz1 = m1*vz1
lmx2 = m2*vx2
lmy2 = m2*vy2
lmz2 = m2*vz2
# Total xyz-momentum before and after merging
tot_lmx_bcol = m1[:len(R_t2)] *vx1[:len(R_t2)] + m2*vx2
tot_lmy_bcol = m1[:len(R_t2)] *vy1[:len(R_t2)] + m2*vy2
tot_lmz_bcol = m1[:len(R_t2)] *vz1[:len(R_t2)] + m2*vz2

tot_lmx_acol = m1[len(R_t2):]* vx1[len(R_t2):]
tot_lmy_acol = m1[len(R_t2):]* vy1[len(R_t2):]
tot_lmz_acol = m1[len(R_t2):]* vz1[len(R_t2):]

print("xBefore:","{:.2e}".format(lmx1[len(R_t2)-1]))
print("xBefore:","{:.2e}".format(lmx2[len(R_t2)-1]))
print("xAfter:","{:.2e}".format(lmx1[len(R_t2)]))
print("{:.2e}".format(tot_lmx_acol[0]))

print("yBefore:","{:.2e}".format(lmy1[len(R_t2)-1]))
print("yBefore:","{:.2e}".format(lmy2[len(R_t2)-1]))
print("yAfter:","{:.2e}".format(lmy1[len(R_t2)]))
print("{:.2e}".format(tot_lmy_acol[0]))

print("zBefore:","{:.2e}".format(lmz1[len(R_t2)-1]))
print("zBefore:","{:.2e}".format(lmz2[len(R_t2)-1]))
print("zAfter:","{:.2e}".format(lmz1[len(R_t2)]))
print("{:.2e}".format(tot_lmz_acol[0]))

plt.show()
