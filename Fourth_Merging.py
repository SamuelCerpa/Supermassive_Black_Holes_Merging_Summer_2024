import pynbody
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

file = "/mnt/data0/jillian/samuel/newmerge/testj.BlackHoles"
s = pynbody.load("/mnt/data0/jillian/samuel/newmerge/testj.000012")

#file = "/mnt/data0/jillian/samuel/momcons/testj_momcons.BlackHoles"
#s = pynbody.load("/mnt/data0/jillian/samuel/momcons/testj_momcons.000012")

#file = "/mnt/data0/jillian/samuel/newmergerecoil/testj_recoil.BlackHoles"
#s = pynbody.load("/mnt/data0/jillian/samuel/newmergerecoil/testj_recoil.000012")

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

x1 = np.array(df["x[kpc]"][BH1])
y1 = np.array(df["y[kpc]"][BH1])
z1 = np.array(df["z[kpc]"][BH1])
vx1 = np.array(df["vx[km/s]"][BH1])
vy1 = np.array(df["vy[km/s]"][BH1]) 
vz1 = np.array(df["vz[km/s]"][BH1]) 
t1 = np.array(df["time[Gyr]"][BH1])
m1 = np.array(df["mass[Msun]"][BH1])
x2 = np.array(df["x[kpc]"][BH2])[1:]
y2 = np.array(df["y[kpc]"][BH2])[1:]
z2 = np.array(df["z[kpc]"][BH2])[1:]
vx2 = np.array(df["vx[km/s]"][BH2])[1:]
vy2 = np.array(df["vy[km/s]"][BH2])[1:] 
vz2 = np.array(df["vz[km/s]"][BH2])[1:]
t2 = np.array(df["time[Gyr]"][BH2])[1:]
m2 = np.array(df["mass[Msun]"][BH2])[1:]

# Finding the Time of Merging of BH1 and BH2
# Special Data
R_t1 = (t1-7.69864571)*(10**3) # 670 elements
R_t2 = (t2-7.69864571)*(10**3) # 195 elements #To work in the same number array [:1]
t1_bcol = R_t1[:len(R_t2)] # 194 elements
t1_acol = R_t1[len(R_t2):] # 476 elements
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

# Plot of BH1 & BH2 velocity before merging
v1_mag = np.sqrt((vx1**2)+(vy1**2)+(vz1**2))
v2_mag = np.sqrt((vx2**2)+(vy2**2)+(vz2**2))

# Determining the xyz-momentum of BH1 & BH2
lmx1 = m1*vx1
lmy1 = m1*vy1
lmz1 = m1*vz1
lmx2 = m2*vx2
lmy2 = m2*vy2
lmz2 = m2*vz2

# Total xyz-momentum before and after merging
tot_lmx_bcol = (m1[:len(R_t2)] * vx1[:len(R_t2)]) + (m2 * vx2)
tot_lmy_bcol = (m1[:len(R_t2)] * vy1[:len(R_t2)]) + (m2 * vy2)
tot_lmz_bcol = (m1[:len(R_t2)] * vz1[:len(R_t2)]) + (m2 * vz2)

tot_lmx_acol = m1[len(R_t2):] * vx1[len(R_t2):]
tot_lmy_acol = m1[len(R_t2):] * vy1[len(R_t2):]
tot_lmz_acol = m1[len(R_t2):] * vz1[len(R_t2):]

# Kinetic Energy Before and After Merging
KE1 = 0.5 * m1 * (v1_mag**2)
KE2 = 0.5 * m2 * (v2_mag**2)

KE1_bcol = KE1[:len(R_t2)]
KE1_acol = KE1[len(R_t2):]

tot_KE_bcol = KE2 + KE1_bcol
'''
# Velocity of BH1, BH2 & SBH
fig = plt.figure( figsize = (10,6))
plt.plot(R_t2 , v1_mag[:len(R_t2)], label = 'BH1 velocity', color = 'red')
plt.plot(R_t2 , v2_mag , label = 'BH2 velocity', color = 'blue')
plt.plot(t1_acol, v1_mag[len(R_t2):] , label = 'SBH velocity', color = 'green')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
plt.xlim(0 , R_t1[len(R_t1)-1])
plt.xlabel("Time [Megayears]")
plt.ylabel("Velocity [km/s]")
plt.title("BH1,BH2 & SBH Velocity")
plt.legend()
'''

#xyz momentum 


#xyz momentum 
fig = plt.figure(figsize = (14,8))
ax1 = fig.add_subplot(2, 3, 1)
ax1.plot(t1_bcol, tot_lmx_bcol, label = 'Total X-momentum before merging' , color = 'purple', linewidth = 2.0)
ax1.plot(t1_acol, tot_lmx_acol, label =  'Total X-momentum after merging' , color = 'green', linewidth = 2.0)
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax1.set_xlim( 0 , R_t1[len(R_t2)+200])
plt.xlabel("Time [Megayears]")
plt.ylabel("Momentum [Msun*km/s]")
ax1.set_title('Total X-momentum before and after Merging')
ax1.legend(loc = 'lower right',fontsize='small')

ax4 = fig.add_subplot(2, 3, 4)
ax4.plot(R_t2, lmx1[:len(R_t2)], label = 'X-momentum BH1' , color = 'red')
ax4.plot(R_t2, lmx2, label = 'Y-momentum BH2' , color = 'blue')
#plt.plot(t1_bcol, tot_lmx_bcol, label = 'Total X-momentum before merging' , color = 'purple')
plt.plot(t1_acol, tot_lmx_acol, label =  'Total X-momentum after merging' , color = 'green', linewidth = 2.0)
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax4.set_xlim(0,R_t1[len(R_t2)+200])
plt.xlabel("Time [Megayears]")
plt.ylabel("Momentum [Msun*km/s]")
ax4.set_title('BH1, BH2 & SBH X-momentum')
#plt.legend(loc = 'upper right')
ax4.legend(loc = 'lower right', fontsize='small')

ax2 = fig.add_subplot(2, 3, 2)
ax2.plot(t1_bcol, tot_lmy_bcol, label = 'Total Y-momentum before merging' , color = 'purple', linewidth = 2.0)
ax2.plot(t1_acol, tot_lmy_acol, label =  'Total Y-momentum after merging' , color = 'green', linewidth = 2.0)
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax2.set_xlim( 0 , R_t1[len(R_t2)+200])
plt.xlabel("Time [Megayears]")
plt.ylabel("Momentum [Msun*km/s]")
ax2.set_title('Total Y-momentum before and after Merging')
ax2.legend(loc = 'upper left',fontsize='small')

ax5 = fig.add_subplot(2, 3, 5)
ax5.plot(R_t2, lmy1[:len(R_t2)], label = 'Y-momentum BH1' , color = 'red')
ax5.plot(R_t2, lmy2, label = 'Y-momentum BH2' , color = 'blue')
#plt.plot(t1_bcol, tot_lmx_bcol, label = 'Total X-momentum before merging' , color = 'purple')
ax5.plot(t1_acol, tot_lmy_acol, label =  'Total Y-momentum after merging' , color = 'green', linewidth = 2.0)
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax5.set_xlim(0,R_t1[len(R_t2)+200])
plt.xlabel("Time [Megayears]")
plt.ylabel("Momentum [Msun*km/s]")
ax5.set_title('BH1, BH2 & SBH Y-momentum')
#plt.legend(loc = 'upper right')
ax5.legend(loc = 'upper left',fontsize='small')

ax3 = fig.add_subplot(2, 3, 3)
ax3.plot(t1_bcol, tot_lmz_bcol, label = 'Total Z-momentum before merging' , color = 'purple', linewidth = 2.0)
ax3.plot(t1_acol, tot_lmz_acol, label =  'Total Z-momentum after merging' , color = 'green', linewidth = 2.0)
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax3.set_xlim( 0 , R_t1[len(R_t2)+200])
plt.xlabel("Time [Megayears]")
plt.ylabel("Momentum [Msun*km/s]")
ax3.set_title('Total Z-momentum before and after Merging')
ax3.legend(loc = 'upper left',fontsize='small')

ax6 = fig.add_subplot(2, 3, 6)
ax6.plot(R_t2, lmz1[:len(R_t2)], label = 'Z-momentum BH1' , color = 'red')
ax6.plot(R_t2, lmz2, label = 'Z-momentum BH2' , color = 'blue')
#plt.plot(t1_bcol, tot_lmx_bcol, label = 'Total X-momentum before merging' , color = 'purple')
ax6.plot(t1_acol, tot_lmz_acol, label =  'Total Z-momentum after merging' , color = 'green', linewidth = 2.0)
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
ax6.set_xlim(0,R_t1[len(R_t2)+200])
plt.xlabel("Time [Megayears]")
plt.ylabel("Momentum [Msun*km/s]")
ax6.set_title('BH1, BH2 & SBH Z-momentum')
#plt.legend(loc = 'upper right')
ax6.legend(loc = 'upper left',fontsize='small')
plt.tight_layout()

'''
# Plot for BH1 & BH2 xyz-velocities before merging
fig = plt.figure(figsize = (9,8))
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

fig = plt.figure(figsize = (10,6))
plt.plot(R_t2, v1_mag[:len(R_t2)], label = 'BH1 velocity', color = 'red')
plt.plot(R_t2, v2_mag, label = 'BH2 velocity', color = 'blue')
plt.plot(t1_acol, v1_mag[len(R_t2):], label = 'SBH velocity', color = 'purple')
plt.axvline(R_t2[len(R_t2)-1] , label = 'Time of merging', linestyle = '--', color ='black')
plt.xlim(0,R_t1[len(R_t1)-1])
plt.xlabel("Time [Megayears]")
plt.ylabel("Velocity [km/s]")
plt.title('BH1 & BH2 Velocity around the time')
plt.legend()
'''
plt.show()