from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

# whisker order: D3, D2, D1, C3, C2, C1, B3, B2, B1
PSTHInVivoAvg = [0.089,0.217,0.356,0.217,0.333,0.339,0.100,0.200,0.156]
PSTHModelAvg = [0.099,0.328,0.374,0.323,0.472,0.440,0.113,0.116,0.264]
ModelPSTH_9_D3 = [0.115,0.120,0.290,0.260,0.640,0.445,0.155,0.140,0.385]
ModelPSTH_4_C1 = [0.055,0.385,0.535,0.430,0.365,0.775,0.165,0.215,0.485]
ModelPSTH_3_B3 = [0.200,0.620,0.505,0.350,0.730,0.460,0.055,0.050,0.235]

L2RF = [0.0086,0.0143,0.0243,0.0271,0.0229,0.0086,0.0057,0.0086,0.0214]
L3RF = [0.0001,0.0029,0.0286,0.0229,0.1314,0.0001,0.0029,0.0001,0.0143]
L4pyRF = [0.0250,0.0250,0.0750,0.0250,0.0250,0.0750,0.0750,0.0750,0.0001]
L4spRF = [0.0063,0.0001,0.0088,0.0175,0.1338,0.1525,0.0125,0.0088,0.0001]
L5stRF = [0.0046,0.0262,0.0238,0.0269,0.0254,0.0085,0.0215,0.0377,0.0377]
L5ttRF = [0.089,0.217,0.356,0.217,0.333,0.339,0.100,0.200,0.156]
L6ccRF = [0.0167,0.0167,0.0001,0.1167,0.4250,0.3583,0.0333,0.0333,0.0417]
L6ctRF = [0.0600,0.0001,0.0001,0.0001,0.0400,0.0100,0.0001,0.0100,0.0001]
L6invRF = [0.0001,0.0001,0.0500,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001]
VPMRF = [0.0001,0.0240,0.0001,0.1040,0.5449,0.0960,0.0001,0.0650,0.0480]

L34spinact = [0.067,0.2,0.205,0.083,0.157,0.239,0.077,0.07,0.176]
L5ttinact = [0.053,0.165,0.178,0.164,0.257,0.209,0.064,0.055,0.194]
L6ccinact = [0.019,0.04,0.092,0.042,0.127,0.125,0.016,0.021,0.047]
VPMinact = [0.089,0.276,0.252,0.249,0.053,0.203,0.11,0.077,0.234]

LFPL3RF = [1.6,3,1.8,1.8,7.25,3.4,1.8,3,1.7]
LFPL4RF = [3.8,4,3.1,5.7,8,3.6,3.4,4,2.4]
LFPL5RF = [3.3,3.3,2.2,5.1,6.2,3.6,3.8,3,1.8]
LFPL6RF = [3.4,3.4,2.3,4.8,4.1,3,3.8,2.4,3.1]

# reversed order to correlate syn. input with spiking output
VPMActiveSyn = [35232,62749,92588,122953,616327,222940,74005,164506,82682]
L6ccActiveSyn = [1150253,1101847,898218,1885572,2685525,1896130,1576878,2008002,1562160]

xpos, ypos = np.meshgrid(np.array(range(3)) + 0.3, np.array(range(3)) + 0.3)
xpos = xpos.flatten()
ypos = ypos.flatten()
zpos = np.zeros_like(xpos)
dx = 0.4 * np.ones_like(zpos)
dy = dx.copy()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.bar3d(xpos, ypos, zpos, dx, dy, PSTHModelAvg, color='b', edgecolor='k', zsort='average')
#faceColors = ['lightgrey','lightgrey','lightgrey','lightgrey',[0.2,0.2,0.2],'lightgrey','lightgrey','lightgrey','lightgrey']
dg = [0.2,0.2,0.2]
lg = [0.95,0.95,0.95]
faceColors = [lg,lg,lg,lg,dg,lg,lg,lg,lg]
#dg = 6*[[0.2,0.2,0.2]]
#lg = 6*[[0.95,0.95,0.95]]
#faceColors = 4*lg + dg + 4*lg
ax.bar3d(xpos, ypos, zpos, dx, dy, PSTHModelAvg, color=faceColors, edgecolor='k', alpha=1.0, linewidth=0.5, zsort='average')
ax.set_zlim3d([0,0.5]) # for IC types
#ax.set_zlim3d([0,0.55]) # for TC
ax.set_xlim3d([0.1,2.9])
ax.set_ylim3d([0.1,2.9])
#ztickLocations = 1.0*np.array([0,2,4,6,8])
ztickLocations = 0.1*np.array([0,1,2,3,4,5])
ax.set_zticks(ztickLocations)
#tickLocations = [0.5,1.5,2.5]
#ax.set_xticks(tickLocations)
#ax.set_xticklabels(['+1', 'PC', '-1'])
#ax.set_yticks(tickLocations)
#ax.set_yticklabels(['+1', 'PC', '-1'])
ax.set_xticks([])
ax.set_yticks([])
plt.show()

