import rk_public as plume_calc
import ice_ocean_public as ice_ocean
import numpy as np
from matplotlib import pylab as plt



##line plume tidewater glacier


# +++++++++ calculate glacier variables +++++++++      

sina_scal = 1 #slope of glacier,1=tidewater
Zx0 = -1000 #depth of glacier
values =5000 #regulates stepsize
X = np.linspace(0.0,(abs(Zx0)/sina_scal),values)
sina = 0*X+sina_scal 

#++++++++++++++++fjord parameters:

Ta_start = 1
Z = (Zx0 + X*sina)  # depth of glacier
Ta = X*0 + Ta_start #fjord temperature					
Sa = X*0 + 24.65 #fjord salinity
E0 = 0.1 #entrainment coefficient

#+++++++++++ plume:

T0 = 0.  #initial  plume temperature
S0 = 1e-06  #initial plume salinity


#line plume:

q = 0.1 #subglacial discharge per glacier width for line plume
U0 = ice_ocean.balance_velocity_lp(Sa[0],S0,Ta[0],T0,sina[0],E0,q) #balance velocity


y0 = np.zeros(4)  #set intial plume values:
D = q/U0
y0[0] = q
y0[1] = y0[0]*U0
y0[2] = y0[0]*T0
y0[3] = y0[0]*S0


melt_1d,y_1d,Sb_1d=plume_calc.plumes(X,Z,Ta,Sa,E0,sina,y0,'line') #calculates line plume

####cone plume

Q = 500 #[m^3/s] total discharge Q=(Pi/2)*D^2 U
S0 = 1e-06
T0 = ice_ocean.t_freeze(S0,Zx0)

U0 = ice_ocean.balance_velocity_cp(Sa[0],S0,Ta[0],T0,sina[0],E0,Q) #balance velocity


y0 = np.zeros(4)  #set intial plume values:
y0[0] = Q*2/np.pi
y0[1] = y0[0]*U0
y0[2] = y0[0]*T0
y0[3] = y0[0]*S0

melt_cone,y_1d,Sb_1d=plume_calc.plumes(X,Z,Ta,Sa,E0,sina,y0,'cone')   # calcutales cone plume



#+++++++++++++plotting
plt.plot(melt_1d *3600*24,Z, label='line')
plt.plot(melt_cone *3600*24,Z, label='cone')

plt.xlabel('melt [m/d]')
plt.ylabel('Z [m]')

plt.legend()
plt.show()

