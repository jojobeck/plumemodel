import numpy as np

Cd      =  2.5e-3     # Drag coefficient
Cd12GammaTS =  5.9e-4   # Stanton number for temperature
Cd12GammaT  =  1.1e-3# thermal number for temperature
Cd12GammaS  =  3.1e-5 # haline Stanton number for temperature
lambda1     = -5.73e-02 # [degC] seawater freezing point: slope
lambda2     =  8.32e-02 # [degC] seawater freezing point: offset
lambda3     =  7.61e-04 # [degC/m] seawater freezing point: depth dependence
beta_S      =  7.86e-04 # haline density coefficient 
beta_T      =  3.87e-05 # [1/degC]thermal density coefficient
g       =  9.81     # [m/s^2] gravity acceleration on earth
c       =  3.974e+03    # [J/kg/degC] specific heat capacity of water
ci      =  2.009e+03    # [J/kg/degC] specific heat capacity of ice
L       =  3.35e+05 # [J/kg] latent heat of ice
Ti       =  -15.   # [degC] temperature of ice
sina        =  1.   #defalut slope of glacier base
Si      =  0.0      # [psu] salinity of ice/glacier
Pi=np.pi

def line_plume(Ta,Z,Sb,Sa,E0,sina,y,m_rk):
    Sb_rk=Sb
    return np.array([((E0*sina)*y[1]/y[0]+m_rk),
            ((g*sina)*y[0]**2/y[1]*(beta_S*(Sa-y[3]/y[0])-beta_T*(Ta-y[2]/y[0]))-Cd*y[1]*y[1]/(y[0]*y[0])),
            ((E0*sina*y[1]/y[0])*Ta+m_rk*(lambda1*Sb_rk+lambda2+lambda3*Z)-Cd12GammaT*y[1]/y[0]*(y[2]/y[0]-(lambda1*Sb_rk+lambda2+lambda3*Z))),
            ((E0*sina*y[1]/y[0])*Sa+m_rk*Sb_rk-Cd12GammaS*y[1]/y[0]*(y[3]/y[0]-Sb_rk) )])

def calculate_Sb(T,S,Z):
    A=(c*Cd12GammaT-ci*Cd12GammaS)*lambda1
    B=Cd12GammaS*ci*(Ti+lambda1*S-lambda2-lambda3*Z)-Cd12GammaT*c*(T+lambda1*Si-lambda2-lambda3*Z)-Cd12GammaS*L
    C=(Cd12GammaS*ci*(lambda2+lambda3*Z-Ti)+Cd12GammaS*L)*S-Cd12GammaT*(lambda2+lambda3*Z-T)*Si
    p=B/A
    q=C/A
    if((p/2.)**2-q<0.):
        print(Ti,Si,T,S)
        Sb=-100. # artifical error:set sb <0 to stop melt loop

    else:
        Sb=-p/2+np.abs(np.sqrt((p/2)**2-q))
    return Sb

def melt(Sb,y2):
    m=Cd12GammaS*y2[1]/y2[0]*(y2[3]/y2[0]-Sb)/(Sb-Si)
    return m

def balance_velocity_lp(Sa,S0,Ta,T0,sina,E0,DU0):
    rho=beta_S*(Sa-S0)-beta_T*(Ta-T0)
    U0jenkins=(sina*9.81*rho*DU0/(E0*sina+Cd))**(1./3)
    return(U0jenkins)

def balance_velocity_cp(Sa,S0,Ta,T0,sina,E0,DU0):
    rho=beta_S*(Sa-S0)-beta_T*(Ta-T0)
    U0jenkins=abs((9.81*sina*rho*DU0/((DU0*2*Pi)**(1/2.)*(E0*sina+2*Cd/Pi)))**(2/5.))
    return U0jenkins

def cone_plume(Ta,Z,Sb,Sa,E0,sina,y,m_rk):
    Sb_rk=Sb
    the_array= np.array([((2*E0*sina)*abs(y[1]**(1./2))+m_rk*4./Pi*y[0]/abs(y[1]**(1./2))),
        ((g*sina)*(y[0]**2)/y[1]*(beta_S*(Sa-y[3]/y[0])-beta_T*(Ta-y[2]/y[0]))-4*Cd/Pi*abs(y[1]**(3/2.))/y[0]),
        (2*E0*sina*abs(y[1]**(1./2))*Ta+4./Pi*m_rk*(lambda1*Sb_rk+lambda2+lambda3*Z)*y[0]/abs(y[1]**(1./2))- \
            4/Pi*Cd12GammaT*abs(y[1]**(1./2))*(y[2]/y[0]-(lambda1*Sb_rk+lambda2+lambda3*Z))),
        (2*E0*sina*abs(y[1]**(1./2))*Sa+4./Pi*m_rk*Sb_rk*y[0]/abs(y[1]**(1./2))-\
            4/Pi*Cd12GammaS*abs(y[1]**(1./2))*(y[3]/y[0]-Sb_rk) )])
    return the_array
                        
def t_freeze(S,Z):
    Tf=lambda1*S+lambda2+lambda3*Z
    return Tf

