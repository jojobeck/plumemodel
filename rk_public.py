import numpy as np
import ice_ocean_public as ice_ocean


def plumes(X,Z,Ta,Sa,E0,sina,y_0,name):
    '''
      calculates plume properties
      --------------

      Parameters
      -----------
      1d arrays:

      X: distance under shelf/ glacier front
      Z: glacier depth
      Ta: fjord temperature profile
      Sa: fjord salinity profile
      E0: entrainement paramter
      sina: glacier slope, 1= tidewater

      array (4,1):
      y_0 = intial plume properties

      str:
      name : plume kind, line or cone plume

      return: 
      ------------

      1d arrays:
      m : melt rate along glacier
      y: plume properties along glacier (q,qU,qT,qS)
      Sb: salinity at boundary layer

      I use the simple (non -adaptive step size) Runge Kutta method, and instead introduce Ufac (#HACK) to control of too big grid steps.'''

    j=0
    values=len(X)
    y=np.zeros((values,4))
    ynew=y_0
    y[0,:]=y_0
    Sb=np.zeros(values)
    U=np.zeros(values)
    Ufac=np.zeros(values)
    m=np.zeros(values)


    commands = {'line': ice_ocean.line_plume,
            'cone':ice_ocean.cone_plume}
    call_it=commands[name]



    while (j < len(X)-1) :
        h=X[j+1]-X[j]
        if(np.divide(y[j,1],y[j,0])<0):
            print('negative velocity at '+str(Z[j]))
            m[j:]=0
            break  
        if np.divide(y[j,1], y[j,0]) < 0:
            raise ValueError('should never happen')
        T=np.divide(y[j,2],y[j,0])
        S=np.divide(y[j,3],y[j,0])
        U[j]=y[j,1]/y[j,0]
        Sb[j]=ice_ocean.calculate_Sb(T,S,Z[j])
        
        if(j>2):#HACK!
            Ufac[j-1]=U[j]/U[j-1]
            Ufac[j-2]=U[j-1]/U[j-2]
            if(Ufac[j-1]/Ufac[j-2]>=2):
                print('to big step j',j)
                U[j]=0.
                break
        if(Sb[j]<=0.):
        	m[j:]=0.
        	print ('sb smaller zero')
        	break
        m[j]=ice_ocean.melt(Sb[j],y[j,:])  
        if(y[j,0]<=0):
        	m[j:]=0.
        	break
        if(name=='jenkins_area'):
        	k1=call_it(Ta[j],Z[j],Sb[j],Sa[j],E0,sina[j],y[j,:],m[j])*h
        	if((y[j,1]+k1[1]/2)<0):break
        	k2=call_it(Ta[j],Z[j],Sb[j],Sa[j],E0,sina[j],y[j,:]+k1/2,m[j])*h
        	if((y[j,1]+k2[1]/2)<0):break
        	k3=call_it(Ta[j],Z[j],Sb[j],Sa[j],E0,sina[j],y[j,:]+k2/2,m[j])*h
        	if((y[j,1]+k3[1])<0):break
        	k4=call_it(Ta[j],Z[j],Sb[j],Sa[j],E0,sina[j],y[j,:]+k3,m[j])*h
        	ynew=1./6*(k1+2*k2+2*k3+k4)
        	y[j+1,:]=y[j,:]+ynew
        	j=j+1
        else:
        	k1=call_it(Ta[j],Z[j],Sb[j],Sa[j],E0,sina[j],y[j,:],m[j])*h
        	k2=call_it(Ta[j],Z[j],Sb[j],Sa[j],E0,sina[j],y[j,:]+k1/2,m[j])*h
        	k3=call_it(Ta[j],Z[j],Sb[j],Sa[j],E0,sina[j],y[j,:]+k2/2,m[j])*h
        	k4=call_it(Ta[j],Z[j],Sb[j],Sa[j],E0,sina[j],y[j,:]+k3,m[j])*h
        	ynew=1./6*(k1+2*k2+2*k3+k4)
        	y[j+1,:]=y[j,:]+ynew
        	j=j+1
    return(m,y,Sb)


