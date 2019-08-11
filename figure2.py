import scipy.optimize as optimize
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from pylab import rcParams
import matplotlib.patches as patches


from scipy.optimize import fsolve

# include the elastic energy implictly
def f(P):
    global E1 , E2 , T1 , a1 , b , Cph , T0 , a , gamma
    return ( 0.5 * a * P**2.0 + 0.25 * b * P**4.0 - P * E1 )

def equation1(p):
    global E1 , E2 , T1 , a1 , b , Cph , T0 , a , gamma , EE
    PE = (p)
    return( a * PE + b *  PE**3.0  - EE  )

#%  %  %  initial temperature


# Define the figure size
rcParams['figure.figsize'] = 6, 6
rcParams['axes.linewidth'] = 1.8


#%  %  %  strength of the coupling between defect and normal polarization
gamma = 3.0
#%  %  %  parameters of P^2
a1 = 1.0 
#%  %  %  curie temperature
T0 = 1.0
#TT1 = np.linspace(0.1,2.0,20)
TT1 = np.array([0.8])
print TT1
#%  %  %  initial temperature
#T1 = 0.8
#%  %  %  parameters of P^4
b = 1.0 / 3.0
#%  %  %  parameters of P^6
c = 1.0 / 3.0
#%  %  %  heat capacity
Cph = 15.0
#######defect dipole
#Pd = 0.01
#Pd = 1.0 / 3.0 * 0.1
Pd = 0.0
#%  %  %  critical field
Ecp = 6.0 * b**2.0 *  np.sqrt( 3.0 * abs(b) / 10.0 / c ) / 25.0 / c


## electric field
E = np.linspace( -0.1 , 0.1 , 1000)
## applied electric field 
E_a = np.linspace( 0.1 , -0.05  , 1000)
#E = np.array( [0.0,0.4] )
# initialize the polarization
P1 = 0.0 * E
P2 = 0.0 * E
# initialize the polarization (applied)
P1_a = 0.0 * E_a
P2_a = 0.0 * E_a

j = 0








for T1 in TT1 :
  fig = plt.figure()
  ax1 = fig.add_subplot(111)
  print T1
  a = a1 * (T1 - T0) 
  PF = - np.sqrt(-a/(3*b))          # negative polarization 
  PD = np.sqrt(-a/(3*b))            # positve polarization 
  EF = a * PF + b * PF**3.0    # positive field 
  EB = EF
  ED = a * PD + b * PD**3.0    # negative field
  EE = ED
  PE = fsolve ( equation1, -1 ) 
  PB = -PE   
  # find the minimum value
  i = 0
  for E1 in E:
    P1[i] = optimize.fmin( f , [  1.0 ] )
    P2[i] = optimize.fmin( f , [  -1.0 ] )
    i = i + 1 
  # find the minimum value
  i = 0  
  for E1 in E_a:
    P1_a[i] = optimize.fmin( f , [  1.0 ] )
    P2_a[i] = optimize.fmin( f , [  -1.0 ] )
    i = i + 1
  delta_E =  E_a[1:]-E_a[:-1]
  delta_P1 = P1_a[1:]-P1_a[:-1]
  delta_P2 = P2_a[1:]-P2_a[:-1]    
  delta_P1_n = delta_P1 / np.sqrt(delta_E**2.0 + delta_P1**2.0)
  delta_P2_n = delta_P2 / np.sqrt(delta_E**2.0 + delta_P2**2.0)
  print delta_E.shape, delta_P1.shape, delta_P2.shape
  if (P1_a[1]*P2_a[1]<0 ) :
    delta_E_n = delta_E / np.sqrt(delta_E**2.0 + delta_P1**2.0)
    ax1.quiver(E_a[::100], P1_a[::100], delta_E_n[::100], delta_P1_n[::100], color='red', angles='xy', linewidths=(0.1,), edgecolors=('k'), headaxislength=5 , zorder = 2 , scale = 20.0)
    ax1.plot( max(E_a) , max(P1_a) , marker = 's' , markersize = 7 , color = 'black' , zorder = 2 )
  else :
    delta_E_n = delta_E / np.sqrt(delta_E**2.0 + delta_P1**2.0)
    ax1.quiver(E_a[::100], P1_a[::100], delta_E_n[::100], delta_P1_n[::100], color='red', angles='xy', linewidths=(0.1,), edgecolors=('k'), headaxislength=5 , zorder = 2 , scale = 20.0)
    ax1.plot( max(E_a) , max(P1_a) , marker = 's' , markersize = 7 , color = 'black' , zorder = 2 )
  ## plot the figure  
  line1 = plt.plot( E , P1 , color = 'blue' ,linewidth=3.0 , linestyle='--', zorder = 1)
  print 'P1=',P1 
  # line of reversible
  ax1.plot( Pd*gamma , 0.0 , 0 , 0 , marker = '+' , markersize = 20 , color = 'red' , linestyle='--', zorder = 2 )
  line2 = plt.plot( E , P2 , color = 'blue' , linewidth=3.0 , linestyle='--', zorder = 1)
  j = j + 1    
  # +++++++++++++++++++++++++++++++++++++++
  # reversible polarization
  P_ir = np.linspace( PD , PB, 200 )
  Pr1 = ( a * P_ir + b* P_ir**3.0 ) / ( a + 3.0 * b * PB**2.0 )
  print a + 3.0 * b * PB**2.0, PB
  # reversible electric field
  Er1 = a * Pr1 + 3.0 * b * PB**2.0 * Pr1 
  # plot
  liner = plt.plot( Er1 , Pr1 ,  color = 'black' , linestyle='--',linewidth=2.0 , zorder = 2 )
  yy2 = [0.0 , 0.0 , min(Pr1) , min(Pr1) ]  
  
  # +++++++++++++++++++++++++++++++++++++++
  # reversible polarization
  Pr = ( a * P1_a + b* P1_a**3.0 ) / ( a + 3.0 * b * PB**2.0 )
  print a + 3.0 * b * PB**2.0, PB
  # reversible electric field
  Er = a * Pr + 3.0 * b * PB**2.0 * Pr 
  print 'PE=',PE,'max(P1_a',max(P1_a)

  # the whole area is shadowed by grey
# +++++++++++++++++++++++++++++++++++++++
# temperature below TA1
# negative total work
  ax1.fill_between( E_a , P1_a , max(PB) , where = ( (E_a>=0.0) & (E_a<EB) ), interpolate=True , facecolor='gray')
  print min(P1_a)  
  ax1.fill_between( E_a , P1_a , min(P1_a) , where = ( (E_a<=0.0) & (P1_a>=min(P1_a)) ), interpolate=True , facecolor='gray')
# negative reversible work
  ax1.fill_between( Er1 , Pr1, max(Pr1), where= ( (Er1>=0.0)  )  , facecolor='green' , interpolate=True )   
# positive reversible work  
  ax1.fill_between( Er , Pr, min(Pr), where= ( (Er<=0.0) )  , facecolor='green' , interpolate=True )   
# positive total work
  #pcr = -0.1
    # plot figure
  #ax1.text( 0.03 , -0.8 , 'T='+str(T1) , fontsize = 20)
  #ax1.text( -0.1 , 0.9 , '(a)' , fontsize = 20)
  ax1.text( 0.095 , 0.82 , 'A' , fontsize = 20)
  ax1.text( 0.06 , 0.75 , 'B' , fontsize = 20),
  ax1.text( 0.06 , -0.03 , 'B$^\prime$' , fontsize = 20)
  ax1.plot( EB-0.001 , PB ,  marker = 's' , markersize = 7 , color = 'black' , zorder = 2 )  
  ax1.plot( EB-0.001 , 0.1 ,  marker = 's' , markersize = 7 , color = 'black' , zorder = 2 )  
  ax1.text( 0.008 , 0.7 , 'C' , fontsize = 20)
  ax1.plot( 0.0 , 0.76 ,  marker = 's' , markersize = 7 , color = 'black' , zorder = 2 )
  ax1.text( -0.073 , 0.5 , 'D' , fontsize = 20)
  ax1.text( 0.073 , -0.5 , 'F' , fontsize = 20)  
  #ax1.plot( EE+0.003 , PD+0.1 ,  marker = 's' , markersize = 7 , color = 'black' , zorder = 2 )
  ax1.plot( EE , PD ,  marker = 's' , markersize = 7 , color = 'black' , zorder = 2 )
  ax1.plot( -EE , -PD ,  marker = 's' , markersize = 7 , color = 'black' , zorder = 2 )
  ax1.text( -0.07 , -0.85 , 'E' , fontsize = 20)
  ax1.plot( EE+0.002 , PE ,  marker = 's' , markersize = 7 , color = 'black' , zorder = 2 )
  #ax1.text( -0.055 , 0.75 , 'M' , fontsize = 20)
  ax1.plot( min(Er) , min(Pr) ,  marker = 'o' , markersize = 7 , color = 'red' , zorder = 2 )
  ax1.text( -0.05 , 0.45 , 'N' , fontsize = 20)
  ax1.text( -0.05 , -0.2 , 'N$^\prime$' , fontsize = 20)
  ax1.plot( min(E_a) , min(P1_a) ,  marker = 'o' , markersize = 7 , color = 'red' , zorder = 2 )
  
  ax1.set_xlabel('Electric field',fontsize=20)
  ax1.set_ylabel('Polarization',fontsize=20)
  
    # all five areas
  ax1.text( -0.025 , 0.6 , '3' , fontsize = 20)
  ax1.text( 0.02 , 0.8 , '1' , fontsize = 20)  
  ax1.text( 0.01 , 0.03 , '2' , fontsize = 20)
  ax1.text( -0.03 , -0.1 , '4' , fontsize = 20) 

  
  # size of the ticklabels
  plt.setp( ax1.get_xticklabels(), fontsize=20)
  plt.setp( ax1.get_yticklabels(), fontsize=20)
  plt.xlim([-0.12 , 0.12 ])
  plt.ylim([-1.0,1.1])
  # the ticklabels show
  #plt.setp( ax1 , xticks = [-0.3,0.0,0.3,0.6,0.9] , yticks = [-1.6,-0.8,0.0,0.8,1.6] , xlim=[-0.45,1.0])
  # tick thickness
  ax1.tick_params(width=3,size=10)
  # S-type shape
  P_im = np.linspace( -max(P1_a) , max(P1_a) , 500 )
  E_im = a * P_im + b * P_im**3.0
  ax1.plot( E_im , P_im , color='cyan' , linewidth=3.0 , zorder=1 )  
  # save the figure
  plt.savefig('figure2.pdf',bbox_inches='tight', format='pdf')  



