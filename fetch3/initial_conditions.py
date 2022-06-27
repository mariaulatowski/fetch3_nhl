"""
##################
Initial conditions
##################
"""
import numpy as np
#import Osmoregulation
#from Osmoregulation import psi_s
eff=0.97 #salt filtration efficiency 
c=479.1 #mol/m^3 of nacl salt conc (28ppt for mangroves port F, LA)
Tw=296.5 #water temperature in kelvins
iv = 2. # van't hoff coefficient for NaCl
def psi_s(E,c,Tw):
    iv=2 #van't hoff coeff
    R=8.314 #universal gas constant
    pi_s=-c*iv*R*Tw
    return pi_s
pi_s=psi_s(eff,c,Tw)

#######################################################################
#INITIAL CONDITIONS
#######################################################################
#soil initial conditions as described in the paper [VERMA et al., 2014]
 #osmotic potential in the soil due to salt conc

def initial_conditions(cfg, q_rain, zind):
    dz = cfg.dz

    initial_H=np.zeros(shape=zind.nz)

    factor_soil=(cfg.H_init_soilbottom-(cfg.H_init_soilmid))/(int((cfg.clay_d-cfg.cte_clay)/dz)) #factor for interpolation

    #soil
    for i in np.arange(0,len(zind.z_soil),1):
        if  0.0<=zind.z_soil[i]<=cfg.cte_clay :
            initial_H[i]=cfg.H_init_soilbottom
        if cfg.cte_clay<zind.z_soil[i]<=zind.z[zind.nz_clay]:
            initial_H[i]=initial_H[i-1]-factor_soil #factor for interpolation
        if cfg.clay_d<zind.z_soil[i]<= zind.z[zind.nz_r-1]:
            initial_H[i]=cfg.H_init_soilmid

    initial_H[zind.nz_s-1]=cfg.H_init_soilmid


    factor_xylem=(cfg.H_init_canopytop-(cfg.H_init_soilbottom))/((zind.z[-1]-zind.z[zind.nz_s])/dz)

    #roots and xylem
    initial_H[zind.nz_s]=cfg.H_init_soilbottom
    for i in np.arange(zind.nz_s+1,zind.nz,1):
        initial_H[i]=initial_H[i-1]+factor_xylem #meters


    #putting initial condition in Pascal
    H_initial=initial_H*cfg.g*cfg.Rho  #Pascals
    for i in np.arange(zind.nz_s+1,zind.nz,1):
      H_initial[i]=H_initial[i]+pi_s #[pascals]

    #osmotic potential in the soil due to salt conc
    #E=0.95 # % salt filtration efficiency 
    #c=200 #mm molar mass salt conc
    #iv=2 #van't hoff coeff for nacl
    #R=8.314 #universal gas constant j/mol K
    #Tw=293 #water temperature in kelvinss    psi_o=(E*c*iv*R*Tw)  #calculate osmotic potential at t=0 in the soil with c=200mm
    #psi_o=(c*iv*R*Tw)
    #osmo=np.zeros(shape=zind.nz) #create zero array of osmotic potential att=0
    #osmo=np.full((osmo),psi_o) #fill array with the initial osmotic potential
    #initial_H[zind.nz_s]=cfg.H_init_soilbottom
    #for i in np.arange(zind.nz_s+1,zind.nz,1):
      # H_initial[i]=H_initial[i-1]+psi_o #[pascals]



    ###########################################################################
    #BOTTOM BOUNDARY CONDITION FOR THE SOIL
    #The model contains different options, therefore this variable is created but
    #only used if you choose a  Dirichlet BC
    ######################################################################
    soil_bottom=np.zeros(shape=len(q_rain))
    for i in np.arange(0,len(q_rain),1):
        soil_bottom[i]=28      #0.28 m3/m3 fixed moisture according to VERMA ET AL., 2014

    #clay - van genuchten
    Head_bottom=((((cfg.theta_R1-cfg.theta_S1)/(cfg.theta_R1-(soil_bottom/100)))**(1/cfg.m_1)-1)**(1/cfg.n_1))/cfg.alpha_1
    Head_bottom_H=-Head_bottom*cfg.g*cfg.Rho  #Pa
    Head_bottom_H=np.flipud(Head_bottom_H) #model starts the simulation at the BOTTOM of the soil

    ############## inital condition #######################
    #setting profile for initial condition
    if cfg.BottomBC==0:
        H_initial[0]=Head_bottom_H[0]

    return H_initial, Head_bottom_H
