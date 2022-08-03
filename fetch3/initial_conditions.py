"""
##################
Initial conditions
##################
Calculates initial conditions based on specified soil moisture and assuming hydrostatic conditions in the plant

Initial conditions in the soil layers
- initial soil moisture conditions [m3 m-3] for each soil layer are specified in the configuration file
- corresponding water potential [Pa] is calculated using the van genuchten equation

Initial conditions in the plant:
- potential at bottom of roots equals the soil potential at that depth
- potential at height z = potential at bottom of roots + rho*g*z, where z=0 is the bottom of the roots
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
    """
    Calculate initial water potential conditions

    Parameters
    ----------
    cfg : dataclass
        model configuration
    q_rain : np.ndarray
        array of rain data
    zind : dataclass
        z index dataclass

    Returns
    -------
    H_initial: np.ndarray
        initial values for water potential [Pa] over the concatenated z domain (soil, roots, xylem)
    Head_bottom_H: np.ndarray
        water potential [Pa] for the bottom boundary. size is len(number of timesteps)

    """

    # soil
    H_initial_soil = np.piecewise(
        zind.z_soil,
        [zind.z_soil <= cfg.clay_d, zind.z_soil > cfg.clay_d],
        [
            calc_potential_vangenuchten(
                cfg.initial_swc_clay,
                cfg.theta_R1,
                cfg.theta_S1,
                cfg.alpha_1,
                cfg.m_1,
                cfg.n_1,
                cfg.Rho,
                cfg.g,
            ),
            calc_potential_vangenuchten(
                cfg.initial_swc_sand,
                cfg.theta_R2,
                cfg.theta_S2,
                cfg.alpha_2,
                cfg.m_2,
                cfg.n_2,
                cfg.Rho,
                cfg.g,
            ),
        ],
    )

    # roots

    # z index where roots begin (round to get rid of floating point precision error so it matches the z array)
    z_root_start = np.round(cfg.Soil_depth - cfg.Root_depth, decimals=5)
    H_initial_root_bottom = H_initial_soil[zind.z_soil == z_root_start]
    H_initial_root = H_initial_root_bottom - (zind.z_root - z_root_start) * cfg.Rho * cfg.g

    # xylem
    H_initial_xylem = H_initial_root_bottom - (zind.z_upper - z_root_start) * cfg.Rho * cfg.g

    # concatenated array for z domain
    H_initial = np.concatenate((H_initial_soil, H_initial_root, H_initial_xylem))

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


    # set bottom boundary for initial condition
    if cfg.BottomBC == 0:
        H_initial[0] = Head_bottom_H[0]

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
