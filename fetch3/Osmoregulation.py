E=0.95 # % salt filtration efficiency 
c=200 #mm molar mass salt conc
Tw=293 #water temperature in kelvins
iv = 2. # van't hoff coefficient for NaCl
#E = 0.95 #filtration efficiency
def psi_s(E,c,Tw):
    iv=2 #van't hoff coeff
    R=8.314 #universal gas constant
    pi_s=c*iv*R*Tw
    return pi_s
