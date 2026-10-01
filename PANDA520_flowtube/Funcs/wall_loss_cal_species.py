import numpy as np
from molmass import Formula
from kinetics.diff_coef import diff_coef


def MeanFreePathBaron4_6(T,press):
    #This function calculates the mean free path (units: [nm]) of air molecules as a function of temperature and pressure.
    #Source: Baron&Willeke, 2001, eq. 4-6, p.65
    #mg, 17.08.2005
    # T            //temperature, [°C]
    #press   //pressure, [mbar=hectopascal]
    
    T=T+273.15         # //convert [°C] to [K]
    Lref=66.4                           #reference lambda at NTP (20 °C, 1013 mbar)#//note: TSI uses 66.5 at NTP
    
    S=110.4           #Sutherland constant for air
    Lambda=Lref*(1013/press)*(T/293.15)*(1+S/293.15)/(1+S/T)   
    return Lambda

def CunninghamSlipCorrBaron4_8(Lambda,Diam):
    #This function calculates the Cunningham slip correction factor for solid particles as a function of the mean free path of the air molecules.
    #slightly different parameter for oil droplets
    #Source: Baron, 2001, eq. 4-8, p. 66
    #Note: That is the parametrisation that I use for the SMPS inversion.
    #martin.gysel@psi.ch, 17.08.2005
    # Lambda              //mean free path of the air molecules, [nm]
    # Diam                    //particle diameter, [nm]
    
    Diam=Diam*1e-9#                     //convert from [nm] to [m]
    Lambda=Lambda*1e-9# //convert from [nm] to [m]
    
    Kn=2*Lambda/Diam
    a=1.142
    b=0.558
    g=0.999
    
    Cc=1+Kn*(a+b*np.exp(-g/Kn))
               
    return Cc

def ViscosityOfAirBaron4_10(TC):
    #This function calculates the dynamic viscosity of air (units: [kg/m/s]=[Pa*s]) as a function of air temperature.
    #Source: Baron&Willeke, 2001, eq. 4-10, p. 66
    #mg, 17.08.2005, tested against Ernest Weingartners "VOODOO"
    #variable TC         //temperature, [°C]
    S=110.4             # //Sutherland constant for air
    TC=TC+273.15 #       //convert [°C] -> [K]
    Vref=1.8325e-5  #            //reference viscosity at NTP (20 °C, 1013 mbar)#//note: TSI uses 1.8203e-5 at NTP
    Tref=293.15
        
    visc= Vref*(Tref+S)/(TC+S)*(TC/Tref)**(3/2)
              
    return visc

def MobilityOfParticle_TP(Diam, TC, press):
    #This function calculates the mobility B of a particle.
    #B = Cc/(3*pi*eta*Diam)                          //4-14, p. 67 in Baron, P.A., and K. Willeke, Wiley, New York, 2nd ed, 2001
    #return value: mobility B in [m/s/N]=[s/kg]
    #martin.gysel@psi.ch;08.05.2006; tested against TSI-aerosol-calc => within few %
    # Diam                    //particle diameter
    #//units: [nm]
    # TC                         //temperature
    #//units: [°C]      
    # press   //pressure
    #//units: [mbar=hectopascal]
    
    #convert diameter from nm to m
    DiamMet=Diam/1e9
    #Cunningham and Viscosity
    Lambda=MeanFreePathBaron4_6(TC,press)                   #units: [nm]
    Cc=CunninghamSlipCorrBaron4_8(Lambda,Diam)                   #    //units: [-]
    eta=ViscosityOfAirBaron4_10(TC)     #    //units: [kg/m/s]=[Pa*s]
    #mobility
    
    B=Cc/(3*np.pi*eta*DiamMet)
    return B

def DiffusionCoefficientParticle_TP(Diam, TC, press):
    ##
    #This function calculates the diffusion constant Ddiff of a particle.
    #Ddiff = kB*TC*Cc/(3*pi*eta*Diam)=kB*TK*B                              4-13, p.67 in Baron, P.A., and K. Willeke, Wiley, New York, 2nd ed, 2001
    #return value: diffusion coefficient Ddiff
    ## units: [m²/s]
    #martin.gysel@psi.ch; 08.05.2006; tested against TSI-aerosol-calc => within few  
    #units: [°C]      
    #units: [mbar=hectopascal]
    kB=1.3806e-23             #Boltzmann constant "k"           [J/K]
    #temperature in Kelvin
    TK=TC+273.16#[K]
    #mobility
    B=MobilityOfParticle_TP(Diam, TC, press)
    #diffusion coefficient
    
    Ddiff=kB*TK*B
    #finish procedure
               
    return Ddiff

def mobility_diameter(species,T, p,modelparams):
    """
    Return Dp in nanometers, given:
      Mw  = molecular weight [g/mol]
      rho = density [g/cm^3]
    """
    NA = 6.02214076e23  
    Mw = np.float32(diff_coef(species, T, p, 'air', modelparams)[1])
    # Mw = Formula(species_m).mass 
    rho = 1.4 # assumed density [g/cm³]
    factor = (6.0 * Mw) / (np.pi * rho * NA)      # [cm³ per molecule]
    dp_cm  = factor**(1/3)                       # [cm]
    Diam  = dp_cm * 1e7                         # 1 cm = 1e7 nm
    return Diam

def calculate_k_wall(species, TC,press,modelparams):
    #calculates wall losses in the CLOUD chamber
    # Input:
    #   Dp = particle mobility diameter (nanometer)
    #   T = temperature (Celsius)
    #   fanspeed = fan rotation speed in percent units (0-100)
    #   press   //pressure %//units: [mbar=hectopascal]
                                                                                
    #  Output:
    #    k_wall = wall loss coefficient (s^-1)
    #Diam = 0.636+0.3 # nm for H2SO4 particles, 0.3 nm is the correction for the mobility size
    Diam = mobility_diameter(species, TC,press,modelparams) # convert to nm
    k_wall=0.83; #this is actually the C factor in the wall loss equation. unit is m-1 s-0.5
    kwall=k_wall*np.sqrt(DiffusionCoefficientParticle_TP(Diam, TC, press))
    num_O = diff_coef(species, TC,press, 'air', modelparams)[2]
    if num_O > 3:
        kwall = 1e-3  # Adjust for oxygen content
    return kwall

def get_wall_loss(TC, press, modelparams):
    modelparams.wall_loss = []
    for i in modelparams.comp_namelist:
        if i in ['OH','HO2','H2O2','HSO3','H2SO4','SO3']:
            modelparams.wall_loss.append(2e-3)
        elif i not in modelparams.const_comp:
            kwall = calculate_k_wall(i, TC, press, modelparams)
            modelparams.wall_loss.append(kwall)
        else:
            modelparams.wall_loss.append(0.0)  # No wall loss for constant compounds
    
    return modelparams.wall_loss
    