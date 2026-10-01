import numpy as np

def TROE(T,M,A0,B0,C0,A1,B1,C1,CF):

    K0 = (A0)*np.exp((B0)/T)*(T/300.0)**(C0)
    K1 = (A1)*np.exp((B1)/T)*(T/300.0)**(C1)
    K0 = K0*M
    KR = K0/K1
    NC = 0.75-1.27*(np.log10((CF)))
    F  = 10.0**(np.log10((CF))/(1+(np.log10(KR)/NC)**2))
    K = K0*K1*F/(K0+K1)
    
    return K

def ALK(T,M,A0,B0,C0,n,X0,Y0):
    K0 = 2.0E-22 * np.exp((n))
    K1 = 4.3E-1*(T/298.0)**(-8)
    K0 = K0*M
    K1 = K0/K1
    K2 = (K0/(1.0+K1))*(4.1E-1)**(1.0/(1.0+(np.log10(K1))**2))
    K3 = (C0)/(K2+(C0))
    K4 = (A0)*((X0)-T*(Y0))
    K = K4*np.exp((B0)/T)*K3

    return K

def TUN(T,A0,B0,C0):
    K = (A0)*np.exp(-(B0)/T)*np.exp((C0)/T**3)

    return K

def NIT(T,M,A0,B0,C0,n,X0,Y0):
        
    K0 = 2.0E-22 * np.exp((n))
    K1 = 4.3E-1*(T/298.0)**(-8)
    K0 = K0*M
    K1 = K0/K1
    K2 = (K0/(1.0+K1))*(4.1E-1)**(1.0/(1.0+(np.log10(K1))**2))
    K3 = K2/(K2+(C0))
    K4 = (A0)*((X0)-T*(Y0))
    K = K4*np.exp((B0)/T)*K3
    
    return K 

def EPO(T,M,A1,E1,M1):
    K1 = 1.0/((M1)*M+1.0)
    K = (A1)*np.exp((E1)/T)*K1
    return K


def ISO1(T,A0,A1,B0,C0,C1,D0,D1):

    K0 = (C0)*np.exp((C1)/T)*np.exp(1E8/T**3)
    K1 = (D0)*np.exp((D1)/T)
    K2 = (B0)*K0/(K0+K1)
    K = (A0) * np.exp((A1)/T)*(1.0 - K2)
    return K

def ISO2(T,A0,A1,B0,C0,C1,D0,D1):
    K0 = (C0)*np.exp((C1)/T)*np.exp(1E8/T**3)
    K1 = (D0)*np.exp((D1)/T)
    K2 = (B0)*K0/(K0+K1)
    K = (A0) * np.exp((A1)/T)*K2
    return K

def KCO(T,M,A1,M1):
    K = (A1) * (1.0 + (M / (M1)))
    return K