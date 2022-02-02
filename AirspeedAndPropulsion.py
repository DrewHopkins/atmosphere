import numpy as np
from Atmosphere import seaLevel, constants, propertiesByAltitude
from scipy.optimize import fsolve
from sympy import Symbol

#make sure constants and seaLevel are called no matter what
gamma, R, g, S, keyAltitudes = constants()
t_0, p_0, d_0, a_0, v_0, beta = seaLevel()

def determineAltitude(d):
    h = layer0(d)
    if  h > 36089:
        h = layer1(d)
        if h > 65616:
            h = layer2(d)
            if h > 104987:
                print('This densisty is out of range')
    return h

# 0 <= h <= 36,089ft 
def layer0(d):
    h = round(((t_0/beta[0])*((d/d_0)**(1/(-g/(beta[2]*R)+1))-1)),2)
    if h == 0:
        h = abs(h)
    return h

# 36,089 < h <= 65,616ft
def layer1(d):  
    t_1, p_1, d_1, a_1, v_1 = propertiesByAltitude(keyAltitudes[1])
    h = -(t_1*R/g)*np.log(d/d_1)+keyAltitudes[1]
    return h

# 65,616 < h < 104,987ft
def layer2(d):  
    t_2, p_2, d_2, a_2, v_2 = propertiesByAltitude(keyAltitudes[2])
    h = (t_2/beta[2])*((d/d_2)**(1/(-g/(beta[2]*R)+1))-1)+keyAltitudes[2]
    return h

def KnotToFtPerSec(Vc):
    Vc = Vc * 1.68781 # ft/s
    return Vc
def FtPerSecToKnot(V):
    V = V / 1.68781 #knots
    return V

def TrueA_irspeed(Vc,sigma):
    V = Vc/sigma
    return V

def Sa_in_tVenan_t(M):  
    return (1+0.2*M**2)**3.5-1  

def LordRayleigh(M):
    return 166.92158*M**7/(7*M**2-1)**2.5-1

def globalVar(ratio):
    global ratioP 
    ratioP = ratio  
    
def Sa_in_tVenan_tFun_ction(mach):
    return (1+0.2*mach**2)**3.5-1 - ratioP

def LordRayleighFun_ction(mach):
    return 166.92158*mach**7/(7*mach**2-1)**2.5-1 - ratioP

def trueAirspeed(VcKnot, h):
    Vc = KnotToFtPerSec(VcKnot)
    # at sea level
    sigma = 1**2
    V = TrueA_irspeed(Vc, sigma)
    M1 = V / a_0
    if M1 < 1 :
        ratioP = Sa_in_tVenan_t(M1)
    if M1 > 1 :
        ratioP = LordRayleigh(M1)
    deltaP = ratioP*p_0
    # find the pressure at a certa_in alt
    temps, p, densities, speedOfSound, viscocity = propertiesByAltitude(h)
    ratioP = deltaP/p
    globalVar(ratioP)
    if ratioP < 0.893 :
        M2 = fsolve(Sa_in_tVenan_tFun_ction, 1)
        M2= M2[0]         
    if ratioP > 0.893 :
        M2 = fsolve(LordRayleighFun_ction, 1) 
        M2 = M2[0]    
    vFt = M2 * speedOfSound
    vKnot = FtPerSecToKnot(vFt)
    return M1, a_0, Vc, M2, speedOfSound ,vFt, vKnot


def ThrustSpecificConsumption(altitude, m_a):    
    t_a, p_a, d_a, a_a, visc_a, v, Cp = flightConditions(altitude,m_a)
    t_a = 504#############need to come off
    p_a = 1862#############need to come off
    d_a = p_a/(R*t_a)#############need to come off
    v = m_a * (gamma*R*t_a)**(.5)
    A_i, pr_C, n_C, gamma_C, pr_F, n_F, gamma_F, Beta, pr_B, n_B, t_04, Q_R, Cp_3, gamma_3, Cp_4, gamma_4, n_T, gamma_T, pr_AB, n_AB, t_06, Cp_5, gamma_5, Cp_6, gamma_6, n_N, gamma_N = engineParameters()
    t_0a, p_0a = freeStreamStagCond(t_a, m_a, p_a)
    t_02, p_02 = compressorFace(t_0a, p_0a, m_a)
    t_03, p_03 = burnerInlet(t_02, n_C, pr_C, gamma_C, p_02)
    t_08, p_08 = fanExit(t_02, n_F, pr_F, gamma_F, p_02)
    t_04, p_04 = burnerExit(t_04, p_03, pr_B)
    f = fuel(Cp_4, t_04, Cp_3, t_03, n_B, Q_R)
    t_05, p_05 = turbineExit(p_04, Cp, t_03, t_02, Beta, t_08, n_T, f, Cp_4, t_04, gamma_T)
    t_06, p_06, f_AB = afterburnerExit(t_06, p_05, pr_AB, Beta, Cp_6, Cp, t_08, Cp_5, t_05, n_AB, Q_R)
    p_e, v_e = nozzleExit(p_a, gamma_N, t_06, n_N, p_06)
    mdot_a, T, mdot_f, TSFC = perfromance(d_a, v, A_i, f, f_AB, v_e, a_a, m_a)
    return mdot_a, T, mdot_f, TSFC 

def flightConditions(alt, m_a):
    t_a, p_a, d_a, a_a, visc_a = propertiesByAltitude(alt)
    v = m_a * a_a
    Cp = 0.24
    return t_a, p_a, d_a, a_a, visc_a, v, Cp

def engineParameters():
    # inlet
    A_i = 6.6 # ft^2
    # compressor
    pr_C = 8
    n_C = 0.88
    gamma_C = gamma
    # fan
    pr_F = 1.5
    n_F = 0.92
    gamma_F = gamma
    Beta = 0.36
    # burner
    pr_B = 0.98
    n_B = 0.98
    T04 = 2850
    Q_R = 19000
    Cp_3 = 0.25
    gamma_3 = 1.38
    Cp_4 = 0.27
    gamma_4 = 1.33
    # turbine
    n_T = 0.95
    gamma_T = 1.35 
    # afterburner
    pr_AB = 0.9
    n_AB = 0.99
    t_06 = 3785
    Q_R = Q_R
    Cp_5 = 0.26
    gamma_5 = 1.35
    Cp_6 = 0.3
    gamma_6 = 1.3 #gamma?
    # nozzle 
    n_N = 0.99
    gamma_N = 1.35
    return A_i, pr_C, n_C, gamma_C, pr_F, n_F, gamma_F, Beta, pr_B, n_B, T04, Q_R, Cp_3, gamma_3, Cp_4, gamma_4, n_T, gamma_T, pr_AB, n_AB, t_06, Cp_5, gamma_5, Cp_6, gamma_6, n_N, gamma_N

def freeStreamStagCond(t_a, m_a, p_a):
    t_0a = t_a * (1+(gamma-1)/2*m_a**2)
    p_0a = p_a *(1+(gamma-1)/2*m_a**2)**(gamma/(gamma-1))
    return t_0a, p_0a

def compressorFace(t_0a, p_0a, m_a):
# also the diffuser exit
    t_02 = t_0a
    p_02 = p_0a*(1-0.5*(m_a/3)**3)
    return t_02, p_02

def burnerInlet(t_02, n_C, pr_C, gamma_C, p_02):
# also the compressor exit
    t_03 = t_02 * (1+1/n_C*(pr_C**((gamma_C-1)/gamma_C)-1))
    p_03 = p_02 * pr_C
    return t_03, p_03

def fanExit(t_02, n_F, pr_F, gamma_F, p_02):
    t_08 = t_02 * (1+1/n_F*(pr_F**((gamma_F-1)/(gamma_F))-1)) 
    p_08 = p_02 * pr_F
    return t_08, p_08

def burnerExit(t_04, p_03, pr_B):
# also the turbine inlet
    t_04 = t_04
    p_04 = p_03 * pr_B
    return t_04, p_04

def fuel(Cp_4, t_04, Cp_3, t_03, n_B, Q_R):
    f = (Cp_4 * t_04 - Cp_3 * t_03)/(n_B * Q_R - Cp_3 * t_03)
    return f 

def turbineExit(p_04, Cp, t_03, t_02, Beta, t_08, n_T, f, Cp_4, t_04, gammaT):
# also the afterburner inlet   
    p_05 = p_04 * (1 - (Cp*(t_03-t_02)+Beta*Cp*(t_08-t_02))/(n_T*(1+f)*Cp_4*t_04))**(gammaT/(gammaT-1))
    t_05 = t_04 * (1-n_T*(1-(p_05/p_04)**((gammaT-1)/(gammaT)))) 
    return t_05, p_05

def afterburnerExit(t_06, p_05, pr_AB, Beta, Cp_6, Cp, t_08, Cp_5, t_05, n_AB, Q_R):
# also the nozzle inlet
    t_06 = t_06
    p_06 = p_05 * pr_AB
    f_AB = ((1+Beta)*Cp_6*t_06-Beta*Cp*t_08-Cp_5*t_05)/(n_AB*Q_R-Cp_6*t_06)
    return t_06, p_06, f_AB

def nozzleExit(p_a, gamma_N, t_06, n_N, p_06):
    p_e = p_a
    v_e = (2*gamma_N/(gamma_N -1)*R*t_06*n_N*(1-(p_a/p_06)**((gamma_N-1)/gamma_N)))**(.5)
    return p_e, v_e

def perfromance(d_a, v, A_i, f, f_AB, v_e, a_a, m_a):
    mdot_a = d_a * a_a * A_i *(0.4136 + 0.03246* m_a + 0.1638*m_a**2+0.06432*m_a**3+0.001676*m_a**4)
    T = mdot_a*((1+f+f_AB)*v_e-v)
    mdot_f = mdot_a*(f+f_AB)* 32.174*3600
    TSFC = mdot_f / T
    return mdot_a, T, mdot_f, TSFC

####################################################################################
####################################################################################
####################################################################################

# d is density
# call any density between 0.000170836 and 0.0023769

def problem1(d):
    h = determineAltitude(d)
    print("--------- Problem 1 ---------")
    print("For a given densisty {} (slug/ft^3)".format(d))
    print("the altitude is {} (ft)".format(h))
    print("-----------------------------")

# VcKnot = 600 #knots
# h = 40000 # ft
def problem2(VcKnot, h):
    M1, a_0, Vc, M2, speedOfSound ,vFt, vKnot = trueAirspeed(VcKnot, h)
    print("--------- Problem 2 ---------")
    print("At Sea Level")    
    print("Mach: {}".format(round(M1,4)))
    print("Speed of Sound: {}".format(round(a_0,2)))
    print("Equivalent Airpeed: {} (ft/s) || {} (Knots)".format(round(Vc,2),round(VcKnot,2)))
    print("Moved to an altitude of {} (ft)".format(h))    
    print("Mach: {}".format(round(M2,4)))
    print("Speed of Sound: {}".format(round(speedOfSound,2)))
    print("True Airspeed: {} (ft/s) || {} (Knots)".format(round(vFt,2),round(vKnot,2)))
    print("-----------------------------")

# altitude = 0
# machs = [0,1]
def problem3(altitude, machs):

    print("--------- Problem 3 ---------")
    for m_a in machs:
        mdot_a, T, mdot_f, TSFC = ThrustSpecificConsumption(altitude, m_a)
        print("Altitude: {} (ft)".format(altitude))
        print("Mach: {}".format(m_a))
        print("mdot_a: {} (lbm/s)".format(mdot_a * 32.174))
        print("thrust: {} (lbf)".format(T))
        print("mdot_f: {} (lbm/hr)".format(mdot_f))
        print("Thrust Specific Fuel Consumption: {} (lbm/hr/lbf)".format(TSFC))
        print("-----------------------------")

####################################################################################
####################################################################################
####################################################################################

problem1(0.0023769)

problem2(600,40000)

problem3(0,[0,1])