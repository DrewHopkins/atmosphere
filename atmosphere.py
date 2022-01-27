import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

# constants
gamma = 1.4
R = 1716.56 # ft*lbf/slug*R
g = 32.174 # ft/s^2
S = 198.72 # R
keyAltitudes = [0, 36089.2, 65616.8, 104987, 154199.5, 167323, 232940]

def seaLevel():
    # set sea level condidtons
    t_0 = 518.67 # rankine
    p_0 = 2116.22 # lbf/ft^2
    d_0 = 0.0023769 # lb/ft^3
    a_0 = np.sqrt(gamma * R * t_0)
    v_0 = 3.7373E-7
    return t_0, p_0, d_0, a_0, v_0

def linearLappseRate(tn, pn, dn, vn, h, hn, bn):
    t = tn + bn*(h-hn)
    p = pn*(t/tn)**(-g/(bn*R))   
    d = dn*(t/tn)**(-g/(bn*R)-1)
    a = np.sqrt(gamma * R * t)
    v = vn*(t/tn)**(1.5)*((tn+S)/(t+S))
    return t, p, d, a, v


def isothermal(tn, pn, dn, vn, h, hn):
    t = tn
    p = pn*np.e**(-g*((h-hn)/(R*tn)))
    d = p*dn/pn
    a = np.sqrt(gamma * R * t)
    v = vn*(t/tn)**(1.5)*((tn+S)/(t+S))
    return t, p, d, a, v


def atmosphere(endAlt):
    temps = []
    pressures = []
    densities = []
    speedOfSound = []
    viscocity = []
    t_0, p_0, d_0, a_0, v_0 = seaLevel()
    t_1 = p_1 = d_1 = t_2 = p_2 = d_2 = t_3 = p_3 = d_3 =t_4 = p_4 = d_4 =t_5 = p_5 = d_5 = 0
    beta = [-0.00356616, 0, 0.00054864, 0.00153619, 0, -0.00153619, 0.00109728]
    points = 1000000 #accuarcy of graphs
    y = np.linspace(keyAltitudes[0], endAlt, points)
    for h in y:
        if h <= keyAltitudes[1]:
            hn = 0
            t, p, d, a, v = linearLappseRate(t_0,p_0,d_0,v_0,h,hn,beta[0])
            temps.append(t)
            pressures.append(p)
            densities.append(d)
            speedOfSound.append(a)
            viscocity.append(v)

        elif h <= keyAltitudes[2]: 
            if t_1 == p_1 == d_1 == 0:
                t_1 = temps[-1]
                p_1 = pressures[-1]
                d_1 = densities[-1]
                v_1 = viscocity[-1]
            hn =  keyAltitudes[1]
            t, p, d, a, v = isothermal(t_1,p_1,d_1,v_1,h,hn)
            temps.append(t)
            pressures.append(p)
            densities.append(d)
            speedOfSound.append(a)
            viscocity.append(v)

        elif h <= keyAltitudes[3]:
            if t_2 == p_2 == d_2 == 0:
                t_2= temps[-1]
                p_2 = pressures[-1]
                d_2 = densities[-1]
                v_2 = viscocity[-1]
            hn =  keyAltitudes[2]
            t, p, d, a, v = linearLappseRate(t_2,p_2,d_2,v_2,h,hn,beta[2])
            temps.append(t)
            pressures.append(p)
            densities.append(d)
            speedOfSound.append(a)
            viscocity.append(v)

        elif h <= keyAltitudes[4]: 
            if t_3 == p_3 == d_3== 0:
                t_3 = temps[-1]
                p_3 = pressures[-1]
                d_3 = densities[-1]
                v_3 = viscocity[-1]
            hn =  keyAltitudes[3]
            t, p, d, a, v = linearLappseRate(t_3,p_3,d_3,v_3,h,hn,beta[3])
            temps.append(t)
            pressures.append(p)
            densities.append(d)
            speedOfSound.append(a)
            viscocity.append(v)

        elif h <= keyAltitudes[5]: 
            if t_4 == p_4 == d_4 == 0:
                t_4 = temps[-1]
                p_4 = pressures[-1]
                d_4 = densities[-1] 
                v_4 = viscocity[-1]
            hn =  keyAltitudes[4]
            t, p, d, a, v = isothermal(t_4,p_4,d_4,v_4,h,hn)
            temps.append(t)
            pressures.append(p)
            densities.append(d)
            speedOfSound.append(a)
            viscocity.append(v)            

        elif h <= keyAltitudes[6]:
            if t_5 == p_5 == d_5 == 0:
                t_5 = temps[-1]
                p_5 = pressures[-1]
                d_5 = densities[-1] 
                v_5 = viscocity[-1]
            hn =  keyAltitudes[5]
            t, p, d, a, v = linearLappseRate(t_5,p_5,d_5,v_5,h,hn,beta[5])
            temps.append(t)
            pressures.append(p)
            densities.append(d)
            speedOfSound.append(a)
            viscocity.append(v)
    
    return y, temps, pressures, densities, speedOfSound, viscocity
    
def plots():
    y, temps, pressures, densities, speedOfSound, viscocity = atmosphere(endAlt = keyAltitudes[-1])
    fig, axis = plt.subplots(3, 2)
    axis[0, 0].plot(temps, y)
    axis[0, 0].set_title("Atmosphere Temperature")
    axis[0, 0].grid()
    axis[0, 0].set_xlabel('temperature (Rankine)')
    plt.setp(axis[0, 0].get_xticklabels(), rotation=30, horizontalalignment='right')

    axis[0, 1].plot(pressures, y)
    axis[0, 1].set_title("Atmosphere Pressure")
    axis[0, 1].grid()
    axis[0, 1].set_xlabel('pressure (lbf/ft^2)')
    plt.setp(axis[0, 1].get_xticklabels(), rotation=30, horizontalalignment='right')

    axis[1, 0].plot(densities, y)
    axis[1, 0].set_title("Atmosphere Density")
    axis[1, 0].grid()
    axis[1, 0].set_xlabel('density (lb/ft^3)')
    plt.setp(axis[1, 0].get_xticklabels(), rotation=30, horizontalalignment='right')

    axis[1, 1].plot(speedOfSound, y)
    axis[1, 1].set_title("Atmosphere Speed of Sound")
    axis[1, 1].grid()
    axis[1, 1].set_xlabel('speed of sound')
    plt.setp(axis[1, 1].get_xticklabels(), rotation=30, horizontalalignment='right')

    axis[2, 0].plot(viscocity, y)
    axis[2, 0].set_title("Atmospere Dynamic Viscocity")
    axis[2, 0].grid()
    axis[2, 0].set_xlabel('dynamic viscocity')
    plt.setp(axis[2, 0].get_xticklabels(), rotation=30, horizontalalignment='right')

    axis[2][1].axis('off')
    fig.canvas.set_window_title("Earth's Atmospheric Properties")
    plt.tight_layout()
    plt.show()

def propertiesByAltitude(currentAlt):
    y, temps, pressures, densities, speedOfSound, viscocity = atmosphere(currentAlt)
    print('The temperature is: ', temps[-1])
    print('The pressure is: ', pressures[-1])
    print('The density is: ', densities[-1])
    print('The speed of sound is: ', speedOfSound[-1])
    print('The dynamic viscocity is: ', viscocity[-1])

plots()