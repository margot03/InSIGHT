# pv modeling, 1/31/21
import math
import numpy as np
import matplotlib.pyplot as plt
# k = boltzmann's constnat = 1.3806503 * 10 ** (-23) J/K

# calculates light generated current
# Only necessary when Rs and Rp are not known
def pv_current(i_nominal, ki, temp_nominal, temp_actual, irradiation_nominal, irradiation_actual):
    delta_temp = temp_nominal - temp_actual
    i_current = (i_nominal + ki*delta_temp) * irradiation_actual / irradiation_nominal
    return i_current

# calculates the diode saturation current
def saturation_current(i_sc_nom, v_oc_nom, v_t_nom, a, temp_nom, temp_act, bandgap_energy):
    # equation 5
    # v_t_nom is thermal voltage of Ns cells connected in series
    q = 1.602 * math.pow(10, -19) # charge on electron
    k = 1.3806503 * math.pow(10, -23) # boltzmann constant
    i0_n = i_sc_nom / (math.exp(v_oc_nom / (a * v_t_nom)) - 1)
    temp_energy_charge = (q * bandgap_energy / (a * k)) * ((1/temp_nom) - (1/temp_act))
    I0 = i0_n * (math.pow(temp_nom / temp_act, 3) * math.exp(temp_energy_charge))
    return I0

# calculates total current leaving pv
def current(Ipv, I0, V, T, a, Rp=None, Rs=None, I_i=None):
    k = 1.3806503 * math.pow(10, -23) # boltzmann constant
    q = 1.602 * math.pow(10, -19) # charge on electron
    I = 0
    if Rp == None and Rs == None:
        # ideal equation
        I = Ipv - I0 * (math.exp(q * V / (a*k*T)) - 1) # take out 1000 later!! Just for math range err
    else:
        #non ideal equation
        I = Ipv - I0*(math.exp((V+Rs*I_i)/(Vt*a)-1)) - (V+Rs*I_i)/Rp
    return I

# plots different voltage current pairs on an IV curve
def iv_curve(Vocn, Ipv, I0, Tn, a):
    v = np.arange(0, Vocn, 0.1)
    i = []
    for V in v:
        i.append(current(Ipv, I0, V, Tn, a))
    figure = plt.figure(1)
    plt.plot(v, i)
    #plt.show()
    return figure

def main():
    v_oc_n = 0
    iv_curve(v_oc_n)

if __name__ == '__main__':
    main()
