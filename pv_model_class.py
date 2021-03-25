import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

class PvModel:

    def __init__(self, Iscn, Vocn, Imp, Vmp, Kv, Ki, Ns, Gn, G, Tn, T, Egap, err):
        # Initializes all standard values
        self.Iscn = Iscn
        self.Vocn = Vocn
        self.Imp = Imp
        self.Vmp = Vmp
        self.Kv = Kv
        self.Ki = Ki
        self.Ns = Ns
        self.Gn = Gn
        self.G = G
        self.Tn = Tn
        self.T = T
        self.Egap = Egap
        self.err = err
        self.q = 1.602 * 10 ** -19
        self.k = 1.3806503 * 10 ** -23
    
    def newton_raphson(self, v, Ipv, I0, Rs, Rp, Vt, Ns, a, err, i):
        """
        Given values needed to solve for the current using the Newton-Raphson
        method. Returns the current array.
        """
        _i = i[0]
        for idx in range(len(v)):
            _v = v[idx]
            _g = Ipv - I0 * (math.exp((_v + _i*Rs)/Vt/Ns/a)-1) - (_v + _i*Rs)/Rp - _i
            while abs(_g) > err:
                _g = Ipv - I0 * (math.exp((_v + _i*Rs)/Vt/Ns/a)-1) - (_v + _i*Rs)/Rp - _i
                _glin = -I0 * Rs/Vt/Ns/a * math.exp((_v + _i*Rs)/Vt/Ns/a) - Rs/Rp - 1
                _i = _i - _g/_glin
                i[idx] = _i
        return i

    
    def power(self, Rs, Rp, a):
        """
        Given values for the series and parallel resistances
        and the coefficient `a`, returns the voltage, current,
        and power arrays.
        """
        
        # set constants, etc
        q = self.q
        k = self.k
        Vocn = self.Vocn
        Iscn = self.Iscn
        Ns = self.Ns
        Egap = self.Egap
        Vmp = self.Vmp
        Imp = self.Imp
        Tn = self.Tn
        T = self.T
        Gn = self.Gn
        G = self.G
        err = self.err
        Ki = self.Ki
        Kv = self.Kv
        Pmax_e = Vmp * Imp

        Vtn = k * Tn / q
        Vt = k * T / q
        deltaT = T - Tn
        
        """
        The difference between the I0 equations is the
        distribution of the resistors term! One distributes it
        to only the leakage current and one distributes it to
        both the leakage current and the deltaT term. New method
        also extends the voltage to the open circuit voltage of
        differing quantities, not just ending at the nominal Voc.
        
        original method: I0 = (Iscn*(Rs+Rp)/Rp + Ki*deltaT) * G/Gn
        new method: I0 = (Iscn + Ki*deltaT) * (Rs+Rp)/Rp * G/Gn
        """
        # current equations
        Ipvn = Iscn * (Rs+Rp)/Rp 
        Ipv = (Ipvn + Ki*deltaT) * G/Gn
        Isc = (Iscn + Ki*deltaT) * G/Gn
        
        ## Calculating the leakage current I0
        # original method:
        # I0n = (Ipv - Vocn/Rp) / (math.exp(Vocn/Vt/a/Ns)-1)
        # updated I0 from original method:
        # I0 = (((Iscn + Ki*deltaT)*(Rs+Rp)/Rp)-((Vocn + Kv*deltaT)/Rp)) / (math.exp((Vocn + Kv*deltaT)/a/Vt/Ns) - 1)
    
        # new method:
        Isc_ = (  Iscn + Ki*deltaT  ) * G/Gn
        Voc_ = (  Vocn + Kv*deltaT  )
        Ipv_ = (Rs+Rp)/Rp * Isc_
        I0 = (Ipv_ - Voc_/Rp)/(math.exp(Voc_/Vt/a/Ns)-1)

        v = np.arange(0.1, Voc_, 0.1) #Vocn
        i = np.ones_like(v) * Iscn

        # solve for the current for all voltages
        _i = i[0]
        i = self.newton_raphson(v, Ipv, I0, Rs, Rp, Vt, Ns, a, err, i)
        
        # calculate power
        P = v * i
        
        return v, i, P
        
    def _plot(self, x, y, x_label, y_label, title):
        """
        Given data to plot, along with labels for the axes and title,
        creates a simple single figure plot of the data.
        Does not return the plot but does call plt.show()
        so plot shows up inline.
        """
        fig = plt.figure(1)
        plt.plot(x, y)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.show()
        #return fig
        
    def find_resistors(self):
        """
        Iteratively solve for the values of the series
        and parallel resistances (Rs, Rp) and the
        coefficient `a`.
        This should be done to initialize a model.
        """
        
        # set constants, etc
        q = self.q
        k = self.k
        Vocn = self.Vocn
        Iscn = self.Iscn
        Ns = self.Ns
        Egap = self.Egap
        Vmp = self.Vmp
        Imp = self.Imp
        Tn = self.Tn
        T = self.T
        Gn = self.Gn
        G = self.G
        err = self.err
        Ki = self.Ki
        Kv = self.Kv
        Pmax_e = Vmp * Imp

        #T = Tn
        #G = Gn
        Vtn = k * Tn / q
        Vt = k * T / q
        deltaT = T - Tn
        a = 1
        
        # set initial values of Rs and Rp
        Rs_max = (Vocn - Vmp) / Imp
        Rp_min = Vmp / (Iscn - Imp) - Rs_max
        Rp = Rp_min
        Rs_vals = np.arange(0, Rs_max, 0.001)
        
        # initialize voltage, current, etc
        itr = 0
        v = np.arange(0.1, Vocn, 0.1)
        i = np.zeros_like(v)
        p_dif = Pmax_e
        
        # iteratively solve for Rs, Rp, `a`
        while p_dif >= err and itr < len(v):
            Rs = Rs_vals[itr]

            Ipvn = Iscn * (Rs+Rp)/Rp
            Ipv = (Ipvn + Ki*deltaT) * G/Gn
            Isc = (Iscn + Ki*deltaT) * G/Gn
            #I0n = (Ipv - Vocn/Rp) / (math.exp(Vocn/Vt/a/Ns)-1)
            #I0 = I0n
            I0n = Iscn / (math.exp(Vocn/a/Vtn/Ns) - 1)
            I0 = (Iscn + Ki * deltaT) / (math.exp((Vocn + Kv * deltaT)/a/Vt/Ns) - 1)
            
            a = (Kv - Vocn/Tn) / (Ns * Vtn * (Ki/Ipvn - 3./Tn - Egap/(k*Tn**2)))

            Rp = Vmp * (Vmp + Imp*Rs) / (Vmp * Ipv - Vmp * I0 * math.exp((Vmp+Imp*Rs)/Vt/Ns/a)+Vmp*I0-Pmax_e)
            
            i = self.newton_raphson(v, Ipv, I0, Rs, Rp, Vt, Ns, a, err, i)
            
            P = v * i
            P_idx = np.argmax(P)
            Pmax_m = P[P_idx]

            p_dif = Pmax_m - Pmax_e
            itr += 1
        
        return Rs, Rp, a

def main():
    # assuming data is held in "pv_data.py"
    import pv_data as pv
    
    pv_mod = PvModel(pv.Iscn, pv.Vocn, pv.Imp, \
                    pv.Vmp, pv.Kv, pv.Ki, pv.Ns, \
                    pv.Gn, pv.Tn, pv.Egap, pv.err)
    Rs, Rp, a = pv_mod.find_resistors()
    v, i, P = pv_mod.power(Rs, Rp, a)
    Pmax = P[np.argmax(P)]
    print(Pmax)


if __name__ == '__main__':
    main()