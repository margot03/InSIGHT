import math
import numpy as np
import matplotlib.pyplot as plt

class PvModel:

    def __init__(self, Iscn, Vocn, Imp, Vmp, Kv, Ki, Ns, Gn, Tn, Egap, err):
        # Initializes by finding values of Rs, Rp, a
        self.Iscn = Iscn
        self.Vocn = Vocn
        self.Imp = Imp
        self.Vmp = Vmp
        self.Kv = Kv
        self.Ki = Ki
        self.Ns = Ns
        self.Gn = Gn
        self.Tn = Tn
        self.Egap = Egap
        self.err = err
        self.q = 1.602 * 10 ** -19
        self.k = 1.3806503 * 10 ** -23

    def find_resistors(self):
        q = self.q # charge on electron
        k = self.k # boltzmann constant
        Vocn = self.Vocn
        Iscn = self.Iscn
        Ns = self.Ns
        Egap = self.Egap
        Vmp = self.Vmp
        Imp = self.Imp
        Tn = self.Tn
        Gn = self.Gn
        err = self.err
        Ki = self.Ki
        Kv = self.Kv
        Pmax_e = Vmp * Imp

        T = Tn
        G = Gn
        Vtn = k * Tn / q
        Vt = k * T / q
        deltaT = T - Tn
        a = 1

        Rs_max = (Vocn - Vmp) / Imp
        Rp_min = Vmp / (Iscn - Imp) - Rs_max

        Rp = Rp_min
        Rs_vals = np.arange(0, Rs_max, 0.001)
        itr = 0
        v = np.arange(0.1, Vocn, 0.1)
        #i = np.arange(Iscn, 0, -0.1)
        i = np.zeros_like(v)

        p_dif = Pmax_e

        while p_dif >= err and itr < len(v):
            #print("outer")
            Rs = Rs_vals[itr]

            Ipvn = Iscn * (Rs+Rp)/Rp
            Ipv = (Ipvn + Ki*deltaT) * G/Gn
            Isc = (Iscn + Ki*deltaT) * G/Gn
            I0n = (Ipv - Vocn/Rp) / (math.exp(Vocn/Vt/a/Ns)-1)
            I0 = I0n

            a = (Kv - Vocn/Tn) / (Ns * Vtn * (Ki/Ipvn - 3./Tn - Egap/(k*(Tn**2))))

            Rp = Vmp * (Vmp + Imp*Rs) / (Vmp * Ipv - Vmp* I0 * math.exp((Vmp+Imp*Rs)/Vt/Ns/a)+Vmp*I0-Pmax_e)

            _i = i[0]
            for idx in range(len(v)):
                _v = v[idx]
                _g = Ipv - I0 * (math.exp((_v + _i*Rs)/Vt/Ns/a)-1) - (_v + _i*Rs)/Rp - _i
                while abs(_g) > err:
                    _g = Ipv - I0 * (math.exp((_v + _i*Rs)/Vt/Ns/a)-1) - (_v + _i*Rs)/Rp - _i
                    _glin = -I0 * Rs/Vt/Ns/a * math.exp((_v + _i*Rs)/Vt/Ns/a) - Rs/Rp - 1
                    _i = _i - _g/_glin
                    i[idx] = _i

            P = np.zeros_like(v)
            for idx1 in range(len(v)):
                _v = v[idx1]
                _i = i[idx1]
                _p = (Ipv - I0 * (math.exp((_v + _i*Rs)/Vt/Ns/a)-1)-(_v + _i*Rs)/Rp)*_v
                P[idx1] = _p
            P_idx = np.argmax(P)
            Pmax_m = P[P_idx]

            p_dif = Pmax_m - Pmax_e
            itr += 1

        return Rs, Rp, a

def main():
    Iscn = 8.7 # nominal short circuit current
    Vocn = 37.7 # nominal open circuit voltage
    Imp = 8.2 # array current at MPP
    Vmp = 30.1 # array voltage at MPP
    Pmax_e = Vmp * Imp # experimental array max power output
    Kv = -0.32/100 * Vocn # voltage temperature coefficient
    Ki = -0.032/100 * Iscn # current temperature coefficient
    Ns = 60. # number of series cells
    Gn = 1000. # nominal irradiance
    Tn = 25. + 273.15 # nominal operating temperature

    Egap = 1.8 * 10 ** -19 # Bandgap of silicon (silicio cristalino)

    q = 1.602 * 10 ** -19 # charge on electron
    k = 1.3806503 * 10 ** -23 # boltzmann constant
    err = 0.001

    pv_mod = PvModel(Iscn, Vocn, Imp, Vmp, Kv, Ki, Ns, Gn, Tn, Egap, err)
    print(pv_mod.find_resistors())


if __name__ == '__main__':
    main()
