import numpy as np
import scipy as sp
from scipy import constants
import matplotlib.pyplot as plt


class PvModel:
    def __init__(self, Iscn, Vocn, Imp, Vmp, Kv, Ki, Ns, Gn, G, Tn, T, Egap, err, inverter, array_dim):
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
        self.inverter = inverter
        # array with [number of strings, number in series (per string)]
        self.array_dim = array_dim

        self.q = constants.e
        self.k = constants.Boltzmann

    def calculate_currents(self, Rs, Rp, a, Iscn, Vocn, Ki, Kv, Vt, Ns, G, Gn, deltaT):
        """Calculates and returns the different current values of interest. Takes
        different PV module specifications and weather data as input.
        """
        # current equations
        Ipvn = Iscn * (Rs + Rp) / Rp
        Ipv = (Ipvn + Ki * deltaT) * G / Gn
        Isc = (Iscn + Ki * deltaT) * G / Gn

        # new method:
        Isc_ = Iscn + Ki * deltaT
        Voc_ = Vocn + Kv * deltaT
        Ipv_ = (Rs + Rp) / Rp * Isc_
        I0 = (Ipv_ - Voc_ / Rp) / (np.exp(Voc_ / Vt / a / Ns) - 1)

        return Ipvn, Ipv, Isc, I0

    def newton_raphson(self, v, Ipv, I0, Rs, Rp, Vt, Ns, a, err, i):
        """Given values needed to solve for the current using the Newton-Raphson method.
        Returns the current array.
        """
        _i = i[0]
        for idx in range(len(v)):
            _v = v[idx]
            _g = (
                Ipv
                - I0 * (np.exp((_v + _i * Rs) / Vt / Ns / a) - 1)
                - (_v + _i * Rs) / Rp
                - _i
            )
            while abs(_g) > err:
                _g = (
                    Ipv
                    - I0 * (np.exp((_v + _i * Rs) / Vt / Ns / a) - 1)
                    - (_v + _i * Rs) / Rp
                    - _i
                )
                _glin = (
                    -I0 * Rs / Vt / Ns / a * np.exp((_v + _i * Rs) / Vt / Ns / a)
                    - Rs / Rp
                    - 1
                )
                _i = _i - _g / _glin
                i[idx] = _i

        return i
    
    def find_irradiance(self, DNI, DHI, GHI, rho, n, hour, lat, az_c, sigma):
        """ Given data for a particular location, calculates the value
        of the irradiance to be used in PV model calculations. This should be
        done for every time increment.
        TEXTBOOK EQUATIONS.
        """
        
        # declination
        decl = 23.45 * np.sin(360/365 * (n - 81)* np.pi/180)
        # hour angle
        H = 15 * (12 - hour)
        
        lat = lat * np.pi/180
        H = H * np.pi/180
        decl = decl * np.pi/180
        sigma = sigma * np.pi/180
        az_c = az_c * np.pi/180
        
        # solar altitude
        beta = np.arcsin(
            np.cos(lat) 
            * np.cos(decl) 
            * np.cos(H) 
            + np.sin(lat) 
            * np.sin(decl)
        )
        # solar azimuth
        az_s = np.arcsin(
            np.cos(decl) 
            * np.sin(H) 
            / np.cos(beta)
        )
        
        
        cos_theta = (
            np.cos(beta) 
            * np.cos(az_s - az_c) 
            * np.sin(sigma) 
            + np.sin(beta) 
            * np.cos(sigma)
        )
        B = DNI * cos_theta
        D = DHI * (1 + np.cos(sigma)) / 2
        R = GHI * rho * (1 - np.cos(sigma)) / 2
        G = B + D + R
        return G

    def find_resistors(self):
        """Iteratively solve for the values of the series and parallel resistances
        (Rs, Rp) and the diode ideality constant (a). This should be done to
        initialize a model.
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
        T = Tn
        Gn = self.Gn
        G = Gn
        err = self.err
        Ki = self.Ki
        Kv = self.Kv
        Pmax_e = Vmp * Imp

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
        v = np.arange(0.1, (Vocn + Kv * deltaT), 0.1)
        i = np.zeros_like(v)
        p_dif = Pmax_e

        # iteratively solve for Rs, Rp, `a`
        while p_dif >= err and itr < len(v):
            Rs = Rs_vals[itr]

            Ipvn, Ipv, Isc, I0 = self.calculate_currents(
                Rs, Rp, a, Iscn, Vocn, Ki, Kv, Vt, Ns, G, Gn, deltaT
            )

            a = (Kv - Vocn / Tn) / (
                Ns * Vtn * (Ki / Ipvn - 3.0 / Tn - Egap / (k * Tn ** 2))
            )

            Rp = (
                Vmp
                * (Vmp + Imp * Rs)
                / (
                    Vmp * Ipv
                    - Vmp * I0 * np.exp((Vmp + Imp * Rs) / Vt / Ns / a)
                    + Vmp * I0
                    - Pmax_e
                )
            )

            i = self.newton_raphson(v, Ipv, I0, Rs, Rp, Vt, Ns, a, err, i)

            P = v * i
            P_idx = np.argmax(P)
            Pmax_m = P[P_idx]

            p_dif = Pmax_m - Pmax_e
            itr += 1

        return Rs, Rp, a

    def calculate_power(self, Rs, Rp, a):
        """Given values for the series and parallel resistances and the diode ideality
        constant (a). Returns the voltage, current, and power arrays.
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

        # current equations
        Ipvn, Ipv, Isc, I0 = self.calculate_currents(
            Rs, Rp, a, Iscn, Vocn, Ki, Kv, Vt, Ns, G, Gn, deltaT
        )

        v = np.arange(0.1, (Vocn + Kv * deltaT), 0.1)
        i = np.ones_like(v) * Iscn

        # solve for the current for all voltages
        i = self.newton_raphson(v, Ipv, I0, Rs, Rp, Vt, Ns, a, err, i)

        # account for array dimensions
        i *= self.array_dim[0]
        v *= self.array_dim[1]
        # calculate power
        P = v * i
        
        # account for inverter efficiency
        P *= self.inverter
        
        return v, i, P

    def _plot(self, x, y, x_label, y_label, title):
        """Given data to plot, along with labels for the axes and title, creates a
        simple single figure plot of the data. Does not return the plot but does call
        plt.show() so plot shows up inline.
        """

        fig = plt.figure(1)
        plt.plot(x, y)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.axis((None, None, 0, None))
        plt.show()


def main():
    # assuming data is held in "pv_data.py"
    import pv_data as pv

    pv_mod = PvModel(
        pv.Iscn,
        pv.Vocn,
        pv.Imp,
        pv.Vmp,
        pv.Kv,
        pv.Ki,
        pv.Ns,
        pv.Gn,
        pv.Tn,
        pv.Egap,
        pv.err,
    )
    Rs, Rp, a = pv_mod.find_resistors()
    v, i, P = pv_mod.calculate_power(Rs, Rp, a)
    Pmax = P[np.argmax(P)]
    print(Pmax)


if __name__ == "__main__":
    main()
