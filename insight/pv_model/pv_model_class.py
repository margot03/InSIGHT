import numpy as np
import scipy as sp
from scipy import constants
import matplotlib.pyplot as plt


class PvModel:
    def __init__(self, Iscn, Vocn, Imp, Vmp, Kv, Ki, Ns, Gn, G, Tn, T, Egap, err=0.0001, inverter_eff=0.95, Rs_increment=0.001, array_dim=[1,1], increment=0.1):
        """ Initializes all standard values used for the model.
        
        :param float Iscn: The nominal rated short circuit current for the PV module.
        :param float Vocn: The nominal rated open circuit voltage for the PV module.
        :param float Imp: The expected current at the maximum power point.
        :param float Vmp: The expected voltage at the maximum power point.
        :param float Kv: The rated voltage temperature coefficient.
        :param float Ki: The rated current temperature coefficient.
        :param int Ns: The number of PV cells in series for the PV module.
        :param float Gn: The nominal irradiance of the environment surrounding the PV module.
        :param float G: The actual irradiance of the environment surrounding the PV module.
        :param float Tn: The nominal temperature of the environment surrounding the PV module.
        :param float T: The actual temperature of the environment surrounding the PV module.
        :param float Egap: The bandgap energy of the semiconductor material that creates the
            PV module, usually for silicon.
        :param float err: The tolerance for how exact to make calculations to expected values.
            Default value of 0.0001 to maintain the accuracy of the model.
        :param float inverter_eff: The efficiency of the inverter. Default value of 0.95.
        :param float Rs_increment: The increment for generating possible values for the
            equivalent series resistance, which is then used to iteratively calculate all the
            equivalent resistance values for the PV module. Default value of 0.001 to maintan
            granularity.
        :param array array_dim: The dimensions of the PV array. Index 0 represents the number of
            modules in parallel and Index 1 represents the number of modules in series.
            Default value of [1,1].
        :param array increment: The increment used to advance quantities such as voltage.
            Default value of 0.1 due to it's granularity and common use.
        """
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
        self.inverter_eff = inverter_eff
        self.Rs_increment = Rs_increment
        # array with [number of strings, number in series (per string)]
        self.array_dim = array_dim
        self.increment = increment
        self.q = constants.e
        self.k = constants.Boltzmann

    def calculate_currents(self, Rs, Rp, a, Iscn, Vocn, Ki, Kv, Vt, Ns, G, Gn, deltaT):
        """Calculates and returns the different current values of interest. Takes
        different PV module specifications and weather data as input.
        
        :param float Rs: The equivalent series resistance of the PV module.
        :param float Rp: The equivalent parallel resistance of the PV module.
        :param float a: The diode ideality constant.
        :param float Iscn: The nominal short circuit current of the PV module.
        :param float Vocn: The nominal open circuit voltage of the PV module.
        :param float Ki: The rated current temperature coefficient.
        :param flaot Kv: The rated voltage temperature coefficient.
        :param float Vt: The thermal voltage of the PV module for Ns cells in series.
        :param int Ns: The number of PV cells in series for the PV module.
        :param float G: The irradiance of the environment surrounding the PV module.
        :param float Gn: The nominal irradiance of the environment surrounding the PV module,
            at 25 C degrees.
        :param float deltaT: The difference between the nominal operating temperature and the
            actual temeprature of the environment surrounding the PV module.
        :return: Lists for the nominal and realistic output currents (Ipvn and Ipv, respectively),
            the short circuit current impacted by environmental factors (Isc), and the 
            resulting saturation current (I0).
        """
        # current equations
        Ipvn = Iscn * (Rs + Rp) / Rp
        Ipv = (Ipvn + Ki * deltaT) * G / Gn
        Isc = (Iscn + Ki * deltaT) * G / Gn

        # Caclulate the saturation current, I0, without including the irradiance values
        # When G/Gn is included in the Isc used for I0 it causes the current values in
        # the irradiance demo file to go crazy!!
        Isc_noG = (Iscn + Ki * deltaT)
        Voc = Vocn + Kv * deltaT
        Ipv_ = (Rs + Rp) / Rp * Isc_noG # change to Isc to test with G/Gn
        I0 = (Ipv_ - Voc / Rp) / (np.exp(Voc / Vt / a / Ns) - 1)

        return Ipvn, Ipv, Isc, I0

    def newton_raphson(self, v, Ipv, I0, Rs, Rp, Vt, Ns, a, err, i):
        """Given values needed to solve for the current using the Newton-Raphson method.
        Returns the current array.
        
        :param array v: The array of voltages to evaluate for current values. Usually ranges
            from 0 to the nominal open circuit voltage.
        :param float Ipv: The current generated by the PV module before the current leaves
            the mdoule (excluding the diode saturation current).
        :param float I0: The diode saturation current which represents the amount of Ipv that
            is lost in the PV module (through recombination) before leaving the module.
        :param float Rs: The equivalent series resistance of the PV module.
        :param float Rp: The equivalent parallel resistance of the PV module.
        :param float Vt: The thermal voltage of the PV module for Ns cells in series.
        :param int Ns: The number of PV cells in series for the PV module.
        :param float a: The diode ideality constant.
        :param float err: The tolerance for how exact to make the calculated current.
            Usually very small (0.0001) to maintain the accuracy of the model.
        :param array i: The array of current values to iterate through and update for 
            the Newton-Raphson method.
        :return: The updated array of current values.
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
        of the irradiance to be used in the PV model calculations. This should be
        done for every time increment. These equations are from Box 4.1 and Box 4.2 
        from Masters' textbook, Renewable and Efficient Power Systems, 2e.
        
        :param float DNI: The Direct Normal Irradiance, TMY data, the measured amount of solar radiation that is
            perpendicular to the surface of the PV module (direct path from the sun).
        :param float DHI: Diffuse Horizontal Irradiance, TMY data, the measured amount of solar radiation that
            strikes the surface of the PV module without taking a direct path from the sun.
        :param float GHI: Global Horizontal Irradiance, TMY data, the measured amount of shortwave radiation that
            strikes the horizontal surface of the PV module.
        :param float rho: Ground reflectance coefficient.
        :param int n: The day of the year the TMY data was recorded (ex. July 3rd has n = 184)
        :param int hour: The hour of the day the TMY data was recorded.
        :param float lat: The latitude of the location where the TMY data was recorded.
        :param int az_c: The azimuth of the collector (the module used to collect the TMY data).
        :param float sigma: The tilt angle of the collector used to collect the TMY data.
        :return: The total irradiance received by the PV module at a single point in time, G.
        """
        
        # solar declination (depends on the day) - angle between the sun and the equator
        decl = 23.45 * np.sin(360/365 * (n - 81)* np.pi/180)
        # hour angle - angluar distance between the sun's meridian and the observer's meridian
        H = 15 * (12 - hour)
        
        # convert all variables to radians
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
        
        # horizontal component of the incidence angle between the sun and the collector 
        # for fixed orientation of collector
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
        and the diode ideality constant. This should be done to
        initialize a model.
        
        :return: The calculated equivalent parallel resistance (Rp)
            and the corresponding equivalent series resistance (Rs), 
            along with the diode ideality constant (a), all of which
            contribute to calculating more realistic current and power
            characteristics of the PV module.
        """

        # set constants and variables
        q = self.q
        k = self.k
        increment = self.increment
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
        Rs_increment = self.Rs_increment
        Vtn = k * Tn / q
        Vt = k * T / q
        deltaT = T - Tn
        
        # set initial values of Rs, Rp, and a
        Rs_max = (Vocn - Vmp) / Imp
        Rp_min = Vmp / (Iscn - Imp) - Rs_max
        Rp = Rp_min
        Rs_vals = np.arange(0, Rs_max, Rs_increment)
        a = 1
        
        # initialize voltage, current, etc
        itr = 0
        v = np.arange(0, (Vocn + Kv * deltaT), increment)
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
        :param float Rs: The equivalent series resistance.
        :param float Rp: The equivalent parallel resistance.
        :param float a: The diode ideality constant.
        :return: The voltage values (array v) and associated current (array i) and 
            power (array P) values that create the I-V and P-V characteristics of
            the PV module model.
        """

        # set constants, etc
        q = self.q
        k = self.k
        increment = self.increment
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

        v = np.arange(0, (Vocn + Kv * deltaT), increment)
        i = np.ones_like(v) * Iscn

        # solve for the current for all voltages
        i = self.newton_raphson(v, Ipv, I0, Rs, Rp, Vt, Ns, a, err, i)

        # account for array dimensions
        i *= self.array_dim[0]
        v *= self.array_dim[1]
        # calculate power
        P = v * i
        
        # account for inverter efficiency
        P *= self.inverter_eff
        
        return v, i, P

    def _plot(self, x, y, x_label, y_label, title):
        """Given data to plot, along with labels for the axes and title, creates a
        simple single figure plot of the data. Does not return the plot but does call
        plt.show() so plot shows up inline.
        :param array x: Data to plot on the x axis.
        :param array y: Data to plot on the y axis.
        :param string x_label: The label for the x axis.
        :param string y_label: The label for the y axs.
        :param string title: Title for the figure.
        """

        plt.plot(x, y)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.axis((0, None, 0, None))
        plt.show()
