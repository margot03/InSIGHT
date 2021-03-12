## Constants
q = 1.602 * 10 ** -19 # charge on electron
k = 1.3806503 * 10 ** -23 # boltzmann constant
Egap = 1.8 * 10 ** -19 # Bandgap of silicon (silicio cristalino)

## Data

Iscn = 8.21 #8.7 # nominal short circuit current
Vocn = 32.9 #37.7 # nominal open circuit voltage
Imp = 7.61 #8.2 # array current at MPP
Vmp = 26.3 #30.1 # array voltage at MPP
Pmax_e = Vmp * Imp # experimental array max power output
Kv = -0.32/100 * Vocn # voltage temperature coefficient
Ki = -0.032/100 * Iscn # current temperature coefficient

Ns = 60. # number of series cells
Gn = 1000. # nominal irradiance
Tn = 25. + 273.15 # nominal operating temperature

err = 0.001 # error for power difference, etc
