# Calculation of solution pH in equilibrium with atmospheric carbon dioxide
# using the AcidGas class
#
from pHcalc import IonAq
from pHcalc import AcidAq
from pHcalc import AcidGasEq
from pHcalc import System

#
# CO2 / carbonic acid properties
#

H_s_CO2 = 0.03429
# Henry solubility [M atm-1] of CO2 in water. The present value is for
# pure water at 298K. Taken from ref. [2], Table II, converted from 
# molality to molarity, i.e. 0.997 * 10**-1.463

P_CO2 = 415e-6
P_CO2_1972 = 372.46e-6
# Partial pressure [atm] of carbon dioxide in the atmosphere. The
# present value is the global average atmospheric carbon dioxide in 2021
# as reported by NOAA’s Global Monitoring Lab. This value is rising over
# the years.
# https://www.climate.gov/news-features/understanding-climate/climate-change-atmospheric-carbon-dioxide

pKa_H2CO3 = [6.35, 10.33]
# The two pKa values of the carbonic acid system. These default values
# are for water at 298K. Refs [1] and [2].

# For now, all constants refer to STP (T=298K, P=1atm). 
#
# References:
# [1] H. S. Harned, S. R. Scholes Jr. 
#     "The Ionization Constant of HC03- from 0 to 50°." 
#     J. Am. Chem. Soc. 1941, 63,1706
# [2] H. S. Harned, R. Davis.
#     "The Ionization Constant of Carbonic Acid in Water and the Solubility
#     of Carbon Dioxide in Water and Aqueous Salt Solutions from 0 to 50°."
#     J. Am. Chem. Soc. 1943, 65 , 2030



#
# 1. Simple example: pure water
#
sys0 = System()
sys0.pHsolve()
print('pure water, closed system                     : pH = ', sys0.pH)

ambient = AcidGasEq(charge = 0,
                  pKa = pKa_H2CO3,
                  Hs = H_s_CO2,
                  Pgas = P_CO2)
syst = System(ambient)
syst.pHsolve()
print('pure water, air-equilibrated system (in 2021) : pH = ', syst.pH)

ambient1972 = AcidGasEq(charge = 0,
                      pKa = pKa_H2CO3,
                      Hs = H_s_CO2,
                      Pgas = P_CO2_1972)
syst1972 = System(ambient1972)
syst1972.pHsolve()
print('pure water, air-equilibrated system (in 1972) : pH = ', syst1972.pH)

print()


#
# 2. Benchmark/example calculations
#
# As benchmark calculations we used the on-line version of aqion which uses
# PHREEQC as a back-end for calculations.
#   https://www.aqion.onl/
#   https://www.usgs.gov/software/phreeqc-version-3
# The on-line results are identical to those obtained with the Windows version
# of aqion.
# Reference temperature: 298 K
#
# At higher ionic strengths, deviations appear between the values calculated by
# pHcalc and those from aqion/PHREEQC. This is due to the fact that PHREEQC
# uses activities and ionic strength-dependent equilibrium constants.
# pHcalc uses simple mass-action law with constant equilibrium constants and
# actual concentrations


# CO2 value used by aqion (watch out: aqion uses pCO2, not partial pressure)
aqionP_CO2 = 10**-3.408 # conversion of pCO2 into P_CO2 [atm]
print('Using partial CO2 pressure: {0:5.3e} atm (aqion value)'.\
      format(aqionP_CO2))

# Test cases
paramlist = [
    {'descr':       '0.1 M HCl(aq) in eq. with atmosphere',
     'ioncharge':  -1,
     'ionconc':     0.1,
     'aqion_pH':    1.08},
    {'descr':       '1 mM HCl(aq) in eq. with atmosphere',
     'ioncharge':  -1,
     'ionconc':     0.001,
     'aqion_pH':    3.02},
    {'descr':       'pure water in eq. with atmosphere',
     'ioncharge':   None,
     'ionconc':     None,
     'aqion_pH':    5.61},
    {'descr':       '1 mM NaOH(aq) in eq. with atmosphere',
     'ioncharge':   1,
     'ionconc':     0.001,
     'aqion_pH':    8.20},
    {'descr':       '10 mM NaOH(aq) in eq. with atmosphere',
     'ioncharge':   1,
     'ionconc':     0.01,
     'aqion_pH':    9.11},
    {'descr':       '0.1 M NaOH(aq) in eq. with atmosphere',
     'ioncharge':   1,
     'ionconc':     0.1,
     'aqion_pH':    9.71},
    ]

for param in paramlist:
    ambient = AcidGasEq(charge = 0,
                      pKa = pKa_H2CO3,
                      Hs = H_s_CO2,
                      Pgas = aqionP_CO2)
    if param['ionconc'] is None:
        syst = System(ambient)
    else:
        ion = IonAq(charge = param['ioncharge'],
                    conc = param['ionconc'])
        syst = System(ion, ambient)
    syst.pHsolve()
    DIC = ambient.conc # the total concentration of dissolved carbonic gas
    print('{0:40s} | DIC = {1:7.3e} M | pH = {2:5.2f} [aqion: {3:5.2f}]'.\
          format(param['descr'], DIC, syst.pH, param['aqion_pH']))
