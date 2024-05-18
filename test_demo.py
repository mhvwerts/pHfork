import numpy as np
from pHfork import AcidAq, IonAq, AcidGasEq, System

print('Demo and test of pHfork.')
print()
print('The graph windows that appear on screen need to be closed one after')
print('the other for the script to continue towards the end.')
print()


# NaCl
na = IonAq(charge=+1, conc=0.1)
cl = IonAq(charge=-1, conc=0.1)
nacl = System(na, cl)
nacl.pHsolve()
print('NaCl 0.1 M pH = ', nacl.pH)
print('            I = ', nacl.I, 'M')
print()


# NaOH, just need to define the amount of Na+, solver takes care of the
# rest.
a = IonAq(charge=+1, conc=0.1)
s = System(a)
s.pHsolve()
print('NaOH 0.1 M pH = ', s.pH)
print('            I = ', s.I, 'M')
print()


# NaOH, 0.1M but now in equilibrium with atmospheric CO2
# (formation of carbonates)
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

CO2atm = AcidGasEq(charge = 0,
                   pKa = pKa_H2CO3,
                   Hs = H_s_CO2,
                   Pgas = P_CO2)
s_CO2atm = System(a, CO2atm)
s_CO2atm.pHsolve()
print('NaOH 0.1 M, in equilibrium with atmospheric CO2, pH = ',
      s_CO2atm.pH)
print()

w = System()
w.pHsolve()
print('pure water, pH =', w.pH)
# print('             I =', w.I, 'M')
print()

w_CO2atm = System(CO2atm)
w_CO2atm.pHsolve()
print('pure water, in equil. with atmospheric CO2, pH =',
      w_CO2atm.pH)
# print('                                             I =',
#       w_CO2atm.I)
print()


# [HCl] = 1.0 x 10**-8   take into account autoprotolysis of water Kw!
# You just need to define the amount of Cl-. The solver will find the
# correct H3O+ concentration
b = IonAq(charge=-1, conc=1e-8)
s = System(b)
s.pHsolve()
print('HCl 1e-8 M pH = ', s.pH)
print()


# (NH4)3PO4
# H3PO4 input as the fully acidic species
a = AcidAq(pKa=[2.148, 7.198, 12.375], charge=0, conc=1.e-3)
# NH4+ again input as fully acidic species
# The concentration is 3x greater than the phosphoric acid
b = AcidAq(pKa=9.498, charge=1, conc=3.e-3)
#k = neutral(charge=1, conc=1.e-4)
#k = neutral(charge=-1, conc=3.e-3)
s = System(a, b)
print('(NH4)3PO4 1e-3 M')
# print report before
print(s)
print()
s.pHsolve()
print('(NH4)3PO4 1e-3 M pH = ', s.pH)
print()
# print report after
print(s)
print()


try:
    import matplotlib.pyplot as plt
except:
    print('Matplotlib not installed. Some examples not run.')
else:
    # Distribution diagram H3PO4
    a = AcidAq(pKa=[2.148, 7.198, 12.375], charge=0, conc=1.e-3)
    pH = np.linspace(0, 14, 1000)
    plt.figure(1)
    plt.clf()
    plt.title('Distribution diagram H3PO4(aq)')
    plt.plot(pH, a.alpha(pH))
    plt.xlabel('pH')
    plt.ylabel('ion fraction')
    plt.show()
    
    # Initial guess pH
    #
    # Here we use a "graphical" way of finding the pH at which the solution
    # composition is closest to electroneutrality.
    # 
    s = System(a)
    print('1mM H3PO4(aq)')
    diffs = s._diff_pos_neg(pH)
    plt.figure(2)
    plt.clf()
    plt.title('1mM H3PO4(aq) - initial guess pH')
    plt.semilogy(pH, np.abs(diffs))
    plt.xlabel('pH')
    plt.ylabel('charge imbalance')
    ixmin = np.argmin(np.abs(diffs))
    plt.plot(pH[ixmin], np.abs(diffs[ixmin]), 'o')
    print('"Graphical" guess pH = ',pH[ixmin])
    plt.show()
    s.pHsolve()
    print('pHsolve pH = ', s.pH)
    print()
    
    # Phosphoric Acid Titration Curve
    # First create a list of sodium hydroxide concentrations (titrant)
    Na_concs = np.linspace(1.e-8, 5.e-3, 500)
    # Here's our Acid
    H3PO4 = AcidAq(pKa=[2.148, 7.198, 12.375], charge=0, conc=1.e-3)
    phs = []
    for conc in Na_concs:
        # Create an inert Na+ with the concentration of the sodium
        # hydroxide titrant added
        Na = IonAq(charge=1, conc=conc)
        # Define the system and solve for the pH
        s = System(H3PO4, Na)
        s.pHsolve()
        phs.append(s.pH)
    plt.figure(3)
    plt.clf()
    plt.title('Phosphoric Acid Titration Curve')
    plt.plot(Na_concs, phs)
    plt.xlabel('total added NaOH concentration (M)')
    plt.ylabel('pH')
    plt.show()
