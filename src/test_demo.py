import numpy as np
from pHcalc import AcidAq, IonAq, System

print('Demo and test of pHcalc.')
print()
print('The graph windows that appear on screen need to be closed one after')
print('the other for the script to continue towards the end.')
print()

# KOH, just need to define the amount of K+, solver takes care of the
# rest.
a = IonAq(charge=+1, conc=0.1)
s = System(a)
s.pHsolve()
print('NaOH 0.1 M pH = ', s.pH)
print()

# [HCl] = 1.0 x 10**-8   aka the undergrad nightmare
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
s.pHsolve()
print('(NH4)3PO4 1e-3 M pH = ', s.pH)
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
    # Here we use a graphical way of finding the pH at which the solution
    # composition is closest to electroneutrality.
    # 
    s = System(a)
    diffs = s._diff_pos_neg(pH)
    plt.figure(2)
    plt.clf()
    plt.title('1mM H3PO4(aq) - initial guess pH')
    plt.semilogy(pH, np.abs(diffs))
    plt.xlabel('pH')
    plt.ylabel('charge imbalance')
    ixmin = np.argmin(np.abs(diffs))
    plt.plot(pH[ixmin], np.abs(diffs[ixmin]), 'o')
    print('Graphical guess pH = ',pH[ixmin])
    plt.show()
    s.pHsolve()
    print('pHsolve pH = ', s.pH)

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
