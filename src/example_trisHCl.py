# Example calculations:
# Tris-HCl buffer system
#
#

from pHcalc import IonAq
from pHcalc import AcidAq
from pHcalc import AcidGasEq
from pHcalc import System


# Tris
# 
# pKa = 8.07, T=298K
# pKa = 8.20, T=293K
#
# Tris.HCl      157.60 g/mol, solubility 667 g/l
# Tris base     121.14 g/mol, solubility 666 g/l (?)
#

TrisH = {'charge' : +1,
           'pKa'      : 8.20}

Cl = {'charge' : -1}


print('100 mol% Tris-HCl')
comp1 = AcidAq(**TrisH,
               conc = 0.100)

comp2 = IonAq(**Cl,
              conc = 0.100)

buff = System(comp1, comp2)

buff.pHsolve()
print('pH = ',buff.pH)
print()


print('85 mol% Tris-HCl + 15 mol% Tris free base')
comp1 = AcidAq(**TrisH,
               conc = 0.100)

comp2 = IonAq(**Cl,
              conc = 0.085)

buff = System(comp1, comp2)

buff.pHsolve()
print('pH = ',buff.pH)
print()

print('for 0.1M stock (85 mol% Tris-HCl + 15 mol% Tris free base):')
print(0.085 * 157.6, 'g/l Tris-HCl')
print(0.015 * 121.14, 'g/l Tris base')

