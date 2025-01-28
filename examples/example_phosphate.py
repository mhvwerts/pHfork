# Example calculations:
# phosphate buffer system
#
#

from pHfork import IonAq
from pHfork import AcidAq
from pHfork import AcidGasEq
from pHfork import System

# phosphate pH 7.4
# H2PO4(-) <-> HPO4(2-) + H+
# 
# pKa = 7.20,  room temp, infinite dilution
# pKa = 6.81, total phosphate = 0.1M
# 
#  see https://www.chem.fsu.edu/chemlab/Mastering/PhosphateBuffers.htm
#
# NaH2PO4.2H2O    156.01 g/mol
# Na2HPO4.12H2O   358.14 g/mol
#

H2PO4 = {'charge' : -1,
         'pKa'    : 6.9} # since we dilute below 0.1M,
                         #  apparent pKa chosen to be slightly higher
                         # TODO: be more precise about this!

Na = {'charge' : +1}



print('100 mol% NaH2PO4')
comp1 = AcidAq(**H2PO4,
               conc = 0.100)

comp2 = IonAq(**Na,
              conc = 0.100)

buff = System(comp1, comp2)

buff.pHsolve()
print('pH = ',buff.pH)
print()



print('25 mol% NaH2PO4 + 75 mol% Na2HPO4')
comp1 = AcidAq(**H2PO4,
               conc = 0.100)

comp2 = IonAq(**Na,
              conc = 0.175)

buff = System(comp1, comp2)

buff.pHsolve()
print('pH = ',buff.pH)
print()



print('for 0.1M stock (25 mol% NaH2PO4 + 75 mol% Na2HPO4): ')
print(0.025 * 156.01, 'g/l NaH2PO4.2H2O')
print(0.075 * 268.07, 'g/l Na2HPO4.7H2O')
