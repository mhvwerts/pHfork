# Example calculations 
# concerning air-equilibrated carbonate/hydroxide solutions
#
#

from pHcalc import Inert as IonAq
from pHcalc import Acid as AcidAq
from pHcalc import AcidGas as AcidGasEq
from pHcalc import System

# Sodium carbonate

H2CO3 = {'charge': 0,
         'pKa'   : [6.35, 10.33]}

H2CO3_aireq = H2CO3 | {'Hs'    : 0.03429,
                       'Pgas'  : 415e-6}

carb = AcidAq(**H2CO3,
              conc = 0.1)

Na = IonAq(charge = +1,
           conc = 0.2)

carbCO2 = AcidGasEq(**H2CO3_aireq)
# Interestingly, to calculate air-equilibrated carbonate buffer, we do a 
# calculation that is equivalent, even identical to air-equilibrated NaOH(aq)!
# This concerns 0.2M NaOH, since 0.1M Na2CO3 contains 0.2M Na+

print()
print("Na2CO3 0.1M closed vs air-equilibrated")

sys0 = System(Na, carb)
sys0.pHsolve()
print(sys0.pH)

sys1 = System(Na, carbCO2)
sys1.pHsolve()
print(sys1.pH)



print()
print("'Air-insensitive' carbonate buffer 0.1M")

carb = AcidAq(**H2CO3,
              conc = 0.1)

Na = IonAq(charge = +1,
           conc = 0.1331)

carbCO2 = AcidGasEq(**H2CO3_aireq)

sys0 = System(Na, carb)
sys0.pHsolve()
print(sys0.pH)

sys1 = System(Na, carbCO2)
sys1.pHsolve()
print(sys1.pH)


print()
print("'Air-insensitive' carbonate buffer with pH = pKa? approximately")

carb = AcidAq(**H2CO3,
              conc = 0.30)

Na = IonAq(charge = +1,
           conc = 0.45)

carbCO2 = AcidGasEq(**H2CO3_aireq)

sys0 = System(Na, carb)
sys0.pHsolve()
print(sys0.pH)

sys1 = System(Na, carbCO2)
sys1.pHsolve()
print(sys1.pH)

# generate full report (system_print)
print()
print(sys1)
