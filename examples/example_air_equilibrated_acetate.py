# Example calculations:
# sodium acetate system with excess NaOH
#  equilibrium with ambient CO2

from pHfork import IonAq, AcidAq, System
from pHfork import AcidGasEq


# Chemical information on buffer components
#
# For the weak acid, charge is the charge of the fully protonated acid.
# The charge of the deprotonated species are taken into account accordingly.

acetate = {'charge' : 0,
           'pKa'    : 4.76,
           'name'   : 'HAc/Ac-'}

# chloride = {'charge' : -1}

sodium = {'charge' : +1,
          'name'   : 'Na+' }

H2CO3 = {'charge': 0,
         'pKa'   : [6.35, 10.33]}

H2CO3_aireq = H2CO3 | {'Hs'    : 0.03429,
                       'Pgas'  : 415e-6,
                       'name'  : 'carbonic'}


# Description of buffer preparation
#
sHAc  = 2e-4 # [M] stock solution acetic acid
sNaAc = 1e-4 # [M] stock solution sodium acetate
sNaOH = 2.2e-3 # [M] stock solution sodium hydroxide
 
VHAc  = 1.0 # [l] volume of stock used 
VNaAc = 0.0 # [l] volume of stock used 
VNaOH = 1.0 # [l]

Vtot = VHAc + VNaAc + VNaOH

cHAc   = (VHAc*sHAc) / Vtot   # [M] added concentration
cNaAc  = (VNaAc*sNaAc) / Vtot # [M] added concentration
cNaOH  = (VNaOH*sNaOH) / Vtot

# Final concentrations of relevant components, without H+ and OH-
# Each component concentration is the sum of all that component's species
#
cAc = cHAc + cNaAc   # total concentration of all acetate species
cNa = cNaAc + cNaOH  # total concentration of sodium ions


# Definition of systems
#
ac = AcidAq(**acetate, conc = cAc)
na = IonAq(**sodium, conc = cNa)
buf = System(ac, na)

carbCO2 = AcidGasEq(**H2CO3_aireq)
bufCO2 = System(ac, na, carbCO2)

# Calculate pH
buf.pHsolve()
bufCO2.pHsolve()

# Print result
print('sodium acetate system')
print()
print(buf)
print()
print(f'pH = {buf.pH:.3f}')
print(f' I = {buf.I:.3e} M')

print()
print(60*'*')
print()

# Print result
print('sodium acetate system in eq. w/ ambient CO2')
print()
print(bufCO2)
print()
print(f'pH = {bufCO2.pH:.3f}')
print(f' I = {bufCO2.I:.3e} M')

