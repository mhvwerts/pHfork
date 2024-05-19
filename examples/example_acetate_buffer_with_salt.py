# Example calculations:
# sodium acetate buffer system
# with sodium chloride
#

from pHfork import IonAq, AcidAq, System


# Chemical information on buffer components
#
# For the weak acid, charge is the charge of the fully protonated acid.
# The charge of the deprotonated species are taken into account accordingly.
acetate = {'charge' : 0,
           'pKa'    : 4.76}
chloride = {'charge' : -1}
sodium = {'charge' : +1}


# Description of buffer preparation
sHAc  = 0.1 # [M] stock solution acetic acid
sNaAc = 0.1 # [M] stock solution sodium acetate
sNaCl = 1.0 # [M] stock solution sodium chloride

VHAc  = 36e-3 # [l] volume of stock used 
VNaAc = 64e-3 # [l] volume of stock used 
VNaCl = 0.0   # [l] volume of stock used 

Vtot = VHAc + VNaAc + VNaCl

cHAc   = (VHAc*sHAc) / Vtot   # [M] added concentration
cNaAc  = (VNaAc*sNaAc) / Vtot # [M] added concentration
cNaCl  = (VNaCl*sNaCl) / Vtot # [M] added concentration

cAc = cHAc + cNaAc   # total concentration of all acetate species
cNa = cNaAc + cNaCl  # total concentration of sodium ions
cCl = cNaCl          # total concentration of chloride ions


# Definition of system
ac = AcidAq(**acetate, conc = cAc)
na = IonAq(**sodium, conc = cNa)
cl = IonAq(**chloride, conc = cCl)
buf = System(ac, na, cl)


# Calculate pH
buf.pHsolve()


# Print result
print('sodium acetate buffer with perhaps some added sodium chloride')
print()
print(buf)
print()
print(f'pH = {buf.pH:.3f}')
print(f' I = {buf.I:.3e} M')

