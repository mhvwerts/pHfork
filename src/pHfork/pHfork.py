import numpy as np
import scipy.optimize as spo


class IonAq:
    """An intert ion class.

    This class defines things like K+ and Cl-, which contribute to the
    overall charge balance, but do not have any inherent reactivity with
    water. Adding ions without adding the corresponding counter-ions will
    incite pHfork to generate OH- or H+ as counter-ions. This is the behavior
    corresponding to 'strong' bases (e.g. KOH) and acids (e.g. HCl)
    
    Parameters
    ----------
    charge : int
        The formal charge of the ion.

    conc : float
        The concentration of this species in solution.

    Attributes
    ----------
    charge : int
        The formal charge of the ion.

    conc : float
        The concentration of this species in solution.

    """
    def __init__(self, charge=None, conc=None, name=None):
        if charge == None:
            raise ValueError(
                "The charge for this ion must be defined.")

        self.charge = charge 
        self.conc = conc
        self.name = name

    def alpha(self, pH):
        '''Return the fraction of each species at a given pH.

        Parameters
        ----------
        pH : int, float, or Numpy Array
            These are the pH value(s) over which the fraction should be
            returned.

        Returns
        -------
        Numpy NDArray
            Because this is a non-reactive ion class, this function will
            always return a Numpy array containing just 1.0's for all pH
            values.

        '''
        if isinstance(pH, (int, float)):
            length = 1
        else:
            length = len(pH)
        ones = np.ones(length).reshape(-1,1)
        return ones


        
class AcidAq:
    """An acidic species class.

    This class describes species that exist in different protonation states
    in equilibrium. It is used to calculate a number of parameters related to 
    a 'weak' acid in an aqueous solution. 
    
    Parameters
    ----------
    Ka : None (default), float, list, Numpy Array
        This defines the Ka values for all acidic protons in this species. It
        can be a single Ka value (float), a list of floats, or a Numpy array
        of floats. Either this value or pKa needs to be defined. The other
        will then be calculated from the given values.

    pKa : None (default), float, list, Numpy Array
        The pKa value(s) for all the acidic protons in this species.  This
        follows the same rules as Ka (See Ka description for more details),
        and either this value or Ka must be defined.

    charge : None (default), int
        This is the charge of the fully protonated form of this acid. This
        must be defined.

    conc : None (default), float
        The formal concentration of this acid in solution. This value must be
        defined.

    Note
    ----
    There is no corresponding Base object. To define a base, you must use a
    combination of an AcidAq and IonAq object. See the documentation for
    examples.

    """
    def __init__(self, Ka=None, pKa=None, charge=None, conc=None, name=None):
        # Do a couple quick checks to make sure that everything has been
        # defined.
        if Ka == None and pKa == None:
            raise ValueError(
                "You must define either Ka or pKa values.")
        elif charge == None:
            raise ValueError(
                "The maximum charge for this acid must be defined.")

        # Make sure both Ka and pKa are calculated. For lists of values, be
        # sure to sort them to ensure that the most acidic species is defined
        # first.
        elif Ka == None:
            if isinstance(pKa, (int, float)):
                self.pKa = np.array( [pKa,], dtype=float)
            else:
                self.pKa = np.array(pKa, dtype=float)
                self.pKa.sort()
            self.Ka = 10**(-self.pKa) 
        elif pKa == None:
            if isinstance(Ka, (int, float)):
                self.Ka = np.array( [Ka,], dtype=float)
            else:
                self.Ka = np.array(Ka, dtype=float)
                # Ka values must be in reverse sort order
                self.Ka.sort()
                self.Ka = self.Ka[::-1]
            self.pKa = -np.log10(self.Ka)
        # This temporary Ka array will be used to calculate alpha values. It
        # starts with an underscore so that it won't be confusing for others.
        self._Ka_temp = np.append(1., self.Ka)
        
        # Make a list of charges for each species defined by the Ka values.
        self.charge = np.arange(charge, charge - len(self.Ka) - 1, -1)
        # Make sure the concentrations are accessible to the object instance.
        self.conc = conc 
        self.name = name

    def alpha(self, pH):
        """Return the fraction of each species at a given pH.

        Parameters
        ----------
        pH : int, float, or Numpy Array
            These are the pH value(s) over which the fraction should be
            returned.

        Returns
        -------
        Numpy NDArray
            These are the fractional concentrations at any given pH. They are
            sorted from most acidic species to least acidic species. If a
            NDArray of pH values is provided, then a 2D array will be
            returned. In this case, each row represents the speciation for
            each given pH.
        """
        # If the given pH is not a list/array, be sure to convert it to one
        # for future calcs.
        if isinstance(pH, (int, float)):
            pH = [pH,]
        pH = np.array(pH, dtype=float)

        # Calculate the concentration of H3O+. If multiple pH values are
        # given, then it is best to construct a two dimensional array of
        # concentrations.
        h3o = 10.**(-pH)
        if len(h3o) > 1:
            h3o = np.repeat( h3o.reshape(-1, 1), len(self._Ka_temp), axis=1)

        # These are the powers that the H3O+ concentrations will be raised.
        power = np.arange(len(self._Ka_temp))
        # Calculate the H3O+ concentrations raised to the powers calculated
        # above (in reverse order).
        h3o_pow = h3o**( power[::-1] )
        # Calculate a cumulative product of the Ka values. The first value
        # must be 1.0, which is why _Ka_temp is used instead of Ka.
        Ka_prod = np.cumprod(self._Ka_temp)
        # Multiply the H3O**power values times the cumulative Ka product.
        h3o_Ka = h3o_pow*Ka_prod

        # Return the alpha values. The return signature will differ is the
        # shape of the H3O array was 2-dimensional. 
        if len(h3o.shape) > 1:
            den = h3o_Ka.sum(axis=1)
            return h3o_Ka/den.reshape(-1,1)
        else:
            den = h3o_Ka.sum()
            return h3o_Ka/den


class AcidGasEq(AcidAq):
    """
    An acidic species with additional solution-gas equilibrium.
    
    This class, a subclass of AcidAq, treats those acidic species that 
    exist in equilibrium with a gas-phase resevoir, which keeps the amount
    of dissolved gas (the neutral form) constant. 
    
    The main application of this species is to model the effect of atmospheric
    CO2 on the acid-base chemistry of the solution. Dissolved CO2(aq) exists 
    partially as carbonic     acid (H2CO3)* which then further dissociates into
    bicarbonate (HCO3-) and carbonate (CO3{2-}). This is favoured at high pH.
    
    *) CO2(aq) and H2CO3(aq) are treated as a single species and simply
       referred to as CO2(aq), whose concentration is kept constant, and which
       is in equilibrium with H+(aq) + HCO3-(aq)

    Parameters
    ----------
    Ka, pKa, charge, conc : 
        See AcidAq.
    Hs : float
        Henry solubility of the species (in [M atm-1] or [M Pa-1], 
        depending on the units Pgas). For example, the value for CO2 in pure
        water at 298K is 0.03429 M atm-1, taken from ref. [1], Table II, 
        converted from molality to molarity, i.e. 0.997 * 10**-1.463
    Pgas : float
        The partial pressure of the species in the gas phase in [atm] or [Pa],
        in consistence with Hs. The partial pressure is considered constant,
        i.e. the gas phase is a reservoir of the species. For example, the
        average Pgas for CO2 in the ambient atmosphere is currently 417e-6 atm
        (it was 372e-6 atm in 1972). This value is rising over the years:
        https://www.climate.gov/news-features/understanding-climate/climate-change-atmospheric-carbon-dioxide    
    Additional keyword arguments : 
        See AcidAq.
    

    [1] H. S. Harned, R. Davis.
        "The Ionization Constant of Carbonic Acid in Water and the Solubility
        of Carbon Dioxide in Water and Aqueous Salt Solutions from 0 to 50°."
        J. Am. Chem. Soc. 1943, 65 , 2030

    """
    def __init__(self, *args, Hs, Pgas, **kwargs):
        super().__init__(*args, **kwargs)

        self.gas_conc = Hs*Pgas
        self.conc = self.gas_conc
        self.gas_idx = np.argmin( np.abs( self.charge ) )


    def alpha(self, pH):
        alphas = super().alpha(pH)

        self.conc = self.gas_conc/alphas[self.gas_idx]

        return alphas
        

class System:
    """An object used to define an a system of acid and neutral species.

    This object accepts an arbitrary number of acid and neutral species
    objects and uses these to calculate the pH of the system. Be sure to
    include all of the species that completely define the contents of a
    particular solution.

    Parameters
    ----------
    *species 
        These are any number of AcidAq and IonAq objects that you'd like to
        use to define your system.

    Kw : float (default 1.01e-14)
        The autoionization constant for water. This may vary, e.g. as a
        function of temperature. The default value is for water at
        298 K and 1 atm.

    Attibutes
    ---------
    species : list
        This is a list containing all of the species that you input.

    Kw : float
        The autoionization of water set using the Kw keyword argument.

    pHsolution 
        This is the full minimization output, which is defined by the function
        scipy.optimize.minimize. This is only available after running the
        pHsolve method.

    pH : float
        The pH of this particular system. This is only calculated after
        running the pHsolve method.
    """
    
    
    #####################################################
    # print strings for system_print as class variables #
    #####################################################

    no_ph = '''### THE CONCENTRATIONS OF THIS SYSTEM ARE NOT AT EQUILIBRIUM ###
To determine the equilibrium species distribution use System.pHsolve()\n\n'''
    
    has_ph = '''### THESE ARE THE EQUILIBRIUM SYSTEM CONCENTRATIONS ###
system pH:      {0.pH:.3f}
ionic strength: {0.I:.3e} M\n\n'''
    
    sep_line1 = '='*65 + '\n'
    
    sep_line2 = '-'*65 + '\n'
    
    header_line = f"{'Species':15}{'Charge':10}{'Ka':15}{'pKa':10}{'Conc. / M':15}\n"
    
    acid_line = '{0:15}{1:<+10d}{2:<15.3e}{3:<10.2f}{4:<15.4e}\n'
    
    ion_line = '{0:15}{1:<+35d}{2:<15.4e}\n'
    
    #####################################################
   
    def __init__(self, *species, Kw=1.01e-14):
        self.species = species
        self.Kw = Kw

    def __str__(self, ):
        """Return a string representing the composition of the system.
        """
        # The ultimate string to be returned from this method
        prt_str = ''

        # Print a header based on the status of the system
        if not hasattr(self, 'pH'):
            prt_str += self.no_ph
        else:
            prt_str += self.has_ph.format(self)
        prt_str += self.header_line
        prt_str += self.sep_line1

        # Print the species information
        # Need counters for acid/ion names and total charge concentration
        acid_num = 1
        ion_num = 1
        charge_conc = 0.
        for cpd in self.species:
            # For acids, print information for all possible species
            # The final species will not have a Ka/pKa value
            if isinstance(cpd, AcidAq):
                name = cpd.name if cpd.name else f'AcidAq{acid_num}'
                acid_num += 1

                kas = list(cpd.Ka) + [np.nan,]
                pkas = list(cpd.pKa) + [np.nan,]

                if hasattr(self, 'pH'):
                    concs = cpd.alpha(self.pH) * cpd.conc
                else:
                    concs = [cpd.conc,] + [0,]*len(cpd.Ka)
                    concs = np.array(concs)

                props = zip(cpd.charge, kas, pkas, concs)
                for charge, ka, pka, conc in props:
                    if ka is not np.nan:
                        prt_str += self.acid_line.format(name, charge, ka, 
                                            pka, conc)
                    else:
                        prt_str += self.ion_line.format(name, charge, conc)

                charge_conc += (cpd.charge*concs).sum()

            # Ions do not have Ka/pKa values either, so they are simpler
            else:
                name = cpd.name if cpd.name else f'Ion{ion_num}'
                ion_num += 1

                prt_str += self.ion_line.format(name, cpd.charge, cpd.conc) 
                charge_conc += cpd.charge*cpd.conc
            # Separate species with a different line
            prt_str += self.sep_line2

        # Hydronium/hydroxide concentrations
        if hasattr(self, 'pH'):
            h3o = 10**-self.pH
            oh = self.Kw/h3o
        else:
            # Solve the following equation
            # [H3O] + charge_conc - [OH] = 0
            # [H3O]^2 + [H3O]*charge_conc - Kw = 0
            possible_h3o = np.roots([1., charge_conc, -self.Kw])
            h3o = possible_h3o.max() # Must be positive num
            oh = self.Kw/h3o
        prt_str += self.ion_line.format('H3O+', 1, h3o)
        prt_str += self.ion_line.format('OH-', -1, oh)

        return prt_str

    def __repr__(self, ):
        """The representation of this object will be the same as printing for
        interactive terminals.
        """
        repr_str = self.__str__()

        return repr_str

    def _diff_pos_neg(self, pH):
        """Calculate the charge balance difference.

        Parameters
        ----------
        pH : int, float, or Numpy Array
            The pH value(s) used to calculate the different distributions of
            positive and negative species.

        Returns
        -------
        float or Numpy Array
            The absolute value of the difference in concentration between the
            positive and negatively charged species in the system. A float is
            returned if an int or float is input as the pH: a Numpy array is
            returned if an array of pH values is used as the input.
        """
        twoD = True
        if isinstance(pH, (int, float)) or pH.shape[0] == 1:
            twoD = False
        else:
            pH = np.array(pH, dtype=float)
        # Calculate the h3o and oh concentrations and sum them up.
        h3o = 10.**(-pH)
        oh = (self.Kw)/h3o
        x = (h3o - oh)

        # Go through all the species that were given, and sum up their
        # charge*concentration values into our total sum.
        for s in self.species:
            alphas = s.alpha(pH)
            if twoD == False:
                x += (s.conc*s.charge*alphas).sum()
            else:
                x += (s.conc*s.charge*alphas).sum(axis=1)
        return x
    
    def _get_ionic_strength(self):
        """Calculate the ionic strength from the equilibrium composition
        
        Can only be called after an equilibrium pH has been calculated.
        
        Uses
        
        $$I = \\tfrac{1}{2} \\sum_{i=1}^{n} c_i z_i^{2}$$        
        """
        
        sum_cz2 = 0.0
        # ionic species, except H+ and OH-
        for s in self.species:
            alphas = s.alpha(self.pH)
            sum_cz2 += (s.conc * s.charge**2 * alphas).sum()
        # H+, OH-
        h3o = pow(10, -self.pH)
        oh = self.Kw/h3o
        sum_cz2 += h3o * (+1.0)**2
        sum_cz2 +=  oh * (-1.0)**2
        return 0.5*sum_cz2

    def pHsolve(self, bracket=(-1, 15), method='brentq', xtol=1e-5, 
                options=None, guess=None, guess_est=None, est_num=None,
                tol=None):
        """Solve the pH of the system.

        The pH solving is done using a root-search algorithm which
        finds a root in the difference of the total positive and negative ion
        concentrations in the system. The root-search algorithm can be
        adjusted using the `method` keyword argument. The available methods
        can be found in the documentation for the scipy.optimize.root_scalar
        function.

        Parameters
        ----------

        bracket : a sequence of 2 floats (default (-1, 15))
            An interval of pH values bracketing a root.

        method : str (default 'brentq')
            The type of solver used to find the pH. The possible values
            for this variable are defined in the documentation for the
            scipy.optimize.root_scalar function.

        xtol : float (default 1e-5)
            The tolerance used to determine convergence of the root-searching
            function.

        options : dict, optional (default None)
            A dictionary of solver options.
        """
        if guess is not None:
            print('Warning: option "guess" will be ignored.')
        if guess_est is not None:
            print('Warning: option "guess_est" will be ignored.')
        if est_num is not None:
            print('Warning: option "est_num" will be ignored.')
        if tol is not None:
            print('Warning: option "tol" will be ignored.')

        options = options or dict()
        self.pHsolution = spo.root_scalar(self._diff_pos_neg,
                                          bracket=bracket,
                                          method=method,
                                          xtol=xtol,
                                          **options)

        if not self.pHsolution.converged:
            print('Warning: Unsuccessful pH optimization!')
            print(self.pHsolution.flag)
            self.pH = None
        else:
            self.pH = self.pHsolution.root
            self.I = self._get_ionic_strength()
