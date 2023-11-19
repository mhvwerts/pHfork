import numpy as np
import scipy.optimize as spo


class IonAq:
    """An intert ion class.

    This class defines things like K+ and Cl-, which contribute to the
    overall charge balance, but do not have any inherent reactivity with
    water. Adding ions without adding the corresponding counter-ions will
    incite pHcalc to generate OH- or H+ as counter-ions. This is the behavior
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
    '''An acidic species class.

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

    '''
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
        '''Return the fraction of each species at a given pH.

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
        '''
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
        Ka_prod = np.cumproduct(self._Ka_temp)
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
    '''An object used to define an a system of acid and neutral species.

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
    '''
    
    
    #####################################################
    # print strings for system_print as class variables #
    #####################################################

    no_ph = '''### THE CONCENTRATIONS OF THIS SYSTEM ARE NOT AT EQUILIBRIUM ###
    To determine the equilibrium species distribution use System.pHsolve\n\n'''
    
    has_ph = '''### THESE ARE THE EQUILIBRIUM SYSTEM CONCENTRATIONS ###
    
    SYSTEM pH: {0.pH:.3f}\n\n'''
    
    sep_line1 = '='*65 + '\n'
    
    sep_line2 = '-'*65 + '\n'
    
    header_line = f"{'Species':15}{'Charge':10}{'Ka':15}{'pKa':10}{'Conc':15}\n"
    
    acid_line = '{0:15}{1:<+10d}{2:<15.3e}{3:<10.2f}{4:<15.4e}\n'
    
    ion_line = '{0:15}{1:<+35d}{2:<15.4e}\n'
    
    #####################################################

    
    def __init__(self, *species, Kw=1.01e-14):
        self.species = species
        self.Kw = Kw

    def __str__(self, ):
        '''Return a string representing the composition of the system.

        The unformatted text strings used in this method are in a separate
        module `print_string.py`.
        '''
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
                    concs = cpd.alpha(self.pH)
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
        '''The representation of this object will be the same as printing for
        interactive terminals.
        '''
        repr_str = self.__str__()

        return repr_str

    def _diff_pos_neg(self, pH):
        '''Calculate the charge balance difference.

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
        '''
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

        # Return the absolute value so it never goes below zero.
        return np.abs(x)
        

    def pHsolve(self, guess=7.0, guess_est=False, est_num=1500, 
                method='Nelder-Mead', tol=1e-5):
        '''Solve the pH of the system.

        The pH solving is done using a simple minimization algorithm which
        minimizes the difference in the total positive and negative ion
        concentrations in the system. The minimization algorithm can be
        adjusted using the `method` keyword argument. The available methods
        can be found in the documentation for the scipy.optimize.minimize
        function.
        
        A good initial guess may help the minimization. It can be set manually
        using the `guess` keyword, which defaults to 7.0. There is an
        automated method that can be run as well if you set the `guess_est`
        argument. This will override whatever you pass is for `guess`. The
        `est_num` keyword sets the number of data points that you'd like to
        use for finding the guess estimate. Too few points might start you
        pretty far from the actual minimum; too many points is probably
        overkill and won't help much. This may or may not speed things up.

        Parameters
        ----------

        guess : float (default 7.0)
            This is used as the initial guess of the pH for the system. 

        guess_est : bool (default False)
            Run a simple algorithm to determine a best guess for the initial
            pH of the solution. This may or may not slow down the calculation
            of the pH.

        est_num : int (default 1500)
            The number of data points to use in the pH guess estimation.
            Ignored unless `guess_est=True`.

        method : str (default 'Nelder-Mead')
            The minimization method used to find the pH. The possible values
            for this variable are defined in the documentation for the
            scipy.optimize.minimize function.

        tol : float (default 1e-5)
            The tolerance used to determine convergence of the minimization
            function.
        '''
        if guess_est == True:
            phs = np.linspace(0, 14, est_num)
            guesses = self._diff_pos_neg(phs)
            guess_idx = guesses.argmin()
            guess = phs[guess_idx]
            
        self.pHsolution = spo.minimize(self._diff_pos_neg, guess, 
                method=method, tol=tol)
        
        if self.pHsolution.success == False:
            print('Warning: Unsuccessful pH optimization!')
            print(self.pHsolution.message)
            
        if len(self.pHsolution.x) == 1:
            self.pH = self.pHsolution.x[0]



if __name__ == '__main__':
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
        # This visualizes the method used internallly by the pHsolve function 
        # if you use the guess_est=True flag
        # It is based on (graphically) finding the pH at which the solution
        # composition is closest to electroneutrality.
        s = System(a)
        diffs = s._diff_pos_neg(pH)
        plt.figure(2)
        plt.clf()
        plt.title('1mM H3PO4(aq) - initial guess pH')
        plt.semilogy(pH, diffs)
        plt.xlabel('pH')
        plt.ylabel('charge imbalance')
        plt.show()
        print('Initial guess pH = ',pH[np.argmin(diffs)])
        s.pHsolve()
        print('pHsolve refined pH = ', s.pH)

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
            s.pHsolve(guess_est=True)
            phs.append(s.pH)
        plt.figure(3)
        plt.clf()
        plt.title('Phosphoric Acid Titration Curve')
        plt.plot(Na_concs, phs)
        plt.xlabel('total added NaOH concentration (M)')
        plt.ylabel('pH')
        plt.show()
