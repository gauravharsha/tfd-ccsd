import sys
sys.path.append('../fort_src/')

import h5py, time, numpy as np

from scipy.integrate import ode
from scipy.misc import comb
from multiprocessing import Pool

from iofuncs import *
from inttran import *

from ThermalCCSD import *
from ThermalCISD import *
from ExpVals import *

import pdb

#
# GLOBAL VARIABLES
#

t2_symm_tol = 5e-8

alpha_step_0g = +5e-2
len_t1 = 0
len_t2 = 0

#
# Compress / DeCompress functions
#

def DecompressT2(T2_Compressed,NSO):
    """
    The Idea is to use the antisymmetric property of the T2 tensor:
        T2[p,q,r,s] = -T2[q,p,r,s]
                    = -T2[p,q,s,r]
        i.e. the first 2 and the last 2 indices are anti-symmetric

    This allows a compressed storage of T2's but in passing the T2 amps
    to the FORTRAN SubRoutine, we need to reconstruct the full T2
    which is done by this function
    """

    if np.size(T2_Compressed) != int(comb(NSO,2)**2):
        raise ValueError('Invalid Size of the compressed T2 array',T2_Compressed)

    t2_out = np.zeros((NSO,NSO,NSO,NSO))
    m = 0

    for i in range(NSO):

        for j in range(i+1,NSO):

            for k in range(NSO):

                for l in range(k+1,NSO):

                    val = T2_Compressed[m]
                    t2_out[i,j,k,l] = val
                    t2_out[j,i,k,l] = -val
                    t2_out[i,j,l,k] = -val
                    t2_out[j,i,l,k] = val

                    m += 1
    return t2_out


def CompressT2(T2):
    """
    Here, we also CHECK that the input tensor has the right symmetry
    And proceed directly to COMPRESS in the following way:

        (I,J,K,L) --> I<J and K<L

        Arranged in lexicological order, for instance

            for first pair of indices:  (1,2) comes before (1,3) before (2,3) and so on.

            the actual order is:        (0,1,0,1) -> (0,1,0,2) -> ... -> (0,1,0,NSO-1) ->
                                        (0,1,1,2) -> .... and so on...
    """
    NSO = np.size(T2,axis=0)
    m = 0
    T2_compressed = np.zeros( int( comb(NSO,2)**2 ) )

    for i in range(NSO):
        for j in range(i+1,NSO):
            for k in range(NSO):
                for l in range(k+1,NSO):

                    val = T2[i,j,k,l]

                    chk1 = np.isclose(T2[j,i,k,l],-val,rtol=t2_symm_tol)
                    chk2 = np.isclose(T2[i,j,l,k],-val,rtol=t2_symm_tol)
                    chk3 = np.isclose(T2[j,i,l,k],val,rtol=t2_symm_tol)

                    if chk1 and chk2 and chk3:
                        T2_compressed[m] = val
                        m += 1
                    else:
                        print('T2[{},{},{},{}] = {}'.format(i,j,k,l,T2[i,j,k,l]))
                        print('T2[{},{},{},{}] = {}'.format(j,i,k,l,T2[j,i,k,l]))
                        print('T2[{},{},{},{}] = {}'.format(i,j,l,k,T2[i,j,l,k]))
                        print('T2[{},{},{},{}] = {}'.format(j,i,l,k,T2[j,i,l,k]))
                        raise ValueError('Incorrect Symmetry for the input Tensor',T2)
    
    return T2_compressed

#
# CI - to - CC amplitude transformation
# 

def ci_to_cc_transform(CI_amps, CC_amps, Nso):
    """
    Function that returns the Z-amplitudes from Z-CI and T amplitudes 
    Inputs:     
                CC_amps     ::  T - amplitudes stacked together as follows
                                T0      ::  T0 amplitude
                                T1      ::  T1 amplitude vectors :: (Nso*Nso) elements 
                                T2      ::  T2 amplitude vectors :: (Nso-choose-2)^2 elements
                CI_amps     ::  CI-amplitudes for the left hand state
                Nso         ::  Number of Spin orbitals
    Returns:
                Z_out_1     ::  Nso x Nso matrix for Z1
                Z_out_2     ::  Nso x Nso x Nso x Nso matrix for Z2
    """

    T0 = CC_amps[0]
    T1 = np.reshape(CC_amps[1:Nso*Nso+1], (Nso,Nso))
    T2 = DecompressT2(CC_amps[Nso*Nso+1:], Nso)

    Z0 = CI_amps[0]
    Z1 = np.reshape(CI_amps[1:Nso*Nso+1], (Nso,Nso))
    Z2 = DecompressT2(CI_amps[Nso*Nso+1:], Nso)

    Z_out_0 = 1 + np.einsum('ab,ab',T1,Z1) + 0.25*np.einsum('abcd,abcd',T2,Z2) + 0.5*np.einsum('ac,bd,abcd',T1,T1,Z2)
    Z_out_1 = np.einsum('cd,acbd->ab',T1,Z2)
    Z_out_1 += np.reshape(Z1,(Nso,Nso))
    Z_out_2 = Z2/4

    return Z_out_1/Z_out_0, Z_out_2/Z_out_0


#
# Evolution Driver Functions
# 

def cc_beta_evolve(Beta,Amps,Alpha,Fug,OneH,ERI):
    """
    Driver function for the Beta evolution of Thermal CCSD
    Inputs:
            Beta    ::  Current Beta Value
            Amps    ::  Collection of CCSD-amplitudes
            Alpha   ::  Current Alpha Value
            Fug     ::  (fugacity?) indicating initial exp(Alpha*Beta)
            OneH    ::  1D array of one-body energy eigenvalues
            ERI     ::  2-electron integral
    Returns:
            Concatenated array of unique R0,R1(pq) and R2(pqrs)
    """

    # Number of Spin Orbitals
    NSO = np.size(OneH,axis=0)

    # Thermal HFB parameters
    U = 1/np.sqrt( 1 + np.exp(-Beta*OneH + Alpha)*Fug )
    V = np.sqrt( 1 - U**2 )

    # Extract the amplitudes
    T1 = np.reshape( Amps[1:1+len_t1], (NSO,NSO) )
    T2 = DecompressT2( Amps[1+len_t1:], NSO )

    # Get the residuals
    r0, r1, r2 = betacc(OneH,ERI,T1,T2,U,V)

    # Antisymmetrize r2
    s2 = (
            r2 \
            - np.einsum('baij->abij',r2) \
            - np.einsum('abji->abij',r2) \
            + np.einsum('baji->abij',r2) \
        )/4.0

    # Compress and concatenate
    dt0_dbeta = r0
    dt1_dbeta = np.reshape( r1, int(NSO**2) )
    dt2_dbeta = CompressT2( s2 )

    # return concatenated
    out = np.concatenate(( [dt0_dbeta], dt1_dbeta, dt2_dbeta ))

    return out


def ci_beta_evolve(Beta,Amps,Alpha,Fug,OneH,ERI):
    """
    Driver function for the Beta evolution of Thermal CISD
    Inputs:
            Beta    ::  Current Beta Value
            Amps    ::  Collection of CCSD-amplitudes
            Alpha   ::  Current Alpha
            Fug     ::  (fugacity?) indicating initial exp(Alpha*Beta)
            OneH    ::  1D array of one-body energy eigenvalues
            ERI     ::  2-electron integral
    Returns:
            Concatenated array of unique R0,R1(pq) and R2(pqrs)
    """

    # Number of Spin Orbitals
    NSO = np.size(OneH,axis=0)

    # Thermal HFB parameters
    U = 1/np.sqrt( 1 + np.exp(-Beta*OneH + Alpha)*Fug )
    V = np.sqrt( 1 - U**2 )

    # Extract the amplitudes
    T1 = np.reshape( Amps[1:1+len_t1], (NSO,NSO) )
    T2 = DecompressT2( Amps[1+len_t1:], NSO )

    # Get the residuals
    r0, r1, r2 = betaci(OneH,ERI,T1,T2,U,V)

    # Antisymmetrize r2
    s2 = (
            r2 \
            - np.einsum('baij->abij',r2) \
            - np.einsum('abji->abij',r2) \
            + np.einsum('baji->abij',r2) \
        )/4.0

    # Compress and concatenate
    dt0_dbeta = r0
    dt1_dbeta = np.reshape( r1, int(NSO**2) )
    dt2_dbeta = CompressT2( s2 )

    # return concatenated
    return np.concatenate(( [dt0_dbeta], dt1_dbeta, dt2_dbeta ))


def cc_alpha_evolve(Alpha,Amps,Beta,Fug,OneH):
    """
    Driver function for the Alpha evolution of Thermal CCSD
    Inputs:
            Alpha   ::  Current Alpha
            Amps    ::  Collection of CCSD-amplitudes
            Beta    ::  Current Temperature
            Fug     ::  (fugacity?) indicating initial exp(Alpha*Beta)
            OneH    ::  1D array of one-body energy eigenvalues
    Returns:
            Concatenated array of unique R0,R1(pq) and R2(pqrs)
    """

    # Number of Spin Orbitals
    NSO = len(OneH)

    # Thermal HFB parameters
    U = 1/np.sqrt( 1 + np.exp(-Beta*OneH + Alpha)*Fug )
    V = np.sqrt( 1 - U**2 )

    # Extract the amplitudes
    T1 = np.reshape( Amps[1:1+len_t1], (NSO,NSO) )
    T2 = DecompressT2( Amps[1+len_t1:], NSO )

    # Get the residuals
    r0, r1, r2 = numbercc(T1,T2,U,V)

    # Antisymmetrize r2
    s2 = (
            r2 \
            - np.einsum('baij->abij',r2) \
            - np.einsum('abji->abij',r2) \
            + np.einsum('baji->abij',r2) \
        )/4.0

    # Compress and concatenate
    dt0_dbeta = r0
    dt1_dbeta = np.reshape( r1, int(NSO**2) )
    dt2_dbeta = CompressT2( s2 )

    # return concatenated
    return np.concatenate(( [dt0_dbeta], dt1_dbeta, dt2_dbeta ))


def ci_alpha_evolve(Alpha,Amps,Beta,Fug,OneH):
    """
    Driver function for the Alpha evolution of Thermal CISD
    Inputs:
            Alpha   ::  Current Alpha
            Amps    ::  Collection of CCSD-amplitudes
            Beta    ::  Current Temperature
            Fug     ::  (fugacity?) indicating initial exp(Alpha*Beta)
            OneH    ::  1D array of one-body energy eigenvalues
    Returns:
            Concatenated array of unique R0,R1(pq) and R2(pqrs)
    """

    # Number of Spin Orbitals
    NSO = len(OneH)

    # Thermal HFB parameters
    U = 1/np.sqrt( 1 + np.exp(-Beta*OneH + Alpha)*Fug )
    V = np.sqrt( 1 - U**2 )

    # Extract the amplitudes
    T1 = np.reshape( Amps[1:1+len_t1], (NSO,NSO) )
    T2 = DecompressT2( Amps[1+len_t1:], NSO )

    # Get the residuals
    r0, r1, r2 = numberci(T1,T2,U,V)

    # Antisymmetrize r2
    s2 = (
            r2 \
            - np.einsum('baij->abij',r2) \
            - np.einsum('abji->abij',r2) \
            + np.einsum('baji->abij',r2) \
        )/4.0

    # Compress and concatenate
    dt0_dbeta = r0
    dt1_dbeta = np.reshape( r1, int(NSO**2) )
    dt2_dbeta = CompressT2( s2 )

    # return concatenated
    return np.concatenate(( [dt0_dbeta], dt1_dbeta, dt2_dbeta ))


#
# Energy and Number Eval Functions
# 

def eval_number(CC_amps, CI_amps, X, Y):
    """
        Function to evaluate the number expectation value
    """

    Nso = len(X)

    T0 = CC_amps[0]
    T1 = np.reshape(CC_amps[1:Nso*Nso+1], (Nso,Nso))
    T2 = DecompressT2(CC_amps[Nso*Nso+1:], Nso)

    # CI amplitudes to CC amplitudes
    Z1, Z2 = ci_to_cc_transform(CI_amps,CC_amps,Nso)

    # Number Exp value
    num = evalnumber(T1,T2,Z1,Z2,X,Y)

    return num

def eval_energy(OneH, ERI, CC_amps, CI_amps, X, Y):
    """
        Function to evaluate the number expectation value
    """

    Nso = len(X)

    T0 = CC_amps[0]
    T1 = np.reshape(CC_amps[1:Nso*Nso+1], (Nso,Nso))
    T2 = DecompressT2(CC_amps[Nso*Nso+1:], Nso)

    # CI amplitudes to CC amplitudes
    Z1, Z2 = ci_to_cc_transform(CI_amps,CC_amps,Nso)

    # Number Exp value
    en = evalenergy(OneH, ERI, T1, T2, Z1, Z2, X, Y)

    return en

#
# ODE integration functions
# 


def DoIntegration(integrator, x_final):
    """
    Intermediate function to perform integration from x_initial to x_final.
    The x_initial and all other information is inherently included in the integrator.
    """
    yout = integrator.integrate(x_final)

    return yout


def _do_beta_integration(integrators, amps, betalpha, beta_step, fug, h1, eri):
    """
        Perform the Beta integration - use parallel pool of processes
    """

    # Extract info
    cc_integrator = integrators[0]
    ci_integrator = integrators[1]

    cc_amps = amps[0]
    ci_amps = amps[1]

    beta_in = betalpha[0]
    alpha_in = betalpha[1]

    # Set the initial condition
    cc_integrator.set_initial_value(cc_amps, beta_in).set_f_params(alpha_in, fug, h1, eri)
    ci_integrator.set_initial_value(ci_amps, beta_in).set_f_params(alpha_in, fug, h1, eri)

    cc_amps = DoIntegration(cc_integrator,beta_in + beta_step)
    ci_amps = DoIntegration(ci_integrator,beta_in + beta_step)

    return cc_amps, ci_amps


def _do_alpha_integration(integrators, amps, betalpha, fug, h1, n_elec, ntol):
    """
        Bisection to find Chemical Pot / Alpha value
        and then do integration
    """
    
    # Extract info
    cc_integrator = integrators[0]
    ci_integrator = integrators[1]

    cc_amps = amps[0]
    ci_amps = amps[1]

    beta_in = betalpha[0]
    alpha_in = betalpha[1]

    # Make local copies of the variables
    nso = len(h1)

    # Alpha step
    global alpha_step_0g

    alpha_step = alpha_step_0g

    global len_t1
    global len_t2

    # HFB coefficients
    x = 1/np.sqrt( 1 + np.exp( -beta_in*h1 + alpha_in )*fug )
    y = np.sqrt( 1 - x**2 )

    num = eval_number(cc_amps, ci_amps, x, y)

    # Record the differenc and sign
    ndiff_sgn = np.sign( num - n_elec )
    ndiff_mag = np.abs( num - n_elec )

    print('Number difference: ',ndiff_sgn,num)
    exit()

    # if the number is already converged, then there is no need to do any of the following
    #   and hence we keep an 'if' statement; if the condition evaluates to FALSE, then
    #   the outputs will be yf, mu_f

    if ndiff_mag > ntol:

        # Reset the alpha step
        alpha_step = alpha_step_0g

        # Obtain the bracket to perform the bisection
        mu1 = alpha_in
        mu2 = alpha_in

        cc1 = cc_amps*1
        ci1 = ci_amps*1

        cc2 = cc_amps*1
        ci2 = ci_amps*1

        count = 0
        sp_count = 0

        while np.sign( num - n_elec ) == ndiff_sgn:

            count += 1
            if count > 300:
                print('Could not find the bracket after ',count,' steps')
                count = 0
                exit()
                break

            # Set up for next iteration
            mu1 = mu2
            mu2 += alpha_step

            # update solver initial conditions
            cc1 = cc2*1.0
            ci1 = ci2*1.0

            # Set initial condition for the ODE solvers
            cc_integrator.set_initial_value(cc1, mu1)
            cc_integrator.set_f_params(beta_in, fug, h1)

            ci_integrator.set_initial_value(ci1, mu1)
            ci_integrator.set_f_params(beta_in, fug, h1)

            # Do Evolution
            cc2 = DoIntegration(cc_integrator,mu2)
            ci2 = DoIntegration(ci_integrator,mu2)

            # Update HFB coefficients and check the number expectation
            x = 1/np.sqrt( 1 + np.exp( -beta_in*h1 + mu2 )*fug)
            y = np.sqrt( 1 - x**2 )

            num = eval_number(cc2, ci2, x, y)

            # Exit if converged
            if np.abs(num - n_elec) <= ntol:
                break

            # Check if we are evolving in the right direction or do we need to switch sign of step
            val = np.abs(num - n_elec) - ndiff_mag
            if (val>0):
                if val<1e-1:
                    sp_count += 1
                else:
                    alpha_step_0g *= -1
                    alpha_step *= -1
                
                if sp_count >= 10:
                    alpha_step_0g *= -1
                    alpha_step *= -1
                    sp_count = 0

            ndiff_mag = np.abs(num - n_elec)

            # # Do some printing
            print('\t\t\tStart value of Mu = {}'.format(mu1-alpha_step))
            print('\t\t\tEnd value of Mu = {}'.format(mu1))
            print('\t\t\tNumber of particles after evolution = {}'.format(num))
            print('\t\t\t----------------------------------------------\n')
            


        # Printing disabled
        print('Bracket found betwee mu = {} and mu = {}'.format(mu1-alpha_step,mu1))

        # Now do the Bisection
        mu_bisect = [mu1,mu2]

        mu_mid = mu_bisect[1]

        ndiff_sgn_right = np.sign(num - n_elec)
        ndiff_sgn_left = -ndiff_sgn_right

        count = 0
        while np.abs(num - n_elec)>ntol:

            # Set the initial condition for solver
            cc_integrator.set_initial_value(cc1, mu_bisect[0])
            cc_integrator.set_f_params(beta_in, fug, h1)
            ci_integrator.set_initial_value(ci1, mu_bisect[0])
            ci_integrator.set_f_params(beta_in, fug, h1)

            # get mu_mid
            mu_mid = np.mean(mu_bisect)

            # Do Evolution
            cc2 = DoIntegration(cc_integrator,mu_mid)
            ci2 = DoIntegration(ci_integrator,mu_mid)

            # Update HFB Coefficients and compute N expectation
            x = 1/np.sqrt( 1 + np.exp( - beta_in * h1 + mu_mid )*fug )
            y = np.sqrt( 1 - x**2 )

            # Check the current number of particles
            num = eval_number(cc2, ci2, x, y)

            # Adjust bisection bracket
            if np.sign( num - n_elec ) == ndiff_sgn_left:
                mu_bisect[0] = mu_mid
                cc1 = cc2
                ci1 = ci2
            else:
                mu_bisect[1] = mu_mid

        # Printing disabled
        # print('\t\tBisection Converged to mu = ',mu_mid)

        # Final one-shot evolution
        cc_integrator.set_initial_value(cc_amps,alpha_in)
        cc_integrator.set_f_params(beta_in, fug, h1)
        ci_integrator.set_initial_value(ci_amps,alpha_in)
        ci_integrator.set_f_params(beta_in, fug, h1)
        
        # Do Evolution
        cc_amps_out = DoIntegration(cc_integrator,mu_mid)
        ci_amps_out = DoIntegration(ci_integrator,mu_mid)

    else:

        cc_amps_out = cc_amps
        ci_amps_out = ci_amps
        mu_mid = alpha_in

    return cc_amps_out, ci_amps_out, mu_mid


#
# ODE SOLVER FUNCTIONS FOR CC - BETA AND MU EVOLUTIONS
#

class Evolution(IOps):
    """
    class to handle all the imaginary-time evolution equations.
    """

    # Length of T1 and T2 - unique values
    global len_t1
    global len_t2

    ## Class attributes
    # Integration step size
    alpha_step_0 = 1e-1
    alpha_step = 1e-1
    # alpha, beta initial / current values
    alpha_in = 0.0
    beta_in = 0.0

    # Fugacity -- Initial number is fixed
    fug = 0.0


    def __init__(self, inp_file='Input', alpha_step=None):

        # Initialize super to read Input file and data
        super().__init__(inp_file=inp_file)

        # Set Up the Integrals
        self.setUpInts()

        # Step Sizes
        if alpha_step is not None:
            self.alpha_step = alpha_step
            global alpha_step_0g
            alpha_step_0g = self.alpha_step

        # Set the global parameters of length of t1 and t2 arrays
        global len_t1
        global len_t2
        len_t1 = self.nso**2
        len_t2 = int( comb(self.nso,2)**2 )

        # Initialize the amplitudes
        self.cc_amps = np.zeros(1+len_t1+len_t2)
        self.ci_amps = np.zeros(1+len_t1+len_t2)

        # ODE integrators
        self.cc_beta_integrator = ode(cc_beta_evolve).set_integrator('vode', method='bdf',rtol=self.deqtol)
        self.ci_beta_integrator = ode(ci_beta_evolve).set_integrator('vode', method='bdf',rtol=self.deqtol)
        self.cc_alpha_integrator = ode(cc_alpha_evolve).set_integrator('vode', method='bdf',rtol=self.deqtol)
        self.ci_alpha_integrator = ode(ci_alpha_evolve).set_integrator('vode', method='bdf',rtol=self.deqtol)


    def setUpInts(self):

        # Read in the integrals
        oneh, eri_in, attrs = self.loadHDF()
        self.attrs = attrs

        # Do integral transformation
        self.h1, self.evecs = IntTran2(oneh)
        self.eri = IntTran4(eri_in,self.evecs)

        # Define fugacity
        self.fug = self.n_elec / (self.nso - self.n_elec)

        # set beta_step
        self.beta_step = self.beta_f/(self.beta_pts - 1)

    def setAmps(self,cc_amps,ci_amps):

        # Check the length of the incoming amplitudes
        len_req = int(1 + self.nso**2 + comb(self.nso,2)**2)
        if (len(cc_amps) != len_req) or (len(ci_amps) != len_req):
            raise ValueError('Invalid length of {} or {} amplitudes'.format(cc_amps,ci_amps))

        # Initial cc, ci amplitudes for the current step
        self.cc_amps = cc_amps
        self.ci_amps = ci_amps

    def DoBetaIntegration(self):

        #
        # class level Beta integration function
        # 

        cc_amps, ci_amps = _do_beta_integration(
            [self.cc_beta_integrator, self.ci_beta_integrator], 
            [self.cc_amps, self.ci_amps],
            [self.beta_in, self.alpha_in], 
            self.beta_step, self.fug, self.h1, self.eri
        )

        self.beta_in += self.beta_step
        self.cc_amps = cc_amps
        self.ci_amps = ci_amps

    def BisectionAndAlphaIntegrate(self):

        #
        # class level Beta integration function
        # 

        intgrs = [self.cc_alpha_integrator, self.ci_alpha_integrator]
        amps = [self.cc_amps, self.ci_amps]
        be_al = [self.beta_in, self.alpha_in]

        cc_amps, ci_amps, al_mid = _do_alpha_integration(
            intgrs,amps,be_al, self.fug, self.h1, self.n_elec, self.ntol
        )

        self.alpha_in = al_mid
        self.cc_amps = cc_amps
        self.ci_amps = ci_amps

