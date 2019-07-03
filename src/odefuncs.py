import sys
sys.path.append('../fort_src/')

import h5py, time, numpy as np

from scipy.integrate import ode
from scipy.misc import comb
from multiprocessing import Pool

from ThermalCCSD import *
from ThermalCISD import *

import pdb

#
# GLOBAL VARIABLES
#

t2_symm_tol = 5e-8

mu_step_0 = +5e-2
len_t1 = 0
len_t2 = 0


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
# ODE SOLVER FUNCTIONS FOR CC - BETA AND MU EVOLUTIONS
#

# Here, we want to write all the ODE solving function which will repeat and clutter the Main
# function otherwise.

def DoIntegration(integrator, x_final):
    """
    Intermediate function to perform integration from x_initial to x_final.
    The x_initial and all other information is inherently included in the integrator.
    """
    yout = integrator.integrate(x_final)
    return yout

def cc_and_ci_beta_integrate(cc_integrator, ci_integrator, bspan, y_cc, y_ci, mu_val, alpha, h1, eri):
    """
    Function to perform the task of integrating the BETA evolution equations for zeroth, first
    and second rank CC amplitudes..

    Inputs:
        cc_integrator   ::  The integrator, an instance of scipy.integrate.ode
        ci_integrator   ::  The integrator, an instance of scipy.integrate.ode
        bspan           ::  Array of length 2 containing initial and final beta value
        y_cc            ::  Initial condition at input BETA i.e. first bspan value
        y_ci            ::  Initial condition at input BETA i.e. first bspan value
        alpha           ::  initial value of the ALPHAA / MU parameter
        h1              ::  One electron integrals
        eri             ::  Two electron integrals

    Outputs:
        fin_cc          ::  Final solution at the respective input order
        fin_ci          ::  Final solution at the respective input order
    """

    ### Check for any possible errors in the input
    if not isinstance(cc_integrator,ode):
        raise ValueError('The argument',cc_integrator,'  needs to be an instance of ',ode)
    if not isinstance(ci_integrator,ode):
        raise ValueError('The argument',ci_integrator,'  needs to be an instance of ',ode)
    if len(bspan) != 2:
        raise ValueError('Inappropriate length of the input array ',bspan)

    # Set the initial conditions
    cc_integrator.set_initial_value(y_cc, bspan[0]).set_f_params(mu_val, alpha, h1, eri)
    ci_integrator.set_initial_value(y_ci, bspan[0]).set_f_params(mu_val, alpha, h1, eri)

    # Create the pool of processes
    p1 = Pool( processes=1 )
    p1_async = p1.apply_async(DoIntegration, args=(cc_integrator, bspan[1]))
    p1.close()

    p2 = Pool( processes=1 )
    p2_async = p2.apply_async(DoIntegration, args=(ci_integrator, bspan[1]))
    p2.close()

    p1.join()
    p2.join()

    fin_cc = p1_async.get(timeout=1)
    fin_ci = p2_async.get(timeout=1)

    return fin_cc, fin_ci


def mu_find_and_integrate(cc_integrator, ci_integrator, mu_in, y_cc, y_ci, nelec, beta, alpha, h1):
    """
    Function that will first findt the chemical potential bracket and then, integrate to the
    correct mu and return the final TAMPS.

    """

    ### Check for any possible errors in the input
    if not isinstance(cc_integrator,ode):
        raise ValueError('The argument',cc_integrator,'  needs to be an instance of ',ode)
    if not isinstance(ci_integrator,ode):
        raise ValueError('The argument',ci_integrator,'  needs to be an instance of ',ode)

    # Global variables
    global mu_step_0
    global ntol

    nso = len(h1)

    # HFB coefficients
    x = 1/np.sqrt( 1 + np.exp( -beta*h1 + mu_in )*alpha )
    y = np.exp( -( beta * h1 - mu_in )/ 2 ) * x * np.sqrt( alpha )

    # find the initial expectation value, i.e. <N>
    num = LinRespNumber(
        y_cc, y_ci, x, y
    )

    print('\t\tNumber of particles after the beta evolution = {}'.format(num))

    ndiff_sgn = np.sign( num - nelec )
    ndiff_mag = np.abs( num - nelec )

    yf_cc = y_cc
    yf_ci = y_ci
    mu_f = mu_in
    
    # if the number is already converged, then there is no need to do any of the following
    #   and hence we keep an 'if' statement; if the condition evaluates to FALSE, then
    #   the outputs will be yf, mu_f

    mu_0 = mu_in

    if ndiff_mag > ntol:

        mu_step = mu_step_0

        # Obtain the bracket to perform BISECTION

        mu_1 = mu_0
        mu_2 = mu_0 + mu_step

        yf_cc1 = y_cc
        yf_cc2 = y_cc

        yf_ci1 = y_ci
        yf_ci2 = y_ci

        count = 0
        sp_count = 0
        
        while np.sign(num - nelec) == ndiff_sgn:
            count += 1
            if count > 300:
                print('Could not find the bracket after 300 steps')
                count = 0
                exit()
                break

            mu_span = [mu_1, mu_2]
            
            # Define the solver for CC and CI

            yf_cc1 = yf_cc2
            yf_ci1 = yf_ci2

            cc_integrator.set_initial_value(yf_cc1, mu_span[0])
            cc_integrator.set_f_params(beta, alpha, h1)
            ci_integrator.set_initial_value(yf_ci1, mu_span[0])
            ci_integrator.set_f_params(beta, alpha, h1)

            # Evolve to the mu_f -- DO Integration
            p1 = Pool( processes = 1)
            p1_async = p1.apply_async(DoIntegration, args=(cc_integrator, mu_span[1]))
            p1.close()

            p2 = Pool( processes = 1)
            p2_async = p2.apply_async(DoIntegration, args=(ci_integrator, mu_span[1]))
            p2.close()

            p1.join() #---Wait for the pool process to finish
            p2.join() #---Wait for the pool process to finish

            yf_cc2 = p1_async.get(timeout=1)
            yf_ci2 = p2_async.get(timeout=1)

            # HFB coefficients
            x = 1/np.sqrt( 1 + np.exp( -beta*h1 + mu_span[1] )*alpha )
            y = np.exp( -( beta * h1 - mu_span[1] )/ 2 ) * x * np.sqrt( alpha )

            # Evaluate the Number Expectation
            num = LinRespNumber(
                yf_cc2, yf_ci2, x, y
            )

            # Finer grid if we are closer to nelec
            val = np.abs(num - nelec) - ndiff_mag
            if (val>0):
                if val<1e-1:
                    sp_count += 1
                else:
                    mu_step_0 *= -1
                    mu_step *= -1
                
                if sp_count >= 10:
                    mu_step_0 *= -1
                    mu_step *= -1
                    sp_count = 0

            ndiff_mag = np.abs(num - nelec)

            # Set up for next iteration
            mu_1 = mu_2
            mu_2 = mu_2 + mu_step

            if np.abs(num - nelec) <= ntol:
                break

            print('\t\t\tStart value of Mu = {}'.format(mu_span[0]))
            print('\t\t\tEnd value of Mu = {}'.format(mu_span[1]))
            print('\t\t\tNumber of particles after evolution = {}'.format(num))
            print('\t\t\t----------------------------------------------\n')

        print('Bracket found between mu = {} and mu = {}'.format(mu_span[0],mu_span[1]))

        # Bisection bracket
        mu_bisect = mu_span
        mu_mid = mu_bisect[1]

        ndiff_sgn2 = np.sign(num - nelec)
        ndiff_sgn1 = -ndiff_sgn2

        while np.abs(num - nelec)>ntol:
            # Set the initial condition for the mu_solver
            cc_integrator.set_initial_value(yf_cc1, mu_bisect[0]).set_f_params(beta, alpha, h1)
            ci_integrator.set_initial_value(yf_ci1, mu_bisect[0]).set_f_params(beta, alpha, h1)

            # Evolve the ODE to mid point of the bracket
            mu_mid = np.mean(mu_bisect)

            p1 = Pool( processes=1 )
            p1_async = p1.apply_async(DoIntegration, args=(cc_integrator, mu_mid))
            p1.close()

            p2 = Pool( processes=1 )
            p2_async = p2.apply_async(DoIntegration, args=(ci_integrator, mu_mid))
            p2.close()

            p1.join()
            p2.join()

            yf_cc_mid = p1_async.get(timeout=1)
            yf_ci_mid = p2_async.get(timeout=1)

            # HFB coefficients
            x = 1/np.sqrt( 1 + np.exp( -beta*h1 + mu_mid )*alpha )
            y = np.exp( -( beta * h1 - mu_mid )/ 2 ) * x * np.sqrt( alpha )

            # Compute the number and update the bracket
            num = LinRespNumber(
                yf_cc_mid, yf_ci_mid, x, y
            )

            if np.sign( num - nelec ) == ndiff_sgn1:
                mu_bisect[0] = mu_mid
                yf_cc1 = yf_cc_mid
                yf_ci1 = yf_ci_mid
            else:
                mu_bisect[1] = mu_mid
                yf_cc2 = yf_cc_mid
                yf_ci2 = yf_ci_mid

        print('Bisection converges to mu = {}'.format(mu_mid))

        # Now that we have found the MU_MID, we can use one-shot evolution to avoid any added errors
        mu_f = mu_mid
        cc_integrator.set_initial_value(yf_cc, mu_0).set_f_params(beta, alpha, h1)
        ci_integrator.set_initial_value(yf_ci, mu_0).set_f_params(beta, alpha, h1)

        p1 = Pool( processes = 1 )
        p1_async = p1.apply_async(DoIntegration, args=(cc_integrator, mu_f))
        p1.close()

        p2 = Pool( processes = 1 )
        p2_async = p2.apply_async(DoIntegration, args=(ci_integrator, mu_f))
        p2.close()

        p1.join()
        p2.join()

        yf_cc = p1_async.get(timeout=1)
        yf_ci = p2_async.get(timeout=1)

    return yf_cc, yf_ci, mu_f


#
# MAIN
# 

def main():
    start_time = time.time()
    
    #################################################################
    #                       READ INPUT                              #
    #       Get the data first - including the attributes           #
    #################################################################

    print('\n\n----------------------------------------------\n')
    print('Reading Input')
    print('----------------------------------------------\n')

    fn, n_elec, beta_f, beta_pts, e_nuc = ParseInput(enuc=True)

    h1_in, eri_in, attrs  = loadHDF(fname=fn)
    nso = np.size(h1_in,axis=0)
    
    input_time = time.time()

    print('Number of Spin Orbitals:',nso)
    print('-------------------------------------------------------\n')

    #################################################################
    #                   INTEGRAL TRANSFORMATION                     #
    #       Transforming the One and Two Electron Integral          #
    #       so that the OneH or h1 is diagonal                      #
    #################################################################

    # Note that h1 is a 1D array of the eigenvalues
    h1, evecs = IntTran2(h1_in)
    eri = IntTran4(eri_in,evecs)

    #################################################################
    #               INITIAL HFB and OTHER  PARAMETERS               #
    #       Initializing the important x and y HFB parameters       #
    #       and other arrays for a better idea further              #
    #################################################################

    # HFB parameters in the initial thermal vacuum reference at Beta = 0
    # alpha = exp( beta * mu ) 
    # when beta -> 0 and mu -> inf
    alpha = n_elec/(nso-n_elec)

    x = np.ones(nso)/np.sqrt(1+alpha)
    y = np.ones(nso)*np.sqrt(alpha)/np.sqrt(1+alpha)
    
    # Initializing the amplitude vectors - only the unique matrix
    # elements in the t1 and t2 matrices
    global len_t1
    global len_t2

    len_t1 = int(nso**2)
    len_t2 = int(comb(nso,2)**2)

    t0 = 0.0
    t1_vec = np.zeros(len_t1)
    t2_vec = np.zeros(len_t2)

    z0 = 0.0
    z1_vec = np.zeros(len_t1)
    z2_vec = np.zeros(len_t2)

    # These are the actual t1 and t2 matrices
    t1 = np.reshape(t1_vec,(nso,nso))
    t2 = T2_Decompress(t2_vec,nso)

    # These are the actual z1 and z2 matrices
    z1 = np.reshape(z1_vec,(nso,nso))
    z2 = T2_Decompress(z2_vec,nso)

    # Hartree Fock Energy at BETA = 0
    E_hf = LinRespEnergy(
        h1, eri, np.concatenate(
            ([t0],t1_vec,t2_vec)
        ), np.concatenate(
            ([z0],z1_vec,z2_vec)
        ), x, y
    )

    # Beta Grid
    beta_0 = 0
    beta_step = beta_f / (beta_pts-1)
    beta_grid = np.array([beta_0])

    # Alpha Grid
    mu_0 = np.log(alpha)
    mu_step_0 = +5e-2
    mu_f = mu_step_0

    n_data = 1

    #################################################################
    #               SETTING UP THE INITIAL VALUE PROB               #
    #    Compute the Hartree fock energy and define the y-stacks    #
    #    which will eventually be involved in computation           #
    #################################################################

    #TODO:  First initialize the various variables needed for linear 
    #       response and then go along in the Main() and then in the
    #       sub-functions if needed

    # initial condition for Tamps at beta = 0 = mu
    y0_cc = np.concatenate(([t0],t1_vec,t2_vec))
    y0_ci = np.concatenate(([z0],z1_vec,z2_vec))

    # Data to be output at BETA GRID POINT
    e_cc = np.zeros(n_data)

    mu_cc = np.zeros(n_data)
    chem_cc = 0.0 # Dummy for mu_cc

    n_cc = np.zeros(n_data)

    t0_rms = np.zeros(n_data)
    t1_rms = np.zeros(n_data)
    t2_rms = np.zeros(n_data)

    z0_rms = np.zeros(n_data)
    z1_rms = np.zeros(n_data)
    z2_rms = np.zeros(n_data)


    e_cc[0] = E_hf + e_nuc

    # Time stamp to record that all the initialization process is completed
    init_time = time.time()

    # Print the first things
    num = LinRespNumber(
        y0_cc, y0_ci, x, y
    )
    n_cc[0] = num

    # Confirm and Print the Number of particles at Beta = 0
    print('Thermal Hartree Fock Energy at T = Inf is :',E_hf+e_nuc)
    print('Number of Particles = {}'.format(num))
    print('-------------------------------------------------------\n')

    
    #################################################################
    #                   FILE OUTPUT HANDLE CREATOR                  #
    #   Creating the h5py file headers and update the files after   #
    #   each iteration - so that if the code blows up we know       #
    #   where it went wrong.                                        #
    #################################################################
    fout = fn[0:-7] + 'thermal_ccsd_cov.h5'
    
    print('Writing output to {}'.format(fout))

    fp1 = File(fout,'w')
    fp2 = File(fout[:-3]+'_amps.h5','w')

    # Data set to output the evecs -- sometimes even if we have diagonal, python goes crazy and shuffles the order
    # of the degenerate orbitals.
    fp1.create_dataset('evecs',data=evecs)

    dsets = ['beta','e_cc','mu_cc','n_cc','t0rms','t1rms','t2rms','z0rms','z1rms','z2rms']

    # Create all but the ystack data sets

    dset_list = createh5(fp1, dsets, beta_pts, attrs)
    
    # Create h5 for t1, t2, z1 and z2
    t1_dset = np.zeros((beta_pts,nso**2))
    t2_dset = np.zeros((beta_pts,int(comb(nso,2)**2)))
    z1_dset = np.zeros((beta_pts,nso**2))
    z2_dset = np.zeros((beta_pts,int(comb(nso,2)**2)))

    # Updating the first data point
    vals = [
        0, e_cc, mu_cc, num, t0_rms, t1_rms, t2_rms, z0_rms, z1_rms, z2_rms
    ]

    updateh5(dset_list, vals, 0)

    t1_dset[0,:] = y0_cc[1:nso*nso+1]
    t2_dset[0,:] = y0_cc[1+nso*nso:]
    z1_dset[0,:] = y0_ci[1:nso*nso+1]
    z2_dset[0,:] = y0_ci[1+nso*nso:]


    #################################################################
    #                   TIME AND CHEMPOT EVOLUTION                  #
    #   The actual evolution takes place in this section of code.   #
    #   First evolve in Beta Direction to a grid point and then     #
    #   evolve in the Mu Direction to fix the Number of Particles   #
    #################################################################

    cc_beta_solver = ode(cc_beta_evolve).set_integrator('vode',method='bdf',rtol=deqtol)
    cc_mu_solver = ode(cc_mu_evolve).set_integrator('vode',method='bdf',rtol=deqtol)

    ci_beta_solver = ode(ci_beta_evolve).set_integrator('vode',method='bdf',rtol=deqtol)
    ci_mu_solver = ode(ci_mu_evolve).set_integrator('vode',method='bdf',rtol=deqtol)

    j = 1

    while beta_grid[j-1]<beta_f:

        beta_2 = beta_grid[j-1] + beta_step
        b_span = [beta_grid[j-1], beta_2]

        ############################### 
        # 1: Beta evolution 
        ############################### 

        print('\t\tBeta-Evolution for CCSD')

        yf_cc, yf_ci = cc_and_ci_beta_integrate(
            cc_beta_solver, ci_beta_solver, b_span, y0_cc, y0_ci, chem_cc, alpha, h1, eri
        )

        y0_cc = yf_cc
        y0_ci = yf_ci

        print('\t\tNew Beta = {}'.format(beta_2))

        ############################### 
        # 2: Mu or ChemPot evolution 
        ############################### 
        
        print('\t\tMu-Evolution for CCSD')
        yf_cc, yf_ci, chem_cc = mu_find_and_integrate(
            cc_mu_solver, ci_mu_solver, chem_cc, y0_cc, y0_ci, n_elec, b_span[1], alpha, h1
        )

        ###############################################
        # 3: Update H5 and set up for the next loop
        ###############################################

        # Setting up for next loop
        y0_cc = yf_cc
        y0_ci = yf_ci

        # Computing the quantities of interest
        beta_grid = np.append(beta_grid,b_span[1])

        # Update the HFB Coefficients for CIS
        x = 1/np.sqrt(1 + np.exp( -b_span[1]*h1 + chem_cc )*alpha )
        y = np.exp( ( -b_span[1]*h1 + chem_cc )/2 )*x*np.sqrt(alpha)

        # amplitudes, energy and chem pot
        t0 = y0_cc[0]
        t1 = np.reshape( y0_cc[1:1+nso**2], (nso, nso) )
        t2 = T2_Decompress( y0_cc[1+nso**2:], nso)
        z0 = y0_ci[0]
        z1 = np.reshape( y0_ci[1:1+nso**2], (nso, nso) )
        z2 = T2_Decompress( y0_ci[1+nso**2:], nso)

        en = LinRespEnergy(
            h1, eri, y0_cc, y0_ci, x, y 
        )

        num = LinRespNumber(
            y0_cc, y0_ci, x, y
        )

        e_cc = np.append( e_cc, en + e_nuc)
        n_cc = np.append( n_cc, num )

        t0_rms = np.append( t0_rms, t0 )
        t1_rms = np.append( t1_rms, np.sqrt( np.mean( t1**2 ) ) )
        t2_rms = np.append( t2_rms, np.sqrt( np.mean( t2**2 ) ) )

        z0_rms = np.append( z0_rms, z0 )
        z1_rms = np.append( z1_rms, np.sqrt( np.mean( z1**2 ) ) )
        z2_rms = np.append( z2_rms, np.sqrt( np.mean( z2**2 ) ) )

        mu_cc = np.append( mu_cc, chem_cc)

        # collecting everything together
        vals = [
            beta_grid[-1], e_cc[-1], mu_cc[-1], n_cc[-1], t0_rms[-1], t1_rms[-1], t2_rms[-1], z0_rms[-1], z1_rms[-1], z2_rms[-1], 
        ]

        # update
        updateh5( dset_list, vals, j)

        t1_dset[j,:] = y0_cc[1:1+nso*nso]
        t2_dset[j,:] = y0_cc[1+nso*nso:]
        z1_dset[j,:] = y0_ci[1:1+nso*nso]
        z2_dset[j,:] = y0_ci[1+nso*nso:]

        print('At beta = {}'.format(b_span[-1]))
        print('mu = {}'.format(chem_cc))
        print('N_elec = {}'.format(num))
        print('ECC = {}'.format(e_cc[-1]))
        print('t1rms = {}'.format(t1_rms[-1]))
        print('t2rms = {}'.format(t2_rms[-1]))
        print('Max(t1) = {}'.format(np.max(t1)))
        print('Max(t2) = {}'.format(np.max(t2)))
        print('Min(X) = {}'.format(np.min(x)))
        print('Min(Y) = {}'.format(np.min(y)))
        print('----------------------------------------------\n')

        j += 1

        if j == beta_pts:
            break

    #################################################################
    #                                                               #
    #       CLOSING THE CODE AND PRINT TIMING                       #
    #                                                               #
    #################################################################

    ode_time = time.time()
    fp1.close()

    t1_amps = fp2.create_dataset('t1_amps',data=t1_dset)
    t2_amps = fp2.create_dataset('t2_amps',data=t2_dset)
    z1_amps = fp2.create_dataset('z1_amps',data=z1_dset)
    z2_amps = fp2.create_dataset('z2_amps',data=z2_dset)

    fp2.close()
    end_time = time.time()

    print('----------------------------------------------')
    print('End of Program and Time Profile')
    print('----------------------------------------------')
    print('----------------------------------------------')
    print('Process              Time')
    print('----------------------------------------------')
    print('Reading Input        {}'.format(input_time-start_time))
    print('Init Variables       {}'.format(init_time-input_time))
    print('Solve ODE for grid   {}'.format(ode_time-init_time))
    print('Writing the Output   {}'.format(end_time-ode_time))
    print('Total Time           {}'.format(end_time-start_time))
    print('----------------------------------------------')
    print('----------------------------------------------\n')


if __name__ == '__main__':
    main()

