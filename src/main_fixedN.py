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

import pdb

#
# GLOBAL VARIABLES
#

mu_step_0 = +5e-2
len_t1 = 0
len_t2 = 0


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

    # Initialize the I/O module
    iops = IOps(inp_file='Input')

    # Read integrals
    h1_in, eri_in, attrs  = iops.loadHDF()

    # Input Parameters
    nso = iops.nso
    n_elec = iops.n_elec
    beta_f = iops.beta_f
    beta_pts = iops.beta_pts
    ntol = iops.ntol
    deqtol = iops.deqtol
    e_nuc = iops.e_nuc
    
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

