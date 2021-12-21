import numpy as np
import time
import h5py
from scipy.special import comb
from tfdccsd.odefuncs import Evolution, eval_energy, eval_number

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

    print('\n==============================================================')
    print('\tReading Input')
    print('==============================================================')

    # Initialize the Evolution module
    evol = Evolution(inp_file='Input', alpha_step=0.1)

    # Extract the parameters
    nso = evol.nso
    beta_pts = evol.beta_pts
    e_nuc = evol.e_nuc
    fug = evol.fug

    # Integrals
    eigs = evol.eigs
    h1 = evol.h1
    eri = evol.eri

    input_time = time.time()

    print('Number of Spin Orbitals:', nso)
    print('--------------------------------------------------------------\n')

    #################################################################
    #               INITIALIZE VARIABLES / PARAMETERS               #
    #       Initializing the arrays for computing various values    #
    #################################################################

    global len_t1
    global len_t2

    len_t1 = int(nso**2)
    len_t2 = int(comb(nso, 2)**2)

    # CC and CI amps
    cc_amps = np.zeros(1 + len_t1 + len_t2)
    ci_amps = np.zeros(1 + len_t1 + len_t2)

    # HFB parameters
    x = 1/np.sqrt(1 + np.exp(-evol.beta_in*eigs + evol.alpha_in)*fug)
    y = np.sqrt(1 - x**2)

    # Make arrays for energy, amplitudes, etc.
    e_cc = np.zeros(evol.beta_pts)
    n_cc = np.zeros(evol.beta_pts)

    #################################################################
    #                   SETUP OUTPUT H5PY FILES                     #
    #################################################################

    output_fn = evol.fn.replace('input', 'output')
    output_fn = output_fn.replace('_data.h5', '_tfd_ccsd.h5')

    fout = h5py.File(output_fn, 'w')

    output_dsets = [
        'beta', 'alpha', 'e_cc', 'n_cc'
    ]

    evol.createh5(fout, output_dsets, beta_pts, evol.attrs)

    # the amplitudes datasets will have to be handled separately
    fout.create_dataset('cc_amps', (beta_pts, 1+len_t1+len_t2))
    fout.create_dataset('ci_amps', (beta_pts, 1+len_t1+len_t2))

    #################################################################
    #                   INITIAL VALUE CHECK                         #
    #################################################################

    # Beta point index
    i_beta = 0

    # Hartree Fock Energy at Beta = 0, Alpha = 0
    e_cc[i_beta] = e_nuc + eval_energy(h1, eri, cc_amps, ci_amps, x, y)
    n_cc[i_beta] = eval_number(cc_amps, ci_amps, x, y)

    # Update data_sets
    vals = [
        evol.beta_in, evol.alpha_in, e_cc[i_beta], n_cc[i_beta]
    ]
    evol.updateh5(vals, i_beta)

    # The CC and CI amplitudes at initial value are already ZEROS

    print('Input File Name = ', evol.fn)
    print('Output File Name = ', output_fn)
    print('--------------------------------------------------------------\n')
    print('Initial Beta value = ', evol.beta_in)
    print('Initial Alpha value = ', evol.alpha_in)
    print('Internal energy at inf Temperature = ', e_cc[0])
    print('Number of particles at inf Temperature = ', n_cc[0])
    print('--------------------------------------------------------------\n')

    init_time = time.time()
    #################################################################
    #                   PERFORM THE ODE EVOLUTION                   #
    #    Integrate the ODE's to preserve the number of particles    #
    #################################################################

    while i_beta < evol.beta_pts-1:

        # Beta point index
        i_beta += 1

        # Do Beta Integration
        evol.DoBetaIntegration()

        print('Beta = ', evol.beta_in)

        # Do Alpha search and integration
        evol.BisectionAndAlphaIntegrate()

        # New HFB parameters
        x = 1/np.sqrt(1 + np.exp(-evol.beta_in*eigs + evol.alpha_in)*fug)
        y = np.sqrt(1 - x**2)

        # Extract the amplitudes
        cc_amps = evol.cc_amps
        ci_amps = evol.ci_amps

        # Find Energy and Number
        e_cc[i_beta] = e_nuc + eval_energy(h1, eri, cc_amps, ci_amps, x, y)
        n_cc[i_beta] = eval_number(cc_amps, ci_amps, x, y)

        # Check by printing
        print('Beta = ', evol.beta_in)
        print('Alpha = ', evol.alpha_in)
        print('Internal energy at new step = ', e_cc[i_beta])
        print('Number of particles at new step = ', n_cc[i_beta])
        print('------------------------------------------------------------\n')

        # Write data to the output files
        vals = [evol.beta_in, evol.alpha_in, e_cc[i_beta], n_cc[i_beta]]
        evol.updateh5(vals, i_beta)
        fout['cc_amps'][i_beta] = cc_amps
        fout['ci_amps'][i_beta] = ci_amps

    ode_time = time.time()

    #################################################################
    #                       CLOSE EVERYTHING                        #
    #################################################################
    fout.close()

    end_time = ode_time

    # Time Profile
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

    exit()


if __name__ == '__main__':
    main()
