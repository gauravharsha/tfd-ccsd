import numpy as np


def IntTran4(ERI, Evecs):
    """Function to transform the Two Body Hamlitonian Matrix
    using the transformation Evecs."""

    ERI2 = np.zeros(np.shape(ERI))
    ERI3 = np.zeros(np.shape(ERI))

    # number of spin orbitals
    nso = np.size(Evecs, axis=0)

    # Get the last index elements
    for p in range(nso):
        for q in range(nso):
            for r in range(nso):
                for a in range(nso):
                    ERI2[p, q, r, a] = np.dot(ERI[p, q, r, :], Evecs[:, a])

    # Get the 3rd index transformed
    for p in range(nso):
        for q in range(nso):
            for a in range(nso):
                for s in range(nso):
                    ERI3[p, q, a, s] = np.dot(ERI2[p, q, :, s], Evecs[:, a])

    # Get the 2nd index transformed
    for p in range(nso):
        for a in range(nso):
            for r in range(nso):
                for s in range(nso):
                    ERI2[p, a, r, s] = np.dot(ERI3[p, :, r, s], Evecs[:, a])

    # Get the 1st index transformed
    for a in range(nso):
        for q in range(nso):
            for r in range(nso):
                for s in range(nso):
                    ERI3[a, q, r, s] = np.dot(ERI2[:, q, r, s], Evecs[:, a])

    # all done and return
    return ERI3


def IntTran2(OneH):
    """Diagonalize the One Body Hamiltonian in the MO basis
    and return both the diagonalized Matrix and the Evecs."""

    evecs = np.zeros(np.shape(OneH))

    # diagonalize h1
    fock, evecs = np.linalg.eigh(OneH)

    return fock, evecs
