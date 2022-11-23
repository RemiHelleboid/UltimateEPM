# Partially from https://gist.github.com/wmedlar/

import scipy.constants as cst
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

import itertools

import BZ_Sampling as bzS

cst.epsilon_0

KINETIC_CONSTANT = cst.hbar**2 / (2 * cst.m_e * cst.e)

def generate_basis_vector(N_basis: int) -> np.ndarray:
    """Generate the basis vector.

    Args:
        N_basis (int): Number of basis vectors.

    Returns:
        np.ndarray: The basis vector.
    """
    b1 = np.array([-1, 1, 1])
    b2 = np.array([1, -1, 1])
    b3 = np.array([1, 1, -1])

    basis_G = []
    basis_Gnorm = []

    for i in range(-N_basis, N_basis+1):
        for j in range(-N_basis, N_basis+1):
            for k in range(-N_basis, N_basis+1):
                G = i * b1 + j * b2 + k * b3
                norme_square = np.linalg.norm(G)**2
                if (norme_square <= 24):
                    basis_G.append(G)
                    basis_Gnorm.append(norme_square)
    # Sort the basis vector by the norm of the vector.
    basis_G = np.array(basis_G)
    basis_Gnorm = np.array(basis_Gnorm)
    basis_G = basis_G[np.argsort(basis_Gnorm)]
    return basis_G



def get_basis_index_non_zero(basis_G: np.ndarray) -> np.ndarray:
    list_couples_non_zero = []
    for i in range(len(basis_G)):
        for j in range(len(basis_G)):
            G1 = basis_G[i]
            G2 = basis_G[j]
            Gdiff = G1 - G2
            norm_square = np.linalg.norm(Gdiff)
            if (np.linalg.norm(Gdiff) < 1e-6):
                list_couples_non_zero.append((i, j))
            

    
def PseudoPotential(G: np.ndarray) -> float:
    Ry_to_eV = 13.605698066
    V3S = -0.2241 * Ry_to_eV
    V8S = 0.0551 * Ry_to_eV
    V11S = 0.0724 * Ry_to_eV
    tau = np.array([1, 1, 1]) / 8.0
    Gtau = 2.0 * np.pi * np.dot(G, tau)
    Gsquare = np.linalg.norm(G)**2
    if np.abs(Gsquare - 3) < 1e-6:
        return V3S * np.cos(Gtau)
    elif np.abs(Gsquare - 8) < 1e-6:
        return V8S * np.cos(Gtau)
    elif np.abs(Gsquare - 11) < 1e-6:
        return V11S * np.cos(Gtau)
    else:
        return 0.0


def non_diag_hamiltonian(basis_G: np.ndarray) -> np.ndarray:
    """Generate the Hamiltonian matrix.

    Args:
        k (np.ndarray): The wave vector.
        basis_G (np.ndarray): The basis vector.

    Returns:
        np.ndarray: The Hamiltonian matrix.
    """
    a_0 = 5.431 * 1e-10
    BasisSize = len(basis_G)
    Hamiltonian = np.zeros((BasisSize, BasisSize))
    for i in range(BasisSize):
        for j in range(i+1):
            G1 = basis_G[i]
            G2 = basis_G[j]
            Gdiff = G1 - G2
            Hamiltonian[i, j] = PseudoPotential(Gdiff)
    return Hamiltonian


EPM_BASIS = generate_basis_vector(10)
EPM_BASIS_NORM = np.linalg.norm(EPM_BASIS, axis=1)
GlobHamiltonianNonDiag = non_diag_hamiltonian(EPM_BASIS)


def hamiltonian(k: np.ndarray, basis_G: np.ndarray) -> np.ndarray:
    """Generate the Hamiltonian matrix.

    Args:
        k (np.ndarray): The wave vector.
        basis_G (np.ndarray): The basis vector.

    Returns:
        np.ndarray: The Hamiltonian matrix.
    """
    a_0 = 5.431 * 1e-10
    BasisSize = len(basis_G)
    kinetic_c = (2.0 * np.pi / a_0)**2
    Hamiltonian = np.copy(GlobHamiltonianNonDiag)
    for i in range(BasisSize):
        G1 = basis_G[i]
        Hamiltonian[i, i] = KINETIC_CONSTANT * kinetic_c * np.linalg.norm(G1 + k)**2
    return Hamiltonian

def eigen_states(k: np.ndarray, basis_G: np.ndarray) -> np.ndarray:
    """Generate the eigen states.

    Args:
        k (np.ndarray): The wave vector.
        basis_G (np.ndarray): The basis vector.

    Returns:
        np.ndarray: The eigen states.
    """
    HamiltonianMatrix = hamiltonian(k, basis_G)
    eigen_values = la.eigvalsh(HamiltonianMatrix, lower=True)
    return eigen_values

# symmetry points in the Brillouin zone
G = np.array([0, 0, 0])
L = np.array([1/2, 1/2, 1/2])
K = np.array([3/4, 3/4, 0])
X = np.array([0, 0, 1])
W = np.array([1, 1/2, 0])
U = np.array([1/4, 1/4, 1])

def linpath(a, b, n=50, endpoint=True):
    '''
    Creates an array of n equally spaced points along the path a -> b, inclusive.

    args:
        a: An iterable of numbers that represents the starting position.
        b: An iterable of numbers that represents the ending position.
        n: The integer number of sample points to calculate. Defaults to 50.
        
    returns:
        A numpy array of shape (n, k) where k is the shortest length of either
        iterable -- a or b.
    '''
    # list of n linear spacings between the start and end of each corresponding point
    spacings = [np.linspace(start, end, num=n, endpoint=endpoint) for start, end in zip(a, b)]
    
    # stacks along their last axis, transforming a list of spacings into an array of points of len n
    return np.stack(spacings, axis=-1)

def k_path(n):
    lambd = linpath(L, G, n, endpoint=False)
    delta = linpath(G, X, n, endpoint=False)
    x_uk = linpath(X, U, n // 4, endpoint=False)
    sigma = linpath(K, G, n, endpoint=True)
    path=[lambd, delta, x_uk, sigma]
    return path


def band_structure(path):
    bands = []
    # vstack concatenates our list of paths into one nice array
    for idx, k in enumerate(np.vstack(path)):
        print(f"\r{idx} values over {len(np.vstack(path))}", end="")
        E = eigen_states(k, EPM_BASIS)
        # picks out the lowest eight eigenvalues
        bands.append(E[:8])
    
    return np.stack(bands, axis=-1)



def dielectric_function_mc(energy, q_vect, n_valence, n_conduction, Nk):
    list_colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    first_condcution_idx = 4
    list_kpoints = bzS.RandomKPointsIrreducibleWedge(Nk)
    list_k_plus_q = np.array([k + q_vect for k in list_kpoints])
    basis_vectors = EPM_BASIS
    basis_size = len(basis_vectors)
    sum_total = 0.0
    list_eps = []
    list_k_contrib = []
    for index_k in range(Nk):
        # print(f"\r{index_k} values over {Nk}", end="")
        k = list_kpoints[index_k]
        k_plus_q = list_k_plus_q[index_k]
        H_k = hamiltonian(k, basis_vectors)
        H_k_plus_q = hamiltonian(k_plus_q, basis_vectors)
        eigval_k, eigvec_k = la.eigh(H_k)
        eigval_k_plus_q, eigvec_k_plus_q = la.eigh(H_k_plus_q)
        sum_k = 0.0
        for idx_cond in range(first_condcution_idx, first_condcution_idx + n_conduction):
            # plt.plot(EPM_BASIS_NORM**2, eigvec_k_plus_q[:, idx_cond], '.', c=list_colors[idx_cond], alpha=0.51, markersize=5)
            # plt.pause(0.1)
            for idx_val in range(n_valence):
                sum_k += np.abs(np.dot(eigvec_k[:, idx_val].T, eigvec_k_plus_q[:, idx_cond]))**2
            # print(sum_k)
        factor = 1.0 / (eigval_k_plus_q[idx_cond] - eigval_k[idx_val] - energy) + 1.0 / (eigval_k_plus_q[idx_cond] - eigval_k[idx_val] + energy)
        sum_k *= factor
        list_k_contrib.append(sum_k)
        total_contrib = sum(list_k_contrib) / len(list_k_contrib)
        epsilon_tot = 1.0 + ((4.0*np.pi) / (np.linalg.norm(q_vect)**2)) * total_contrib
        list_eps.append(epsilon_tot)
    # fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    #     print(f"\r{index_k} values over {Nk} ----> epsilon_real = {list_eps[-1]:.2e}", end="", flush=True)
    # axs.plot(list_eps, '-o', color='blue')
    # plt.show()
    return list_eps[-1]

def dielectric_function_Monkhorst_Pack(energy, q_vect, n_valence, n_conduction, Nxyz):
    first_condcution_idx = 4
    list_kpoints = bzS.Monkhorst_Pack(Nxyz, Nxyz, Nxyz)
    list_k_plus_q = np.array([k + q_vect for k in list_kpoints])
    basis_vectors = EPM_BASIS
    basis_size = len(basis_vectors)
    sum_total = 0.0
    list_eps = []
    list_k_contrib = []
    # print(len(list_kpoints))
    for index_k in range(len(list_kpoints)):
        # print(f"\r{index_k} values over {Nk}", end="")
        k = list_kpoints[index_k]
        k_plus_q = list_k_plus_q[index_k]
        H_k = hamiltonian(k, basis_vectors)
        H_k_plus_q = hamiltonian(k_plus_q, basis_vectors)
        eigval_k, eigvec_k = la.eigh(H_k)
        eigval_k_plus_q, eigvec_k_plus_q = la.eigh(H_k_plus_q)
        sum_k = 0.0
        for idx_cond in range(first_condcution_idx, first_condcution_idx + n_conduction):
            for idx_val in range(n_valence):
                mat_elem = np.abs(np.dot(eigvec_k[:, idx_val].T, eigvec_k_plus_q[:, idx_cond]))**2
                factor = 1.0 / (eigval_k_plus_q[idx_cond] - eigval_k[idx_val] - energy) + 1.0 / (eigval_k_plus_q[idx_cond] - eigval_k[idx_val] + energy)
                mat_elem *= factor
                sum_k += mat_elem
        list_k_contrib.append(sum_k)
        total_contrib = sum(list_k_contrib) / len(list_k_contrib)
        epsilon_tot = 1.0 + ((4.0*np.pi) / (np.linalg.norm(q_vect)**2)) * total_contrib
        list_eps.append(epsilon_tot)
    return list_eps[-1]

def convergence_Monkhorst_Pack(energy, q_vect, n_valence, n_conduction, maxNxyz):
    list_eps = []
    list_Nxyz = []
    for Nxyz in range(4, maxNxyz):
        print(f"Nxyz = {Nxyz}")
        list_eps.append(dielectric_function_Monkhorst_Pack(energy, q_vect, n_valence, n_conduction, Nxyz))
        list_Nxyz.append(Nxyz)
    fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    axs.plot(list_Nxyz, list_eps, '-o', color='blue')
    plt.show()
        
def main_epsilon_mc(N_k, energy, q_vect, n_valence, n_conduction):
    list_energies = np.linspace(0.0, 3.5, 2, endpoint=True)
    list_eps = []
    for e in list_energies:
        print(f"\r{e} values over {list_energies[-1]}", end="")
        eps = dielectric_function_mc(e, q_vect, n_valence, n_conduction, N_k)
        list_eps.append(eps)
    fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    axs.plot(list_energies, list_eps, '-o', color='blue')
    axs.set_xlabel('Energy (eV)')
    axs.set_ylabel('$\epsilon_r$')
    plt.show()
    
    for idx, eps in enumerate(list_eps):
        print(f"Energy = {list_energies[idx]} eV ----> epsilon_real = {eps:.2e}")
        
def main_epsilon(Nxyz, q_vect, n_valence, n_conduction):
    list_energies = np.linspace(0.0, 10, 41, endpoint=True)
    list_eps = []
    for e in list_energies:
        eps = dielectric_function_Monkhorst_Pack(e, q_vect, n_valence, n_conduction, Nxyz)
        list_eps.append(eps)
        print(f"E = {e} eV ----> epsilon_real = {eps:.2e}")
    fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    axs.plot(list_energies, list_eps, '-o', color='blue')
    axs.set_xlabel('Energy (eV)')
    axs.set_ylabel('$\epsilon_r$')
    plt.show()
    
    for idx, eps in enumerate(list_eps):
        print(f"Energy = {list_energies[idx]} eV ----> epsilon_real = {eps:.2e}")
    
    
def main_band_structure():
    N_basis = 20
    basis_G = generate_basis_vector(N_basis)
    n = 50
    k = k_path(n)
    bands = band_structure(k)
    bands -= max(bands[3])

    plt.figure(figsize=(15, 9))

    ax = plt.subplot(111)

    # remove plot borders
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # limit plot area to data
    plt.xlim(0, len(bands))
    plt.ylim(min(bands[0]) - 1, max(bands[7]) + 1)

    # custom tick names for k-points
    xticks = n * np.array([0, 0.5, 1, 1.5, 2, 2.25, 2.75, 3.25])
    plt.xticks(xticks, ('$L$', '$\Lambda$', '$\Gamma$', '$\Delta$', '$X$', '$U,K$', '$\Sigma$', '$\Gamma$'), fontsize=18)
    plt.yticks(fontsize=18)

    # horizontal guide lines every 2.5 eV
    for y in np.arange(-25, 25, 2.5):
        plt.axhline(y, ls='--', lw=0.3, color='black', alpha=0.3)

    # hide ticks, unnecessary with gridlines
    plt.tick_params(axis='both', which='both',
                    top='off', bottom='off', left='off', right='off',
                    labelbottom='on', labelleft='on', pad=5)

    plt.xlabel('k-Path', fontsize=20)
    plt.ylabel('E(k) (eV)', fontsize=20)

    plt.text(135, -18, 'Fig. 1. Band structure of Si.', fontsize=12)

    # tableau 10 in fractional (r, g, b)
    colors = 1 / 255 * np.array([
        [31, 119, 180],
        [255, 127, 14],
        [44, 160, 44],
        [214, 39, 40],
        [148, 103, 189],
        [140, 86, 75],
        [227, 119, 194],
        [127, 127, 127],
        [188, 189, 34],
        [23, 190, 207]
    ])

    for band, color in zip(bands, colors):
        plt.plot(band, lw=2.0, color=color)

    plt.show()

    
if __name__ == "__main__":
    energy = 0.0
    q_vect = np.array([1.0, 1.0, 1.0]) * 1.0e-6
    n_valence = 4
    n_conduction = 8
    Nk = 500
    # main_band_structure()
    # main_epsilon(N_k=Nk, energy=energy, q_vect=q_vect, n_valence=n_valence, n_conduction=n_conduction)
    # convergence_Monkhorst_Pack(energy, q_vect, n_valence, n_conduction, 25)
    # eps_mc = dielectric_function_mc(energy, q_vect, n_valence, n_conduction, Nk)
    # print(f"epsilon_real MC = {eps_mc:.2e}")
    
    # Nxyz = 40
    # eps = dielectric_function_Monkhorst_Pack(energy, q_vect, n_valence, n_conduction, Nxyz)
    # print(f"epsilon_real = {eps:.2e}")
    
    main_epsilon(Nxyz=20, q_vect=q_vect, n_valence=n_valence, n_conduction=n_conduction)
    
    