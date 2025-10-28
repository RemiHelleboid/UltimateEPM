# Partially from https://gist.github.com/wmedlar/

import scipy.constants as cst
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib.lines import Line2D
import matplotlib as mpl
import os, sys

import ParametersMat as pMat

try:
    plt.style.use(['science', 'high-vis', 'grid'])
except Exception:
    None
# plt.style.use(['seaborn-paper'])

mpl.rcParams['figure.figsize'] = [3.5, 2.8]


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
    print(f"Number of basis vectors: {len(basis_G)}")
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
            

    
def PseudoPotential(G: np.ndarray, mat_params: pMat.MaterialParameterEPM) -> float:
    Ry_to_eV = 13.605698066
    V3S = mat_params.V3S * Ry_to_eV
    V8S = mat_params.V8S * Ry_to_eV
    V11S = mat_params.V11S * Ry_to_eV
    V3A = mat_params.V3A * Ry_to_eV
    V8A = mat_params.V4A * Ry_to_eV
    V11A = mat_params.V11A * Ry_to_eV
    tau = np.array([1, 1, 1]) / 8.0
    Gtau = 2.0 * np.pi * np.dot(G, tau)
    Gsquare = np.linalg.norm(G)**2

    if np.abs(Gsquare - 3) < 1e-6:
        return V3S * np.cos(Gtau) + V3A * np.sin(Gtau)
    elif np.abs(Gsquare - 8) < 1e-6:
        return V8S * np.cos(Gtau) + V8A * np.sin(Gtau)
    elif np.abs(Gsquare - 11) < 1e-6:
        return V11S * np.cos(Gtau) + V11A * np.sin(Gtau)
    else:
        return 0.0


def non_diag_hamiltonian(mat_params: pMat.MaterialParameterEPM, basis_G: np.ndarray) -> np.ndarray:
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
    nnz = 0
    for i in range(BasisSize):
        for j in range(i+1):
            G1 = basis_G[i]
            G2 = basis_G[j]
            Gdiff = G1 - G2
            Hamiltonian[i, j] = PseudoPotential(Gdiff, mat_params)
            if Hamiltonian[i, j] != 0:
                nnz += 1
    # print("Ratio of non zero elements in the Hamiltonian matrix: ", nnz / (BasisSize**2))
    return Hamiltonian


EPM_BASIS = generate_basis_vector(10)
EPM_BASIS_NORM = np.linalg.norm(EPM_BASIS, axis=1)


def hamiltonian(mat_params: pMat.MaterialParameterEPM, k: np.ndarray, basis_G: np.ndarray) -> np.ndarray:
    """Generate the Hamiltonian matrix.

    Args:
        k (np.ndarray): The wave vector.
        basis_G (np.ndarray): The basis vector.

    Returns:
        np.ndarray: The Hamiltonian matrix.
    """
    a_0 = mat_params.lattice_constant * 1e-10
    BasisSize = len(basis_G)
    kinetic_c = (2.0 * np.pi / a_0)**2
    Hamiltonian = non_diag_hamiltonian(mat_params, basis_G)
    for i in range(BasisSize):
        G1 = basis_G[i]
        Hamiltonian[i, i] = KINETIC_CONSTANT * kinetic_c * np.linalg.norm(G1 + k)**2
    return Hamiltonian


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


def band_structure(mat_params: pMat.MaterialParameterEPM, path):
    bands = []
    # vstack concatenates our list of paths into one nice array
    for idx, k in enumerate(np.vstack(path)):
        print(f"\r{idx} values over {len(np.vstack(path))}", end="")
        E, _ = eigen_states(mat_params, k, EPM_BASIS)
        # picks out the lowest eigt eigenvalues
        bands.append(E[:8])

    
    return np.stack(bands, axis=-1)


def eigen_states(mat_params: pMat.MaterialParameterEPM, k: np.ndarray, basis_G: np.ndarray=None) -> np.ndarray:
    """Generate the eigen states.

    Args:
        k (np.ndarray): The wave vector.
        basis_G (np.ndarray): The basis vector.

    Returns:
        np.ndarray: The eigen states.
    """
    if basis_G is None:
        basis_G = EPM_BASIS
    HamiltonianMatrix = hamiltonian(mat_params, k, basis_G)
    # Get eigenvalues and eigenvectors of the Hamiltonian matrix
    eigen_values, eigen_vectors = la.eigh(HamiltonianMatrix)
    # Sort the eigenvalues and eigenvectors
    idx = eigen_values.argsort()
    eigen_values = eigen_values[idx]
    eigen_vectors = eigen_vectors[:, idx]
    return eigen_values, eigen_vectors

def plot_eigen_states(k: np.ndarray, basis_G: np.ndarray=None, n_states: int=8) -> None:
    """Plot the eigen states.

    Args:
        k (np.ndarray): The wave vector.
        basis_G (np.ndarray): The basis vector.
        n_states (int): The number of states to plot.
    """
    if basis_G is None:
        basis_G = EPM_BASIS
    eigen_values, eigen_vectors = eigen_states(k, basis_G)
    norms = np.linalg.norm(eigen_vectors, axis=0)
    N = len(basis_G)
    Ns = 10
    # Plot the fourier factors colors is shade of blue 
    cmap = plt.cm.get_cmap('Blues')
    colors = [cmap(i) for i in np.linspace(0, 1, Ns)]
    
    fig, ax = plt.subplots(nrows=1, ncols=Ns, figsize=(40, 5), sharey=True)
    for ix in range(Ns):
        i = ix + 0
        ax[ix].scatter(np.arange(N), eigen_vectors[:, i], label=f"{i}", color='b', s=4)
        # Bar plot of the fourier factors
        ax[ix].bar(np.arange(N), eigen_vectors[:, i], color='k', alpha=1, width=1)
        ax[ix].set_title(f"State {i}")
        ax[ix].set_xlim(-1.5, N+1.5)
        ax[ix].set_ylim(-1.5, 1.5)
    fig.tight_layout()

    # ax.set_title(f"Fourier factors for k = {k}")

    fig.tight_layout()
    fig.savefig(f"fourier_factors_k_{k}.svg")
    plt.show()

def get_wave_function_real_space(mat_params: pMat.MaterialParameterEPM, k: np.ndarray, r_min, r_max, n_points, basis_G: np.ndarray=None) -> np.ndarray:
    """Generate the wave function in real space.

    Args:
        k (np.ndarray): The wave vector.
        r_min (float): The minimum value of the real space.
        r_max (float): The maximum value of the real space.
        n_points (int): The number of points in the real space.
        basis_G (np.ndarray): The basis vector.

    Returns:
        np.ndarray: The wave function in real space.
    """
    if basis_G is None:
        basis_G = EPM_BASIS
    eigen_values, eigen_vectors = eigen_states(mat_params, k, basis_G)
    NVal = 4
    eigen_vectors = eigen_vectors[:, :NVal]
    print(f"Shape of eigen_vectors: {eigen_vectors.shape}")
    # Get the wave function in real space
    xyz = np.linspace(r_min, r_max, n_points)
    X, Y, Z = np.meshgrid(xyz, xyz, xyz)
    wave_function = np.zeros((n_points, n_points, n_points), dtype=complex)
    for idx_x in range(n_points):
        for idx_y in range(n_points):
            for idx_z in range(n_points):
                r_point = np.array([X[idx_x, idx_y, idx_z], Y[idx_x, idx_y, idx_z], Z[idx_x, idx_y, idx_z]])
                for i in range(NVal):
                    G_vect = basis_G[i]
                    for idx_G in range(len(G_vect)):
                        wave_function[idx_x, idx_y, idx_z] +=  np.exp(1j * np.dot(G_vect, r_point)) * eigen_vectors[i, idx_G]
                wave_function[idx_x, idx_y, idx_z] *= np.exp(1j * np.dot(k, r_point))
                # wave_function[idx_x, idx_y, idx_z] = np.abs(wave_function[idx_x, idx_y, idx_z])**2
    return wave_function



def plot_wave_function_real_space(mat_params: pMat.MaterialParameterEPM, k: np.ndarray, r_min, r_max, n_points, basis_G: np.ndarray=None):
    """Plot the wave function in real space.

    Args:
        k (np.ndarray): The wave vector.
        r_min (float): The minimum value of the real space.
        r_max (float): The maximum value of the real space.
        n_points (int): The number of points in the real space.
        basis_G (np.ndarray): The basis vector.
    """
    if basis_G is None:
        basis_G = EPM_BASIS
    wave_function = get_wave_function_real_space(mat_params, k, r_min, r_max, n_points, basis_G)
    # Export in a file with the form X,Y,Z,Re(WaveFunction),Im(WaveFunction),Abs(WaveFunction)
    xyz = np.linspace(r_min, r_max, n_points)
    X, Y, Z = np.meshgrid(xyz, xyz, xyz)
    with open("wave_function_real_space.csv", "w") as f:
        f.write("X,Y,Z,Re(WaveFunction),Im(WaveFunction),Abs(WaveFunction)\n")
        for idx_x in range(n_points):
            for idx_y in range(n_points):
                for idx_z in range(n_points):
                    f.write(f"{X[idx_x, idx_y, idx_z]},{Y[idx_x, idx_y, idx_z]},{Z[idx_x, idx_y, idx_z]},{wave_function[idx_x, idx_y, idx_z].real},{wave_function[idx_x, idx_y, idx_z].imag},{np.abs(wave_function[idx_x, idx_y, idx_z])}\n")
    


    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    im = ax.imshow(np.abs(wave_function[:, :, 0])**2, cmap='jet', interpolation='gaussian', extent=[r_min, r_max, r_min, r_max], origin='lower')
    # # Isocontours
    cf = ax.contour(np.abs(wave_function[:, :, 0])**2, c='k', extent=[r_min, r_max, r_min, r_max], levels=10)
    ax.clabel(cf, inline=1, fontsize=10)
    

    ax.set_title(f"z=0")
    fig.colorbar(im, ax=ax)

    fig.tight_layout()

    plt.pause(4)
    plt.show()




def dielectric_function_mc(mat_params: pMat.MaterialParameterEPM, energy, q_vect, n_valence, n_conduction, Nk):
    list_colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    first_condcution_idx = 4
    list_kpoints = bzS.RandomKPointsIrreducibleWedge(Nk)
    list_k_plus_q = np.array([k + q_vect for k in list_kpoints])
    basis_vectors = EPM_BASIS
    basis_size = len(basis_vectors)
    sum_total = 0.0
    list_eps = []
    list_k_contrib = []
    for index_k in range(len(list_kpoints)):
        # print(f"\r{index_k} values over {Nk}", end="")
        k = list_kpoints[index_k]
        k_plus_q = list_k_plus_q[index_k]
        H_k = hamiltonian(mat_params, k, basis_vectors)
        H_k_plus_q = hamiltonian(mat_params, k_plus_q, basis_vectors)
        # print(f"Norm Sq: {np.linalg.norm(H_k)}")
        eigval_k, eigvec_k = la.eigh(H_k)
        eigval_k_plus_q, eigvec_k_plus_q = la.eigh(H_k_plus_q)
        sum_k = 0.0
        for idx_cond in range(first_condcution_idx, first_condcution_idx + n_conduction):
            for idx_val in range(n_valence):
                mat_elem = np.abs(np.dot(eigvec_k[:, idx_val], eigvec_k_plus_q[:, idx_cond]))**2
                print("Mean k state" , np.mean(np.abs(eigvec_k[:, idx_val])**2))
                # print("Mean k+q state" , np.linalg.norm(eigvec_k_plus_q[:, idx_cond]))
                # print(f"mat_elem = {mat_elem:.2e}")
                factor = 1.0 / (eigval_k_plus_q[idx_cond] - eigval_k[idx_val] - energy) + 1.0 / (eigval_k_plus_q[idx_cond] - eigval_k[idx_val] + energy)
                mat_elem *= factor
                # print("Factor = ", factor)
                sum_k += mat_elem
        list_k_contrib.append(sum_k)
        total_contrib = sum(list_k_contrib) / len(list_k_contrib)
        epsilon_tot = 1.0 + ((4.0*np.pi) / (np.linalg.norm(q_vect)**2)) * total_contrib
        list_eps.append(epsilon_tot)
    # fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    #     print(f"\r{index_k} values over {Nk} ----> epsilon_real = {list_eps[-1]:.2e}", end="", flush=True)
    # axs.plot(list_eps, '-o', color='blue')
    # plt.show()
    return list_eps[-1]

def dielectric_function_Monkhorst_Pack(mat_params: pMat.MaterialParameterEPM, energy, q_vect, n_valence, n_conduction, Nxyz, basis_G=None):
    first_condcution_idx = 4
    irreducible_wedge =False
    list_kpoints = bzS.Monkhorst_Pack(Nxyz, Nxyz, Nxyz, irreducible_wedge)
    list_k_plus_q = np.array([k + q_vect for k in list_kpoints])
    basis_vectors = EPM_BASIS
    basis_size = len(basis_vectors)
    if basis_G is None:
        basis_G = EPM_BASIS
    sum_total = 0.0
    list_eps = []
    list_k_contrib = []
    # print(len(list_kpoints))
    for index_k in range(len(list_kpoints)):
        # print(f"\r{index_k} values over {Nk}", end="")
        k = list_kpoints[index_k]
        k_plus_q = list_k_plus_q[index_k]
        H_k = hamiltonian(mat_params, k, basis_G)
        H_k_plus_q = hamiltonian(mat_params, k_plus_q, basis_G)
        eigval_k, eigvec_k = la.eigh(H_k)
        eigval_k_plus_q, eigvec_k_plus_q = la.eigh(H_k_plus_q)
        sum_k = 0.0
        for idx_cond in range(first_condcution_idx, first_condcution_idx + n_conduction):
            for idx_val in range(n_valence):
                mat_elem = np.abs(np.dot(eigvec_k[:, idx_val].T, eigvec_k_plus_q[:, idx_cond]))**2
                # print(f"mat_elem = {mat_elem:.2e}")
                factor = 1.0 / (eigval_k_plus_q[idx_cond] - eigval_k[idx_val] - energy) + 1.0 / (eigval_k_plus_q[idx_cond] - eigval_k[idx_val] + energy)
                mat_elem *= factor
                # print("Factor = ", factor)
                sum_k += mat_elem
        list_k_contrib.append(sum_k)
        total_contrib = sum(list_k_contrib) / len(list_k_contrib)
        epsilon_tot = 1.0 + ((4.0*np.pi) / (np.linalg.norm(q_vect)**2)) * total_contrib
        list_eps.append(epsilon_tot)
    return list_eps[-1]

def convergence_Monkhorst_Pack(mat_params: pMat.MaterialParameterEPM, energy, q_vect, n_valence, n_conduction, maxNxyz):
    list_eps = []
    list_Nxyz = []
    for Nxyz in range(4, maxNxyz):
        print(f"Nxyz = {Nxyz}")
        list_eps.append(dielectric_function_Monkhorst_Pack(energy, q_vect, n_valence, n_conduction, Nxyz))
        list_Nxyz.append(Nxyz)
    fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    axs.plot(list_Nxyz, list_eps, '-o', color='blue')
    plt.show()
        
def main_epsilon_mc(N_k, q_vect, n_valence, n_conduction):
    list_energies = np.linspace(0.0, 10, 11, endpoint=True)
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
        
def main_epsilon(mat_params: pMat.MaterialParameterEPM, Nxyz, q_vect, n_valence, n_conduction):
    print("Calculating epsilon...")
    list_energies = np.linspace(0.0, 10, 41, endpoint=True)
    list_eps = []
    for e in list_energies:
        eps = dielectric_function_Monkhorst_Pack(mat_params, e, q_vect, n_valence, n_conduction, Nxyz)
        list_eps.append(eps)
        print(f"E = {e} eV ----> epsilon_real = {eps:.2e}")
    fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    axs.plot(list_energies, list_eps, '-o', color='blue')
    axs.set_xlabel('Energy (eV)')
    axs.set_ylabel('$\epsilon_r$')
    plt.show()
    
    for idx, eps in enumerate(list_eps):
        print(f"Energy = {list_energies[idx]} eV ----> epsilon_real = {eps:.2e}")
    
    
def main_band_structure(mat_params: pMat.MaterialParameterEPM, n_points):
    N_basis = 20
    basis_G = generate_basis_vector(N_basis)
    n = n_points
    k = k_path(n)
    bands = band_structure(mat_params, k)
    bands -= max(bands[3])

    fig, ax = plt.subplots(figsize=(10, 10))

    # limit plot area to data
    plt.xlim(0, len(bands))
    plt.ylim(min(bands[0]) - 1, max(bands[7]) + 1)

    # custom tick names for k-points
    xticks = n * np.array([0, 0.5, 1, 1.5, 2, 2.25, 2.75, 3.25])
    plt.xticks(xticks, ('$L$', '$\Lambda$', '$\Gamma$', '$\Delta$', '$X$', '$U,K$', '$\Sigma$', '$\Gamma$'), fontsize=18)
    plt.yticks(fontsize=18)

    for band in bands:
        ax.plot(band, lw=2.0)
        print(f"Band min = {min(band):.2f} eV")

    plt.show()

    
if __name__ == "__main__":
    material = sys.argv[1]
    # Get the material parameters
    MyMaterial = pMat.get_material_by_symbol(material)

    energy = 0.0
    q_vect = np.array([1.0, 1.0, 1.0]) * 1.0e-9
    n_valence = 4
    n_conduction = 8
    Nk = 1000
    main_band_structure(MyMaterial, 150)
    k_gamma = np.array([1.0, 1.0, 0.0]) * 1.0e-12
    rmin = -1.0 * np.pi
    rmax = 1.0 * np.pi
    # plot_eigen_states(k_gamma)
    # plot_wave_function_real_space(MyMaterial, k_gamma, rmin, rmax, Nxyz)

    # Nxyz = 20
    # main_epsilon(MyMaterial, Nxyz, q_vect, n_valence, n_conduction)
    
    