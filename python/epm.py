import numpy as np
import matplotlib.pyplot as plt
from scipy import constants



def pseudopotential(G):
    pass


def hamiltonian_diagonal(k, G):
    # h_bar ^ 2 / 2m
    constant_factor = constants.hbar**2 / (2 * constants.m_e * constants.e)
    return constant_factor * np.dot(k + G, k + G)



def hamiltonian_off_diagonal(k, GA, GB, tau):
    G = GA - GB
    return pseudopotential(G) * np.cos(np.dot(GA - GB, tau))