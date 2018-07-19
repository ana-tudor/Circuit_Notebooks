import glob
import h5py
import numpy as np
import os
import pylab
import scipy.linalg
import scipy as sp
import sys
from os.path import dirname
sys.path.append(dirname("C:\Users\hp\Documents\GitHub\Circuit_Notebooks\lib\ProjectQ-develop\projectq"))

from numpy import array, concatenate, zeros
from numpy.random import randn
from scipy.optimize import minimize

from functools import reduce

from openfermion.config import *
from openfermionprojectq import *

from openfermion.hamiltonians import MolecularData
from openfermion.transforms import get_fermion_operator, jordan_wigner, bravyi_kitaev, get_sparse_operator
from openfermion.utils import uccsd_singlet_paramsize
from openfermion.ops import FermionOperator

from openfermionpsi4 import run_psi4

from projectq.ops import X, All, Measure
from projectq.backends import CommandPrinter, CircuitDrawer
from projectq.ops import (BasicGate,
                          H,
                          X,
                          CNOT,
                          Measure,
                          Z,
                          Swap,
                          C,
                          Rx,
                          Rz)
from projectq.backends._circuits import to_latex
from projectq.backends import _printer


hartree_to_kcal = 627.50947415

original_electrons = 16
occupied_indices = range(8 - 2 // 2)
active_indices = range(8 - 2 // 2, 8 -2 // 2 + 2)
active_electrons = original_electrons - len(occupied_indices) * 2
active_qubits = len(active_indices) * 2

# Define a qubit number operator
number_operator = jordan_wigner(
    sum([FermionOperator( ((i, 1), (i, 0)), 1.0) for i in range(active_qubits)], FermionOperator()))
number_operator.compress()
spinz_operator = jordan_wigner(sum([FermionOperator( ((2*i, 1), (2*i, 0)), 0.5)
                                    for i in range(active_qubits // 2)], FermionOperator()) +
                               sum([FermionOperator( ((2*i+1, 1), (2*i+1, 0)), -0.5)
                                    for i in range(active_qubits // 2)], FermionOperator()))
spinz_operator.compress()

n_amplitudes = int(uccsd_singlet_paramsize(active_qubits,
                                       active_electrons))
#         print("Running CAS({},{}) with {} coupled cluster amplitudes".format(n_electrons, n_orbitals,n_amplitudes))
current_amplitudes = [0.0] * n_amplitudes + 0.001 * randn(n_amplitudes)



molecule = MolecularData(filename='C:\Users\hp\Documents\GitHub\Circuit_Notebooks\data\H4-C2_DZP_singlet_ethylene_1.57080.hdf5')
# print("Number of spatial basis functions: {}".format(molecule.n_orbitals))
# print("Number of electrons: {}".format(molecule.n_electrons))


# Extract active space integrals
hamiltonian = (molecule.
               get_molecular_hamiltonian(
                   occupied_indices=occupied_indices,
                   active_indices=active_indices))
fermion_hamiltonian = get_fermion_operator(hamiltonian)


# Use a Jordan-Wigner encoding, and compress to remove 0 imaginary components
qubit_hamiltonian = jordan_wigner(fermion_hamiltonian)

qubit_hamiltonian.compress()

# Set standard UCCSD Compiler engine
compiler_engine = uccsd_trotter_engine()

n_amplitudes = int(uccsd_singlet_paramsize(molecule.n_qubits, molecule.n_electrons))
#         print("Running CAS({},{}) with {} coupled cluster amplitudes".format(n_electrons, n_orbitals,n_amplitudes))
current_amplitudes = [0.0] * n_amplitudes + 0.001 * randn(n_amplitudes)


input_amplitudes = np.array([float(sys.argv[1]), float(sys.argv[2])])
print(input_amplitudes)

from projectq.ops.noise_traits import scale
print(scale)

def energy_objective(packed_amplitudes):

    """Evaluate the energy of a UCCSD singlet wavefunction with packed_amplitudes
    Args:
        packed_amplitudes(ndarray): Compact array that stores the unique
            amplitudes for a UCCSD singlet wavefunction.

    Returns:
        energy(float): Energy corresponding to the given amplitudes
    """

    #Variables for compatibility with imfil.m
    ifail = 0

    # Set Jordan-Wigner initial state with correct number of electrons
    wavefunction = compiler_engine.allocate_qureg(active_qubits)


    # Set some of the qubits to |1>, or occupied
    for i in range(active_electrons):
        X | wavefunction[i]

    # print(active_electrons)

    # Build the circuit and act it on the wavefunction
    evolution_operator = uccsd_singlet_evolution(packed_amplitudes,
                                                 active_qubits,
                                                 active_electrons)
    evolution_operator | wavefunction

    compiler_engine.flush()

    try:

        # Evaluate the energy and reset wavefunction
        energy = compiler_engine.backend.get_expectation_value(qubit_hamiltonian, wavefunction)
        number = compiler_engine.backend.get_expectation_value(number_operator, wavefunction)
        spinz = compiler_engine.backend.get_expectation_value(spinz_operator, wavefunction)
        # print("Energy: {}\t Number:{}\t Sz:{}".format(energy, number, spinz))
        All(Measure) | wavefunction
        compiler_engine.flush()

    except:
        ifail = 1
        energy = float('NaN')

    return (energy, ifail)



name = "C:\Users\hp\Documents\GitHub\Circuit_Notebooks\imfil\output.txt"
energy, ifail = energy_objective(input_amplitudes)
results = [energy,ifail,1]
with open(name, 'wb') as f:
    np.savetxt(f, results)
