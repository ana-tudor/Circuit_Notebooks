import itertools
import os
import time

import numpy
import matplotlib.pyplot as plt
import sys

import cirq
import openfermion
import openfermioncirq as ofc
from openfermioncirq.optimization import OptimizationParams, ScipyOptimizationAlgorithm
from openfermionpyscf import run_pyscf



def h2_geometry(bond_length):
    return [
        ('H', (0.0, 0.0, 0.0)),
        ('H', (0.0, 0.0, bond_length))
    ]


def n2_geometry(bond_length):
    return [
        ('N', (0.0, 0.0, 0.0)),
        ('N', (0.0, 0.0, bond_length))
    ]


def h2o_geometry(bond_length):
    bond_angle = numpy.pi * 104.5 / 180  # 104.5 degrees
    a = bond_length * numpy.cos(bond_angle / 2.0)
    b = bond_length * numpy.sin(bond_angle / 2.0)
    return [
        ('H', (a, b, 0.0)),
        ('H', (a, -b, 0.0)),
        ('O', (0.0, 0.0, 0.0))
    ]


def generate_molecular_hamiltonian(geometry,
                                   n_active_electrons,
                                   n_active_orbitals,
                                   basis='cc-pvdz',
                                   multiplicity=1):

    # Run electronic structure calculations
    molecule = run_pyscf(
        openfermion.hamiltonians.MolecularData(
            geometry,
            basis,
            multiplicity
        )
    )

    # Freeze core orbitals and truncate to active space
    n_core_orbitals = (molecule.n_electrons - n_active_electrons) // 2
    occupied_indices = list(range(n_core_orbitals))
    active_indices = list(range(n_core_orbitals,
                                n_core_orbitals + n_active_orbitals))

    return molecule.get_molecular_hamiltonian(occupied_indices, active_indices)

def cost_with_noise(noise, variance):
    if noise == 0.0:
        return None
    return variance/noise**2 * 2/numpy.pi

def eval_with_noise(study, noise, params):
    cost = cost_with_noise(noise, study.objective.variance_bound)
    return study.value_of(params) + study.objective.noise(cost)

# Create or load a variational study
# ----------------------------------

# Set directory to save studies in
STUDIES_DIR = 'studies'


# Set Hamiltonian parameters
geometry_factory = h2_geometry
bond_length = 1.4
n_active_electrons = 2
n_active_orbitals = 2
hamiltonian_name = 'H2_cc-pvdz_singlet_1.4_2-2'


# Set ansatz parameters
ansatz_class = ofc.SwapNetworkTrotterAnsatz
iterations = 1
ansatz_kwargs = {'include_all_xxyy': True}


# Generate Hamiltonian
hamiltonian = generate_molecular_hamiltonian(
    geometry_factory(bond_length),
    n_active_electrons,
    n_active_orbitals
)


# Create or load study
study_name = '{}_{}_iterations{}'.format(
    hamiltonian_name, ansatz_class.__name__, iterations)

if os.path.isfile(os.path.join(STUDIES_DIR, '{}.study'.format(study_name))):
    # Load study
    study = ofc.VariationalStudy.load(
        study_name,
        datadir=STUDIES_DIR)
    print("LOADED a variational study with {} qubits and {} parameters.".format(
        len(study.ansatz.qubits), study.num_params))
else:
    # Create study
    # Generate ansatz and objective
    hamiltonian_ferm_op = openfermion.get_fermion_operator(hamiltonian)
    ansatz_hamiltonian = openfermion.get_diagonal_coulomb_hamiltonian(
        hamiltonian_ferm_op,
        ignore_incompatible_terms=True)
    ansatz = ansatz_class(
        ansatz_hamiltonian,
        iterations=iterations,
        **ansatz_kwargs)
    objective = ofc.HamiltonianObjective(hamiltonian)

    # Use preparation circuit for mean-field state
    preparation_circuit = cirq.Circuit.from_ops(
        ofc.prepare_gaussian_state(
            ansatz.qubits,
            openfermion.QuadraticHamiltonian(ansatz_hamiltonian.one_body),
            occupied_orbitals=range(n_active_electrons)))

    study = ofc.VariationalStudy(
        study_name,
        ansatz,
        objective,
        preparation_circuit=preparation_circuit,
        datadir=STUDIES_DIR)
    print("CREATED a variational study with {} qubits and {} parameters.".format(
        len(study.ansatz.qubits), study.num_params))


hamiltonian_sparse = openfermion.get_sparse_operator(study.objective.hamiltonian)
energy, _ = openfermion.get_ground_state(hamiltonian_sparse)

print(energy)

# print(ansatz.param_bounds())

params = numpy.array([float(sys.argv[1]), float(sys.argv[2]), 
                    float(sys.argv[3]), float(sys.argv[4]),
                    float(sys.argv[5]), float(sys.argv[6]),
                    float(sys.argv[7]), float(sys.argv[8]),
                    float(sys.argv[9]), float(sys.argv[10]),
                    float(sys.argv[11]), float(sys.argv[12]),
                    float(sys.argv[13]), float(sys.argv[14]),
                    float(sys.argv[15]), float(sys.argv[16]),
                    ])

noise = float(sys.argv[17])

''' Run study using blackbox, evaluate at new x'''
print(params)
print(noise)
result=eval_with_noise(study,noise,params)
print(result)

study.save()

name = 'out_param_energies.txt'

with open(name, 'wb') as f:
    numpy.savetxt(f, [result, 0, 1])
