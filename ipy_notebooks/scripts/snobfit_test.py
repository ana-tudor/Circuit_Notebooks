import numpy as np
import pylab
from obj_function import energy_objective, uccsd_trotter_engine, get_fermion_operator, jordan_wigner, FermionOperator, uccsd_singlet_paramsize, MolecularData
from snobfit.snobfit import Snobfit

# Load only 0 and pi/2
molecule_filenames = ['data/H4-C2_DZP_singlet_ethylene_0.00000.hdf5','data/H4-C2_DZP_singlet_ethylene_1.57080.hdf5']

results = {}
results['0.00000rad Iterations'] = []
results['0.00000rad Initial Energies'] = []
results['0.00000rad Initial Amplitudes'] = []
results['0.00000rad Optimal Energies'] = []
results['0.00000rad Optimal Amplitudes'] = []
results['1.57080rad Iterations'] = []
results['1.57080rad Initial Energies'] = []
results['1.57080rad Initial Amplitudes'] = []
results['1.57080rad Optimal Energies'] = []
results['1.57080rad Optimal Amplitudes'] = []

uccsd_energies = {}
uccsd_amplitudes = {}
original_electrons = 16

N = 3

for iteration in range(0, N):
    for n_electrons, n_orbitals in [(2,2)]:

        occupied_indices = range(8 - n_electrons // 2)
        active_indices = range(8 - n_electrons // 2, 8 - n_electrons // 2 + n_orbitals)
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
        current_amplitudes = [0.0] * n_amplitudes + 0.001 * np.random.randn(n_amplitudes)
        for file_index, filename in enumerate(molecule_filenames):
            angle = filename.split('_')[-1].rstrip('.hdf5')
            iters = 0

#             print("Running CAS({},{}) on {}".format(n_electrons, n_orbitals, filename))
            molecule = MolecularData(filename=filename)
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

            initial_energy = energy_objective(current_amplitudes)

            # Run SNOBFIT optimization to find new CCSD parameters)
            snob = Snobfit(energy_objective, current_amplitudes,
                [np.array([-5, 5]), np.array([-5, 10])],
                dx=np.array([0.00005,0.00005]), maxiter=30, maxfun=20)
            [opt_amplitudes, opt_energy, opt_result] = snob.solve()
            print(opt_result)
            print(opt_energy)
            print(opt_amplitudes)


            results['{}rad Initial Amplitudes'.format(angle)].append(current_amplitudes)

            # Use previous iteration as guess for next iteration
            current_amplitudes = opt_amplitudes[:]

            #Store for file dumping
            results['{}rad Iterations'.format(angle)].append(iters)
            results['{}rad Initial Energies'.format(angle)].append(initial_energy)
            results['{}rad Optimal Energies'.format(angle)].append(opt_energy)
            results['{}rad Optimal Amplitudes'.format(angle)].append(opt_amplitudes)

#             print("\nInitial Energy: {}".format(initial_energy))
#             print("Optimal UCCSD Singlet Energy: {}".format(opt_energy))
#             print("Optimal UCCSD Singlet Amplitudes: {}".format(opt_amplitudes))
#             print("Exact FCI Energy: {} Hartrees".format(
#                     cas_energies['CAS({},{})'.format(n_electrons, n_orbitals)][file_index]))

print(uccsd_energies)
print(uccsd_amplitudes)

# Set up file for storing repeated results
import pickle
name = "data/ethylene_2_2_ethylene4_results_py"

f = open(name, "wb")
pickle.dump(results, f, -1)
f.close()

f2 = open(name, "rb")
b = pickle.load(f2)
print(b)
f2.close()
