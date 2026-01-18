import numpy as np
import pyperf
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

def read_dataset(filename):
    file = open(filename)
    data = file.read()
    file.close()
    data = data.replace('i32', '').replace('i64','')
    lines = data.splitlines()
    seed = int(lines[0])
    qubit_num = int(lines[1])
    opcodes = eval(lines[2])
    qubits1 = eval(lines[3])
    qubits2 = eval(lines[4])
    operations = list(zip(opcodes, qubits1, qubits2))
    return seed, qubit_num, operations

def create_circuit(seed, qubit_num, operations):
    qc = QuantumCircuit(qubit_num, qubit_num)
    for opcode, qubit1, qubit2 in operations:
        if opcode == 0:
            qc.measure(qubit1, qubit1)
        elif opcode == 1:
            qc.h(qubit1)
        elif opcode == 2:
            qc.s(qubit1)
        elif opcode == 3:
            qc.cx(qubit1, qubit2)
    qc.save_statevector()
    return qc

def run_simulation(backend, qc, seed):
    job = backend.run(qc, seed_simulator=seed)
    result = job.result()
    final_vector = result.get_statevector(0)
    return final_vector

def add_cmdline_args(cmd, args):
    cmd.extend(('--file', args.filename))

def assert_state(name, actual_vector, expected_vector):
    # Compute the Euclidean distance (norm of difference)
    diff = np.linalg.norm(actual_vector - expected_vector)
    if diff < 1e-10:
        print(f"[PASS] {name}")
    else:
        print(f"[FAIL] {name}")
        print(f"  Diff: {diff}")
        print(f"  Actual Indices:   {np.where(np.abs(actual_vector) > 1e-5)}")
        print(f"  Expected Indices: {np.where(np.abs(expected_vector) > 1e-5)}")

def verify_implementation():    
    backend = AerSimulator(method='statevector')

    # Test 1: Hadamard (Opcode 1) |0> -> |+>
    ops_h = [(1, 0, 0)] 
    qc_h = create_circuit(0, 1, ops_h)
    
    vec_h = run_simulation(backend, qc_h, 0)
    
    # Expected: 1/sqrt(2) * (|0> + |1>)
    expected_h = np.zeros(2, dtype=complex)
    expected_h[0] = 1/np.sqrt(2)
    expected_h[1] = 1/np.sqrt(2)
    
    assert_state("Hadamard (Opcode 1)", vec_h, expected_h)


    # Test 2: Phase/S-Gate (Opcode 2) |1> -> i|1>
    ops_s = [(2, 0, 0)]
    qc_s_base = create_circuit(0, 1, ops_s)

    qc_s = QuantumCircuit(1, 1)
    qc_s.x(0)
    qc_s = qc_s.compose(qc_s_base)
    
    vec_s = run_simulation(backend, qc_s, 0)
    
    # Expected: i|1> (Index 1 should be i)
    expected_s = np.zeros(2, dtype=complex)
    expected_s[1] = 1j
    
    assert_state("Phase/S (Opcode 2)", vec_s, expected_s)

    # Test 3: CNOT Forward (Opcode 3) |10> -> |11>
    ops_cx = [(3, 0, 1)]
    qc_cx_base = create_circuit(0, 2, ops_cx)
    
    # Initialize q0=1
    qc_cx = QuantumCircuit(2, 2)
    qc_cx.x(0) 
    qc_cx = qc_cx.compose(qc_cx_base)
    
    vec_cx = run_simulation(backend, qc_cx, 0)
    
    expected_cx = np.zeros(4, dtype=complex)
    expected_cx[3] = 1.0
    
    assert_state("CNOT Forward (c=0, t=1)", vec_cx, expected_cx)

    # Test 4: CNOT Reverse (Opcode 3) |01> -> |11> (Logic: q0=0, q1=1 -> q0=1, q1=1)
    ops_rcx = [(3, 1, 0)]
    qc_rcx_base = create_circuit(0, 2, ops_rcx)
    
    qc_rcx = QuantumCircuit(2, 2)
    qc_rcx.x(1)
    qc_rcx = qc_rcx.compose(qc_rcx_base)
    
    vec_rcx = run_simulation(backend, qc_rcx, 0)
    
    expected_rcx = np.zeros(4, dtype=complex)
    expected_rcx[3] = 1.0
    
    assert_state("CNOT Reverse (c=1, t=0)", vec_rcx, expected_rcx)

    print("Testing Complete.")

if __name__ == "__main__":
    runner = pyperf.Runner(add_cmdline_args=add_cmdline_args)
    runner.argparser.add_argument("--file", dest="filename", required=False, help="Path to the dataset file")

    args = runner.parse_args()

    if not args.filename:
        print("Running tests")
        verify_implementation()
    
    else:

        seed, qubit_num, operations = read_dataset(args.filename)

        qc = create_circuit(seed, qubit_num, operations)
        
        backend = AerSimulator()
        if 'GPU' in backend.available_devices():
            print("Running GPU")
            backend = AerSimulator(method='statevector', device='GPU')
        else:
            print("Running CPU")
            backend = AerSimulator(method='statevector', device='CPU')

        qc = transpile(qc, backend)

        runner.bench_func(
            args.filename,      
            run_simulation,    
            backend, qc, seed  
        )





    