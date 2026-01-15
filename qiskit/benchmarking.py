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

if __name__ == "__main__":
    runner = pyperf.Runner(add_cmdline_args=add_cmdline_args)
    runner.argparser.add_argument("--file", dest="filename", required=True, help="Path to the dataset file")
    args = runner.parse_args()

    seed, qubit_num, operations = read_dataset(args.filename)

    qc = create_circuit(seed, qubit_num, operations)

    backend = AerSimulator(method='statevector')

    qc = transpile(qc, backend)

    runner.bench_func(
        args.filename,      
        run_simulation,    
        backend, qc, seed  
    )





    