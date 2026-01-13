module Main where

import HQP
import qualified HQP.QOp.MatrixSemantics as MS
import System.Environment (getArgs)
import System.Random (mkStdGen, randoms, randomR)

-- An single Clifford operation is defined as a tuple of 3 integers: (OperationCode, Qubit Number, Qubit Number 2)
-- Opcodes are:
-- 0: Measurement, 1: Hadamard, 2: Phase, 3: CNOT
type OperationInput = (Int, Int, Int)

-- qubit_num starts indexing from 1
interpretOperation :: Int -> OperationInput -> Step
interpretOperation qubit_num (opcode, qubit1, qubit2)
    | opcode == 0 = Measure [qubit1]
    | opcode == 1 = Unitary (Id qubit1 ⊗ H ⊗ Id (qubit_num - 1 - qubit1))
    | opcode == 2 = Unitary (Id qubit1 ⊗ R Z (1/2) ⊗ Id (qubit_num - 1 - qubit1))
    | opcode == 3 = cnotHelper qubit_num qubit1 qubit2
    | otherwise = error "Invalid opcode"

cnotHelper :: Int -> Int -> Int -> Step
cnotHelper qubit_num c t
    | c == t = error "Cannot apply CNOT with both control and target qubit the same"
    | c < t = Unitary (Id c ⊗ C (Id (t-c-1) ⊗ X) ⊗ Id (qubit_num - 1 - t))
    | otherwise =
        Unitary (hh ∘ Id t ⊗ C (Id (c-t-1) ⊗ X) ⊗ Id (qubit_num-c-1) ∘ hh)
    where
        hh = hAt qubit_num c ∘ hAt qubit_num t

hAt :: Int -> Int -> QOp
hAt qubit_num qubit1 = Id qubit1 ⊗ H ⊗ Id (qubit_num-qubit1-1)

interpretOperations :: Int -> [OperationInput] -> Program
interpretOperations qubit_num ops = map (interpretOperation qubit_num) ops







