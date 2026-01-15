module Main where

import HQP
import qualified HQP.QOp.MatrixSemantics as MS
import Criterion.Main
import System.Environment (getArgs, withArgs)
import Numeric.LinearAlgebra (atIndex)
import System.Random (mkStdGen, randoms)
import Data.Complex (Complex)
import Control.DeepSeq (NFData, rnf)

instance NFData Step where
    rnf x = seq x ()

-- Wrapper to prevent DeepSeq from traversing the infinite list
newtype InfiniteRNG = InfiniteRNG [Double]

-- Instance that does NOTHING when asked to evaluate deeply
instance NFData InfiniteRNG where
    rnf (InfiniteRNG x) = seq x ()

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

removeSuffixes :: String -> String
removeSuffixes [] = []
removeSuffixes ('i':'3':'2':xs) = removeSuffixes xs
removeSuffixes ('i':'6':'4':xs) = removeSuffixes xs
removeSuffixes (x:xs) = x : removeSuffixes xs


readDataset :: String -> (Int, Int, [OperationInput])
readDataset file =
    let filteredFile = removeSuffixes file
        (line1 : line2 : line3 : line4 : line5 :_) = lines filteredFile
        seed = read line1 :: Int
        qubit_num = read line2 :: Int
        opcodes = read line3 :: [Int]
        qubit1s = read line4 :: [Int]
        qubit2s = read line5 :: [Int]
        operations = zip3 opcodes qubit1s qubit2s
    in (seed, qubit_num, operations)


runSimulation :: Program -> MS.StateT -> RNG -> Complex Double
runSimulation prog psi rng = 
    let (end_state, _, _) = MS.evalProg prog psi rng
    in end_state `atIndex` (0,0) -- Force evaluation here


setupEnvironment :: FilePath -> IO (Program, MS.StateT, InfiniteRNG)
setupEnvironment filename = do
    file <- readFile filename
    let (seed, qubit_num, operations) = readDataset file
        program = interpretOperations qubit_num operations
        starting_state = MS.ket (replicate qubit_num 0)
        rng  = randoms (mkStdGen seed) :: [Double]
    return (program, starting_state, InfiniteRNG rng)

main :: IO ()
main = do
    args <- getArgs
    case args of
        [] -> putStrLn "Usage: cabal run benchmark-criterion -- <filename>"
        (filename : critArgs) -> 
            withArgs critArgs $ defaultMain [
                bgroup "simulation" [
                    env (setupEnvironment filename) $ \ ~(prog, psi, InfiniteRNG rng) ->
                        bench filename $ whnf (\(p, s, r) -> runSimulation p s r) (prog, psi, rng)
                ]
            ]







