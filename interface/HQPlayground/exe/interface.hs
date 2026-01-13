module Main where

import Control.Monad (forM_)
import Data.Int (Int32)
import qualified Data.Vector.Storable as VS

import HQP.QOp.GPU.Compile (compileProgram, Instr(..))
import HQP.QOp.GPU.RunFuthark (withFuthark, runSimulate, GPUResult(..))
import Programs.RepeaterProtocol (teleport)

numQubitsFromInstrs :: [Instr] -> Int
numQubitsFromInstrs is =
  1 + maximum (0 : [ fromIntegral (a i) | i <- is ] ++ [ fromIntegral (b i) | i <- is ])

validate :: Int -> [Instr] -> IO ()
validate nQ instrs = do
  let badOpcode = [ i | i <- instrs, opcode i < 0 || opcode i > 3 ]
      badA      = [ i | i <- instrs, a i < 0 || a i >= fromIntegral nQ ]
      badB      = [ i | i <- instrs, opcode i == 3, b i < 0 || b i >= fromIntegral nQ ]
  if null badOpcode then pure () else error ("Bad opcode(s): " <> show badOpcode)
  if null badA      then pure () else error ("Bad qubit index in a: " <> show badA)
  if null badB      then pure () else error ("Bad qubit index in b (CNOT): " <> show badB)

runTest :: String -> Int -> [Instr] -> IO ()
runTest name nQ instrs = do
  putStrLn ""
  putStrLn ("=== " <> name <> " ===")
  withFuthark $ \ctx ->
    forM_ ([1..20] :: [Int32]) $ \seed -> do
      -- run the simulation and get both tableau and measurements
      GPUResult { gpuTableau = tab, gpuMeasurements = meas } <-
        runSimulate ctx seed nQ instrs
      -- print both tableau and measurements for debugging
      -- putStrLn (show seed <> " -> tableau: " <> show (VS.toList tab) <> ", measurements: " <> show (VS.toList meas))
      putStrLn (show seed <> " -> measurements: " <> show (VS.toList meas))

testDeterministicZero :: IO ()
testDeterministicZero =
  runTest "Deterministic |0> measurement"
    1
    [ Instr 0 0 0 ]

testHadamardRandom :: IO ()
testHadamardRandom =
  runTest "Single-qubit H then Measure (50/50)"
    1
    [ Instr 1 0 0
    , Instr 0 0 0
    ]

testBellPair :: IO ()
testBellPair =
  runTest "Bell pair correlation"
    2
    [ Instr 1 0 0
    , Instr 3 0 1
    , Instr 0 0 0
    , Instr 0 1 0
    ]

testIndependentQubits :: IO ()
testIndependentQubits =
  runTest "Two independent qubits"
    2
    [ Instr 1 0 0
    , Instr 1 1 0
    , Instr 0 0 0
    , Instr 0 1 0
    ]

-- testTeleport :: IO ()
-- testTeleport = do
--   putStrLn ""
--   putStrLn "=== Teleport program ==="

--   let prog = teleport 3 0 1 2
--   case compileProgram prog of
--     Left err -> print err
--     Right (_nQ_from_HQP, instrs) -> do
--       let nQ = numQubitsFromInstrs instrs
--       validate nQ instrs
--       putStrLn ("instr count = " <> show (length instrs))
--       withFuthark $ \ctx ->
--         forM_ ([2026..2045] :: [Int32]) $ \seed -> do
--           meas <- runSimulateMeasurements ctx seed nQ instrs
--           putStrLn (show seed <> " -> " <> show (VS.toList meas))

main :: IO ()
main = do
  testDeterministicZero
  testHadamardRandom
  testBellPair
  testIndependentQubits
