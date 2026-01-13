module Main where

import Control.Monad (forM_)
import Data.Int (Int32)
import qualified Data.Vector.Storable as VS

import HQP.QOp.GPU.Compile (compileProgram, Instr(..))
import HQP.QOp.GPU.RunFuthark (withFuthark, runSimulate, GPUResult(..))
import Programs.RepeaterProtocol (teleport)

import Data.Int (Int8)

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

printTableauPrefix :: Int -> Int -> Int -> Int -> VS.Vector Int8 -> IO ()
printTableauPrefix rows cols maxR maxC tab = do
  let rLim = min rows maxR
      cLim = min cols maxC
      idx r c = tab VS.! (r * cols + c)
  forM_ [0..rLim-1] $ \r -> do
    let row = [ idx r c | c <- [0..cLim-1] ]
    putStrLn ("Row " <> show r <> ": " <> show row <> if cLim < cols then " ..." else "")

printTableau :: Int -> Int -> VS.Vector Int8 -> IO ()
printTableau rows cols tab = do
  let idx r c = tab VS.! (r * cols + c)
  printTableauPrefix rows cols rows cols tab
  forM_ [0..rows-1] $ \r -> do
    let row = [ idx r c | c <- [0..cols-1] ]
    putStrLn ("Row " <> show r <> ": " <> show row)

runTest :: String -> Int -> [Instr] -> IO ()
runTest name nQ instrs = do
  putStrLn ""
  putStrLn ("=== " <> name <> " ===")
  withFuthark $ \ctx ->
    forM_ ([1..20] :: [Int32]) $ \seed -> do
      res <- runSimulate ctx seed nQ instrs
      putStrLn (show seed <> " (measurement) -> " <> show (VS.toList (gpuMeasurements res)))
      putStrLn ("Tableau shape = " <> show (gpuRows res, gpuCols res))
      printTableauPrefix (gpuRows res) (gpuCols res) 8 16 (gpuTableau res)

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

testTeleport :: IO ()
testTeleport =
  runTest "Quantum Teleportation"
    3
    [ Instr 1 0 0      -- H on qubit 0
    , Instr 3 0 1      -- CNOT qubit 0 -> qubit 1
    , Instr 0 0 0      -- Measure qubit 0
    , Instr 0 1 0      -- Measure qubit 1
    , Instr 3 1 2      -- CNOT qubit 1 -> qubit 2
    , Instr 1 1 0      -- H on qubit 1
    , Instr 0 1 0      -- Measure qubit 1
    , Instr 3 0 2      -- CNOT qubit 0 -> qubit 2
    , Instr 2 0 0      -- S on qubit 0
    ]


popCountOnes :: VS.Vector Int8 -> Int
popCountOnes = VS.foldl' (\acc x -> acc + if x /= 0 then 1 else 0) 0

showPrefix :: Int -> VS.Vector Int8 -> String
showPrefix k v = show (VS.toList (VS.take k v))

runBigTest :: String -> Int -> [Instr] -> IO ()
runBigTest name nQ instrs = do
  putStrLn ""
  putStrLn ("=== " <> name <> " ===")
  withFuthark $ \ctx ->
    forM_ ([1..10] :: [Int32]) $ \seed -> do
      res <- runSimulate ctx seed nQ instrs
      let m = gpuMeasurements res
          ones = popCountOnes m
      putStrLn (show seed <> " -> ones=" <> show ones <> "/" <> show nQ <> ", first32=" <> showPrefix 32 m)

test200AllHadamards :: IO ()
test200AllHadamards =
  runBigTest "200 qubits: H on all, then measure all"
    200
    ( [ Instr 1 (fromIntegral q) 0 | q <- [0..199] ]
   <> [ Instr 0 (fromIntegral q) 0 | q <- [0..199] ]
    )

test200CNOTChain :: IO ()
test200CNOTChain =
  runBigTest "200 qubits: H on 0, CNOT chain, then measure all"
    200
    ( [ Instr 1 0 0 ]
   <> [ Instr 3 (fromIntegral q) (fromIntegral (q+1)) | q <- [0..198] ]
   <> [ Instr 0 (fromIntegral q) 0 | q <- [0..199] ]
    )

testMANYAllHadamards :: Int -> IO ()
testMANYAllHadamards qubits =
  runBigTest (show qubits <> " qubits: H on all, then measure all")
    qubits
    ( [ Instr 1 (fromIntegral q) 0 | q <- [0..qubits-1] ]
   <> [ Instr 0 (fromIntegral q) 0 | q <- [0..qubits-1] ]
    )

testMANYCNOTChain :: Int -> IO ()
testMANYCNOTChain qubits =
  runBigTest (show qubits <> " qubits: H on 0, CNOT chain, then measure all")
    qubits
    ( [ Instr 1 0 0 ]
   <> [ Instr 3 (fromIntegral q) (fromIntegral (q+1)) | q <- [0..qubits-2] ]
   <> [ Instr 0 (fromIntegral q) 0 | q <- [0..qubits-1] ]
    )

main :: IO ()
main = do
  testDeterministicZero
  testHadamardRandom
  testBellPair
  testIndependentQubits
  testTeleport
  test200AllHadamards
  test200CNOTChain
  testMANYAllHadamards 2000
  testMANYCNOTChain 2000