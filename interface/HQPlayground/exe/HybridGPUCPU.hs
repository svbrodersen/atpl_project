{-# LANGUAGE LambdaCase #-}

module Main where

import Control.Monad (forM_)
import Data.Int (Int32, Int8)
import Data.Ratio ((%))
import System.Random (mkStdGen, randomRs)
import qualified Data.Vector.Storable as VS

import HQP.QOp.Syntax (Program, Step(..), QOp(..))
import HQP.QOp.HelperFunctions (prog_qubits)

import HQP.QOp.GPU.Compile (compileProgram)
import HQP.QOp.GPU.RunFuthark (withFuthark, runSimulateMeasurements)

import HQP.QOp.MatrixSemantics (evalProg, ket)

--------------------------------------------------------------------------------
-- Small HQP programs (actual HQP Program values)

-- Deterministic: measure |0>
progDet0 :: Program
progDet0 =
  [ Measure [0]
  ]

-- Random: H then measure
progHMeas :: Program
progHMeas =
  [ Unitary H
  , Measure [0]
  ]

-- Bell pair: (H ⊗ I); (C X) i.e. CNOT(0->1); measure both
progBell :: Program
progBell =
  [ Unitary (Tensor H (Id 1))
  , Unitary (C X)
  , Measure [0,1]
  ]

-- S on |0> then measure (still deterministic 0)
progSOnZero :: Program
progSOnZero =
  [ Unitary (R Z (1 % 2))
  , Measure [0]
  ]

-- Intentionally unsupported by GPU compiler: T = Rz(1/4)
progUnsupported :: Program
progUnsupported =
  [ Unitary (R Z (1 % 4))
  , Measure [0]
  ]

--------------------------------------------------------------------------------
-- RNG for CPU semantics (HQP uses RNG = [Double])

rngFromSeed :: Int32 -> [Double]
rngFromSeed s = randomRs (0.0, 1.0) (mkStdGen (fromIntegral s))

--------------------------------------------------------------------------------
-- GPU runner

runGPU :: Int32 -> Program -> IO (Either String [Int8])
runGPU seed prog =
  case compileProgram prog of
    Left err -> pure (Left (show err))
    Right (nQ, instrs) ->
      withFuthark $ \ctx -> do
        meas <- runSimulateMeasurements ctx seed nQ instrs
        pure (Right (VS.toList meas))

--------------------------------------------------------------------------------
-- CPU runner (MatrixSemantics)

runCPU :: Int32 -> Program -> IO ()
runCPU seed prog = do
  let nQ = prog_qubits prog
      initState = ket (replicate nQ 0)
      rng = rngFromSeed seed
      (_finalState, outcomes, _rng') = evalProg prog initState rng
  putStrLn ("CPU outcomes = " <> show outcomes)

--------------------------------------------------------------------------------
-- Hybrid policy: try GPU, fall back to CPU

runHybrid :: Int32 -> Program -> IO ()
runHybrid seed prog = do
  e <- runGPU seed prog
  case e of
    Right meas ->
      putStrLn ("GPU measurements = " <> show meas)
    Left why -> do
      putStrLn ("GPU not applicable (" <> why <> "), falling back to CPU.")
      runCPU seed prog

--------------------------------------------------------------------------------
-- Main

main :: IO ()
main = do
  let seed = 2026 :: Int32

  putStrLn ""
  putStrLn "=== progDet0 ==="
  runHybrid seed progDet0

  putStrLn ""
  putStrLn "=== progHMeas ==="
  runHybrid seed progHMeas

  putStrLn ""
  putStrLn "=== progBell ==="
  runHybrid seed progBell

  putStrLn ""
  putStrLn "=== progSOnZero ==="
  runHybrid seed progSOnZero

  putStrLn ""
  putStrLn "=== progUnsupported (should fall back) ==="
  runHybrid seed progUnsupported

  putStrLn ""
  putStrLn "=== progHMeas across seeds (GPU if supported) ==="
  forM_ ([1..20] :: [Int32]) $ \s -> do
    putStrLn ("seed " <> show s <> ":")
    runHybrid s progHMeas
