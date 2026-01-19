{-# LANGUAGE LambdaCase #-}


import Control.Monad (forM_)
import Data.Int (Int32, Int8)
import qualified Data.Vector.Storable as VS

import HQP.QOp.Syntax (Program, Step(..), QOp(..), Nat)
import HQP.QOp.GPU.Compile (compileProgram, Instr(..))
import HQP.QOp.GPU.RunFuthark (withFuthark, runSimulate, GPUResult(..))

--------------------------------------------------------------------------------
-- Small helpers

printTableauPrefix :: Int -> Int -> Int -> Int -> VS.Vector Int8 -> IO ()
printTableauPrefix rows cols maxR maxC tab = do
  let rLim = min rows maxR
      cLim = min cols maxC
      idx r c = tab VS.! (r * cols + c)
  forM_ [0..rLim-1] $ \r -> do
    let row = [ idx r c | c <- [0..cLim-1] ]
    putStrLn ("Row " <> show r <> ": " <> show row <> if cLim < cols then " ..." else "")

runHQP :: String -> Program -> IO ()
runHQP name prog = do
  putStrLn ""
  putStrLn ("=== " <> name <> " ===")
  case compileProgram prog of
    Left err -> print err
    Right (nQ, instrs) -> do
      putStrLn ("nQ = " <> show nQ <> ", instrs = " <> show (length instrs))
      withFuthark $ \ctx ->
        forM_ ([1..20] :: [Int32]) $ \seed -> do
          res <- runSimulate ctx seed nQ instrs
          putStrLn (show seed <> " -> meas " <> show (VS.toList (gpuMeasurements res)))
          putStrLn ("tableau shape " <> show (gpuRows res, gpuCols res))
          printTableauPrefix (gpuRows res) (gpuCols res) 4 16 (gpuTableau res)

--------------------------------------------------------------------------------
-- HQP programs that stay inside your supported subset:
-- Supported by compiler: H, S (as R Z (1/2)), CNOT (as C body with exactly one X), Measure, Initialize.

-- Deterministic: measure |0>
progDeterministic0 :: Program
progDeterministic0 =
  [ Measure [0]
  ]

-- 1-qubit random: H then measure
progHThenMeasure :: Program
progHThenMeasure =
  [ Unitary H
  , Measure [0]
  ]

-- Bell pair: H on q0, CNOT 0->1, measure both
-- CNOT 0->1 on 2 qubits is just C X (control is implicit as the extra qubit).
progBell :: Program
progBell =
  [ Unitary (H `Tensor` Id 1)
  , Unitary (C X)
  , Measure [0,1]
  ]

-- Two independent qubits: H on both, measure both
progIndependent :: Program
progIndependent =
  [ Unitary (H `Tensor` H)
  , Measure [0,1]
  ]

-- 200 qubits: H on qubit 0, then measure qubit 0.
-- This exercises the large-tableau path without needing huge measurement output.
prog200 :: Program
prog200 =
  [ Unitary (H `Tensor` Id 199)
  , Measure [0]
  ]

--------------------------------------------------------------------------------

main :: IO ()
main = do
  runHQP "HQP deterministic |0> (expect always 0)" progDeterministic0
  runHQP "HQP H then measure (expect varying 0/1)" progHThenMeasure
  runHQP "HQP Bell pair (expect [0,0] or [1,1])" progBell
  runHQP "HQP two independent qubits (expect all four combos across seeds)" progIndependent
  runHQP "HQP 200 qubits (H on q0, measure q0)" prog200
