-- exe/GPUHQPTests.hs
{-# LANGUAGE GHC2021 #-}

module Main where

import Control.Monad (forM_)
import Data.Int (Int32)
import Data.Ratio ((%))
import qualified Data.Vector.Storable as VS

import HQP.QOp.Syntax (Program, Step(..), QOp(..), (@>), (<@))
import HQP.QOp.GPU.Compile (compileProgram)
import HQP.QOp.GPU.RunFuthark (withFuthark, runSimulate, GPUResult(..))

--------------------------------------------------------------------------------
-- HQP helpers: build full-width unitaries that match the backend requirement
-- that every Unitary step spans the whole circuit width.

on1 :: Int -> Int -> QOp -> QOp
on1 n q op = (q <@ op) @> (n - q - 1)

sGate :: QOp
sGate = R Z (1 % 2)   -- compiler maps this to OpS (Phase kernel)

--------------------------------------------------------------------------------
-- A small test runner

runHQP :: String -> Program -> IO ()
runHQP name prog = do
  putStrLn ""
  putStrLn ("=== " <> name <> " ===")
  case compileProgram prog of
    Left err -> putStrLn ("GPU compile failed: " <> show err)
    Right (nQ, instrs) ->
      withFuthark $ \ctx ->
        forM_ ([1..10] :: [Int32]) $ \seed -> do
          res <- runSimulate ctx seed nQ instrs
          putStrLn (show seed <> " -> meas " <> show (VS.toList (gpuMeasurements res))
                               <> " | tableau " <> show (gpuRows res, gpuCols res))

--------------------------------------------------------------------------------
-- Programs (all within the supported gate subset)

progDeterministic0 :: Program
progDeterministic0 =
  [ Measure [0] ]

progHMeasure :: Program
progHMeasure =
  [ Unitary H
  , Measure [0]
  ]

progSThenMeasure :: Program
progSThenMeasure =
  [ Unitary sGate
  , Measure [0]
  ]

progBell :: Program
progBell =
  [ Unitary (on1 2 0 H)      -- H on qubit 0
  , Unitary (C X)            -- control qubit 0, target qubit 1 (body spans remaining qubits)
  , Measure [0,1]
  ]

progGHZ3 :: Program
progGHZ3 =
  [ Unitary (on1 3 0 H)                -- H on qubit 0
  , Unitary (C (Tensor X (Id 1)))      -- CNOT 0 -> 1  (X at relative index 0 in the body)
  , Unitary (C (Tensor (Id 1) X))      -- CNOT 0 -> 2  (X at relative index 1 in the body)
  , Measure [0,1,2]
  ]

progLarge200 :: Program
progLarge200 =
  [ Unitary (on1 200 0 H)   -- put qubit 0 in |+>
  , Measure [0]             -- measurement should vary with the seed
  ]

main :: IO ()
main = do
  runHQP "Deterministic |0> measurement" progDeterministic0
  runHQP "Single-qubit H then measure (should vary by seed)" progHMeasure
  runHQP "Single-qubit S then measure (still deterministic on |0>)" progSThenMeasure
  runHQP "Bell pair (measurements correlated)" progBell
  runHQP "GHZ on 3 qubits (all bits equal)" progGHZ3
  runHQP "Large n=200 (H on qubit 0 then measure)" progLarge200
