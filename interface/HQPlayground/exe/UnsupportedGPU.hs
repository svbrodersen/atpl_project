-- exe/UnsupportedGPU.hs
{-# LANGUAGE GHC2021 #-}

module Main where

import HQP.QOp.Syntax (Program, Step(..), QOp(..))
import HQP.QOp.GPU.Compile (compileProgram)

-- Z is a stabilizer gate, but your current GPU compiler does not lower it.
progUnsupportedZ :: Program
progUnsupportedZ =
  [ Unitary Z
  , Measure [0]
  ]

-- This is the exact pattern you previously saw: controlled body contains Z.
progUnsupportedControlledZ :: Program
progUnsupportedControlledZ =
  [ Unitary (C (Tensor (Id 1) Z))   -- 2-qubit operator: control on qubit 0, Z somewhere in the body
  , Measure [0,1]
  ]

main :: IO ()
main = do
  putStrLn "=== Unsupported Z ==="
  print (compileProgram progUnsupportedZ)

  putStrLn ""
  putStrLn "=== Unsupported controlled-Z-shaped op ==="
  print (compileProgram progUnsupportedControlledZ)
