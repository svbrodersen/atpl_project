{-# LANGUAGE LambdaCase #-}

module HQP.QOp.GPU.Compile
  ( Opcode(..)
  , Instr(..)
  , CompileError(..)
  , compileProgram
  , compileStep
  , compileUnitary
  ) where

import HQP.QOp.Syntax (Program, Step(..), QOp(..), op_qubits)
import HQP.QOp.Simplify (cleanop)
import HQP.QOp.HelperFunctions (prog_qubits)

import Data.Int (Int64)
import Data.Ratio ((%))

data Opcode
  = OpMeasure  -- 0
  | OpH        -- 1
  | OpS        -- 2
  | OpCNOT     -- 3
  | OpInit     -- 4
  deriving (Show, Eq)

data Instr = Instr
  { opcode :: !Int64
  , a      :: !Int64
  , b      :: !Int64
  } deriving (Show, Eq)

data CompileError
  = UnsupportedQOp QOp
  | DimMismatch { expected :: Int, got :: Int, badop :: QOp }
  | BadInitArity { ksLen :: Int, vsLen :: Int }
  deriving (Show, Eq)

opcodeTag :: Opcode -> Int64
opcodeTag = \case
  OpMeasure -> 0
  OpH       -> 1
  OpS       -> 2
  OpCNOT    -> 3
  OpInit    -> 4

mk :: Opcode -> Int -> Int -> Instr
mk op x y = Instr (opcodeTag op) (fromIntegral x) (fromIntegral y)

compileProgram :: Program -> Either CompileError (Int, [Instr])
compileProgram prog = do
  let n = prog_qubits prog
  chunks <- traverse (compileStep n) prog
  pure (n, concat chunks)

compileStep :: Int -> Step -> Either CompileError [Instr]
compileStep n = \case
  Unitary op ->
    let q = op_qubits op
    in if q == n
       then compileUnitary op
       else Left DimMismatch { expected = n, got = q, badop = op }

  Measure ks ->
    pure [ mk OpMeasure k 0 | k <- ks ]

  Initialize ks vs ->
    if length ks == length vs
    then pure [ mk OpInit k (if v then 1 else 0) | (k,v) <- zip ks vs ]
    else Left BadInitArity { ksLen = length ks, vsLen = length vs }

compileUnitary :: QOp -> Either CompileError [Instr]
compileUnitary op0 = go 0 (cleanop op0)
  where
    go :: Int -> QOp -> Either CompileError [Instr]
    go off = \case
      Id _    -> pure []
      Phase _ -> pure []  -- HQP global phase

      H       -> pure [ mk OpH off 0 ]
      R Z th | th == (1 % 2) -> pure [ mk OpS off 0 ]

      SX      -> pure [ mk OpH off 0, mk OpS off 0, mk OpH off 0 ]

      Tensor a b -> do
        as <- go off a
        bs <- go (off + op_qubits a) b
        pure (as <> bs)

      Compose op1 op2 -> do
        is2 <- go off op2
        is1 <- go off op1
        pure (is2 <> is1)

      C body -> do
        case findUniquePauli 0 body of
          Right (rel, 'X') ->
            pure [ mk OpCNOT off (off + 1 + rel) ]
          Right (rel, 'Z') ->
            let t = off + 1 + rel
            in pure [ mk OpH t 0
                    , mk OpCNOT off t
                    , mk OpH t 0
                    ]
          _ ->
            Left (UnsupportedQOp (C body))


      other -> Left (UnsupportedQOp other)

    findUniquePauli :: Int -> QOp -> Either CompileError (Int, Char)
    findUniquePauli rel op =
      case collect rel op of
        [p] -> Right p
        _   -> Left (UnsupportedQOp (C op))
      where
        collect :: Int -> QOp -> [(Int, Char)]
        collect r = \case
          X          -> [(r,'X')]
          Z          -> [(r,'Z')]
          Id _       -> []
          Phase _    -> []
          Tensor a b -> collect r a <> collect (r + op_qubits a) b
          _          -> []


    collectXPositions :: Int -> QOp -> [Int]
    collectXPositions rel = \case
      X          -> [rel]
      Id _       -> []
      Phase _    -> []
      Tensor a b -> collectXPositions rel a <> collectXPositions (rel + op_qubits a) b
      _          -> []