module HQP.QOp.TensorSemantics
  ( StateT
  , OpT(..)
  , apply
  , evalOp      -- :: QOp -> OpT
  , evalStep    -- :: (StateT, Outcomes, RNG) -> Step -> (StateT, Outcomes, RNG)
  , evalProg    -- :: Program -> StateT -> RNG -> (StateT, Outcomes, RNG)
  , measure1    -- :: (StateT, Outcomes, RNG) -> Int -> (StateT, Outcomes, RNG)
  , applyPauliProduct
  ) where

import Data.Complex

import HQP.QOp.Syntax          
import HQP.QOp.HelperFunctions
import HQP.PrettyPrint.PrettyOp
import Data.Bits(shiftL,xor,testBit, setBit, clearBit)
import Data.Array(accumArray,elems)
--import qualified Torch as T
--import qualified Torch.Functional as F
--import qualified Torch.Tensor as TT
import Data.Bits (testBit, setBit)
import Data.List (foldl')
import Data.Massiv.Array
  ( Array, Comp(Seq)
  , Ix2((:.)), IxN(Ix3,(:>))     -- constructor imports (this is the key fix)
  , Sz(..), Sz2
  , (<!)
  , size, resize', makeArray
  , backpermute'
  , D
  )
import Data.Massiv.Array.Unsafe (unsafeIndex)
import qualified Data.Massiv.Array as A

type StateT   = Array D Ix2 ComplexT

--------------------------------------------------------------------------------
-- Backend-wide conventions
--------------------------------------------------------------------------------

newtype OpT = OpT { runOpT :: StateT -> StateT }
instance Semigroup OpT where
  (OpT f) <> (OpT g) = OpT (f . g)
instance Monoid OpT where
  mempty = OpT id
instance HasQubits StateT where
  n_qubits ψ =
    case size ψ of
      Sz2 d _ -> ilog2 d   -- matrix view [2^n, r]
    

evalOp   :: QOp -> OpT
evalStep :: (StateT, Outcomes, RNG) -> Step -> (StateT, Outcomes, RNG)
evalProg :: Program -> StateT -> RNG -> (StateT, Outcomes, RNG) 



apply :: OpT -> StateT -> StateT
apply = runOpT

--------------------------------------------------------------------------------
-- Pretty infix ops
--------------------------------------------------------------------------------

infixl 6 .+., .-.
infixl 7 .*.

--- BACKEND-SPECIFIC IMPLEMENTATION HERE - MOVE TO BACKEND FILE LATER ---
{-| Torch implementation -- but Torch turns out to be difficult to install. :'-(. Doing a slower Massiv implementation for now.
type StateT = T.Tensor

instance HasQubits StateT where
  n_qubits ψ =
    case T.shape ψ of
    []        -> 0                      -- scalar
    [d]       -> ilog2 d                -- vector view [2^n]
    (d:_)     -> ilog2 d                -- matrix view [2^n, rest] -- possibly do product of dims instead.

  

-- | Scalar multiplication: (c .*. Ψ) = c Ψ.
(.*.) :: ComplexT -> StateT -> StateT
--(.*.) = F.mulScalar
(.*.) c ψ = F.mul (T.asTensor c) ψ


-- Precondition: product ds = numel x.
view :: [Int] -> StateT -> StateT
view dims x =
  let x' = if TT.isContiguous x then x else TT.contiguous x
  in F.view dims x'

-- | 'permuteIx ks x' permutes the axes of x according to the permutation ks.
permuteIx :: [Int] -> StateT -> StateT
permuteIx = F.permute

--------------------------------------------------------------------------------
-- Selector split/stack (DirectSum / control)
--------------------------------------------------------------------------------
-- | 'splitSel k r' splits a (k+1)-qubit matrix-view Ψ : ℂ^{2^(k+1) × r} by the
-- front selector bit, producing (Ψ0, Ψ1) with Ψc : ℂ^{2^k × r}.
--
-- Write a basis index x∈{0,1}^{k+1} as x=(c,r) with c∈{0,1}, r∈{0,1}^k.
-- Then:
--   (Ψ0)_{r,β} = Ψ_{(0,r),β}
--   (Ψ1)_{r,β} = Ψ_{(1,r),β}.
splitFront :: Int -> Int -> StateT -> (StateT, StateT)
splitFront k r psi =
  let psi3 = view [2,2^k, r] psi
  in case F.split 1 (T.Dim 0) psi3 of
       [p0,p1] -> (view [2^k, r] p0, view [2^k, r] p1)
       _       -> error "splitFront: expected leading dim == 2"

-- | 'stackSel k r Ψ0 Ψ1' is the inverse of 'splitSel':
-- given Ψ0,Ψ1 : ℂ^{2^k × r}, it returns Ψ : ℂ^{2^(k+1) × r} defined by
--   Ψ_{(0,r),β} = (Ψ0)_{r,β}
--   Ψ_{(1,r),β} = (Ψ1)_{r,β}.
stackFront :: Int -> Int -> StateT -> StateT -> StateT
stackFront k r psi0 psi1 =
  let p0 = view [1, 2^k, r] psi0
      p1 = view [1, 2^k, r] psi1
  in view [2^(k+1), r] (F.cat (T.Dim 0) [p0,p1])

-}

{-| Massiv implementation -}


--------------------------------------------------------------------------------
-- O(1) reshape/view (resize') for only  the ranks we need
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- reshape/view
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- view (reshape, delayed)
--------------------------------------------------------------------------------

view :: [Int] -> StateT -> StateT
view [m,n] = resize' (Sz2 m n)
view _     = error "view: massiv backend only supports Ix2 here (matrix view)"

--------------------------------------------------------------------------------
-- swap axes 0 and 1 on Ix3 (delayed)
--------------------------------------------------------------------------------

swap01_3 :: Array D (IxN 3) e -> Array D (IxN 3) e
swap01_3 x =
  let Sz (a :> (b :. c)) = size x
      outSz = Sz (b :> (a :. c))
  in backpermute' outSz (\(j :> (i :. k)) -> (i :> (j :. k))) x

--------------------------------------------------------------------------------
-- split/stack on front selector (delayed)
--------------------------------------------------------------------------------

splitFront :: Int -> Int -> StateT -> (StateT, StateT)
splitFront k r psi =
  let psi3 = resize' (Sz (2 :> (2^k :. r))) psi  -- Ix3 view
      p0   = resize' (Sz2 (2^k) r) (psi3 <! 0)    -- outer slice, delayed
      p1   = resize' (Sz2 (2^k) r) (psi3 <! 1)
  in (p0,p1)

stackFront :: Int -> Int -> StateT -> StateT -> StateT
stackFront k r psi0 psi1 =
  let a0   = resize' (Sz2 (2^k) r) psi0
      a1   = resize' (Sz2 (2^k) r) psi1
      out3 = makeArray Seq (Sz (2 :> (2^k :. r))) $ \(c :> (i :. j)) ->
               if c == 0 then unsafeIndex a0 (i :. j)
                         else unsafeIndex a1 (i :. j)
  in resize' (Sz2 (2^(k+1)) r) out3




--------------------------------------------------------------------------------
-- permuteLocal on row-index bits (delayed)
--------------------------------------------------------------------------------

permuteLocal :: [Int] -> Int -> Int -> StateT -> StateT
permuteLocal ks k r psi =
  let nRows = 2^k
      psi2  = resize' (Sz2 nRows r) psi
      outSz = Sz2 nRows r
  in backpermute' outSz (\(x :. β) -> (permBits ks x :. β)) psi2

-- | permBits ks x = Σ_p bit(x,p) · 2^{ks(p)}.
permBits :: [Int] -> Int -> Int
permBits ks x =
  foldl' (\acc (p,q) -> if testBit x p then setBit acc q else acc) 0 (zip [0..] ks)


--------------------------------------------------------------------------------
-- Scalar mul / add / sub (fusible)
--------------------------------------------------------------------------------
(.*.) :: ComplexT -> StateT -> StateT
(.*.) c = A.map (c *)             -- delayed map fuses

(.+.) :: StateT -> StateT -> StateT
(.+.) = A.zipWith (+)

(.-.) :: StateT -> StateT -> StateT
(.-.) = A.zipWith (-)
--------------------------------------------------------------------------------
-- 1-qubit atoms on the *frontmost* qubit (direct tensor ops, no matmul)
--------------------------------------------------------------------------------

-- Each of these acts on Ψ : ℂ^{2^k × B} by exposing Ψ as ℂ^{2 × 2^(k-1) × B}
-- (front qubit axis explicit), applying the 2×2 rule along that axis, and packing back.

-- | applyX1 implements (X ⊗ Id (k-1)) on k qubits.
-- Action on basis:
--   X|0⟩ = |1⟩,  X|1⟩ = |0⟩
-- so it swaps the two slices of the front axis:
--   Ψ'_{0,r,β} = Ψ_{1,r,β},  Ψ'_{1,r,β} = Ψ_{0,r,β}.

-- | applyZ1 implements (Z ⊗ Id (k-1)) on k qubits.
-- Action on basis:
--   Z|0⟩ = |0⟩,  Z|1⟩ = -|1⟩
-- hence:
--   Ψ'_{0,r,β} = Ψ_{0,r,β},  Ψ'_{1,r,β} = -Ψ_{1,r,β}.


-- | applyY1 implements (Y ⊗ Id (k-1)) on k qubits.
-- Action on basis:
--   Y|0⟩ =  i|1⟩,   Y|1⟩ = -i|0⟩
-- hence:
--   Ψ'_{0,r,β} = -i Ψ_{1,r,β}
--   Ψ'_{1,r,β} =  i Ψ_{0,r,β}.


-- | applyH1 implements (H ⊗ Id (k-1)) on k qubits.
-- Action on basis:
--   H|0⟩ = (|0⟩+|1⟩)/√2
--   H|1⟩ = (|0⟩-|1⟩)/√2
-- hence:
--   Ψ'_{0,r,β} = (Ψ_{0,r,β} + Ψ_{1,r,β})/√2
--   Ψ'_{1,r,β} = (Ψ_{0,r,β} - Ψ_{1,r,β})/√2.

-- Assumed to be provided by the backend module:
--   splitFront :: Int -> Int -> StateT -> (StateT, StateT)
--   stackFront :: Int -> Int -> StateT -> StateT -> StateT

-- | (X ⊗ Id) acting on the frontmost qubit
applyX1 :: Int -> Int -> StateT -> StateT
applyX1 k r ψ =
  let (p0, p1) = splitFront (k-1) r ψ
  in  stackFront (k-1) r p1 p0

-- | (Z ⊗ Id) acting on the frontmost qubit
applyZ1 :: Int -> Int -> StateT -> StateT
applyZ1 k r ψ =
  let (p0, p1) = splitFront (k-1) r ψ
  in  stackFront (k-1) r p0 ((-1) .*. p1)

-- | (Y ⊗ Id) acting on the frontmost qubit
applyY1 :: Int -> Int -> StateT -> StateT
applyY1 k r ψ =
  let (p0, p1) = splitFront (k-1) r ψ
      i  = 0 :+ 1
      mi = 0 :+ (-1)
  in  stackFront (k-1) r (mi .*. p1) (i .*. p0)

-- | (H ⊗ Id) acting on the frontmost qubit
applyH1 :: Int -> Int -> StateT -> StateT
applyH1 k r ψ =
  let (p0, p1) = splitFront (k-1) r ψ
      s = (1 / sqrt 2) :+ 0
  in  stackFront (k-1) r
        (s .*. (p0 .+. p1))
        (s .*. (p0 .-. p1))

applySX1 :: Int -> Int -> StateT -> StateT
applySX1 k r ψ =
  let (ψ0, ψ1) = splitFront (k-1) r ψ
      (p, m)   = ((1/2) :+ (1/2) :: ComplexT, (1/2) :+ (-1/2) :: ComplexT)
      sxψ0 = (p .*. ψ0) .+. (m .*. ψ1)
      sxψ1 = (m .*. ψ0) .+. (p .*. ψ1)
  in stackFront (k-1) r sxψ0 sxψ1


-- HERFRA ER ALT BACKEND-UAFHÆNGIGT
--------------------------------------------------------------------------------
-- Structural combinators
--------------------------------------------------------------------------------
{-|onLeftBlock kA kB r opA applies opA = [[A]] to the left block of a state tensor
   viewed as [2^kA, 2^kB, r], i.e. it forms the Torch tensor operation

    Ψ ↦ [[A ⊗ (Id kB) Ψ]] 

    by reshaping the tensor [[\Psi]] : [2^(kA+kB), r] to [2^kA, 2^kB, r],
    contrating [[A]]: [2^kA, 2^kA] with the left index, and reshaping back.
-}
onLeftBlock :: Int -> Int -> Int -> OpT -> OpT
onLeftBlock kA kB r (OpT f) = OpT $ \x ->
  let xL = view [2^kA, 2^kB * r] x
  in  view [2^(kA+kB), r] (f xL)

{-|onRightBlock kA kB r opB applies opB = [[B]] to the right block of a state tensor, i.e.:

    Ψ ↦ [[(Id kA) ⊗ B Ψ]] 

  analogously to onLeftBlock.
-}
onRightBlock :: Int -> Int -> Int -> OpT -> OpT  -- TODO: Rewrite backend-independent using view 
onRightBlock kA kB r (OpT f) = OpT $ \ψ ->
  let t0 = resize' (Sz (2^kA :> (2^kB :. r))) ψ    -- Ix3 view  [2^kA,2^kB,r]
      tS = swap01_3 t0                               -- Ix3 view  [2^kB,2^kA,r]
      xM = resize' (Sz2 (2^kB) (2^kA * r)) tS        -- Ix2 view  [2^kB, 2^kA·r]
      yM = f xM
      yS = resize' (Sz (2^kB :> (2^kA :. r))) yM     -- Ix3 view  [2^kB,2^kA,r]
      tU = swap01_3 yS                               -- Ix3 view  [2^kA,2^kB,r]
  in  resize' (Sz2 (2^(kA+kB)) r) tU                 -- Ix2 view  [2^(kA+kB), r]


{-| onSelector k r f0 f1 implements branching on a selector qubit:
-- given f0,f1 acting on k-qubit states, produce an endomorphism on (k+1) qubits:
--
--   |0⟩⊗ψ  ↦  |0⟩⊗f0(ψ)
--   |1⟩⊗ψ  ↦  |1⟩⊗f1(ψ)
--
-- i.e. the block-diagonal operatora
--   |0⟩⟨0| ⊗ A0  +  |1⟩⟨1| ⊗ A1
-- where A0,A1 are the operators denoted by f0,f1.
--
-- Special case: controlled U is (Id k) ⊕ U, i.e. f0 = Id, f1 = U.
-}
onSelector::
     (StateT -> StateT)   -- f0 : action when selector = 0
  -> (StateT -> StateT)   -- f1 : action when selector = 1
  -> Int      -- k  : number of target qubits
  -> Int      -- r  : batch dimension
  -> StateT   -- input state
  -> StateT   -- output state
onSelector f0 f1 k r ψ =
  let (ψ0, ψ1) = splitFront k r ψ
  in  stackFront k r (f0 ψ0) (f1 ψ1)


--------------------------------------------------------------------------------
-- Fast application of Pauli-product (as OpT + phase)
--------------------------------------------------------------------------------

-- Apply a Pauli product P (up to global phase) to a matrix-view Psi:[2^k,r].
-- Returns (globalPhase, opT) where endomorphism is the actual action on Psi.
--
-- Allowed constructors for P (recommended):
--   Id, Phase, X,Y,Z, Permute, Tensor, Compose, Adjoint
-- (No H, no DirectSum/C, no R).
pauliOp :: QOp -> Int -> Int -> (Complex Double, OpT)
pauliOp op k r = go op
  where
    go = \case
      Id _       -> (1:+0,    mempty)
      Phase q    -> (exp(0:+ (pi * fromRational q)), mempty)

      X          -> (1:+0, OpT (applyX1 k r))
      Y          -> (1:+0, OpT (applyY1 k r))
      Z          -> (1:+0, OpT (applyZ1 k r))

      Permute ks -> (1:+0, OpT (permuteLocal ks k r))

      Compose a b ->
        let (ϕA, opA) = go a
            (ϕB, opB) = go b
        in (ϕA * ϕB, opA <> opB)

      Tensor a b ->
        let kA = op_qubits a
            kB = op_qubits b
            (ϕA, opA') = pauliOp a kA (2^kB * r)
            (ϕB, opB') = pauliOp b kB (2^kA * r)
            opA = onLeftBlock  kA kB r opA'
            opB = onRightBlock kA kB r opB'
        in (ϕA * ϕB, opB <> opA)

      Adjoint a ->
        -- X,Y,Z,H are Hermitian; permutation inverts; Phase negates.
        go (dagger a)

      _ ->
        error "pauliEndo: non-Pauli constructor encountered"

-- apply a Pauli string with phase op=[[ϕ⋅P]] to a concrete Ψ, returning (ϕ, [[PΨ]])).
applyPauliProduct :: QOp -> Int -> Int -> StateT -> (Complex Double, StateT)
applyPauliProduct op k r ψ =
  let (ϕ, opP) = pauliOp op k r
  in  (ϕ, apply opP ψ)

--------------------------------------------------------------------------------
-- QOp semantics on matrix-views Psi:[2^k,r]
--------------------------------------------------------------------------------

-- Internal: matrix-view semantics, closes over (k,r).
{-| evalOpMat op k r applies the operator op to a state tensor viewed as [2^k, r],
    returning an OpT representing the action on that view.

    Note: this function assumes that op is well-typed for k qubits, i.e. op_qubits op == k.
-}
evalOpMat :: QOp -> Int -> Int -> OpT
evalOpMat op k r = go op
  where
    go = \case
      Id _       -> mempty
      Phase q    -> OpT (exp( 0:+ (pi * fromRational q)) .*. )

      X          -> OpT (applyX1 k r)
      Y          -> OpT (applyY1 k r)
      Z          -> OpT (applyZ1 k r)
      H          -> OpT (applyH1 k r)
      SX         -> OpT (applySX1 k r)

      Permute ks -> OpT (permuteLocal ks k r)

      Compose a b -> go a <> go b
      Adjoint a   -> go (dagger a)

      Tensor a b  ->
        let kA  = op_qubits a
            kB  = op_qubits b
            opA = evalOpMat a kA (2^kB * r)
            opB = evalOpMat b kB (2^kA * r)
        in onRightBlock kA kB r opB <> 
           onLeftBlock  kA kB r opA

      DirectSum a b -> OpT $ onSelector (sub a) (sub b) (k-1) r 
      C a           -> OpT $ onSelector id      (sub a) (k-1) r 

      R axis θ | θ == 0 -> mempty
               | otherwise  -> OpT $ \x -> let
                  t         = pi * fromRational θ / 2
                  c         = cos t :+ 0
                  s         = sin t :+ 0
                  i        = 0 :+ 1
                  (ϕ, px) = applyPauliProduct axis k r x
                in  (c .*. x) .+.  (i*s*ϕ) .*. px
    
    sub a = apply $ evalOpMat a (k-1) r
          
--------------------------------------------------------------------------------
-- Public “backend semantics”: QOp -> OpT on full state tensors
--------------------------------------------------------------------------------

-- Convention: the runtime state is a rank-n tensor of shape [2,2,..,2] (n qubits).
-- evalOp uses only structural reshapes; it does not build global dense matrices.
evalOp op =
  let n = op_qubits op
  in OpT $ \psiN ->
      let psiT = view (replicate n 2) psiN
          psiM = view [2^n, 1] psiT
          outM = apply (evalOpMat op n 1) psiM
      in view (replicate n 2) outM

evalStep (st, outs, rng) step = let n = n_qubits st in 
  case step of
    Unitary op | (n_qubits op == n) -> (apply (evalOp op) st, outs, rng)
               | otherwise -> error $ "Dim-mismatch between " ++ showOp op ++ " and n="++show n

    -- outcomes are latest-first, so ks is reversed on input
    Measure ks -> foldl (measure1) (st, outs, rng) (reverse ks)
    Initialize ks vs ->
      let
        (st', os, rng') = evalStep (st, [], rng) (Measure ks)
        -- List of outcomes xor values for each initialized qubit
        corrections     = zipWith xor os vs
        -- Now we build the full list, including unaffected qubits        
        corrFull        = accumArray xor False (0,n-1) (zip ks corrections)
        corrOp          = foldl (⊗) One [ if c then X else I | c <- elems corrFull ]
      in (evalStep) (st', outs, rng') (Unitary corrOp)

  -- | evalProg simply executes a program step by step
evalProg steps psi0 rng0 = foldl evalStep (psi0, [], rng0) steps

measure1 :: (StateT, Outcomes, RNG) -> Int -> (StateT, Outcomes, RNG)
measure1 _ _ = error "Implement the measure1 function for Tensor backend"

--------------------------------------------------------------------------------
-- Dagger (structural)
--------------------------------------------------------------------------------

dagger :: QOp -> QOp
dagger = \case
  Id n          -> Id n
  Phase q       -> Phase (-q)
  X             -> X
  Y             -> Y
  Z             -> Z
  H             -> H
  SX            -> Adjoint SX -- Do we need applySXdagger?
  R a theta     -> R a (-theta)
  C a           -> C (dagger a)
  Permute ks    -> Permute (invertPerm ks)
  Tensor a b    -> Tensor (dagger a) (dagger b)
  DirectSum a b -> DirectSum (dagger a) (dagger b)
  Compose a b   -> Compose (dagger b) (dagger a)
  Adjoint a     -> a
