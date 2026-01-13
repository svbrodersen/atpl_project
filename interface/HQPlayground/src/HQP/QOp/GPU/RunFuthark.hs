{-# LANGUAGE ForeignFunctionInterface #-}

module HQP.QOp.GPU.RunFuthark
  ( withFuthark
  , GPUResult(..)
  , runSimulate
  , runSimulateMeasurements
  , runSimulateDebug
  ) where

import Control.Exception (bracket)
import Control.Monad (when)
import Foreign
import Foreign.C.Types
import Foreign.C.String (CString, peekCString)
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import Data.Int (Int64, Int8, Int32)

import HQP.QOp.GPU.Compile (Instr(..))

data Ctx
data Cfg
data ArrI64_1d
data ArrI8_1d
data ArrI8_2d

data GPUResult = GPUResult
  { gpuMeasurements :: !(VS.Vector Int8)
  , gpuRows         :: !Int
  , gpuCols         :: !Int
  , gpuTableau      :: !(VS.Vector Int8)
  } deriving (Show)

foreign import ccall unsafe "futhark_context_config_new"
  cfg_new :: IO (Ptr Cfg)

foreign import ccall unsafe "futhark_context_config_free"
  cfg_free :: Ptr Cfg -> IO ()

foreign import ccall unsafe "futhark_context_new"
  ctx_new :: Ptr Cfg -> IO (Ptr Ctx)

foreign import ccall unsafe "futhark_context_free"
  ctx_free :: Ptr Ctx -> IO ()

foreign import ccall unsafe "futhark_context_sync"
  ctx_sync :: Ptr Ctx -> IO CInt

foreign import ccall unsafe "futhark_context_get_error"
  ctx_get_error :: Ptr Ctx -> IO CString

foreign import ccall unsafe "futhark_new_i64_1d"
  new_i64_1d :: Ptr Ctx -> Ptr Int64 -> Int64 -> IO (Ptr ArrI64_1d)

foreign import ccall unsafe "futhark_free_i64_1d"
  free_i64_1d :: Ptr Ctx -> Ptr ArrI64_1d -> IO CInt

foreign import ccall unsafe "futhark_entry_simulate"
  entry_simulate
    :: Ptr Ctx
    -> Ptr (Ptr ArrI8_2d)   -- out0 tableau
    -> Ptr (Ptr ArrI8_1d)   -- out1 measurements
    -> Int32                -- seed
    -> Int64                -- num_qubits
    -> Ptr ArrI64_1d        -- gates
    -> Ptr ArrI64_1d        -- cQ
    -> Ptr ArrI64_1d        -- tQ
    -> IO CInt

foreign import ccall unsafe "futhark_free_i8_2d"
  free_i8_2d :: Ptr Ctx -> Ptr ArrI8_2d -> IO CInt

foreign import ccall unsafe "futhark_free_i8_1d"
  free_i8_1d :: Ptr Ctx -> Ptr ArrI8_1d -> IO CInt

foreign import ccall unsafe "futhark_shape_i8_1d"
  shape_i8_1d :: Ptr Ctx -> Ptr ArrI8_1d -> IO (Ptr Int64)

foreign import ccall unsafe "futhark_values_i8_1d"
  values_i8_1d :: Ptr Ctx -> Ptr ArrI8_1d -> Ptr Int8 -> IO CInt

foreign import ccall unsafe "futhark_shape_i8_2d"
  shape_i8_2d :: Ptr Ctx -> Ptr ArrI8_2d -> IO (Ptr Int64)

foreign import ccall unsafe "futhark_values_i8_2d"
  values_i8_2d :: Ptr Ctx -> Ptr ArrI8_2d -> Ptr Int8 -> IO CInt

withFuthark :: (Ptr Ctx -> IO a) -> IO a
withFuthark k = bracket cfg_new cfg_free $ \cfg ->
  bracket (ctx_new cfg) ctx_free k

packInstrs :: [Instr] -> (VS.Vector Int64, VS.Vector Int64, VS.Vector Int64)
packInstrs instrs =
  ( VS.fromList [ opcode i | i <- instrs ]
  , VS.fromList [ a i      | i <- instrs ]
  , VS.fromList [ b i      | i <- instrs ]
  )

dieRC :: Ptr Ctx -> String -> CInt -> IO a
dieRC ctx where_ rc = do
  cstr <- ctx_get_error ctx
  msg <- peekCString cstr
  error (where_ <> " failed rc=" <> show rc <> " err=" <> msg)

syncOrDie :: Ptr Ctx -> IO ()
syncOrDie ctx = do
  rc <- ctx_sync ctx
  when (rc /= 0) $ dieRC ctx "futhark_context_sync" rc

readI8_1d_values :: Ptr Ctx -> Ptr ArrI8_1d -> IO (VS.Vector Int8)
readI8_1d_values ctx arr = do
  shp <- shape_i8_1d ctx arr
  len <- fromIntegral <$> peekElemOff shp 0
  mv <- VSM.new len
  VSM.unsafeWith mv $ \pOut -> do
    rc <- values_i8_1d ctx arr pOut
    when (rc /= 0) $ dieRC ctx "futhark_values_i8_1d" rc
  syncOrDie ctx
  VS.unsafeFreeze mv

readI8_2d_values :: Ptr Ctx -> Ptr ArrI8_2d -> IO (Int, Int, VS.Vector Int8)
readI8_2d_values ctx arr = do
  shp <- shape_i8_2d ctx arr
  r64 <- peekElemOff shp 0
  c64 <- peekElemOff shp 1
  let r = fromIntegral r64
      c = fromIntegral c64
      len = r * c
  mv <- VSM.new len
  VSM.unsafeWith mv $ \pOut -> do
    rc <- values_i8_2d ctx arr pOut
    when (rc /= 0) $ dieRC ctx "futhark_values_i8_2d" rc
  syncOrDie ctx
  v <- VS.unsafeFreeze mv
  pure (r, c, v)

callSimulate
  :: Ptr Ctx -> Int32 -> Int -> [Instr]
  -> IO (Ptr ArrI8_2d, Ptr ArrI8_1d)
callSimulate ctx seed numQubits instrs = do
  let (gates, cQ, tQ) = packInstrs instrs
      n64 = fromIntegral (VS.length gates) :: Int64
  VS.unsafeWith gates $ \pG ->
    VS.unsafeWith cQ $ \pC ->
      VS.unsafeWith tQ $ \pT -> do
        gArr <- new_i64_1d ctx pG n64
        cArr <- new_i64_1d ctx pC n64
        tArr <- new_i64_1d ctx pT n64
        alloca $ \outTabPtr ->
          alloca $ \outMeasPtr -> do
            rc <- entry_simulate ctx outTabPtr outMeasPtr seed (fromIntegral numQubits) gArr cArr tArr
            _ <- free_i64_1d ctx gArr
            _ <- free_i64_1d ctx cArr
            _ <- free_i64_1d ctx tArr
            when (rc /= 0) $ dieRC ctx "futhark_entry_simulate" rc
            syncOrDie ctx
            tabArr  <- peek outTabPtr
            measArr <- peek outMeasPtr
            pure (tabArr, measArr)

runSimulateMeasurements :: Ptr Ctx -> Int32 -> Int -> [Instr] -> IO (VS.Vector Int8)
runSimulateMeasurements ctx seed numQubits instrs = do
  (tabArr, measArr) <- callSimulate ctx seed numQubits instrs
  meas <- readI8_1d_values ctx measArr
  _ <- free_i8_1d ctx measArr
  _ <- free_i8_2d ctx tabArr
  pure meas

runSimulateDebug
  :: Ptr Ctx -> Int32 -> Int -> [Instr]
  -> IO (VS.Vector Int8, (Int, Int, VS.Vector Int8))
runSimulateDebug ctx seed numQubits instrs = do
  (tabArr, measArr) <- callSimulate ctx seed numQubits instrs
  meas <- readI8_1d_values ctx measArr
  (r, c, tab) <- readI8_2d_values ctx tabArr
  _ <- free_i8_1d ctx measArr
  _ <- free_i8_2d ctx tabArr
  pure (meas, (r, c, tab))

runSimulate :: Ptr Ctx -> Int32 -> Int -> [Instr] -> IO GPUResult
runSimulate ctx seed numQubits instrs = do
  (tabArr, measArr) <- callSimulate ctx seed numQubits instrs
  meas <- readI8_1d_values ctx measArr
  (r, c, tab) <- readI8_2d_values ctx tabArr
  _ <- free_i8_1d ctx measArr
  _ <- free_i8_2d ctx tabArr
  pure GPUResult
    { gpuMeasurements = meas
    , gpuRows = r
    , gpuCols = c
    , gpuTableau = tab
    }
