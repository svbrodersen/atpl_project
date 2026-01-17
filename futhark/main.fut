import "definitions"

entry simulate [n] (seed: i32) (num_qubits: i64) (gates: [n]i64) (cQ: [n]i64) (tQ: [n]i64) : ([][]i8, []i8) =
  let zipped_gates = zip3 gates cQ tQ
  let num_measurements = reduce (+) 0 <| map (\(gate, _, _) -> if gate == 0 then 1 else 0) zipped_gates
  let (tableu, _, _, measurements, _) =
    loop (tableu, i, rng, measurements, measurement_count) = (initial_tableu num_qubits, 0, rng_engine.rng_from_seed [seed], replicate num_measurements 0, 0)
    while i < n do
      let (gate, control, target) = zipped_gates[i]
      in if gate == 0
         then let (new_rng, new_tab, measurement) = Measurement rng tableu control
              in ( new_tab
                 , i + 1
                 , new_rng
                 , measurements with [measurement_count] = measurement
                 , measurement_count + 1
                 )
         else if gate == 1
         then (Hadamard tableu control, i + 1, rng, measurements, measurement_count)
         else if gate == 2
         then (Phase tableu control, i + 1, rng, measurements, measurement_count)
         else if gate == 3
         then (CNOT tableu control target, i + 1, rng, measurements, measurement_count)
         else -- do nothing
              (tableu, i + 1, rng, measurements, measurement_count)
  in (tableu, measurements)

-- Note the first nobench input is the Teleportation mentioned in the paper
-- ==
-- entry: main
-- nobench input {2026i32 1i64 [1i64, 0, 0] [0i64, 0, 0] [0i64, 0, 0]}
-- output {[1i8, 1]}
-- nobench input {1i32 1i64 [1i64, 0, 0] [0i64, 0, 0] [0i64, 0, 0]}
-- output {[0i8, 0]}
-- notest compiled input @ data/200_10.in
-- notest compiled input @ data/200_100.in
-- notest compiled input @ data/200_1000.in
-- notest compiled input @ data/200_10000.in
-- notest compiled input @ data/200_100000.in
-- notest compiled input @ data/400_10.in
-- notest compiled input @ data/400_100.in
-- notest compiled input @ data/400_1000.in
-- notest compiled input @ data/400_10000.in
-- notest compiled input @ data/400_100000.in
-- notest compiled input @ data/800_10.in
-- notest compiled input @ data/800_100.in
-- notest compiled input @ data/800_1000.in
-- notest compiled input @ data/800_10000.in
-- notest compiled input @ data/800_100000.in
-- notest compiled input @ data/1600_10.in
-- notest compiled input @ data/1600_100.in
-- notest compiled input @ data/1600_1000.in
-- notest compiled input @ data/1600_10000.in
-- notest compiled input @ data/1600_100000.in
-- notest compiled input @ data/3200_10.in
-- notest compiled input @ data/3200_100.in
-- notest compiled input @ data/3200_1000.in
-- notest compiled input @ data/3200_10000.in
-- notest compiled input @ data/3200_100000.in
-- notest compiled input @ data/6400_10.in
-- notest compiled input @ data/6400_100.in
-- notest compiled input @ data/6400_1000.in
-- notest compiled input @ data/6400_10000.in
-- notest compiled input @ data/6400_100000.in
-- notest compiled input @ data/12000_10.in
-- notest compiled input @ data/12000_100.in
-- notest compiled input @ data/12000_1000.in
-- notest compiled input @ data/12000_10000.in
-- notest compiled input @ data/12000_100000.in
-- notest compiled input @ data/24000_10.in
-- notest compiled input @ data/24000_100.in
-- notest compiled input @ data/24000_1000.in
-- notest compiled input @ data/24000_10000.in
-- notest compiled input @ data/24000_100000.in
-- notest compiled input @ data/48000_10.in
-- notest compiled input @ data/48000_100.in
-- notest compiled input @ data/48000_1000.in
-- notest compiled input @ data/48000_10000.in
-- notest compiled input @ data/48000_100000.in
-- notest compiled input @ data/96000_10.in
-- notest compiled input @ data/96000_100.in
-- notest compiled input @ data/96000_1000.in
-- notest compiled input @ data/96000_10000.in
-- notest compiled input @ data/96000_100000.in
entry main [n] (seed: i32) (num_qubits: i64) (gates: [n]i64) (cQ: [n]i64) (tQ: [n]i64) : []i8 =
  let (_, measurements) = simulate seed num_qubits gates cQ tQ
  in measurements

-- Simple function without the seed
-- ==
-- entry: simple
-- nobench input {2i64 [0i64, 0] [0i64, 1] [0i64, 0]}
-- output {[0i8, 0]}
-- nobench input {2i64 [3i64, 0, 0] [0i64, 0, 1] [1i64, 0, 0]}
-- output {[0i8, 0]}
-- nobench input {2i64 [3i64, 0, 0] [1i64, 0, 1] [0i64, 0, 0]}
-- output {[0i8, 0]}
-- nobench input {1i64 [1i64, 2, 2, 1, 0] [0i64, 0, 0, 0, 0] [0i64, 0, 0, 0, 0]}
-- output {[1i8]}
-- nobench input {2i64 [1i64, 2, 2, 1, 3, 0] [0i64, 0, 0, 0, 0, 1] [0i64, 0, 0, 0, 1, 0]}
-- output {[1i8]}
entry simple [n] (num_qubits: i64) (gates: [n]i64) (cQ: [n]i64) (tQ: [n]i64) : []i8 =
  main 2026 num_qubits gates cQ tQ

def teleportation (seed: i32) : ([][]i8, []i8) =
  let (tg, tc, tt) =
    ( [1i64, 3, 3, 1, 0, 0, 3, 3, 3, 1, 3, 1]
    , [1i64, 1, 0, 0, 0, 1, 0, 1, 4, 2, 3, 2]
    , [0i64, 2, 1, 0, 0, 0, 3, 4, 2, 0, 2, 0]
    )
  let (gates, cQ, tQ) = ([1] ++ tg ++ [0], [0] ++ tc ++ [4], [0] ++ tt ++ [0])
  in simulate seed 5 gates cQ tQ

def estimate_teleportation (n: i64) : f64 =
  let measurement seed =
    let tmp = teleportation (i32.i64 seed) |> (.1)
    in tmp[2]
  let mean =
    map (\i -> f64.i8 <| measurement i) (iota n) |> reduce (+) 0
    |> (/ (f64.i64 n))
  in mean
