import "definitions"

def simulate [n] (seed: i32) (num_qubits: i64) (gates: [n]i64) (cQ: [n]i64) (tQ: [n]i64) : ([][]i8, []i8) =
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

entry main [n] (seed: i32) (num_qubits: i64) (gates: [n]i64) (cQ: [n]i64) (tQ: [n]i64) : []i8 =
  let (_, measurements) = trace (simulate seed num_qubits gates cQ tQ)
  in measurements
