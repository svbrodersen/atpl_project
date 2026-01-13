import "definitions"

module rng_engine = minstd_rand
module rand_i8 = uniform_int_distribution i8 rng_engine

def seed : i64 = 2026

def simulate [n] [m] [z] (num_qubits: i64) (gates: [n](i64, i64, i64)) : (tab [m], [z]t) =
  let num_measurements = reduce (+) 0 <| map (\(gate, _, _) -> if gate == 0 then 1 else 0) gates
  let res =
    loop (tableu, i, rng, measurements, measurement_count) = (initial_tableu num_qubits, 0, rng_engine.rng_from_seed seed, replicate num_measurements 0, 0)
    while i < n do
      let (gate, control, target) = gates[i]
      in if gate == 0
         then let (new_rng, new_tab, measurement) = Measurement tableu control target
              in (new_tab, i + 1, new_rng, measurements with [measurement_count] = measurement, measurement_count + 1)
         else if gate == 1
         then (Hadamard tableu control, i + 1, rng, measurements, measurement_count + 1)
         else if gate == 2
         then (Phase tableu control, i + 1, rng, measurements, measurement_count + 1)
         else if gate == 3
         then (CNOT tableu control target, i + 1, rng, measurements, measurement_count + 1)
         else -- do nothing
              (tableu, i + 1, rng, measurements, measurement_count + 1)
  in res

entry main [n] [m] (num_qubits: i64) (gates: [n](i64, i64, i64)) : ([m]t) =
  let (_, measurements) = simulate num_qubits gates
  in measurements
