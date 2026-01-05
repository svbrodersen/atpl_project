-- [asdf;asdlfkasjfd;lkj, H, H, H, M, H, H]
--
-- function(comp) -> measurements
-- convert to list
-- possible potimisations
-- create the stuff below

-- |
-- |
-- |
-- |
-- |
-- |
-- |
-- |
--\ /

import Definitions

let rng0 = create Random Engine
let t0 = initial_tableu 10
let t1 = Hadamard t0 6
let t2 = Phase t1 4
let (rng1, t3, v1) = Measurement rng0 t2 6
let (rng2, t4, v2) = Measurement rng1 t1 4
in (t4, (v1, v2))
