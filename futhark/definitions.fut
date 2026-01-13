import "lib/github.com/diku-dk/cpprandom/random"

module rng_engine = minstd_rand
module rand_i8 = uniform_int_distribution i8 u32 rng_engine
type t = i8
type tab [n] = [2 * n + 1][2 * n + 1]t

local def size (n: i64) : i64 = 2 * n + 1

local
def g (x1: t) (x2: t) (z1: t) (z2: t) =
  if x1 == z1 && x1 == 0
  then 0
  else if x1 == z1 && x1 == 1
  then z2 - x2
  else if x1 == 1 && z1 == 0
  then z2 * (2 * x2 - 1)
  else if x1 == 0 && z1 == 1
  then x2 * (1 - 2 * z2)
  else 0

local
def rowsum [n] (tableu: tab [n]) (h: i64) (i: i64) : ([n * 2 + 1](i64, i64), [n * 2 + 1]i8) =
  let max_idx = (size n) - 1
  let tot_sum = map (\j -> g tableu[i][j] tableu[i][j + n] tableu[h][j] tableu[h][j + n]) (iota n) |> reduce (+) 0
  let res = (2 * tableu[h][max_idx] + 2 * tableu[i][max_idx] + tot_sum) % 4
  let v = if res == 0 then 0 else 1
  let is = map (\j -> [(h, j), (h, j + n)]) (iota n) |> flatten
  let vs = map (\j -> [tableu[i][j] ^ tableu[h][j], tableu[i][n + j] ^ tableu[h][n + j]]) (iota n) |> flatten
  in (is ++ [(h, max_idx)], vs ++ [v])

def initial_tableu (n: i64) : tab [n] =
  let tmp = 2 * n + 1
  in tabulate_2d tmp tmp (\i j ->
                            if j == tmp - 1 || i == tmp - 1
                            then -- r_i and last row
                                 0
                            else if i == j
                            then -- diagonal
                                 1
                            else 0)

def CNOT [n] (tableu: *tab [n]) (a: i64) (b: i64) : *tab [n] =
  let tmp = 2 * n
  let indices =
    map (\i ->
           [(i, tmp), (i, b), (i, n + a)])
        (iota tmp)
    |> flatten
  let values =
    map (\i ->
           let zia = tableu[i][n + a]
           let zib = tableu[i][n + b]
           let xia = tableu[i][a]
           let xib = tableu[i][b]
           let ri = tableu[i][tmp]
           in [ri ^ xia * zib * (xib ^ zia ^ 1i8), xib ^ xia, zia ^ zib])
        (iota tmp)
    |> flatten
  in scatter_2d tableu indices values

def Hadamard [n] (tableu: *tab [n]) (a: i64) : *tab [n] =
  let tmp = 2 * n
  let indices = map (\i -> [(i, tmp), (i, a), (i, n + a)]) (iota tmp) |> flatten
  let values =
    map (\i ->
           let ri = tableu[i][tmp]
           let xia = tableu[i][a]
           let zia = tableu[i][a + n]
           in [ri ^ xia * zia, zia, xia])
        -- swap in regards to indices, which are ri, xia, zia
        (iota tmp)
    |> flatten
  in scatter_2d tableu indices values

def Phase [n] (tableu: *tab [n]) (a: i64) : *tab [n] =
  let tmp = 2 * n
  let indices = map (\i -> [(i, tmp), (i, a + n)]) (iota tmp) |> flatten
  let values =
    map (\i ->
           let ri = tableu[i][tmp]
           let xia = tableu[i][a]
           let zia = tableu[i][n + a]
           in [ri ^ xia * zia, zia ^ xia])
        (iota tmp)
    |> flatten
  in scatter_2d tableu indices values

def Measurement [n] (eng: rng_engine.rng) (tableu: *tab [n]) (a: i64) : (rng_engine.rng, *tab [n], t) =
  let (p, xpa) = reduce_comm (\(p1, xp1a) (p2, xp2a) -> 
    if xp1a == 1 && xp2a == 1 then
      if p1 <= p2 then (p1, xp1a)
      else (p2, xp2a)
    else if xp1a == 1 then (p1, xp1a)
    else if xp2a == 1 then (p2, xp2a)
    else (0, 0)
  ) (0, 0) <| map (\p -> (p, tableu[p][a])) (n...2 * n - 1) -- both inclusive 
  in if xpa == 1
     then -- Call rowsum for all i in {1...2n}
          let filtered_is = filter (\i -> i != p && tableu[i][a] == 1) (iota (2 * n))
          let (is1, vs1) =
            map (\i ->
                   rowsum tableu i p)
                (filtered_is)
            |> unzip
          let tmp1 = scatter_2d tableu (is1 |> flatten) (vs1 |> flatten)
          -- set (p - n) row equal to pth row
          let is2 = map (\i -> (p - n, i)) (iota (size n))
          let vs2 = map (\i -> tmp1[p][i]) (iota (size n))
          let tmp2 = scatter_2d tmp1 is2 vs2
          -- set pth row 0 except rp is 0 or 1 50/50 and zpa = 1
          let (eng1, rand_val) = rand_i8.rand (0, 1i8) eng
          let is3 = map (\i -> (p, i)) (iota (size n))
          let vs3 =
            map (\i ->
                   -- zpa
                   if i == (n + a)
                   then 1
                   else -- rp
                   if i == (2 * n)
                   then rand_val
                   else 0)
                (iota (size n))
          let tmp3 = scatter_2d tmp2 is3 vs3
          in (eng1, tmp3, rand_val)
     else -- set (2n + 1) st row to be identically 0
          let max_idx = (size n) - 1
          let is1 = map (\i -> (max_idx, i)) (iota (size n))
          let vs1 = replicate (size (n)) 0
          let tmp1 = scatter_2d tableu is1 vs1
          -- call rowsum 2n+1, i+n
          let filtered_is = filter (\i -> tmp1[i][a] == 1) (iota n)
          let (is2, vs2) = map (\i -> rowsum tmp1 (max_idx) (i + n)) (filtered_is) |> unzip
          let tmp2 = scatter_2d tmp1 (is2 |> flatten) (vs2 |> flatten)
          in (eng, tmp2, tmp2[max_idx][max_idx])
