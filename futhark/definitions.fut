import "lib/github.com/diku-dk/cpprandom/random"

module rng_engine = minstd_rand
module rand_i8 = uniform_int_distribution i8 u32 rng_engine
type t = i8
type tab [n] = [2 * (n + 1)][2 * (n + 1)]t

local def size (n: i64) : i64 = 2 * (n + 1)

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
  let tot_sum = map (\j -> g tableu[i][j] tableu[i][j + n] tableu[h][j] tableu[h][j + n]) (iota n) |> reduce (+) 0
  let res = (2 * tableu[h][size (n) - 1] + 2 * tableu[i][size (n) - 1] + tot_sum) % 4
  let v = if res == 0 then 0 else 1
  let is = map (\j -> [(h, j), (h, j + n)]) (iota n) |> flatten
  let vs = map (\j -> [tableu[i][j] ^ tableu[h][j], tableu[i][n + j] ^ tableu[h][n + j]]) (iota n) |> flatten
  in (is ++ [(h, size (n))], vs ++ [v])

def initial_tableu (n: i64) : tab [n] =
  let tmp = 2 * (n + 1)
  in tabulate_2d tmp tmp (\i j ->
                            if j == n - 1
                            then -- r_i
                                 0
                            else if i == j
                            then -- diagonal
                                 1
                            else 0)

def CNOT [n] (tableu: *tab [n]) (a: i64) (b: i64) : *tab [n] =
  let tmp = 2 * n
  let indices =
    map (\i ->
           [(i, tmp - 1), (i, b), (i, n + a)])
        (iota tmp)
    |> flatten
  let values =
    map (\i ->
           let zia = tableu[i][n + a]
           let zib = tableu[i][n + b]
           let xia = tableu[i][a]
           let xib = tableu[i][b]
           let ri = tableu[i][tmp - 1]
           in [ri ^ xia * zib * (xib ^ zia ^ 1i8), xib ^ xia, zia ^ zib])
        (iota tmp)
    |> flatten
  in scatter_2d tableu indices values

def Hadamard [n] (tableu: *tab [n]) (a: i64) : *tab [n] =
  let tmp = 2 * n
  let indices = map (\i -> [(i, tmp - 1), (i, a), (i, n + a)]) (iota tmp) |> flatten
  let values =
    map (\i ->
           let ri = tableu[i][tmp - 1]
           let xia = tableu[i][a]
           let zia = tableu[i][a + n]
           in [ri ^ xia * zia, zia, xia])
        -- swap in regards to indices, which are ri, xia, zia
        (iota tmp)
    |> flatten
  in scatter_2d tableu indices values

def Phase [n] (tableu: *tab [n]) (a: i64) : *tab [n] =
  let tmp = 2 * n
  let indices = map (\i -> [(i, tmp - 1), (i, a + n)]) (iota tmp) |> flatten
  let values =
    map (\i ->
           let ri = tableu[i][tmp - 1]
           let xia = tableu[i][a]
           let zia = tableu[i][n + a]
           in [ri ^ xia * zia, zia ^ xia])
        (iota tmp)
    |> flatten
  in scatter_2d tableu indices values

def Measurement [n] (eng: rng_engine.rng) (tableu: *tab [n]) (a: i64) : (rng_engine.rng, *tab [n], t) =
  let p =
    reduce_comm (\p1 p2 ->
                   -- 0 used as neutral element, as we check n+1...2n, so if p is 0, we know that none is found that is equal to 0
                   let xp1a = if p1 != 0 then tableu[p1][a] else 0
                   let xp2a = if p2 != 0 then tableu[p2][a] else 0
                   in if xp1a == 1 && xp1a == 1
                      then if p1 <= p2
                           then p1
                           else p2
                      else if xp1a == 1
                      then p1
                      else if xp2a == 1
                      then p2
                      else 0)
                (0)
                (n + 1...2 * n)
  in if p != 0
     then -- Call rowsum for all i in {1...2n}
          let filtered_is = filter (\i -> i != p && tableu[i][a] == 1) (iota (2 * n))
          let (is1, vs1) =
            map (\i ->
                   rowsum tableu i p)
                (filtered_is)
            |> unzip
          let tmp1 = scatter_2d tableu (is1 |> flatten) (vs1 |> flatten)
          -- set (p - n) row equal to pth row
          let is2 = map (\i -> (p - n, i)) (iota (size (n)))
          let vs2 = map (\i -> tmp1[p][i]) (iota (size (n)))
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
     else -- set (2n +1) st row to be identically 0
          let is1 = map (\i -> (size (n), i)) (iota (size (n)))
          let vs1 = replicate (size (n)) 0
          let tmp1 = scatter_2d tableu is1 vs1
          -- set call rowsum 2n+1, i+n
          let filtered_is = filter (\i -> tmp1[i][a] == 1) (iota n)
          let (is2, vs2) = map (\i -> rowsum tmp1 (size n) (i + n)) (filtered_is) |> unzip
          let tmp2 = scatter_2d tmp1 (is2 |> flatten) (vs2 |> flatten)
          in (eng, tmp2, tmp2[size (n) - 1][size (n) - 1])
