module Quantum = {
  type t = i8
  type tab [n] = [2 * (n + 1)][2 * (n + 1)]t

  def initial_tableu (n: i64) : tab [n] =
    let tmp = 2 * (n + 1)
    in tabulate_2d tmp tmp (\i j ->
                              if j == n - 1
                              then 0
                              else -- r_i
                              if i == j
                              then 1
                              else -- diagonal
                                   0)

  -- all other

  def CNOT [n] (table: *tab [n]) (a: i64) (b: i64) : *tab [n] =
    let tmp = 2 * (n + 1)
    let indices =
      map (\i ->
             [(i, tmp - 1), (i, b), (i, n + a)])
          (iota tmp)
      |> flatten
    let values =
      map (\i ->
             let zia = table[i][n + a]
             let zib = table[i][n + b]
             let xia = table[i][a]
             let xib = table[i][b]
             let ri = table[i][tmp - 1]
             in [ri ^ xia * zib * (xib ^ zia ^ 1i8), xib ^ xia, zia ^ zib])
          (iota tmp)
      |> flatten
    in scatter_2d table indices values

  def Hadamard [n] (table: *tab [n]) (a: i64) : *tab [n] =
    let tmp = 2 * (n + 1)
    let indices = map (\i -> [(i, tmp - 1), (i, a), (i, n + a)]) (iota tmp) |> flatten
    let values =
      map (\i ->
             let ri = table[i][tmp - 1]
             let xia = table[i][a]
             let zia = table[i][a + n]
             in [ri ^ xia * zia, zia, xia])
          -- swap in regards to indices, which are ri, xia, zia
          (iota tmp)
      |> flatten
    in scatter_2d table indices values

  def Phase [n] (table: *tab [n]) (a: i64) : *tab [n] =
    let tmp = 2 * (n + 1)
    let indices = map (\i -> [(i, tmp - 1), (i, a + n)]) (iota tmp) |> flatten
    let values =
      map (\i ->
             let ri = table[i][tmp - 1]
             let xia = table[i][a]
             let zia = table[i][n + a]
             in [ri ^ xia * zia, zia ^ xia])
          (iota tmp)
      |> flatten
    in scatter_2d table indices values
}
