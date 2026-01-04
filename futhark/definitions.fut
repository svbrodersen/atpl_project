module Quantum = {
  type t = i8
  type tab [n] = [2 * (n + 1)][2 * (n+1)]t

  def initial_tableu (n: i64) : tab[n] =  
  let tmp = 2 * (n + 1) 
  in tabulate_2d tmp tmp (\i j -> 
    if j == n-1 then 0 -- r_i
    else if i == j then 1 -- diagonal
    else 0 -- all other
  )

  def CNOT [n] (tab: *tab[n]) (a: i64) (b: i64) : *tab[n] = 
    let tmp = 2 * (n + 1)
    let indices = map (\i -> 
      [(i, tmp-1), (i, b), (i, n + a)]
    ) (iota tmp) |> flatten
    let values = map (\i ->
      let zia = tab[i][n + a]
      let zib = tab[i][n + b]
      let xia = tab[i][a]
      let xib = tab[i][b]
      let ri = tab[i][tmp-1]
      in [ri ^ xia * zib * (xib ^ zia ^ 1i8), xib ^ xia, zia ^ zib]
    ) (iota tmp) |> flatten
    in scatter_2d tab indices values

  def Hadamard [n] (tab: *tab[n]) (a: i64) : *tab[n] = 
    let tmp = 2 * (n + 1) 
    let indices = map (\i -> [(i, tmp-1), (i, a), (i, n + a)]) (iota tmp) |> flatten
    let values = map (\i -> 
      let ri = tab[i][tmp-1]
      let xia = tab[i][a]
      let zia = tab[i][a + n]
      in [ri ^ xia * zia, zia, xia] -- swap in regards to indices, which are ri, xia, zia
    ) (iota tmp) |> flatten
    in scatter_2d tab indices values
}
