module Quantum = {
  type t = i8

  def initial_tableu (n: i64) : [2 * (n + 1)][2 * (n+1)]t =  -- is there some way we can define this 2 * (n-1)... for the module itself instead??
  let tmp = 2 * (n + 1) 
  in tabulate_2d tmp tmp (\i j -> 
    if j == n-1 then 0 -- r_i
    else if i == j then 1 -- diagonal
    else 0 -- all other
  )

  def CNOT [n] (tab: *[2 * (n + 1)][2 * (n + 1)]t) (a: i64) (b: i64) : *[2 * (n + 1)][2 * (n + 1)]t = 
    let len = 2 * (n - 1)
    let indices = map (\i -> 
      [(i, len-1), (i, b), (i, n + a)]
    ) (iota len) |> flatten
    let values = map (\i ->
      let zia = tab[i][n + a]
      let zib = tab[i][n + b]
      let xia = tab[i][a]
      let xib = tab[i][b]
      let ri = tab[i][len-1]
      in [ri ^ xia * zib * (xib ^ zia ^ 1i8), xib ^ xia, zia ^ zib]
    ) (iota len) |> flatten
    in scatter_2d tab indices values

}
