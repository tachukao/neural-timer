open Owl
open Printf

let total_time = 3300E-3
let n_steps = 3300
let dt = 1E-3
let tau = 10E-3
let n = 200

let t_ready = 0.1
let ts_set = [| 0.24; 0.28; 0.32; 0.39; 0.45; 0.33; 0.41; 0.49; 0.57; 0.65 |]

let n_ts = Array.length ts_set 

let input_type = Cmdargs.(get_string "-input" |> force ~usage:"-input")


let tonic_input g ts t = 
  let inp = 
    if g = 1.5 then 0.4 
    else if g = 1.0 then 0.3
    else 0.35 in
  if (t >= t_ready) && (t < (ts +. t_ready)) then inp
  else 0.

let transient_input g ts t = 
  let inp = 
    if g = 1.5 then 0.4 
    else if g = 1.0 then 0.3
    else 0.35 in
  if (t >= t_ready) && (t < (ts +. t_ready)) then inp
  else 0.


let set_input t ts = 
  if ( (t >= t_ready) && (t < (t_ready +. 0.02)) ) || 
     ( (t >= (ts +. t_ready) ) && (t < (ts +. t_ready +. 0.02)))
  then 0.4
  else 0.

let target_output a t g ts = 
  let alpha = 1. /. Maths.(log (1. +. 1. /. a)) in
  let tt = ts *. g in
  let target_start = t_ready +. ts in
  if t < target_start then 0. else a *. Maths.( (exp ((t -. target_start) /. alpha /. tt ) ) -. 1.)

let context_input = if input_type = "tonic" then tonic_input else if input_type = "transient" then transient_input else (assert false)

let context_inputs g = Mat.init_2d n_ts n_steps 
    (fun ts_idx t -> 
       let t = (float t) *. dt in
       let ts = ts_set.(ts_idx) in
       context_input t g ts
    )

let target_outputs a g = 
  Mat.init_2d n_ts n_steps  
    (fun ts_idx t -> 
       let t = (float t) *. dt in
       let ts = ts_set.(ts_idx) in
       target_output a t g ts 
    )

let set_inputs = Mat.init_2d n_ts n_steps
    (fun ts_idx t -> 
       let t = (float t) *. dt in
       let ts = ts_set.(ts_idx) in
       set_input t ts
    )


let unpack prms = 
  let x0 = Mat.get_slice [ [ ]; [0] ] prms  in
  let w0 = Mat.get_slice [ [ ]; [1] ] prms  in
  let bc = Mat.get_slice [ [ ]; [2] ] prms in
  let bs = Mat.get_slice [ [ ]; [3] ] prms  in
  let j = Mat.get_slice [ [ ]; [4; Pervasives.(n + 3); ] ] prms  in
  let cx = Mat.get_slice [ [ ]; [Pervasives.(n+4)] ] prms  in
  let cz = Mat.get prms 0 Pervasives.(n + 5) in
  x0, w0, bc, bs, j, cx, cz

let simulate ts g prms = 
  let x0, w0, bc, bs, j, cx, cz = unpack prms in
  let t_target = t_ready +. ts in 
  let tt = g *. ts in
  let uc_noise = Mat.gaussian ~sigma:5E-3 1 n_steps |> Mat.to_array in
  let tau = 10E-3  in
  let dt = 1E-3  in
  let rec iterate x xacc zacc step = 
    let t = (float step) *. 1E-3 in
    if t < (t_target +. tt) then begin
      let us = set_input t ts in
      let uc = context_input t g ts in
      let dx = Mat.( (j *@ (tanh x) - x + (bc *$ Pervasives.(uc +. uc_noise.(step))) + bs *$ us  + cx) /$ tau ) in
      let next_x = Mat.(x + dx *$ dt) in
      let z = Mat.( sum' ( (transpose w0) *@ (tanh x) +$ cz ) ) in
      iterate next_x (next_x::xacc) (z::zacc) (succ step) 
    end
    else 
      let xacc = xacc |> List.rev |> Array.of_list |> Mat.concatenate ~axis:1 in
      let zacc = [| zacc |> List.rev |> Array.of_list |] |> Mat.of_arrays in
      xacc, zacc
  in iterate x0 [x0] [] 1






