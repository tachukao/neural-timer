open Owl
open Printf
open Defaults

(* tau dx/dt = -x + Jr + Bu + cx + px *)
(* r = tanh (x) *)
(* x0 ~ uniform (-1, 1) *)
(* J  ~ gaussian (0, 1/N) *)
(* bc  ~ uniform (-1, 1) *)
(* bs  ~ uniform (-1, 1) *)
(* w0  zero in the beginning *)
(* cz  zero in the beginning *)
 
let _ = Printexc.record_backtrace true
let input_type = Cmdargs.(get_string "-input" |> force ~usage:"-input")
let a = Cmdargs.(get_float "-a" |> force ~usage:"-a")

let in_dir s = Dir.in_dir (sprintf "%s/%s" input_type s)

let initialize_prms () = 
  let x0 = Mat.uniform ~a:(-1.) ~b:1. n 1 in
  let w0 = Mat.zeros n 1 in
  let bc = Mat.uniform ~a:(-1.) ~b:1. n 1 in
  let bs = Mat.uniform ~a:(-1.) ~b:1. n 1 in
  let j = Mat.((gaussian n n) /$ (Maths.sqrt (float n)) ) in
  let cx = Mat.uniform ~a:(-1.) ~b:1. n 1 in
  let cz = Mat.zeros n 1 in
  Mat.concatenate ~axis:1 [| x0; w0; bc; bs; j; cx; cz |]

let n_bc_prms = n
let n_bs_prms = n
let n_j_prms = n * n
let n_x0_prms = n * 1
let n_w0_prms = n * 1
let n_cx_prms = n
let n_cz_prms = n

let n_total_prms = n_bc_prms + n_bs_prms + n_j_prms + n_x0_prms + n_w0_prms + n_cz_prms + n_cx_prms

let forward a ts g x0 w0 bc bs j cx cz = 
  let t_target = t_ready +. ts in 
  let tt = g *. ts in
  (* let uc_noise = Mat.gaussian ~sigma:5E-3 1 n_steps |> Mat.to_array |> Array.map (fun x -> Algodiff.D.(F x)) in *)
  let uc_noise = Mat.zeros 1 n_steps |> Mat.to_array |> Array.map (fun x -> Algodiff.D.(F x)) in 
  let open Algodiff.D in
  let tau = F 10E-3  in
  let dt = F 1E-3  in
  let rec iterate x acc step = 
    let t = (float step) *. 1E-3 in
    if t < (t_target +. tt) then begin
      let us = F (Defaults.set_input t ts) in
      let uc = F (Defaults.context_input t g ts) in
      let dx = Maths.( (j *@ (tanh x) - x + bc * (uc + uc_noise.(step)) + bs * us  + cx) / tau ) in
      let next_x = Maths.(x + dx * dt) in
      let z = Maths.( sum' ( (transpose w0) *@ (tanh x) + cz ) ) in
      let acc = if (t >= t_target) && (t < (t_target +. tt)) then  
          begin 
            let target = F Defaults.(target_output a t g ts) in
            (* printf "%f, %f, %f\n" t (unpack_flt target) (unpack_flt z);
               flush_all ();*)
            Maths.(  acc + (sqr (z - target)) / (F tt) / (F 1000.)  ) 
          end
        else acc in
      iterate next_x acc (succ step) 
    end
    else acc
  in iterate x0 (F 0.) 1

let unpack prms = 
  let open Algodiff.D in
  let x0 = Algodiff.D.Maths.get_slice [ [ ]; [0] ] prms  in
  let w0 = Algodiff.D.Maths.get_slice [ [ ]; [1] ] prms  in
  let bc = Algodiff.D.Maths.get_slice [ [ ]; [2] ] prms in
  let bs = Algodiff.D.Maths.get_slice [ [ ]; [3] ] prms  in
  let j = Algodiff.D.Maths.get_slice [ [ ]; [4; Pervasives.(n + 3); ] ] prms  in
  let cx = Algodiff.D.Maths.get_slice [ [ ]; [Pervasives.(n+4)] ] prms  in
  let cz = Algodiff.D.Maths.get_item prms 0 Pervasives.(n + 5) in
  x0, w0, bc, bs, j, cx, cz


let cost a prms = 
  let x0, w0, bc, bs, j, cx, cz = unpack prms in
  let repetitions = 1 in
  let open Algodiff.D in
  (* let regulariser = Maths.( (l2norm_sqr' w0) + (l2norm_sqr' cx) + (l2norm_sqr' bc) + (l2norm_sqr' bs) ) in *)
  (* printf "regulariser %f\n" (regulariser |> unpack_flt);*)
  let sum_cost_for g = 
    let open Algodiff.D in
    let ts_set = [|0.5; 0.55; 0.6; 0.65; 0.7; 0.75; 0.8; 0.85; 0.9; 0.95; 1.0|] in
    (* let ts_set = Array.concat [ts_set; Array.init 1 (fun _ -> Stats.uniform_rvs ~a:0.5 ~b:1.0 )] in  *)
    Array.map (fun ts -> 
        Array.init repetitions (fun _ -> forward a ts g x0 w0 bc bs j cx cz)
        |> Array.fold_left Maths.(+) (F 0.)
      ) ts_set
    |> Array.fold_left Maths.(+)  (F 0.) in
  Maths.(  (sum_cost_for 1.0) + (sum_cost_for 1.25) + (sum_cost_for 1.5)  ) (* + (F 0.0001 * regulariser) ) *)

let dcost a = Algodiff.D.grad (fun prms -> cost a prms)

let for_prms prms = 
  let open Bigarray in 
  let z = genarray_of_array1 prms in
  (* printf "num total params %i\n" n_total_prms; *)
  genarray_of_array2 (reshape_2 z n (n+6) ) 

let f_df a prms g = 
  let open Bigarray in
  let open Algodiff.D in
  let prms = Arr (for_prms prms) in
  let cost = cost a prms in 
  let dprms = (dcost a) prms in 
  Bigarray.Genarray.blit (unpack_arr (dprms)) (for_prms g);
  unpack_flt cost


let prms = initialize_prms () 
           |> (fun x -> Mat.reshape x [| n_total_prms; 1|])
           |> (fun x -> Bigarray.reshape_1 x n_total_prms) 


let run_optimisation i a =
  let in_dir s = in_dir (sprintf "%s_a_%i" s (int_of_float a) ) in
  let stop max_iter k c =
    printf "input %s | a %f | iter %5i | cost = %.8f%!\n" input_type a k c;
    Gc.minor ();
    Mat.save_txt (for_prms prms) (in_dir "prms");
    flush_all ();
    k >= max_iter in
  let max_iter = 2000 in
  let open Lbfgs in
  Gc.minor () ;
  let stop st = stop max_iter (iter st) (previous_f st) in
  C.min ~print:(No) ~pgtol:0. ~factr:1E1 ~corrections:20 ~stop (f_df a) prms |> ignore


let _ = 
  let a_list = [a]  in
  List.iteri run_optimisation a_list


(* 
(*  vanilla gradient descrent *)
let learn prms0 = 
  printf "start training\n";
  flush_all ();
  let open Algodiff.D in
  let beta = F 1E-2 in
  let rec loop step old_c prms = 
    if step < 200 then
      let c = cost prms
              |> unpack_flt in
      let dprms  = dcost prms in
      let prms = Maths.(prms - beta * dprms) in
      printf "costs %.9f\n" c;
      flush_all ();
      Owl.Mat.save_txt (unpack_arr prms) (in_dir "prms");
      loop (succ step) c prms
    else prms 
  in loop 0 0. prms0
let () = 
  let prms0 = Algodiff.D.(Arr (initialize_prms () ))  in
  let prms = learn prms0 in
  printf "done!\n"
*)

