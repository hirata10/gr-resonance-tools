"""
Resonance-aware evolution of the actions J and angles psi of both bodies (jump prescription).

Both bodies are evolved together, one RK4 step at a time (td_res_exe/J_psi_step), in physical time t:

    dJ/dt   = J_dot_sf(J)      self-force rate of each body (the right-hand side of step 4)
    dpsi/dt = Omega(J)         fundamental frequencies of each body

After each step the resonance function of every mode pair,

    f = N_outer . Omega_outer - N_inner . Omega_inner,      N = (n, k, m),  m_inner = m_outer,

is compared with its value at the start of the step. When a pair has crossed f = 0 the step is repeated so that it ends at
the crossing; the resonant angle is read from the angles there, both bodies get the tidal jump Delta J of the step 6
executable (Eq. (38) of arXiv:2207.07733), and the evolution continues from that time with the kicked actions.

Angles and the resonant angle:
    psi = (psi_r, psi_theta, psi_phi) of a body is measured from the torus origin of CKerr_TorusOrigin (kerrtraj.c):
    psi = 0 is the point of its orbit at pericentre (u_r = 0), crossing the equator with u_theta > 0, at phi = 0.
    The resonant angle of a mode pair is
        Theta = N_outer . psi_outer - N_inner . psi_inner        (Eq. (5) of arXiv:2207.07733)
    and theta_res, the phase that enters Delta J, is Theta at the crossing, where Theta = theta_res + Gamma (t - t_res)^2 / 2.
    Only the frequencies drive the angles: the order-mu corrections of each body's own self-force to dpsi/dt are not included.
    At a crossing the angles are continuous (the angle jump of the resonant Hamiltonian is not included); afterwards they
    advance with the frequencies of the kicked actions.

Usage (from the repository root, after "make td_res" and the Cython build):

    python src/python/J_evolve_resonant.py <label> <t_max> --psi-inner a b c --psi-outer a b c
                                           [--mode-max 2] [--omp-threads 1] [--max-rows N] [--no-apply] [--theta-res X]

    label         system label: first column of action_angle_pairs.txt (label, inner J_r J_theta J_phi, outer J_r J_theta J_phi)
    t_max         stop after the first accepted state with t >= t_max (physical time, units of M)
    --psi-inner   the three angles of the inner body at t0 (radians); --psi-outer the same for the outer body
    --mode-max    every integer of the mode vectors lies within +-mode_max; 0 = no resonances (plain self-force evolution)
    --omp-threads OpenMP threads of the RK4 step and of the Delta J runs
    --max-rows    stop when the trajectory file has this many rows after row 0 (default: n_time_inner - 1 of pars.txt)
    --no-apply    compute and log the jumps but do not add them to J
    --theta-res   use this resonant phase for every crossing instead of the one from the angles (to reproduce a stored row)

What is used:
    td_res_exe/J_psi_step       one RK4 step of (J, psi) of both bodies; its last argument dt_max shortens the step
    td_res_exe/Delta_J_single   Delta J of both bodies at one crossing; pars.txt must have Gamma_source 0 (Gamma supplied)
    gr_wrapper                  (I, r_p, r_a) at a crossing and Omega of a kicked state
The mass ratios, the self-force sum, the tidal sum, the spin and the grids are the ones compiled in from pars.txt.

How a crossing is handled:
    1. Trial step with the time step of step 4 (0.01 min |J/J_dot|). f_0 at its start, f_1 at its end.
    2. A pair has crossed when f ends the step on the other side from its stored sign and on the side it is drifting toward.
       The stored sign is the side the pair started on; it changes only when the pair is kicked. A shift of f caused by a kick
       (the pair's own or another pair's) therefore neither triggers a second kick nor hides a crossing.
    3. First estimate of the crossing, linear in the step: t_c = t + dt f_0 / (f_0 - f_1), for the earliest pair.
    4. The step is repeated from the same start with dt_max = t_c - t. Its end state is the crossing state: J and psi of both
       bodies at t_c. There f = f_c, small but not zero.
    5. Gamma is the slope at t_c of the parabola through (t, f_0), (t_c, f_c), (t + dt, f_1), and t_res = t_c - f_c / Gamma.
       If |t_res - t_c| > 0.01 dt for the earliest pair the step is repeated once more, ending at t_res.
    6. theta_res = Theta(t_c) - f_c^2 / (2 Gamma), reduced to [0, 2 pi): Theta is stationary at the crossing, so the small
       offset t_res - t_c enters only at second order and this term removes it.
    7. Every pair whose crossing lies within 0.01 dt of t_c is kicked at this state. Delta J of both bodies is added to the
       actions; the angles are not changed. The next trial step starts here.
    Logged but not applied: a jump that is not finite, or a crossing with tau = sqrt(2 pi / |Gamma|) > 10 dt.
    A pair that another kick pushed across f = 0 cannot be placed at its own crossing; it is kicked where it is found, with
    Theta as it is (no correction), and flagged.

Output (outputs_data/):
    J_evolve_res_<label>.txt         one row per accepted state, 22 columns:
        0 row   1 kind (0 = end of an ordinary step, 1 = crossing state, before its jumps)   2 t   3 dt
        4-6 J_inner   7-9 psi_inner   10-12 Omega_inner   13-15 J_outer   16-18 psi_outer   19-21 Omega_outer
    J_evolve_res_kicks_<label>.txt   one row per crossing, 57 columns:
        0 row of the crossing state   1 applied (1/0)   2 flag (0, or 1 = Theta used without the correction of point 6)
        3 t_c   4 t_res - t_c   5 Gamma   6 tau/dt   7 f_c (outer - inner frequency at the crossing state)
        8-10 psi_inner   11-13 psi_outer (at t_c, not reduced)   14 Theta(t_c)   15 the correction -f_c^2/(2 Gamma)
        16 theta_res passed to Delta_J_single   17-19, 20-22 change of Omega_inner, Omega_outer by all the jumps of this state
        23-56 the 34 columns of the Delta_J_single row (README, step 6): 23 crossing number, 24 system label,
        25-29 n_in k_in n_out k_out m, 30 theta_res, 31 Gamma, 32-37 Omega_in Omega_out, 38-43 J_in J_out,
        44-46 Delta_J_in, 47-49 Delta_J_in / J_in, 50-52 Delta_J_out, 53-55 Delta_J_out / J_out, 56 Gamma source
    The state after row i is row i plus the applied jumps logged with row i. Running the same command again continues the run.
"""
import argparse
import copy
import itertools
import math
import os
import subprocess
import sys

import numpy as np

import globalpars
from gr_wrapper import ckerr_eql2j, ckerr_j2eql, ckerr_minverse, ckerr_minv2omega

parser = argparse.ArgumentParser()
parser.add_argument("label", type=int)
parser.add_argument("t_max", type=float)
parser.add_argument("--psi-inner", type=float, nargs=3, required=True)
parser.add_argument("--psi-outer", type=float, nargs=3, required=True)
parser.add_argument("--mode-max", type=int, default=2)
parser.add_argument("--omp-threads", type=int, default=1)
parser.add_argument("--max-rows", type=int, default=globalpars.GLOBALPAR_n_time_inner - 1)
parser.add_argument("--no-apply", action="store_true")
parser.add_argument("--theta-res", type=float, default=None)
args = parser.parse_args()

M = globalpars.GLOBALPAR_M
astar = globalpars.GLOBALPAR_astar
mu_inner = globalpars.GLOBALPAR_mu_inner
mu_outer = globalpars.GLOBALPAR_mu_outer
label = args.label
exe_dir = "td_res_exe"
output_dir = "outputs_data"
traj_file = f"{output_dir}/J_evolve_res_{label}.txt"
kick_file = f"{output_dir}/J_evolve_res_kicks_{label}.txt"
TAU_MAX = 10.0   # a crossing that lasts longer than this many steps is logged, not applied
TOL = 0.01       # crossings within this fraction of a step of the crossing state are kicked there
MAX_REPEATS = 4  # most repetitions of a step to end it at a crossing
KICK_COL = 23    # first column of the Delta_J_single row in the kick log
settings = (f"# settings: label {label} mode_max {args.mode_max} apply {int(not args.no_apply)} theta_res {args.theta_res}"
            f" psi_inner {args.psi_inner} psi_outer {args.psi_outer}")
env = dict(os.environ, OMP_NUM_THREADS=str(args.omp_threads))


def fmt(x):
    return "%.17g" % x


def mode_list(mode_max):
    """Mode pairs (n_in, k_in, n_out, k_out, m) with every integer within +-mode_max, by the rules of the resonance finder:
    k_in - k_out even (selection rule, Eq. (34) of arXiv:2207.07733); one of N and -N (first nonzero of m, n_in, n_out, k_in, k_out
    positive); gcd 1, or gcd 2 when (k_in - k_out)/2 is odd (then the vector with gcd 1 is forbidden by the selection rule and this
    is its first allowed harmonic). The higher harmonics of each pair are summed inside the tidal calculation (N_res).
    A pair in which the vector of either body is all zero is left out: its frequency is zero, which the jump formula does not
    cover (1/omega^3 in the amplitudes; the paper sets the omega = 0 terms aside). Its Delta J need not vanish."""
    modes = []
    r = range(-mode_max, mode_max + 1)
    for n_in, k_in, n_out, k_out, m in itertools.product(r, repeat=5):
        if (k_in - k_out) % 2 != 0:
            continue
        if (n_in, k_in, m) == (0, 0, 0) or (n_out, k_out, m) == (0, 0, 0):
            continue
        if next(x for x in (m, n_in, n_out, k_in, k_out) if x != 0) < 0:
            continue
        gd = math.gcd(n_in, k_in, n_out, k_out, m)
        if gd == 1 or (gd == 2 and ((k_in - k_out) // 2) % 2 == 1):
            modes.append((n_in, k_in, n_out, k_out, m))
    return np.array(modes, dtype=int).reshape(-1, 5)


def combine(modes, inner, outer):
    """N_outer . x_outer - N_inner . x_inner of every mode pair, for x = Omega (the resonance function f) or x = psi (Theta)"""
    n_in, k_in, n_out, k_out, m = modes.T
    return (n_out * outer[0] + k_out * outer[1] + m * outer[2]) - (n_in * inner[0] + k_in * inner[1] + m * inner[2])


def omega(J):
    """(Omega_r, Omega_theta, Omega_phi) of the orbit with actions J, or None when there is no bound orbit"""
    flag, Minv = ckerr_minverse([float(x) for x in J], M, astar)
    return ckerr_minv2omega(Minv) if flag == 1 else None


def ancillary(J):
    """(I, r_p, r_a) of the orbit with actions J, or None when there is no bound orbit"""
    flag, EQL = ckerr_j2eql([float(x) for x in J], M, astar)
    return ckerr_eql2j(list(EQL), M, astar)[2] if flag == 1 else None


def numbers(line):
    """The entries of a line if all of them are numbers (kept as strings), otherwise an empty list"""
    parts = line.split()
    try:
        [float(x) for x in parts]
    except ValueError:
        return []
    return parts


def append(path, text):
    with open(path, "a") as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())


class State:
    """t, dt of the step that led here, J and psi of both bodies (y), their frequencies (Omega)"""
    def __init__(self, values):
        v = np.array(values, dtype=float)  # t dt J_in(3) psi_in(3) Omega_in(3) J_out(3) psi_out(3) Omega_out(3)
        self.t, self.dt = v[0], v[1]
        self.y = np.concatenate([v[2:8], v[11:17]])        # J_in, psi_in, J_out, psi_out
        self.Omega = np.concatenate([v[8:11], v[17:20]])   # Omega_in, Omega_out
        self.text = " ".join(fmt(x) for x in v)

    def f(self):
        return combine(modes, self.Omega[0:3], self.Omega[3:6])


def rk4_step(state, dt_max):
    """One RK4 step of both bodies from state, not longer than dt_max if dt_max > 0. Returns the new State, or None and what was printed"""
    cmd = [f"{exe_dir}/J_psi_step", fmt(state.t)] + [fmt(x) for x in state.y] + [fmt(dt_max)]
    out = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, env=env).stdout
    lines = [line for line in out.splitlines() if line.strip()]
    if len(lines) == 1 and len(numbers(lines[0])) == 20:
        return State(lines[0].split()), ""
    return None, " | ".join(lines)


def slope_at_cut(f_0, f_c, f_1, h_c, h):
    """Slope at h_c of the parabola through (0, f_0), (h_c, f_c), (h, f_1); the chord (f_1 - f_0)/h if h_c is at an end of the step"""
    if h_c < 1e-6 * h or h - h_c < 1e-6 * h:
        return (f_1 - f_0) / h
    return f_0 * (h_c - h) / (h_c * h) + f_c * (2.0 * h_c - h) / (h_c * (h_c - h)) + f_1 * h_c / (h * (h - h_c))


def delta_J(crossings, number):
    """Run Delta_J_single for the crossings of one state, at the same time. Fills c["row"] (34 strings, or None on failure)"""
    child_env = dict(os.environ, OMP_NUM_THREADS=str(max(1, args.omp_threads // max(1, len(crossings)))))
    procs = []
    for c in crossings:
        number += 1
        n_in, k_in, n_out, k_out, m = c["mode"]
        J_in, J_out, Om = c["state"].y[0:3], c["state"].y[6:9], c["state"].Omega
        # argument order of Delta_J_single (calling.c): modes, (ra, rp, I) of each body, Gamma, theta_res, mu_outer, labels, Omegas, Js, mu_inner
        cmd = [f"{exe_dir}/Delta_J_single"] + [str(x) for x in (n_in, n_out, k_in, k_out, m, m)] \
            + [fmt(x) for x in (c["anc_in"][2], c["anc_in"][1], c["anc_in"][0], c["anc_out"][2], c["anc_out"][1], c["anc_out"][0], c["Gamma"], c["theta_res"], mu_outer)] \
            + [str(label), str(number)] + [fmt(x) for x in (*Om, *J_in, *J_out, mu_inner)]
        procs.append(subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, env=child_env))
    for c, proc in zip(crossings, procs):
        lines = [line for line in proc.communicate()[0].splitlines() if line.strip()]
        row = numbers(lines[-1]) if lines else []
        if len(row) == 34 and int(row[33]) != 0:
            sys.exit("Delta_J_single was built with Gamma_source != 0: set Gamma_source 0 in pars.txt and run make clean-td_res, make td_res")
        c["row"] = row if len(row) == 34 else None
    return number


# ---------------------------------------------------------------- set-up
os.makedirs(output_dir, exist_ok=True)
J_data = np.loadtxt("action_angle_pairs.txt", delimiter=" ", ndmin=2)
J_initial = J_data[J_data[:, 0] == label]
if len(J_initial) != 1:
    sys.exit(f"action_angle_pairs.txt has {len(J_initial)} rows with label {label}")
J_initial = J_initial[0]
modes = mode_list(args.mode_max)
print(f"{len(modes)} mode pairs with every integer within +-{args.mode_max}", flush=True)

header = ("# angles psi = (psi_r, psi_theta, psi_phi): zero at pericentre, at the equatorial crossing with u_theta > 0, at phi = 0 (CKerr_TorusOrigin)\n"
          "# resonant angle Theta = N_outer.psi_outer - N_inner.psi_inner; theta_res = Theta at the crossing\n")
if not os.path.exists(traj_file):
    Om_in, Om_out = omega(J_initial[1:4]), omega(J_initial[4:7])
    if Om_in is None or Om_out is None:
        sys.exit("The initial actions are not bound orbits")
    row_0 = [globalpars.GLOBALPAR_t0, 0.0, *J_initial[1:4], *args.psi_inner, *Om_in, *J_initial[4:7], *args.psi_outer, *Om_out]
    append(traj_file, settings + "\n" + header
           + "# row kind(0 = end of a step, 1 = crossing state before its jumps) t dt J_inner(3) psi_inner(3) Omega_inner(3) J_outer(3) psi_outer(3) Omega_outer(3)\n"
           + "0 0 " + " ".join(fmt(x) for x in row_0) + "\n")
if not os.path.exists(kick_file):
    append(kick_file, settings + "\n" + header
           + "# row applied flag t_c t_res-t_c Gamma tau/dt f_c psi_inner(3) psi_outer(3) Theta(t_c) correction theta_res dOmega_inner(3) dOmega_outer(3), then the 34 columns of the Delta_J_single row\n")

# Continue from the last row of the trajectory file (row 0 for a new run)
with open(traj_file) as f:
    lines = f.read().splitlines()
if lines[0] != settings:
    sys.exit(f"{traj_file} was written with other settings:\n  {lines[0]}\nthis run:\n  {settings}\nUse the same settings or move the old files away")
if lines[-1].startswith("# stopped"):
    sys.exit(f"{traj_file} is finished: {lines[-1]}")
rows = [numbers(line) for line in lines if not line.startswith("#")]
rows = [row for row in rows if len(row) == 22]
n_rows = int(rows[-1][0])
state = State(rows[-1][2:])

# Stored sign of each pair: its side at the start, then the side it was kicked to (the sign of Gamma) for every logged crossing
side = np.sign(State(rows[0][2:]).f())
index_of = {tuple(mode): i for i, mode in enumerate(modes)}
number = 0  # crossings logged so far
kicked = False
with open(kick_file) as f:
    kick_lines = [line for line in f.read().splitlines() if line.strip()]
kept = [line for line in kick_lines if line.startswith("#") or int(line.split()[0]) <= n_rows]
if len(kept) != len(kick_lines):  # crossings of a state whose row was never written: that step is redone
    with open(kick_file, "w") as f:
        f.write("\n".join(kept) + "\n")
for line in kept:
    if line.startswith("#"):
        continue
    p = line.split()
    number += 1
    mode = tuple(int(x) for x in p[KICK_COL + 2:KICK_COL + 7])
    if mode in index_of:
        side[index_of[mode]] = np.sign(float(p[5]))
    if int(p[1]) == 1 and int(p[0]) == n_rows:
        state.y[0:3] += np.array(p[KICK_COL + 21:KICK_COL + 24], dtype=float)
        state.y[6:9] += np.array(p[KICK_COL + 27:KICK_COL + 30], dtype=float)
        kicked = True
if kicked:
    Om_in, Om_out = omega(state.y[0:3]), omega(state.y[6:9])
    if Om_in is None or Om_out is None:
        sys.exit(f"The kicked state after row {n_rows} is not a bound orbit")
    state.Omega = np.concatenate([Om_in, Om_out])
f_start = state.f()
print(f"Starting from row {n_rows}, t = {fmt(state.t)}, {number} crossings logged so far", flush=True)

# ---------------------------------------------------------------- evolution
while True:
    if state.t >= args.t_max:
        print(f"Reached t_max = {args.t_max}")
        break
    if n_rows >= args.max_rows:
        print(f"Reached max_rows = {args.max_rows}")
        break

    # 1. Trial step with the time step of step 4
    trial, message = rk4_step(state, 0.0)
    if trial is None:
        append(traj_file, f"# stopped after row {n_rows}: {message}\n")
        print(f"Stopped after row {n_rows}: {message}")
        break
    h, f_1 = trial.dt, trial.f()

    # 2. Pairs that crossed f = 0 in the trial step
    crossed = np.nonzero((np.sign(f_1) != side) & (np.sign(f_1) == np.sign(f_1 - f_start)))[0]
    if len(crossed) == 0:
        n_rows += 1
        append(traj_file, f"{n_rows} 0 {trial.text}\n")
        state, f_start = trial, f_1
        print(f"row {n_rows}  t = {state.t:.9g}  dt = {h:.6g}  J_inner = {state.y[0]:.10g} {state.y[1]:.10g} {state.y[2]:.10g}", flush=True)
        continue

    # 3.-5. Repeat the step so that it ends at the earliest crossing
    h_c = h * min(1.0, max(0.0, np.min(f_start[crossed] / (f_start[crossed] - f_1[crossed]))))
    stopped = False
    for repeat in range(MAX_REPEATS):
        if h_c > 1e-9 * h:
            cut, message = rk4_step(state, h_c)
            if cut is None:
                stopped = True
                break
            h_c, f_c = cut.dt, cut.f()
        else:  # the crossing is at the start of the step: the crossing state is the present state
            cut, h_c, f_c = state, 0.0, f_start
        Gamma = np.array([slope_at_cut(f_start[i], f_c[i], f_1[i], h_c, h) for i in crossed])
        offset = -f_c[crossed] / Gamma  # t_res - t_c of each crossed pair
        earliest = np.min(offset)
        if abs(earliest) <= TOL * h or (h_c == 0.0 and earliest < 0.0):
            break
        h_c = min(h, max(0.0, h_c + earliest))
    if stopped:
        append(traj_file, f"# stopped after row {n_rows}: {message}\n")
        print(f"Stopped after row {n_rows}: {message}")
        break

    # 6. The crossings at this state: resonant angle from the angles, orbits for the tidal sum
    psi_in, psi_out = cut.y[3:6], cut.y[9:12]
    crossings = []
    for i, G, off in zip(crossed, Gamma, offset):
        if off > TOL * h:
            continue  # this pair crosses later; a following step finds it
        Theta = combine(modes[i:i + 1], psi_in, psi_out)[0]
        flag = int(abs(off) > TOL * h)  # not at its own crossing (pushed across f = 0 by another kick): no correction
        correction = 0.0 if flag else -f_c[i] ** 2 / (2.0 * G)
        theta_res = args.theta_res if args.theta_res is not None else (Theta + correction) % (2.0 * math.pi)
        crossings.append({"mode": modes[i], "index": i, "Gamma": G, "offset": off, "f_c": f_c[i], "tau_dt": math.sqrt(2.0 * math.pi / abs(G)) / h,
                          "Theta": Theta, "flag": flag, "correction": correction, "theta_res": theta_res, "state": cut,
                          "anc_in": ancillary(cut.y[0:3]), "anc_out": ancillary(cut.y[6:9]), "row": None})
    if len(crossings) == 0:  # no pair within TOL of this state: an ordinary state, the crossing is found from here
        if cut is state:
            sys.exit(f"No crossing could be placed after row {n_rows} (offsets {offset}); stopping")
        n_rows += 1
        append(traj_file, f"{n_rows} 0 {cut.text}\n")
        state, f_start = cut, f_c
        continue

    # 7. Delta J of both bodies (not for a crossing longer than TAU_MAX steps, or without a bound orbit)
    number = delta_J([c for c in crossings if c["tau_dt"] <= TAU_MAX and c["anc_in"] is not None and c["anc_out"] is not None], number)
    Delta_J_inner, Delta_J_outer = np.zeros(3), np.zeros(3)
    for c in crossings:
        values = np.array(c["row"][21:24] + c["row"][27:30], dtype=float) if c["row"] is not None else np.full(6, np.nan)
        c["applied"] = int(np.all(np.isfinite(values)) and not args.no_apply)
        if c["applied"]:
            Delta_J_inner += values[0:3]
            Delta_J_outer += values[3:6]
        if c["row"] is None:  # no Delta J was computed: fill the row with what is known
            number += 1
            c["row"] = [str(number), str(label)] + [str(x) for x in c["mode"]] + [fmt(c["theta_res"]), fmt(c["Gamma"])] + ["nan"] * 24 + ["0"]
        side[c["index"]] = np.sign(c["Gamma"])

    # The kicked state: J of both bodies jumps, the angles are continuous
    new = copy.deepcopy(cut)
    new.y[0:3] += Delta_J_inner
    new.y[6:9] += Delta_J_outer
    if np.any(Delta_J_inner != 0.0) or np.any(Delta_J_outer != 0.0):
        Om_in, Om_out = omega(new.y[0:3]), omega(new.y[6:9])
        if Om_in is None or Om_out is None:
            message = "the kicked state is not a bound orbit"
            append(traj_file, f"# stopped after row {n_rows}: {message}\n")
            print(f"Stopped after row {n_rows}: {message}")
            break
        new.Omega = np.concatenate([Om_in, Om_out])
    dOmega = new.Omega - cut.Omega

    # Log the crossings, then the row of the crossing state (if the step was repeated; otherwise the state is the last row)
    row_index = n_rows + (1 if cut is not state else 0)
    kick_text = ""
    for c in crossings:
        kick_text += " ".join([str(row_index), str(c["applied"]), str(c["flag"])]
                              + [fmt(x) for x in (cut.t, c["offset"], c["Gamma"], c["tau_dt"], c["f_c"], *psi_in, *psi_out, c["Theta"], c["correction"], c["theta_res"], *dOmega)]
                              + c["row"]) + "\n"
    append(kick_file, kick_text)
    if cut is not state:
        n_rows += 1
        append(traj_file, f"{n_rows} 1 {cut.text}\n")
    state, f_start = new, new.f()
    for c in crossings:
        print(f"row {n_rows}  t = {cut.t:.9g}  crossing {tuple(int(x) for x in c['mode'])}  theta_res = {c['theta_res']:.6f}  Gamma = {c['Gamma']:.6e}"
              f"  t_res - t_c = {c['offset']:.3g}  applied {c['applied']}  Delta_J_inner = {' '.join(c['row'][21:24])}", flush=True)
