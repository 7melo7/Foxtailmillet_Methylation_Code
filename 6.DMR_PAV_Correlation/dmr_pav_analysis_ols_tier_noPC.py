#!/usr/bin/env python3
"""
DMR-PAV Linear Regression with Tiered Permutation (OLS, no PC)
================================================================
Linear model per DMR-PAV pair WITHOUT PC adjustment:
  DMR = α + γ·G_p + ε

Closed-form OLS → t-test on γ̂ → p, beta, SE, 95% CI.
1000 label-shuffling permutations → tier-based empirical p → FDR.

Tier definition (p-value gradient):
  [0, 0.001]   → tier 0
  (0.001, 0.01] → tier 1
  (0.01, 0.02]  → tier 2
  (0.02, 0.03]  → tier 3
  (0.03, 0.04]  → tier 4
  (0.04, 0.05]  → tier 5
  > 0.05        → tier 6

Identical to OLS+PC except no PC covariates.
"""

import numpy as np
from scipy.stats import t as t_dist
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
import argparse
import time
import os
import csv
import warnings

warnings.filterwarnings("ignore")

N_PERMUTATIONS = 1000
N_THREADS     = 8
BIN_STEP      = 100
MAX_DIST      = 10000
FDR_THRESHOLD = 0.05


def timer(msg):
    class _T:
        def __enter__(slf):
            slf.t0 = time.time()
            return slf
        def __exit__(slf, *a):
            print(f"[{msg}]  {time.time() - slf.t0:.1f}s", flush=True)
    return _T()


# ════════════════════════════════════════════
#  OLS — closed-form linear regression
# ════════════════════════════════════════════

def ols_fit(y, X):
    """
    Ordinary least squares: β = (X'X)⁻¹X'y.
    Returns beta, se, p_value (t-test for last coefficient).

    For degenerate y (all 0 or all 1): returns NaN.
    """
    n, p = X.shape
    if np.all(y == 0) or np.all(y == 1):
        return np.full(p, np.nan), np.full(p, np.nan), np.nan

    try:
        XtX = X.T @ X
        XtX_inv = np.linalg.solve(XtX, np.eye(p))
        beta = XtX_inv @ X.T @ y
    except np.linalg.LinAlgError:
        return np.full(p, np.nan), np.full(p, np.nan), np.nan

    residuals = y - X @ beta
    df = n - p
    if df <= 0:
        return beta, np.full(p, np.nan), np.nan
    sigma2 = np.sum(residuals ** 2) / df
    var_beta = sigma2 * np.diag(XtX_inv)
    se = np.sqrt(np.maximum(var_beta, 0))

    t_val = beta[-1] / se[-1] if se[-1] > 0 else 0.0
    p_val = 2.0 * t_dist.sf(abs(t_val), df)

    return beta, se, p_val


PAIR_DTYPE = np.dtype([
    ("dmr_idx",         "i4"),
    ("pav_idx",         "i4"),
    ("distance_bp",     "i4"),
    ("abs_distance",    "i4"),
    ("t_stat",          "f4"),
    ("p_value",         "f4"),
    ("tier",            "i1"),
    ("beta",            "f4"),
    ("beta_se",         "f4"),
    ("beta_ci_low",     "f4"),
    ("beta_ci_high",    "f4"),
    ("dmr_mean",        "f4"),
    ("pav_maf",         "f4"),
])


def get_tier(p_val):
    """Map Wald p-value to significance tier (0=most sig, 6=not sig)."""
    if np.isnan(p_val):
        return 99
    if p_val <= 0.001:  return 0
    if p_val <= 0.01:   return 1
    if p_val <= 0.02:   return 2
    if p_val <= 0.03:   return 3
    if p_val <= 0.04:   return 4
    if p_val <= 0.05:   return 5
    return 6


def _compute_pairs_batch(dmr_indices, dmr_matrix, pav_matrix, X_base, dmr_pav_pairs):
    """For each DMR in batch, fit OLS for every nearby PAV. Returns numpy array."""
    max_est = sum(len(dmr_pav_pairs.get(d, [])) for d in dmr_indices)
    max_est = max(max_est, 1)
    buf = np.empty(max_est, dtype=PAIR_DTYPE)
    pos = 0
    nan = np.float32('nan')

    for d_idx in dmr_indices:
        y = dmr_matrix[d_idx, :].astype(np.float64)
        pairs = dmr_pav_pairs.get(d_idx, [])
        if not pairs:
            continue

        dmr_mean_val = float(np.mean(y))

        # degenerate DMR
        if np.all(y == 0) or np.all(y == 1):
            for pi, dist in pairs:
                g = pav_matrix[pi, :].astype(np.float64)
                maf = float(np.minimum(np.mean(g), 1.0 - np.mean(g)))
                buf[pos] = (d_idx, pi, int(dist), abs(dist),
                            nan, nan, 99, nan, nan, nan, nan, dmr_mean_val, maf)
                pos += 1
            continue

        for loc, (pi, dist) in enumerate(pairs):
            g_vec = pav_matrix[pi, :].astype(np.float64)
            X_full = np.column_stack([X_base, g_vec])
            beta, se, p_val = ols_fit(y, X_full)

            if np.any(np.isnan(beta)):
                g = pav_matrix[pi, :].astype(np.float64)
                maf = float(np.minimum(np.mean(g), 1.0 - np.mean(g)))
                buf[pos] = (d_idx, pi, int(dist), abs(dist),
                            nan, nan, 99, nan, nan, nan, nan, dmr_mean_val, maf)
                pos += 1
                continue

            n = len(y)
            p_model = X_full.shape[1]
            df = n - p_model

            if df > 0 and se[-1] > 0:
                t_val = beta[-1] / se[-1]
                p_wald = float(2.0 * t_dist.sf(abs(t_val), df))
                tier    = get_tier(p_wald)
                beta_hat = float(beta[-1])
                se_val   = float(se[-1])
                ci_low   = beta_hat - 2.0 * se_val
                ci_high  = beta_hat + 2.0 * se_val
                t_stat   = float(t_val)
            else:
                t_stat   = nan
                p_wald   = nan
                tier     = 99
                beta_hat = nan
                se_val   = nan
                ci_low   = nan
                ci_high  = nan

            g = pav_matrix[pi, :].astype(np.float64)
            maf = float(np.minimum(np.mean(g), 1.0 - np.mean(g)))
            buf[pos] = (d_idx, pi, int(dist), abs(dist),
                        t_stat, p_wald, tier, beta_hat, se_val, ci_low, ci_high,
                        dmr_mean_val, maf)
            pos += 1

    return buf[:pos]


def compute_all_pairs(dmr_matrix, pav_matrix, X_base, dmr_pav_pairs, n_threads):
    """Real-data: Wald test for all 514K pairs (threaded)."""
    import threading
    all_idx = list(range(dmr_matrix.shape[0]))
    chunk_size = max(1, len(all_idx) // n_threads)
    chunks = [all_idx[i:i + chunk_size] for i in range(0, len(all_idx), chunk_size)]

    total = sum(len(dmr_pav_pairs.get(d, [])) for d in all_idx)
    result = np.empty(total, dtype=PAIR_DTYPE)
    next_pos = [0]
    lock = threading.Lock()

    with ThreadPoolExecutor(max_workers=n_threads) as ex:
        futures = {ex.submit(_compute_pairs_batch, c, dmr_matrix, pav_matrix,
                             X_base, dmr_pav_pairs): c for c in chunks}
        for fut in as_completed(futures):
            try:
                arr = fut.result()
                n = len(arr)
                if n:
                    with lock:
                        s = next_pos[0]; next_pos[0] += n
                    result[s:s + n] = arr
            except Exception as e:
                import traceback
                print(f"  Error: {type(e).__name__}: {e}", flush=True)
                traceback.print_exc()
    return result[:next_pos[0]]


# ════════════════════════════════════════════
#  Permutation
# ════════════════════════════════════════════

def _compute_pairs_perm_batch(dmr_indices, dmr_matrix, pav_matrix, X_base,
                              dmr_pav_pairs, real_tier_lookup, emp_counter):
    """Permuted data: OLS fit, tier-based comparison."""
    for d_idx in dmr_indices:
        y = dmr_matrix[d_idx, :].astype(np.float64)
        pairs = dmr_pav_pairs.get(d_idx, [])
        if not pairs or (np.all(y == 0) or np.all(y == 1)):
            continue

        for loc, (pi, dist) in enumerate(pairs):
            g_vec = pav_matrix[pi, :].astype(np.float64)
            X_full = np.column_stack([X_base, g_vec])
            beta, se, p_val = ols_fit(y, X_full)

            if np.any(np.isnan(beta)) or se[-1] <= 0:
                continue

            t_val = beta[-1] / se[-1]
            n, p_model = len(y), X_full.shape[1]
            df = n - p_model
            if df <= 0:
                continue
            p_perm = float(2.0 * t_dist.sf(abs(t_val), df))
            perm_tier = get_tier(p_perm)

            real_tier = real_tier_lookup.get((d_idx, pi))
            if real_tier is not None and perm_tier < real_tier:
                emp_counter[d_idx][loc] += 1


def run_permutations(dmr_matrix, pav_matrix, X_base, dmr_pav_pairs,
                     real_tier_lookup, n_perms, n_threads):
    """1000 label-shuffles, per-pair tier-based counter update."""
    print(f"\nRunning {n_perms} tier-based permutations ({n_threads} threads)...", flush=True)
    n_dmrs = dmr_matrix.shape[0]

    # per-DMR counter vectors
    emp_counter = {}
    for d_idx in dmr_pav_pairs:
        emp_counter[d_idx] = np.zeros(len(dmr_pav_pairs[d_idx]), dtype=np.int16)

    for perm_i in range(n_perms):
        t0 = time.time()
        dmr_shuffled = dmr_matrix.copy()
        rng = np.random.default_rng(perm_i)
        for i in range(n_dmrs):
            dmr_shuffled[i, :] = rng.permutation(dmr_shuffled[i, :])

        all_idx = list(range(n_dmrs))
        chunk_size = max(1, n_dmrs // n_threads)
        chunks = [all_idx[i:i + chunk_size] for i in range(0, n_dmrs, chunk_size)]

        with ThreadPoolExecutor(max_workers=n_threads) as ex:
            futures = [ex.submit(_compute_pairs_perm_batch, c, dmr_shuffled,
                                 pav_matrix, X_base, dmr_pav_pairs,
                                 real_tier_lookup, emp_counter)
                       for c in chunks]
            for fut in as_completed(futures):
                try:
                    fut.result()
                except Exception as e:
                    import traceback
                    print(f"  Perm error: {e}", flush=True); traceback.print_exc()

        elapsed = time.time() - t0
        if (perm_i + 1) % 10 == 0 or perm_i == 0:
            print(f"  Permutation {perm_i + 1}/{n_perms}  ({elapsed:.1f}s)", flush=True)

    return emp_counter


# ════════════════════════════════════════════
#  Empirical p-values & FDR
# ════════════════════════════════════════════

def benjamini_hochberg(pvals):
    """Returns BH-adjusted q-values."""
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    qvals = np.minimum(1.0, sorted_p * n / np.arange(1, n + 1))
    for i in range(n - 2, -1, -1):
        qvals[i] = min(qvals[i], qvals[i + 1])
    reorder = np.argsort(order)
    return qvals[reorder]


def compute_emp_p_and_fdr(flat_arr, emp_counter, dmr_pav_pairs, n_perms):
    """Add empirical_p_value and fdr_q columns. Returns emp_p, fdr_q arrays."""
    n = len(flat_arr)
    emp_p = np.full(n, np.float32(1.0))

    # Build (d_idx, p_idx) → position mapping
    pav_pos = {}
    for i in range(n):
        d = int(flat_arr["dmr_idx"][i])
        p = int(flat_arr["pav_idx"][i])
        pav_pos[(d, p)] = i

    for d_idx, pairs in dmr_pav_pairs.items():
        counter_vec = emp_counter.get(d_idx)
        if counter_vec is None or len(counter_vec) == 0:
            continue
        for loc, (p_idx, _) in enumerate(pairs):
            pos = pav_pos.get((d_idx, p_idx))
            if pos is not None:
                emp_p[pos] = np.float32((int(counter_vec[loc]) + 1) / (n_perms + 1))

    # FDR on empirical p-values
    fdr_q = benjamini_hochberg(emp_p.astype(np.float64)).astype(np.float32)
    return emp_p, fdr_q


# ════════════════════════════════════════════
#  Distance enrichment
# ════════════════════════════════════════════

def distance_enrichment(flat_arr, emp_p, fdr_q, bin_step=100, max_dist=10000,
                        fdr_thresh=0.05):
    abs_d = flat_arr["abs_distance"]
    bins = np.arange(0, max_dist + bin_step, bin_step)
    results = []
    for i in range(len(bins) - 1):
        lo, hi = bins[i], bins[i + 1]
        mask = (abs_d >= lo) & (abs_d < hi)
        n_total = int(np.sum(mask))
        if n_total == 0:
            continue
        n_sig_emp = int(np.sum(emp_p[mask] < fdr_thresh))
        n_sig_fdr = int(np.sum(fdr_q[mask] < fdr_thresh))
        median_beta = float(np.median(flat_arr["beta"][mask])) if n_total > 0 else np.nan
        median_t = float(np.median(flat_arr["t_stat"][mask])) if n_total > 0 else np.nan
        results.append({
            "distance_bin_start": int(lo),
            "distance_bin_end":   int(hi),
            "n_pairs":            n_total,
            "n_sig_emp":          n_sig_emp,
            "n_sig_fdr":          n_sig_fdr,
            "sig_ratio_emp":      n_sig_emp / n_total,
            "sig_ratio_fdr":      n_sig_fdr / n_total,
            "median_beta":        median_beta,
            "median_t":           median_t,
            "mean_empirical_p":   float(np.mean(emp_p[mask])),
        })
    return results


# ════════════════════════════════════════════
#  Data loading (unchanged)
# ════════════════════════════════════════════

def load_data(data_dir, pairwise_file="DMR_overlap_PAV.20kbp.pairwise"):
    print("Loading PAV genotypes ...", flush=True)
    pav_rows, pav_ids = [], []
    with open(os.path.join(data_dir, "PAV.GenoMatrix.txt"), "r") as f:
        pav_samples = f.readline().strip().split("\t")[1:]
        for line in f:
            parts = line.strip().split("\t")
            pav_rows.append([int(x) for x in parts[1:]])
            pav_ids.append(parts[0])
    pav_matrix = np.array(pav_rows, dtype=np.int8)
    pav_id_to_idx = {sid: i for i, sid in enumerate(pav_ids)}
    print(f"    {len(pav_ids)} PAVs x {pav_matrix.shape[1]} samples", flush=True)

    print("Loading DMR genotypes ...", flush=True)
    dmr_rows, dmr_ids = [], []
    with open(os.path.join(data_dir, "DMR.genotypev1.txt"), "r") as f:
        dmr_samples = f.readline().strip().split("\t")[1:]
        for line in f:
            parts = line.strip().split("\t")
            dmr_rows.append([int(x) for x in parts[1:]])
            dmr_ids.append(parts[0])
    dmr_matrix = np.array(dmr_rows, dtype=np.int8)
    dmr_id_to_idx = {sid: i for i, sid in enumerate(dmr_ids)}
    print(f"    {len(dmr_ids)} DMRs x {dmr_matrix.shape[1]} samples", flush=True)

    canonical = pav_samples
    dmr_reorder = [dmr_samples.index(s) for s in canonical]
    dmr_matrix = dmr_matrix[:, dmr_reorder]
    print(f"    Aligned: {len(canonical)} shared", flush=True)

    X_base = np.ones((len(canonical), 1), dtype=np.float64)

    print("Loading pairwise data ...", flush=True)
    dmr_pav_pairs = defaultdict(list)
    with open(os.path.join(data_dir, pairwise_file), "r") as f:
        for line in f:
            dmr_id, pav_id, d_str = line.strip().split("\t")
            if dmr_id in dmr_id_to_idx and pav_id in pav_id_to_idx:
                dmr_pav_pairs[dmr_id_to_idx[dmr_id]].append(
                    (pav_id_to_idx[pav_id], int(d_str)))
    n_pairs = sum(len(v) for v in dmr_pav_pairs.values())
    print(f"    {n_pairs} pairs", flush=True)

    return (dmr_matrix, pav_matrix, X_base, dmr_ids, pav_ids,
            dmr_id_to_idx, pav_id_to_idx, dmr_pav_pairs)


# ════════════════════════════════════════════
#  Output
# ════════════════════════════════════════════

def write_outputs(output_dir, flat_arr, emp_p, fdr_q, dmr_ids, pav_ids, enrichment):
    with open(os.path.join(output_dir, "dmr_pav_pairwise_results.csv"), "w", newline="") as f:
        fields = ["dmr_id", "pav_id", "distance_bp", "abs_distance_bp",
                  "t_stat", "p_value", "tier", "beta", "beta_se",
                  "beta_ci_low", "beta_ci_high",
                  "dmr_mean", "pav_maf", "empirical_p_value", "fdr_q"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        order = np.argsort(flat_arr["abs_distance"])
        for i in order:
            w.writerow({
                "dmr_id":             dmr_ids[int(flat_arr["dmr_idx"][i])],
                "pav_id":             pav_ids[int(flat_arr["pav_idx"][i])],
                "distance_bp":        int(flat_arr["distance_bp"][i]),
                "abs_distance_bp":    int(flat_arr["abs_distance"][i]),
                "t_stat":              float(flat_arr["t_stat"][i]),
                "p_value":             float(flat_arr["p_value"][i]),
                "tier":                int(flat_arr["tier"][i]),
                "beta":                float(flat_arr["beta"][i]),
                "beta_se":             float(flat_arr["beta_se"][i]),
                "beta_ci_low":         float(flat_arr["beta_ci_low"][i]),
                "beta_ci_high":        float(flat_arr["beta_ci_high"][i]),
                "dmr_mean":           float(flat_arr["dmr_mean"][i]),
                "pav_maf":            float(flat_arr["pav_maf"][i]),
                "empirical_p_value":  float(emp_p[i]),
                "fdr_q":              float(fdr_q[i]),
            })

    # Distance enrichment
    with open(os.path.join(output_dir, "distance_enrichment.csv"), "w", newline="") as f:
        fields = ["distance_bin_start", "distance_bin_end", "n_pairs",
                  "n_sig_emp", "n_sig_fdr", "sig_ratio_emp", "sig_ratio_fdr",
                  "median_beta", "median_t", "mean_empirical_p"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in enrichment:
            w.writerow(row)

    print(f"\nOutput → {output_dir}/", flush=True)
    print(f"  dmr_pav_pairwise_results.csv  ({len(flat_arr):,} rows)", flush=True)
    print(f"  distance_enrichment.csv       ({len(enrichment)} rows)", flush=True)


# ════════════════════════════════════════════
#  Main
# ════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(description="DMR-PAV OLS Tier Permutation (no PC)")
    p.add_argument("--data-dir", default=".")
    p.add_argument("--pairwise-file", default="DMR_overlap_PAV.20kbp.pairwise")
    p.add_argument("--output-dir", default="./output_ols_tier_noPC")
    p.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS)
    p.add_argument("--n-threads", type=int, default=N_THREADS)
    p.add_argument("--skip-permutations", action="store_true")
    p.add_argument("--bin-step", type=int, default=BIN_STEP)
    p.add_argument("--fdr-threshold", type=float, default=FDR_THRESHOLD)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"OLS tier (no PC)  |  Perms: {args.n_permutations}  |  "
          f"Threads: {args.n_threads}")

    # ── Phase 1 ──
    with timer("Phase 1 – Load"):
        (dmr_mat, pav_mat, X_base, dmr_ids, pav_ids,
         _, _, dmr_pav_pairs) = load_data(args.data_dir, args.pairwise_file)

    n_pairs = sum(len(v) for v in dmr_pav_pairs.values())
    print(f"\n  {dmr_mat.shape[0]} DMRs  {pav_mat.shape[0]} PAVs  "
          f"{n_pairs:,} pairs\n")

    # ── Phase 2 ──
    with timer("Phase 2 – OLS fit (all pairs)"):
        flat_arr = compute_all_pairs(dmr_mat, pav_mat, X_base, dmr_pav_pairs,
                                     args.n_threads)
    print(f"  Pairs computed: {len(flat_arr):,}", flush=True)

    # ── Phase 3 ──
    emp_p = np.ones(len(flat_arr), dtype=np.float32)
    fdr_q = np.ones(len(flat_arr), dtype=np.float32)
    if not args.skip_permutations:
        # Build tier lookup for fast perm comparison
        real_tier = {}
        for i in range(len(flat_arr)):
            d = int(flat_arr["dmr_idx"][i])
            p = int(flat_arr["pav_idx"][i])
            t = int(flat_arr["tier"][i])
            if t < 99:
                real_tier[(d, p)] = t

        with timer("Phase 3 – Tier-based permutations"):
            emp_counter = run_permutations(dmr_mat, pav_mat, X_base,
                                           dmr_pav_pairs, real_tier,
                                           args.n_permutations, args.n_threads)

        with timer("Empirical p + FDR"):
            emp_p, fdr_q = compute_emp_p_and_fdr(
                flat_arr, emp_counter, dmr_pav_pairs, args.n_permutations)
    else:
        print("\n[SKIP] Permutation test", flush=True)

    # ── Phase 4 ──
    with timer("Phase 4 – Distance enrichment"):
        enrichment = distance_enrichment(
            flat_arr, emp_p, fdr_q, bin_step=args.bin_step,
            fdr_thresh=args.fdr_threshold)

    # ── Phase 5 ──
    with timer("Phase 5 – Write outputs"):
        write_outputs(args.output_dir, flat_arr, emp_p, fdr_q, dmr_ids, pav_ids,
                      enrichment)

    # ── Summary ──
    n_sig_emp = int(np.sum(emp_p < args.fdr_threshold))
    n_sig_fdr = int(np.sum(fdr_q < args.fdr_threshold))
    print(f"\n{'='*60}")
    print(f"RESULTS  |  {len(flat_arr):,} pairs")
    print(f"  emp_p  < {args.fdr_threshold}: {n_sig_emp:,}")
    print(f"  FDR q  < {args.fdr_threshold}: {n_sig_fdr:,}")
    print(f"  emp_p resolution: 1/{args.n_permutations+1} = {1.0/(args.n_permutations+1):.4f}")
    print(f"{'='*60}")

    if enrichment:
        top = sorted(enrichment, key=lambda r: r["sig_ratio_emp"], reverse=True)[:6]
        print(f"\nTop bins by emp_p sig ratio:")
        for r in top:
            print(f"  {r['distance_bin_start']:>4d}-{r['distance_bin_end']:<4d} bp  "
                  f"n={r['n_pairs']:>6d}  sig_emp={r['n_sig_emp']:>5d}  "
                  f"ratio={r['sig_ratio_emp']:.4f}  median_beta={r['median_beta']:.4f}")

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
