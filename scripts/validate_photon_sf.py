#!/usr/bin/env python3
"""
Type 1 diagnostic validation of photon scale-factor branches in a single
nano2pico "unskimmed" pico ROOT file (output of apply_corrections.exe).

This is a standalone script -- it does NOT modify nano2pico source. It only
reads a pico ROOT file you already produced and writes a PDF report.

Checks performed
-----------------
  1. NaN / Inf counts on w_photon, w_phshape, w_fakephoton, weight
  2. Non-positive (<=0) counts on w_photon, w_phshape
  3. w_fakephoton != 1 fraction (expected ~0 for prompt/signal MC --
     fake-photon SF is meant for non-prompt/background photons)
  4. sys_photon / sys_photon_csev presence + emptiness (flags the known
     apply_corrections.cxx gap where these branches are left unfilled)
  5. Weight-decomposition self-consistency: weight / product(known
     per-event weight components) should be CONSTANT across all events
     (it equals the scalar corr.weight() normalization from merge_corrections,
     which isn't itself stored in the pico tree, so we can't check its
     absolute value -- but we CAN check that the ratio doesn't vary
     event-to-event, which verifies the multiplicative formula is being
     applied consistently)

Plots produced
---------------
  - w_photon, w_phshape, weight, w_fakephoton distributions
  - weight / product(components) ratio distribution (sanity spike)
  - nphoton, photon_pt, photon_eta distributions
  - sys_photon / sys_photon_csev vs w_photon nominal (if populated)

Usage
-----
  python3 validate_photon_sf.py <input.root> [--tree tree] [--out report.pdf]
"""

import argparse
import sys
from datetime import datetime

import numpy as np
import uproot
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Candidate per-event multiplicative weight components that make up
# pico.out_weight() in apply_corrections.cxx (see lines ~139-144), MINUS
# the scalar corr.weight() normalization factor which isn't stored in the
# output tree. Only those actually present in the file are used.
WEIGHT_COMPONENTS = [
    "w_lumi", "w_lep", "w_fs_lep", "w_bhig", "w_jetpuid", "w_trig",
    "w_isr", "w_pu", "w_prefire", "w_photon", "w_phshape",
    "w_fakephoton", "w_nnlo",
]


def safe_get(tree, name):
    """Return a numpy/awkward array for branch `name`, or None if absent."""
    if name not in tree.keys():
        return None
    try:
        return tree[name].array(library="ak")
    except Exception as exc:
        print(f"  [warn] could not read branch '{name}': {exc}")
        return None


def flatten(arr):
    """Flatten a possibly-jagged awkward array to a flat numpy array."""
    if arr is None:
        return np.array([])
    try:
        flat = ak.flatten(arr, axis=None)
        return ak.to_numpy(flat)
    except Exception:
        return ak.to_numpy(arr)


def counts_line(label, n, total):
    pct = 100.0 * n / total if total else 0.0
    return f"{label}: {n} / {total}  ({pct:.3f}%)"


def run_checks(tree, n_entries):
    results = []  # list of (check_name, status, detail_lines)

    # scalar branches (one value per event)
    w_photon = safe_get(tree, "w_photon")
    w_phshape = safe_get(tree, "w_phshape")
    w_fakephoton = safe_get(tree, "w_fakephoton")
    weight = safe_get(tree, "weight")

    # --- Check 1: NaN / Inf -------------------------------------------------
    lines = []
    status = "PASS"
    for name, arr in [("w_photon", w_photon), ("w_phshape", w_phshape),
                       ("w_fakephoton", w_fakephoton), ("weight", weight)]:
        if arr is None:
            lines.append(f"{name}: branch NOT FOUND, skipped")
            continue
        np_arr = flatten(arr) if arr.ndim > 1 else ak.to_numpy(arr)
        n_nan = int(np.sum(np.isnan(np_arr)))
        n_inf = int(np.sum(np.isinf(np_arr)))
        lines.append(counts_line(f"{name} NaN", n_nan, len(np_arr)))
        lines.append(counts_line(f"{name} Inf", n_inf, len(np_arr)))
        if n_nan or n_inf:
            status = "FAIL"
    results.append(("NaN / Inf check", status, lines))

    # --- Check 2: non-positive values ---------------------------------------
    lines = []
    status = "PASS"
    for name, arr in [("w_photon", w_photon), ("w_phshape", w_phshape)]:
        if arr is None:
            lines.append(f"{name}: branch NOT FOUND, skipped")
            continue
        np_arr = ak.to_numpy(arr)
        n_bad = int(np.sum(np_arr <= 0))
        lines.append(counts_line(f"{name} <= 0", n_bad, len(np_arr)))
        if n_bad:
            status = "WARN"
    results.append(("Non-positive weight check", status, lines))

    # --- Check 3: w_fakephoton should be 1 for signal MC ---------------------
    lines = []
    status = "PASS"
    if w_fakephoton is None:
        lines.append("w_fakephoton: branch NOT FOUND, skipped")
        status = "WARN"
    else:
        np_arr = ak.to_numpy(w_fakephoton)
        n_ne1 = int(np.sum(~np.isclose(np_arr, 1.0)))
        lines.append(counts_line("w_fakephoton != 1", n_ne1, len(np_arr)))
        lines.append("(expected ~0 for prompt/signal MC -- fake-photon SF "
                      "targets non-prompt/background photons)")
        if n_ne1 > 0:
            status = "WARN"
    results.append(("w_fakephoton sanity (signal MC)", status, lines))

    # --- Check 4: sys_photon / sys_photon_csev population --------------------
    lines = []
    status = "PASS"
    for name in ["sys_photon", "sys_photon_csev"]:
        arr = safe_get(tree, name)
        if arr is None:
            lines.append(f"{name}: branch NOT FOUND in file")
            status = "FAIL"
            continue
        lengths = ak.num(arr, axis=1) if arr.ndim > 1 else None
        if lengths is None:
            lines.append(f"{name}: not a jagged/vector branch, unexpected "
                          f"structure")
            status = "WARN"
            continue
        n_empty = int(np.sum(ak.to_numpy(lengths) == 0))
        lines.append(counts_line(f"{name} empty (len==0)", n_empty, n_entries))
        if n_empty == n_entries:
            lines.append(f"  -> {name} is EMPTY for every event. This "
                          f"matches the known apply_corrections.cxx gap "
                          f"(lines ~147,151-152 are commented out) -- "
                          f"nominal w_photon/w_phshape are still fine, but "
                          f"this file cannot be used for photon-SF "
                          f"systematic templates until that's fixed.")
            status = "FAIL"
        elif n_empty > 0:
            status = "WARN"
    results.append(("sys_photon / sys_photon_csev population", status, lines))

    # --- Check 5: weight decomposition self-consistency -----------------------
    lines = []
    status = "PASS"
    if weight is None:
        lines.append("weight: branch NOT FOUND, cannot run this check")
        status = "WARN"
        ratio = None
    else:
        weight_np = ak.to_numpy(weight).astype(np.float64)
        product = np.ones_like(weight_np)
        found, missing = [], []
        for comp in WEIGHT_COMPONENTS:
            arr = safe_get(tree, comp)
            if arr is None:
                missing.append(comp)
                continue
            found.append(comp)
            product *= ak.to_numpy(arr).astype(np.float64)
        lines.append(f"components found:   {', '.join(found) if found else '(none)'}")
        lines.append(f"components missing: {', '.join(missing) if missing else '(none)'}")
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(product != 0, weight_np / product, np.nan)
        valid = ratio[np.isfinite(ratio)]
        if len(valid) == 0:
            lines.append("no valid entries to compute ratio")
            status = "FAIL"
        else:
            mean = float(np.mean(valid))
            std = float(np.std(valid))
            rel_spread = std / abs(mean) if mean != 0 else float("nan")
            lines.append(f"weight / product(components):  mean={mean:.6g}  "
                         f"std={std:.3g}  rel.spread={rel_spread:.2e}")
            lines.append("This ratio equals the missing scalar corr.weight() "
                         "normalization factor -- it should be CONSTANT "
                         "across all events in the file (not necessarily "
                         "close to 1). A non-negligible spread means one of "
                         "the components is being read/applied "
                         "inconsistently event-to-event.")
            if rel_spread > 1e-3:
                status = "FAIL"
            elif rel_spread > 1e-5:
                status = "WARN"
    results.append(("Weight decomposition self-consistency", status, ratio, lines) if False else
                    ("Weight decomposition self-consistency", status, lines))

    return results, dict(w_photon=w_photon, w_phshape=w_phshape,
                          w_fakephoton=w_fakephoton, weight=weight,
                          ratio=ratio)


def make_title_page(pdf, in_file, tree_name, n_entries, checks):
    fig = plt.figure(figsize=(8.5, 11))
    fig.suptitle("Photon SF Type 1 Validation Report", fontsize=16, y=0.97)
    y = 0.90
    fig.text(0.08, y, f"Input file:", fontsize=10, fontweight="bold"); y -= 0.025
    fig.text(0.10, y, in_file, fontsize=8, wrap=True); y -= 0.035
    fig.text(0.08, y, f"Tree: {tree_name}     Entries: {n_entries}     "
             f"Generated: {datetime.now().isoformat(timespec='seconds')}",
             fontsize=9); y -= 0.05

    fig.text(0.08, y, "Summary of checks:", fontsize=11, fontweight="bold")
    y -= 0.035
    status_symbol = {"PASS": "[PASS]", "WARN": "[WARN]", "FAIL": "[FAIL]"}
    for name, status, lines in checks:
        fig.text(0.10, y, f"{status_symbol.get(status, status)}  {name}",
                 fontsize=9.5, fontweight="bold")
        y -= 0.022
        for line in lines:
            fig.text(0.13, y, line, fontsize=7.5)
            y -= 0.018
        y -= 0.01
        if y < 0.05:
            break
    plt.axis("off")
    pdf.savefig(fig)
    plt.close(fig)


def hist_page(pdf, title, arrays_labels, bins=60, xrange=None, xlabel="", logy=False):
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    any_data = False
    # When multiple arrays are overlaid and no explicit xrange is given, each
    # array must share the SAME bin edges -- otherwise np.histogram/ax.hist
    # auto-ranges each array independently to its own min/max, giving every
    # series a different bin width and making peak heights not directly
    # comparable (looks like "fewer entries" even when totals are equal).
    finite_arrays = [arr[np.isfinite(arr)] for arr, _ in arrays_labels
                      if arr is not None and len(arr) > 0]
    finite_arrays = [a for a in finite_arrays if len(a) > 0]
    if xrange is None and len(finite_arrays) > 1:
        combined = np.concatenate(finite_arrays)
        if len(combined) > 0:
            xrange = (float(np.min(combined)), float(np.max(combined)))
    for arr, label in arrays_labels:
        if arr is None or len(arr) == 0:
            continue
        arr = arr[np.isfinite(arr)]
        if len(arr) == 0:
            continue
        any_data = True
        ax.hist(arr, bins=bins, range=xrange, histtype="step", label=label, linewidth=1.5)
    if not any_data:
        ax.text(0.5, 0.5, "No data available for this plot", ha="center", va="center")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("events")
    if logy:
        ax.set_yscale("log")
    if any_data:
        ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="path to unskimmed pico ROOT file")
    parser.add_argument("--tree", default="tree", help="TTree name (default: tree)")
    parser.add_argument("--out", default=None, help="output PDF path")
    args = parser.parse_args()

    out_pdf = args.out or (args.input.rsplit(".root", 1)[0] + "_photon_sf_validation.pdf")

    print(f"Opening {args.input} ...")
    f = uproot.open(args.input)
    if args.tree not in [k.split(";")[0] for k in f.keys()]:
        # try to find a likely tree
        candidates = [k.split(";")[0] for k in f.keys()]
        print(f"Tree '{args.tree}' not found. Available keys: {candidates}")
        sys.exit(1)
    tree = f[args.tree]
    n_entries = tree.num_entries
    print(f"Tree '{args.tree}' has {n_entries} entries and "
          f"{len(tree.keys())} branches.")

    checks, arrs = run_checks(tree, n_entries)

    print("\n=== Check summary ===")
    for name, status, lines in checks:
        print(f"[{status}] {name}")
        for line in lines:
            print(f"    {line}")

    photon_pt = safe_get(tree, "photon_pt")
    photon_eta = safe_get(tree, "photon_eta")
    sys_photon = safe_get(tree, "sys_photon")
    sys_photon_csev = safe_get(tree, "sys_photon_csev")

    with PdfPages(out_pdf) as pdf:
        make_title_page(pdf, args.input, args.tree, n_entries, checks)

        hist_page(pdf, "w_photon (nominal ID x CSEV SF product)",
                  [(ak.to_numpy(arrs["w_photon"]) if arrs["w_photon"] is not None else None, "w_photon")],
                  bins=60, xlabel="w_photon")

        hist_page(pdf, "w_phshape (shower-shape DNN reweight)",
                  [(ak.to_numpy(arrs["w_phshape"]) if arrs["w_phshape"] is not None else None, "w_phshape")],
                  bins=60, xlabel="w_phshape")

        hist_page(pdf, "w_fakephoton (expect spike at 1 for signal MC)",
                  [(ak.to_numpy(arrs["w_fakephoton"]) if arrs["w_fakephoton"] is not None else None, "w_fakephoton")],
                  bins=60, xlabel="w_fakephoton")

        hist_page(pdf, "Total event weight",
                  [(ak.to_numpy(arrs["weight"]) if arrs["weight"] is not None else None, "weight")],
                  bins=60, xlabel="weight", logy=True)

        if arrs["ratio"] is not None:
            valid_ratio = arrs["ratio"][np.isfinite(arrs["ratio"])]
            hist_page(pdf, "weight / product(components)  [should be a single spike]",
                      [(valid_ratio, "ratio")], bins=60,
                      xlabel="weight / product(known components)")

        if photon_pt is not None:
            nphoton = ak.to_numpy(ak.num(photon_pt, axis=1))
            hist_page(pdf, "nphoton", [(nphoton, "nphoton")],
                      bins=range(0, int(max(nphoton, default=1)) + 2),
                      xlabel="nphoton")
            hist_page(pdf, "photon_pt (all photons, flattened)",
                      [(flatten(photon_pt), "photon_pt")], bins=60,
                      xlabel="photon pt [GeV]")

        if photon_eta is not None:
            hist_page(pdf, "photon_eta (all photons, flattened)",
                      [(flatten(photon_eta), "photon_eta")], bins=60,
                      xlabel="photon eta")

        if sys_photon is not None and arrs["w_photon"] is not None:
            lengths = ak.num(sys_photon, axis=1)
            if int(np.sum(ak.to_numpy(lengths) == 2)) > 0:
                mask = ak.to_numpy(lengths) == 2
                s0 = ak.to_numpy(sys_photon[mask][:, 0])
                s1 = ak.to_numpy(sys_photon[mask][:, 1])
                nom = ak.to_numpy(arrs["w_photon"])[mask]
                hist_page(pdf, "sys_photon[0] and sys_photon[1] vs nominal w_photon "
                          "(ID SF up/down -- CSEV held at nominal)",
                          [(s0, "sys_photon[0]"), (s1, "sys_photon[1]"), (nom, "w_photon nominal")],
                          bins=60, xlabel="weight")
            else:
                fig = plt.figure(figsize=(8.5, 3))
                fig.text(0.1, 0.5, "sys_photon branch is present but empty for "
                         "every event -- see 'sys_photon / sys_photon_csev "
                         "population' check on the title page.", fontsize=11)
                plt.axis("off")
                pdf.savefig(fig)
                plt.close(fig)

        if sys_photon_csev is not None and arrs["w_photon"] is not None:
            lengths = ak.num(sys_photon_csev, axis=1)
            if int(np.sum(ak.to_numpy(lengths) == 2)) > 0:
                mask = ak.to_numpy(lengths) == 2
                s0 = ak.to_numpy(sys_photon_csev[mask][:, 0])
                s1 = ak.to_numpy(sys_photon_csev[mask][:, 1])
                nom = ak.to_numpy(arrs["w_photon"])[mask]
                hist_page(pdf, "sys_photon_csev[0] and sys_photon_csev[1] vs nominal w_photon "
                          "(CSEV SF up/down -- ID held at nominal)",
                          [(s0, "sys_photon_csev[0]"), (s1, "sys_photon_csev[1]"), (nom, "w_photon nominal")],
                          bins=60, xlabel="weight")
            else:
                fig = plt.figure(figsize=(8.5, 3))
                fig.text(0.1, 0.5, "sys_photon_csev branch is present but empty for "
                         "every event -- see 'sys_photon / sys_photon_csev "
                         "population' check on the title page.", fontsize=11)
                plt.axis("off")
                pdf.savefig(fig)
                plt.close(fig)

    print(f"\nWrote report: {out_pdf}")


if __name__ == "__main__":
    main()
