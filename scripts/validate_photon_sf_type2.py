#!/usr/bin/env python3
"""
Type 2 diagnostic: distributional impact of the photon SF chain on a
single nano2pico "unskimmed" pico ROOT file.

Where Type 1 (validate_photon_sf.py) checked that w_photon/w_phshape/weight
are calculated correctly, Type 2 checks HOW MUCH they change the kinematic
distributions you actually care about, and whether the photon-SF
systematic uncertainty distorts the shape of those distributions (as
opposed to just shifting the overall yield).

This is a standalone script -- it does NOT modify nano2pico source. It only
reads a pico ROOT file you already produced and writes a PDF report.

What it does
------------
  1. Yield summary: unweighted event count, sum(weight), sum(w_photon*w_phshape),
     and how much of the total weight is attributable to the photon SF chain.
  2. Shape comparison (unweighted vs weight-applied, normalized to unit area)
     for leading/subleading photon pt, eta, phi, idmva, plus event-level
     diphoton mass, diphoton pt, and diphoton dR. A ratio panel underneath
     shows local shape distortion, not just an overall normalization shift.
  3. Systematic-band comparison for mgg, ptgg, and leading photon pt:
     overlays the nominal-weighted distribution against the
     sys_photon[0]/[1]-shifted distributions (holding all other weight
     components fixed), using
         weight_shifted = weight / w_photon * sys_photon[i]
     which is the standard way to build up/down systematic templates.

Inclusive vs signal-region-cut view
------------------------------------
By default this looks at the inclusive (no selection) population, which is
best for catching pathological behavior anywhere in phase space. Pass
--photon-sr-cuts to additionally require leading pt > 35 GeV, subleading
pt > 25 GeV, and |eta| < 2.5 excluding the 1.4442-1.566 barrel/endcap gap,
for both photons -- this is NOT a full signal-region definition (no
b-jet/MET cuts), just the photon-object acceptance, so you can see the
distributions that actually matter for the photon-only part of the chain.

Usage
-----
  python3 validate_photon_sf_type2.py <input.root> [--tree tree] [--out report.pdf]
      [--weight-branches w_lumi,w_photon,w_phshape] [--photon-sr-cuts]
      [--mgg-branch NAME] [--ptgg-branch NAME] [--drgg-branch NAME]
      [--idmva-branch NAME] [--phi-branch NAME]

If a branch name guess below doesn't match your tree, override it with the
corresponding --*-branch flag. Run `python3 -c "import uproot; print(list(uproot.open('file.root:tree').keys()))"`
to see your actual branch names if needed.
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

MGG_CANDIDATES = ["photonphoton_m", "photon_photon_m", "diphoton_m", "mgg",
                   "photonphoton_am", "gammagamma_m"]
PTGG_CANDIDATES = ["photonphoton_pt", "photon_photon_pt", "diphoton_pt", "ptgg"]
DRGG_CANDIDATES = ["photonphoton_dr", "photon_photon_dr", "diphoton_dr", "drgg"]
IDMVA_CANDIDATES = ["photon_idmva", "photon_mva", "photon_id_mva", "photon_mvaid"]
PHI_CANDIDATES = ["photon_phi"]


def safe_get(tree, name):
    if name is None or name not in tree.keys():
        return None
    try:
        return tree[name].array(library="ak")
    except Exception as exc:
        print(f"  [warn] could not read branch '{name}': {exc}")
        return None


def find_branch(tree, candidates, override=None):
    if override is not None:
        if override in tree.keys():
            return override
        print(f"  [warn] requested branch '{override}' not found in file")
        return None
    for c in candidates:
        if c in tree.keys():
            return c
    return None


def get_weight(tree, branch_spec):
    """Build the 'weight' array used throughout this report

    branch_spec is a comma-separated list of branch names, multiplied
    together elementwise. Use this to test a partial weight chain (e.g.
    "w_lumi,w_photon,w_phshape") before other components of the full
    production `weight` formula are trusted.
    """
    names = [b.strip() for b in branch_spec.split(",") if b.strip()]
    product = None
    used, missing = [], []
    for name in names:
        arr = safe_get(tree, name)
        if arr is None:
            missing.append(name)
            continue
        vals = ak.to_numpy(arr).astype(np.float64)
        product = vals if product is None else product * vals
        used.append(name)
    return product, used, missing


def full_length_value(arr, index, n_entries):
    """Return a length-n_entries float64 array holding arr[event][index]
    where that index exists (jagged branch), or arr[event] directly if the
    branch is already flat/scalar. NaN where the index doesn't exist -- the
    plotting helper drops NaN rows automatically, so this is how missing
    values and any additional cuts (via apply_cut) get excluded."""
    if arr is None:
        return np.full(n_entries, np.nan)
    try:
        is_jagged = arr.ndim > 1
    except AttributeError:
        is_jagged = False
    if not is_jagged:
        return ak.to_numpy(arr).astype(np.float64)
    lengths = ak.to_numpy(ak.num(arr, axis=1))
    mask = lengths > index
    out = np.full(n_entries, np.nan)
    out[mask] = ak.to_numpy(arr[mask][:, index]).astype(np.float64)
    return out


def apply_cut(vals, keep_mask):
    """Return a copy of vals with entries failing keep_mask set to NaN."""
    out = vals.copy()
    out[~keep_mask] = np.nan
    return out


def eb_ee_fiducial(eta):
    """CMS ECAL barrel/endcap fiducial acceptance: |eta|<2.5 excluding the
    1.4442-1.566 transition gap. NaN-safe (NaN eta -> False)."""
    with np.errstate(invalid="ignore"):
        abs_eta = np.abs(eta)
        ok = (abs_eta < 2.5) & ~((abs_eta > 1.4442) & (abs_eta < 1.566))
    return np.nan_to_num(ok.astype(float), nan=0.0).astype(bool)


def ratio_hist_page(pdf, title, xlabel, datasets, bins=40, xrange=None, normalize=True):
    """datasets: list of (values, weights, label). Draws normalized shape
    overlay on top, ratio-to-first-dataset panel below. NaN entries in
    either values or weights are dropped."""
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8.5, 6.5), sharex=True,
                                    gridspec_kw={"height_ratios": [3, 1]})
    hists = []
    edges = None
    any_data = False
    # Same fix as validate_photon_sf.py's hist_page: without an explicit
    # xrange, np.histogram would auto-range each dataset to its own min/max,
    # giving different bin widths per series and making shapes/heights not
    # directly comparable. Force a shared range across all datasets here.
    if xrange is None:
        finite_vals = []
        for vals, wts, _ in datasets:
            if vals is None or len(vals) == 0:
                continue
            finite = np.isfinite(vals)
            if wts is not None:
                finite &= np.isfinite(wts)
            v = vals[finite]
            if len(v) > 0:
                finite_vals.append(v)
        if finite_vals:
            combined = np.concatenate(finite_vals)
            xrange = (float(np.min(combined)), float(np.max(combined)))
    for vals, wts, label in datasets:
        if vals is None or len(vals) == 0:
            continue
        finite = np.isfinite(vals)
        if wts is not None:
            finite &= np.isfinite(wts)
        vals_f = vals[finite]
        wts_f = wts[finite] if wts is not None else None
        if len(vals_f) == 0:
            continue
        any_data = True
        h, e = np.histogram(vals_f, bins=bins, range=xrange, weights=wts_f)
        if normalize and h.sum() != 0:
            h = h / h.sum()
        hists.append((h, label))
        edges = e
    if not any_data:
        ax0.text(0.5, 0.5, "No data available for this plot", ha="center", va="center")
        ax0.set_title(title)
        plt.axis("off")
        pdf.savefig(fig)
        plt.close(fig)
        return
    centers = 0.5 * (edges[:-1] + edges[1:])
    for h, label in hists:
        ax0.step(centers, h, where="mid", label=label, linewidth=1.5)
    ax0.set_title(title, fontsize=10)
    ax0.set_ylabel("normalized events" if normalize else "events")
    ax0.legend(fontsize=8)
    ax0.grid(alpha=0.3)

    ref_h = hists[0][0]
    for h, label in hists[1:]:
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(ref_h != 0, h / ref_h, np.nan)
        ax1.step(centers, ratio, where="mid", label=f"{label} / {hists[0][1]}", linewidth=1.2)
    ax1.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("ratio")
    ax1.set_ylim(0.5, 1.5)
    ax1.grid(alpha=0.3)
    if len(hists) > 1:
        ax1.legend(fontsize=7)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def missing_page(pdf, message):
    fig = plt.figure(figsize=(8.5, 2.5))
    fig.text(0.1, 0.5, message, fontsize=10, wrap=True)
    plt.axis("off")
    pdf.savefig(fig)
    plt.close(fig)


def make_summary_page(pdf, in_file, tree_name, n_entries, sr_cuts, n_pass_sr, summary_lines):
    fig = plt.figure(figsize=(8.5, 11))
    fig.suptitle("Photon SF Type 2 Validation Report (distributional impact)",
                 fontsize=14, y=0.97)
    y = 0.90
    fig.text(0.08, y, "Input file:", fontsize=10, fontweight="bold"); y -= 0.025
    fig.text(0.10, y, in_file, fontsize=8, wrap=True); y -= 0.035
    fig.text(0.08, y, f"Tree: {tree_name}     Entries: {n_entries}     "
             f"Generated: {datetime.now().isoformat(timespec='seconds')}",
             fontsize=9); y -= 0.035
    if sr_cuts:
        fig.text(0.08, y, f"photon-SR cuts APPLIED: pt_lead>35, pt_sublead>25, "
                 f"|eta|<2.5 (EB/EE gap excluded), both photons. "
                 f"{n_pass_sr}/{n_entries} events pass "
                 f"({100*n_pass_sr/n_entries:.1f}%).", fontsize=8.5)
    else:
        fig.text(0.08, y, "photon-SR cuts NOT applied (inclusive view). "
                 "Pass --photon-sr-cuts to restrict to the photon "
                 "acceptance region.", fontsize=8.5)
    y -= 0.05

    fig.text(0.08, y, "Yield summary:", fontsize=11, fontweight="bold")
    y -= 0.035
    for line in summary_lines:
        fig.text(0.10, y, line, fontsize=9)
        y -= 0.024
    plt.axis("off")
    pdf.savefig(fig)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="path to unskimmed pico ROOT file")
    parser.add_argument("--tree", default="tree", help="TTree name (default: tree)")
    parser.add_argument("--out", default=None, help="output PDF path")
    parser.add_argument("--mgg-branch", default=None)
    parser.add_argument("--ptgg-branch", default=None)
    parser.add_argument("--drgg-branch", default=None)
    parser.add_argument("--idmva-branch", default=None)
    parser.add_argument("--phi-branch", default=None)
    parser.add_argument("--weight-branches", default="weight",
                         help="comma-separated branch name(s) to multiply "
                              "together to form the 'weight' used in this "
                              "report. Default: 'weight' (full production "
                              "weight). E.g. 'w_lumi,w_photon,w_phshape' to "
                              "restrict to just the validated photon SF "
                              "chain.")
    parser.add_argument("--photon-sr-cuts", action="store_true",
                         help="require pt_lead>35, pt_sublead>25, |eta|<2.5 "
                              "(EB/EE gap excluded) for both photons. NOT a "
                              "full signal-region definition (no b-jet/MET "
                              "cuts) -- just photon-object acceptance.")
    args = parser.parse_args()

    out_pdf = args.out or (args.input.rsplit(".root", 1)[0] + "_photon_sf_type2.pdf")

    print(f"Opening {args.input} ...")
    f = uproot.open(args.input)
    if args.tree not in [k.split(";")[0] for k in f.keys()]:
        candidates = [k.split(";")[0] for k in f.keys()]
        print(f"Tree '{args.tree}' not found. Available keys: {candidates}")
        sys.exit(1)
    tree = f[args.tree]
    n_entries = tree.num_entries
    print(f"Tree '{args.tree}' has {n_entries} entries and {len(tree.keys())} branches.")

    photon_pt = safe_get(tree, "photon_pt")
    photon_eta = safe_get(tree, "photon_eta")
    photon_phi_branch = find_branch(tree, PHI_CANDIDATES, args.phi_branch)
    photon_phi = safe_get(tree, photon_phi_branch) if photon_phi_branch else None
    idmva_branch = find_branch(tree, IDMVA_CANDIDATES, args.idmva_branch)
    photon_idmva = safe_get(tree, idmva_branch) if idmva_branch else None
    w_photon = safe_get(tree, "w_photon")
    sys_photon = safe_get(tree, "sys_photon")

    mgg_branch = find_branch(tree, MGG_CANDIDATES, args.mgg_branch)
    mgg = safe_get(tree, mgg_branch) if mgg_branch else None
    ptgg_branch = find_branch(tree, PTGG_CANDIDATES, args.ptgg_branch)
    ptgg = safe_get(tree, ptgg_branch) if ptgg_branch else None
    drgg_branch = find_branch(tree, DRGG_CANDIDATES, args.drgg_branch)
    drgg = safe_get(tree, drgg_branch) if drgg_branch else None

    print(f"diphoton mass branch:  {mgg_branch or 'NOT FOUND (--mgg-branch)'}")
    print(f"diphoton pt branch:    {ptgg_branch or 'NOT FOUND (--ptgg-branch)'}")
    print(f"diphoton dR branch:    {drgg_branch or 'NOT FOUND (--drgg-branch)'}")
    print(f"photon phi branch:     {photon_phi_branch or 'NOT FOUND (--phi-branch)'}")
    print(f"photon idmva branch:   {idmva_branch or 'NOT FOUND (--idmva-branch)'}")

    weight_np, weight_used, weight_missing = get_weight(tree, args.weight_branches)
    weight_label = " * ".join(weight_used) if weight_used else "weight"
    print(f"weight = {weight_label}"
          + (f"   (missing, skipped: {', '.join(weight_missing)})" if weight_missing else ""))
    w_photon_np = ak.to_numpy(w_photon).astype(np.float64) if w_photon is not None else None

    # --- optional photon-SR-like acceptance cuts ---
    sr_mask = np.ones(n_entries, dtype=bool)
    if args.photon_sr_cuts and photon_pt is not None and photon_eta is not None:
        pt0 = full_length_value(photon_pt, 0, n_entries)
        pt1 = full_length_value(photon_pt, 1, n_entries)
        eta0 = full_length_value(photon_eta, 0, n_entries)
        eta1 = full_length_value(photon_eta, 1, n_entries)
        with np.errstate(invalid="ignore"):
            sr_mask = (pt0 > 35) & (pt1 > 25) & eb_ee_fiducial(eta0) & eb_ee_fiducial(eta1)
        sr_mask = np.nan_to_num(sr_mask.astype(float), nan=0.0).astype(bool)
    n_pass_sr = int(sr_mask.sum())

    weight_np_cut = apply_cut(weight_np, sr_mask) if weight_np is not None else None
    weight_none = np.ones(n_entries) if weight_np is not None else None
    weight_none_cut = apply_cut(weight_none, sr_mask) if weight_none is not None else None

    # --- yield summary (computed on the sr_mask-selected population) ---
    summary = []
    summary.append(f"unweighted N events (selected): {n_pass_sr} / {n_entries}")
    summary.append(f"weight = {weight_label}"
                    + (f"  (missing: {', '.join(weight_missing)})" if weight_missing else ""))
    if weight_np_cut is not None:
        wsum = np.nansum(weight_np_cut)
        summary.append(f"sum(weight):                {wsum:.6g}")
        summary.append(f"mean(weight):               {np.nanmean(weight_np_cut):.6g}")
    if w_photon_np is not None:
        w_photon_cut = apply_cut(w_photon_np, sr_mask)
        summary.append(f"mean(w_photon):             {np.nanmean(w_photon_cut):.6g}  "
                        f"(overall normalization shift from photon ID/CSEV SF: "
                        f"{100*(np.nanmean(w_photon_cut)-1):+.2f}%)")

    with PdfPages(out_pdf) as pdf:
        make_summary_page(pdf, args.input, args.tree, n_entries, args.photon_sr_cuts,
                           n_pass_sr, summary)

        # base weight for systematic reweighting (everything except w_photon)
        base_np = None
        if weight_np is not None and w_photon_np is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                base_np = weight_np / np.where(w_photon_np != 0, w_photon_np, np.nan)

        s0_full = s1_full = None
        if sys_photon is not None:
            s0_full = full_length_value(sys_photon, 0, n_entries)
            s1_full = full_length_value(sys_photon, 1, n_entries)

        def shape_page(title, xlabel, values_full, xrange, bins=50):
            vals_cut = apply_cut(values_full, sr_mask)
            ratio_hist_page(
                pdf, f"{title} (weight = {weight_label})", xlabel,
                [(vals_cut, weight_none_cut, "unweighted"),
                 (vals_cut, weight_np_cut, "weighted")],
                bins=bins, xrange=xrange)

        def sysband_page(title, xlabel, values_full, xrange, bins=50):
            if base_np is None or s0_full is None or s1_full is None:
                missing_page(pdf, f"{title}: sys_photon/w_photon/weight not "
                             f"all available -- skipping systematic-band plot.")
                return
            vals_cut = apply_cut(values_full, sr_mask)
            w_up = apply_cut(base_np * s0_full, sr_mask)
            w_down = apply_cut(base_np * s1_full, sr_mask)
            ratio_hist_page(
                pdf, f"{title}: nominal vs sys_photon up/down", xlabel,
                [(vals_cut, weight_np_cut, "nominal"),
                 (vals_cut, w_up, "sys_photon[0] (up)"),
                 (vals_cut, w_down, "sys_photon[1] (down)")],
                bins=bins, xrange=xrange)

        # --- per-photon variables: leading + subleading ---
        per_photon_vars = [
            ("pt", photon_pt, "p_T [GeV]", (0, 200)),
            ("eta", photon_eta, "eta", (-2.5, 2.5)),
            ("phi", photon_phi, "phi", (-3.2, 3.2)),
            ("idmva", photon_idmva, "ID MVA score", None),
        ]
        for name, arr, xlabel, xrange in per_photon_vars:
            if arr is None:
                missing_page(pdf, f"photon {name} branch not found -- skipping.")
                continue
            for idx, tag in [(0, "leading"), (1, "subleading")]:
                vals_full = full_length_value(arr, idx, n_entries)
                shape_page(f"{tag.capitalize()} photon {name}", xlabel, vals_full, xrange)

        # --- event-level diphoton variables ---
        if mgg is not None:
            mgg_full = full_length_value(mgg, 0, n_entries)
            shape_page("Diphoton mass", "m(gamma gamma) [GeV]", mgg_full, (100, 150))
            sysband_page("Diphoton mass", "m(gamma gamma) [GeV]", mgg_full, (100, 150))
        else:
            missing_page(pdf, "Diphoton mass branch not found -- rerun with "
                         "--mgg-branch <name> if your tree uses a different name.")

        if ptgg is not None:
            ptgg_full = full_length_value(ptgg, 0, n_entries)
            shape_page("Diphoton pt", "p_T(gamma gamma) [GeV]", ptgg_full, (0, 200))
            sysband_page("Diphoton pt", "p_T(gamma gamma) [GeV]", ptgg_full, (0, 200))
        else:
            missing_page(pdf, "Diphoton pt branch not found -- rerun with "
                         "--ptgg-branch <name> if your tree uses a different name.")

        if drgg is not None:
            drgg_full = full_length_value(drgg, 0, n_entries)
            shape_page("Diphoton dR", "dR(gamma, gamma)", drgg_full, (0, 4))
        else:
            missing_page(pdf, "Diphoton dR branch not found -- rerun with "
                         "--drgg-branch <name> if your tree uses a different name.")

        # --- leading photon pt systematic band (in addition to mgg/ptgg above) ---
        if photon_pt is not None:
            pt0_full = full_length_value(photon_pt, 0, n_entries)
            sysband_page("Leading photon pt", "leading photon pt [GeV]", pt0_full, (0, 200))

    print(f"\nWrote report: {out_pdf}")
    for line in summary:
        print(" ", line)


if __name__ == "__main__":
    main()
