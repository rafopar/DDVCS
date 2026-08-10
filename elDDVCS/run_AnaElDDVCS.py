#!/usr/bin/env python3
"""
run_AnaElDDVCS.py

Usage:
    python run_AnaElDDVCS.py <dataset>

    dataset: 8p477GeV or 6p395GeV

Scans SkimDDVCS/<dataset>/ for files named elDDVCS_<RUN>.hipo, extracts each
run number, and launches

    bin/AnaElDDVCS.exe <RUN> <dataset>

for every unique run. At most MAX_PARALLEL jobs run concurrently; additional
runs queue up until a slot opens.

When all jobs finish, hadds the per-run ROOT files into
    OutHists/AnaElDDVCS/DDVCS_RGK_<dataset>_All.root

Run this script from the installation directory where SkimDDVCS/<dataset>/
exists.
"""

import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

EXECUTABLE   = "bin/AnaElDDVCS.exe"
MAX_PARALLEL = 18

FILE_PAT = re.compile(r"^elDDVCS_(\d+)\.hipo$")


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------

def run_ana(run, dataset):
    """
    Launch bin/AnaElDDVCS.exe <run> <dataset> and wait for it to finish.
    Returns (run, returncode, elapsed_seconds, stderr_tail).
    """
    t0 = time.time()
    try:
        result = subprocess.run(
            [EXECUTABLE, str(run), dataset],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
        rc  = result.returncode
        err = result.stderr or ""
    except FileNotFoundError:
        return run, -1, 0.0, f"executable not found: {EXECUTABLE}"
    except OSError as exc:
        return run, -1, 0.0, f"launch failed: {exc}"

    elapsed = time.time() - t0
    tail = "\n".join(err.strip().splitlines()[-5:])
    return run, rc, elapsed, tail


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def fmt_seconds(s):
    if s < 60:
        return f"{s:.1f}s"
    m, s = divmod(s, 60)
    if m < 60:
        return f"{int(m)}m{int(s):02d}s"
    h, m = divmod(m, 60)
    return f"{int(h)}h{int(m):02d}m{int(s):02d}s"


def print_table(title, headers, rows):
    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            col_widths[i] = max(col_widths[i], len(str(cell)))
    sep = "+" + "+".join("-" * (w + 2) for w in col_widths) + "+"
    def fmt(cells):
        return "|" + "|".join(
            f" {str(c):<{col_widths[i]}} " for i, c in enumerate(cells)
        ) + "|"
    print(title)
    print(sep)
    print(fmt(headers))
    print(sep)
    for row in rows:
        print(fmt(row))
    print(sep)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print("Usage: python run_AnaElDDVCS.py <dataset>")
        print("  dataset: 8p477GeV or 6p395GeV")
        sys.exit(1)

    dataset  = sys.argv[1]
    skim_dir = f"SkimDDVCS/{dataset}"

    if not os.path.isdir(skim_dir):
        print(f"Error: directory '{skim_dir}' not found.")
        sys.exit(1)
    if not os.path.isfile(EXECUTABLE) or not os.access(EXECUTABLE, os.X_OK):
        print(f"Error: '{EXECUTABLE}' is not an executable file.")
        sys.exit(1)

    runs = set()
    for entry in os.scandir(skim_dir):
        if not entry.is_file():
            continue
        m = FILE_PAT.match(entry.name)
        if m:
            runs.add(int(m.group(1)))

    if not runs:
        print(f"No files matching elDDVCS_<RUN>.hipo found in {skim_dir}/.")
        sys.exit(0)

    runs = sorted(runs)
    n_total = len(runs)

    print(f"Found {n_total} unique run(s) in {skim_dir}/")
    print(f"Launching {EXECUTABLE} <RUN> {dataset}  with up to {MAX_PARALLEL} concurrent job(s).")
    print("=" * 70)

    results  = {}
    t_start  = time.time()
    done_cnt = 0
    width    = len(str(n_total))

    with ThreadPoolExecutor(max_workers=MAX_PARALLEL) as pool:
        futures = {pool.submit(run_ana, r, dataset): r for r in runs}

        for fut in as_completed(futures):
            run, rc, elapsed, tail = fut.result()
            results[run] = (rc, elapsed, tail)
            done_cnt += 1
            status = "OK" if rc == 0 else f"FAIL(rc={rc})"
            print(f"[{done_cnt:>{width}}/{n_total}]  Run {run:>7}  "
                  f"{status:>10}   {fmt_seconds(elapsed)}")

    total_elapsed = time.time() - t_start

    n_ok   = sum(1 for rc, _, _ in results.values() if rc == 0)
    n_fail = n_total - n_ok

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Total runs        : {n_total}")
    print(f"  Successful        : {n_ok}")
    print(f"  Failed            : {n_fail}")
    print(f"  Wall-clock time   : {fmt_seconds(total_elapsed)}")

    failed = [(r, rc, elapsed, tail)
              for r, (rc, elapsed, tail) in results.items() if rc != 0]
    if failed:
        failed.sort(key=lambda x: x[0])
        print()
        print_table(
            "Failed runs:",
            ["Run", "Exit code", "Elapsed"],
            [(r, rc, fmt_seconds(el)) for r, rc, el, _ in failed],
        )
        print()
        print("Stderr tail for failed runs:")
        for r, rc, _, tail in failed:
            print(f"  Run {r}:")
            if tail:
                for line in tail.splitlines():
                    print(f"    {line}")
            else:
                print("    (no stderr output)")

    # ── hadd all per-run ROOT files into a single output ──────────────────────
    ok_runs    = sorted(r for r, (rc, _, _) in results.items() if rc == 0)
    input_root = [f"OutHists/AnaElDDVCS/Hists_{r}.root" for r in ok_runs]
    output_root = f"OutHists/AnaElDDVCS/DDVCS_RGK_{dataset}_All.root"

    print()
    print("=" * 70)
    if not input_root:
        print("No successful runs — skipping hadd.")
    else:
        print(f"Running hadd → {output_root}")
        hadd_result = subprocess.run(
            ["hadd", "-f", output_root] + input_root,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if hadd_result.returncode == 0:
            print("hadd completed successfully.")
        else:
            print(f"hadd FAILED (rc={hadd_result.returncode})")
            if hadd_result.stderr:
                print(hadd_result.stderr.strip())


if __name__ == "__main__":
    main()
