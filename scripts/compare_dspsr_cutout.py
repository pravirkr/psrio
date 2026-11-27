#!/usr/bin/env python3
"""Check `psrio cutout` against dspsr on a search-mode PSRFITS file.

Development tool only (needs dspsr, astropy, numpy, and optionally psrchive's
Python bindings). For each case it

1. picks a window on dspsr's sub-integration grid (one SUBINT row in, so the
   window lies in dspsr's first block),
2. runs `psrio cutout -S ...` and `dspsr -S -T -c -cepoch -b -D [-K] -s -A
   -scloffs` over the same data,
3. compares our single sub-integration with dspsr's matching one:
   DATA, DAT_SCL, DAT_OFFS, DAT_WTS, DAT_FREQ, TSUBINT, PERIOD byte for
   byte, the absolute epoch (STT + OFFS_SUB) to 1 ns, and every header and
   HISTORY value except provenance (DATE, DATE_PRO, PROC_CMD, NSUB, NAXIS2),
4. loads our archive with psrchive and compares its amplitudes with dspsr's.

dspsr is run with HOME pointing at a directory holding an empty
.psrchive.cfg (a missing one aborts psrchive on some macOS builds).

Usage:
    scripts/compare_dspsr_cutout.py DATA.sf [--psrio build/applications/psrio]
        [--dm 413] [--work DIR] [--case nok-1 k-1 nok-8 k-8] [--keep]
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits

DM_CONST = 2.41e-4


@dataclass
class Case:
    name: str
    dedisperse: bool
    nbin: int
    samples_per_bin: int


CASES = {
    "nok-1": Case("nok-1", False, 1024, 1),
    "k-1": Case("k-1", True, 1024, 1),
    "nok-8": Case("nok-8", False, 256, 8),
    "k-8": Case("k-8", True, 256, 8),
}


def band_delays(obsfreq: float, obsbw: float, nchan: int, rate: float, dm: float):
    """dspsr's per-channel sample delays, span centre, and zero delay."""
    base = obsfreq - 0.5 * obsbw + 0.5 * obsbw / nchan
    freqs = [base + (c * obsbw / nchan) for c in range(nchan)]
    disp = dm / DM_CONST

    def delay(f: float) -> int:
        return math.floor(disp * (1.0 / obsfreq**2 - 1.0 / f**2) * rate + 0.5)

    d = [delay(f) for f in freqs]
    mid = (freqs[0] + freqs[-1]) / 2
    span = math.floor(disp * (1.0 / obsfreq**2 - 1.0 / mid**2) * rate + 0.5)
    zero = max(0, max(x - span for x in d))
    applied = [span + zero - x for x in d]
    return d, span, zero, applied


def run(cmd, env=None, cwd=None, quiet=False):
    if not quiet:
        print("+", " ".join(str(c) for c in cmd), flush=True)
    res = subprocess.run(cmd, env=env, cwd=cwd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout[-4000:], res.stderr[-4000:], sep="\n")
        raise SystemExit(f"command failed ({res.returncode})")
    return res


def epoch_seconds(hdr, offs_sub, ref_imjd):
    """Seconds of STT + OFFS_SUB after day ref_imjd, keeping precision."""
    return (
        (hdr["STT_IMJD"] - ref_imjd) * 86400.0
        + hdr["STT_SMJD"]
        + hdr["STT_OFFS"]
        + offs_sub
    )


def compare(ours_path: Path, ref_path: Path, subint: int, home: Path) -> list[str]:
    problems: list[str] = []
    ours = fits.open(ours_path, memmap=True)
    ref = fits.open(ref_path, memmap=True)

    names_o = [h.name for h in ours]
    names_r = [h.name for h in ref]
    if names_o != names_r:
        problems.append(f"HDUs differ: {names_o} vs {names_r}")

    # Primary header values.
    po, pr = ours[0].header, ref[0].header
    skip = {"DATE", "COMMENT", ""}
    same_start = subint == 0
    if not same_start:
        skip |= {"STT_IMJD", "STT_SMJD", "STT_OFFS"}
    for key in pr:
        if key in skip:
            continue
        if key not in po:
            problems.append(f"PRIMARY {key} missing")
        elif po[key] != pr[key]:
            problems.append(f"PRIMARY {key}: {po[key]!r} != {pr[key]!r}")
    for key in po:
        if key not in skip and key not in pr:
            problems.append(f"PRIMARY {key} is extra")

    # SUBINT header values.
    so, sr = ours["SUBINT"].header, ref["SUBINT"].header
    for key in sr:
        if key in {"NAXIS2", "COMMENT", ""} or key.startswith(("TTYPE", "TFORM", "TUNIT", "TDIM")):
            continue
        if key not in so:
            problems.append(f"SUBINT {key} missing")
        elif so[key] != sr[key]:
            problems.append(f"SUBINT {key}: {so[key]!r} != {sr[key]!r}")
    co = [(c.name, c.format, c.unit, c.dim) for c in ours["SUBINT"].columns]
    cr = [(c.name, c.format, c.unit, c.dim) for c in ref["SUBINT"].columns]
    if co != cr:
        problems.append(f"SUBINT columns differ:\n  {co}\n  {cr}")

    # HISTORY: the last row describes the archive.
    ho, hr = ours["HISTORY"].data, ref["HISTORY"].data
    for name in hr.columns.names:
        if name in {"DATE_PRO", "PROC_CMD", "NSUB"}:
            continue
        a, b = ho[name][-1], hr[name][-1]
        if a != b:
            problems.append(f"HISTORY {name}: {a!r} != {b!r}")

    # SUBINT columns, byte for byte.
    do, dr = ours["SUBINT"].data, ref["SUBINT"].data
    if len(do) != 1:
        problems.append(f"ours has {len(do)} sub-integrations")
    for name in ["DAT_FREQ", "DAT_WTS", "DAT_OFFS", "DAT_SCL", "DATA", "TSUBINT", "PERIOD"]:
        a = np.asarray(do[name][0])
        b = np.asarray(dr[name][subint])
        if a.tobytes() != b.tobytes():
            diff = np.count_nonzero(a != b)
            extra = ""
            if name == "DATA":
                extra = f", max |diff| {np.max(np.abs(a.astype(np.int32) - b.astype(np.int32)))}"
            problems.append(f"{name} differs in {diff} of {a.size} values{extra}")
    if same_start:
        if do["OFFS_SUB"][0] != dr["OFFS_SUB"][subint]:
            problems.append(f"OFFS_SUB {do['OFFS_SUB'][0]!r} != {dr['OFFS_SUB'][subint]!r}")
    t_o = epoch_seconds(po, do["OFFS_SUB"][0], pr["STT_IMJD"])
    t_r = epoch_seconds(pr, dr["OFFS_SUB"][subint], pr["STT_IMJD"])
    if abs(t_o - t_r) > 1e-9:
        problems.append(f"epoch differs by {t_o - t_r:.3e} s")
    print(f"  epoch difference {t_o - t_r:.3e} s")

    # psrchive reads our archive and gets dspsr's metadata and amplitudes.
    # (The command-line tools are used: the Python bindings abort on some
    # macOS builds while loading ObsDescription.)
    env = dict(os.environ, HOME=str(home))
    fields = "nbin,nchan,npol,state,dm,dmc,freq,bw,period,name,site,rcvr:name,be:name,length"
    vo = run(["vap", "-n", "-c", fields, str(ours_path)], env=env, quiet=True).stdout.split()[1:]
    vr = run(["vap", "-n", "-c", fields, str(ref_path)], env=env, quiet=True).stdout.split()[1:]
    for key, a, b in zip(fields.split(","), vo, vr):
        if key == "length":
            continue  # whole-archive length; dspsr's holds every subint
        if a != b:
            problems.append(f"vap {key}: {a} != {b}")
    for chan in (0, so["NCHAN"] // 3, so["NCHAN"] - 1):
        to = run(["pdv", "-R", "-t", "-i", "0", "-n", str(chan), str(ours_path)], env=env, quiet=True).stdout
        tr = run(["pdv", "-R", "-t", "-i", str(subint), "-n", str(chan), str(ref_path)], env=env, quiet=True).stdout
        lo = [line for line in to.splitlines() if line[:1].isdigit()]
        lr = [line.split(None, 1)[1] for line in tr.splitlines() if line[:1].isdigit()]
        lo = [line.split(None, 1)[1] for line in lo]
        if lo != lr or not lo:
            problems.append(f"pdv channel {chan} differs")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("data", type=Path, nargs="+", help="search-mode PSRFITS file(s)")
    parser.add_argument("--psrio", type=Path, default=Path("build/applications/psrio"))
    parser.add_argument("--dspsr", default="dspsr")
    parser.add_argument("--dm", type=float, default=413.0)
    parser.add_argument("--row", type=int, default=10, help="SUBINT row dspsr starts reading at")
    parser.add_argument("--skip", type=int, default=3,
                        help="dspsr sub-integrations before the compared one")
    parser.add_argument("--plain", action="store_true",
                        help="leave out --scloffs on both sides")
    parser.add_argument("--case", nargs="+", default=list(CASES), choices=list(CASES))
    parser.add_argument("--work", type=Path, default=None)
    parser.add_argument("--keep", action="store_true", help="keep the archives")
    args = parser.parse_args()

    work = args.work or Path(tempfile.mkdtemp(prefix="psrio_cutout_"))
    work.mkdir(parents=True, exist_ok=True)
    home = work / "home"
    home.mkdir(exist_ok=True)
    (home / ".psrchive.cfg").touch()
    env = dict(os.environ, HOME=str(home))

    data = [p.resolve() for p in args.data]
    with fits.open(data[0], memmap=True) as hdul:
        primary, sub = hdul[0].header, hdul["SUBINT"].header
        tbin, nsblk, nchan = sub["TBIN"], sub["NSBLK"], sub["NCHAN"]
        obsfreq, obsbw = primary["OBSFREQ"], primary["OBSBW"]
        total = sum(fits.getheader(p, "SUBINT")["NAXIS2"] for p in data) * nsblk
    rate = 1.0 / tbin

    failures = 0
    for name in args.case:
        case = CASES[name]
        n = case.nbin * case.samples_per_bin
        _, span, zero, applied = band_delays(obsfreq, obsbw, nchan, rate, args.dm)
        if not case.dedisperse:
            span, zero, applied = 0, 0, [0]
        load = args.row * nsblk
        # dspsr divides output time from STT (+zero/rate): output sample o is
        # in division o // n. The first whole division at or after the row.
        k = -(-load // n) + args.skip
        first = k * n
        sub_index = (first - load) // n
        last_input = first + n + max(applied)
        end = -(-(last_input + 1) // nsblk) * nsblk
        if end > total:
            raise SystemExit(f"{name}: the window needs {end} samples, data has {total}")
        period = n / rate
        seek = load / rate
        # dspsr ends at uint64((S+T)*rate), floored to a row: pad by half a
        # sample so the truncation cannot drop the last row.
        length = (end - load + 0.5) / rate
        # Time stamp of the window start, for -cepoch.
        imjd, smjd, offs = primary["STT_IMJD"], primary["STT_SMJD"], primary["STT_OFFS"]
        cepoch_sec = smjd + offs + (load / rate) + ((zero / rate) if case.dedisperse else 0.0) + (first - load) / rate
        day_frac = f"{cepoch_sec / 86400.0:.15f}"
        cepoch = f"{imjd + int(float(day_frac))}{day_frac[day_frac.index('.'):]}"

        print(f"== {name}: {n} samples, {case.nbin} bins, window at output sample {first}"
              f" (dspsr sub-integration {sub_index})")
        case_dir = work / name
        case_dir.mkdir(exist_ok=True)
        ours = case_dir / "ours.ar"
        ref = case_dir / "dspsr"
        cmd = [str(args.psrio), "cutout", *map(str, data), "-D", repr(args.dm), "-b", str(case.nbin),
               "-S", repr((first + span + zero) / rate), "-T", repr(period), "-O", str(ours)]
        if not args.plain:
            cmd.append("--scloffs")
        if case.dedisperse:
            cmd.append("-K")
        run(cmd)
        cmd = [args.dspsr, "-S", repr(seek), "-T", repr(length), "-c", repr(period),
               "-cepoch", cepoch, "-b", str(case.nbin), "-D", repr(args.dm), "-s", "-A",
               "-U", "16384", "-O", str(ref), *map(str, data)]
        if not args.plain:
            cmd.insert(1, "-scloffs")
        if case.dedisperse:
            cmd.insert(1, "-K")
        run(cmd, env=env, cwd=case_dir)
        problems = compare(ours, case_dir / "dspsr.ar", sub_index, home)
        if problems:
            failures += 1
            print(f"FAIL {name}")
            for p in problems:
                print("   ", p)
        else:
            print(f"PASS {name}")
        if not args.keep:
            shutil.rmtree(case_dir, ignore_errors=True)
    print(f"{len(args.case) - failures}/{len(args.case)} cases match dspsr; work dir {work}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
