# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 10:01:50 2025

@author: 李津莹
"""
import os
import sys
import shutil
import hashlib
import time
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple, Dict

import pandas as pd
from Bio import SeqIO
import streamlit as st
import subprocess

# =============================
# App Config
# =============================
st.set_page_config(page_title="Species-Specific Marker Finder", layout="wide")

# -----------------------------
# Utilities
# -----------------------------
def which(cmd: str) -> str:
    """Return absolute path to an executable or raise RuntimeError."""
    path = shutil.which(cmd)
    if not path:
        raise RuntimeError(f"Executable not found in PATH: {cmd}")
    return path


def run_cmd(cmd: List[str], timeout: int = 1800) -> Tuple[int, str, str]:
    """Run a command safely with timeout; return (returncode, stdout, stderr)."""
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        out, err = proc.communicate(timeout=timeout)
        return proc.returncode, out, err
    except subprocess.TimeoutExpired:
        proc.kill()
        return -9, "", f"Timeout running: {' '.join(cmd)}"


def safe_write_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(text)


def hash_files(paths: List[Path]) -> str:
    h = hashlib.sha256()
    for p in paths:
        h.update(str(p).encode())
        with open(p, 'rb') as f:
            while True:
                chunk = f.read(1024 * 1024)
                if not chunk:
                    break
                h.update(chunk)
    return h.hexdigest()[:16]


# -----------------------------
# FASTA helpers
# -----------------------------
def concat_fastas(fastas: List[Path], out_fa: Path) -> int:
    count = 0
    with open(out_fa, 'w') as out:
        for fa in fastas:
            for rec in SeqIO.parse(str(fa), 'fasta'):
                SeqIO.write(rec, out, 'fasta')
                count += 1
    return count


def sliding_windows(seq: str, win: int, step: int) -> List[Tuple[int, int, str]]:
    arr = []
    n = len(seq)
    if n < win:
        return [(0, n, seq)]
    i = 0
    while i + win <= n:
        arr.append((i, i + win, seq[i:i + win]))
        i += step
    if arr and arr[-1][1] < n:
        arr.append((n - win, n, seq[n - win:n]))
    return arr


def write_slices(fasta: Path, win: int, step: int, out_fa: Path, info_csv: Path) -> int:
    """Slice sequences with sliding window; write FASTA + info CSV; return slice count."""
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    info_csv.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    total = 0
    with open(out_fa, 'w') as fout:
        for rec in SeqIO.parse(str(fasta), 'fasta'):
            windows = sliding_windows(str(rec.seq), win, step)
            for s, e, subseq in windows:
                sid = f"{rec.id}|{s+1}-{e}"
                fout.write(f">{sid}\n{subseq}\n")
                rows.append({
                    'record_id': rec.id,
                    'start': s + 1,
                    'end': e,
                    'length': e - s,
                    'slice_id': sid
                })
                total += 1
    pd.DataFrame(rows).to_csv(info_csv, index=False)
    return total


# -----------------------------
# BLAST helpers
# -----------------------------
def make_blastdb(fasta: Path, db_prefix: Path, dbtype: str = 'nucl'):
    makeblastdb = which('makeblastdb')
    cmd = [makeblastdb, '-in', str(fasta), '-dbtype', dbtype, '-out', str(db_prefix), '-parse_seqids']
    rc, out, err = run_cmd(cmd)
    if rc != 0:
        raise RuntimeError(f"makeblastdb failed: {err}\n{out}")


def blastn_query(query_fa: Path, db_prefix: Path, out_tsv: Path, evalue: float = 1e-5, max_target_seqs: int = 5):
    blastn = which('blastn')
    cmd = [
        blastn, '-query', str(query_fa), '-db', str(db_prefix),
        '-evalue', str(evalue), '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend sstart send bitscore evalue',
        '-max_target_seqs', str(max_target_seqs)
    ]
    rc, out, err = run_cmd(cmd)
    if rc != 0:
        raise RuntimeError(f"blastn failed: {err}\n{out}")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, 'w') as f:
        f.write(out)


def parse_hits(tsv: Path) -> pd.DataFrame:
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'evalue']
    if not tsv.exists() or tsv.stat().st_size == 0:
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(tsv, sep='\t', header=None, names=cols)
    df['qcov'] = df['length'] / df['qlen']
    return df


# -----------------------------
# Specificity filter pipeline
# -----------------------------
def filter_specific_slices(
    slices_fa: Path,
    nontarget_db: Path,
    identity_thr: float,
    coverage_thr: float,
    threads: int,
    chunk_size: int = 1000,
    workdir: Path = Path('.')
) -> Tuple[Path, Path]:
    """Return (specific_fasta, specific_table_csv)."""
    all_records = list(SeqIO.parse(str(slices_fa), 'fasta'))
    if not all_records:
        raise RuntimeError("No slices to process.")
    tmpdir = workdir / 'tmp_chunks'
    tmpdir.mkdir(parents=True, exist_ok=True)

    chunks = [all_records[i:i + chunk_size] for i in range(0, len(all_records), chunk_size)]

    def worker(i_ch: int, recs):
        qfa = tmpdir / f"chunk_{i_ch:05d}.fa"
        with open(qfa, 'w') as f:
            SeqIO.write(recs, f, 'fasta')
        tsv = tmpdir / f"chunk_{i_ch:05d}.blast.tsv"
        blastn_query(qfa, nontarget_db, tsv)
        hits = parse_hits(tsv)
        keep_records, keep_rows = [], []
        grouped = hits.groupby('qseqid') if not hits.empty else {}
        for r in recs:
            dfq = grouped.get_group(r.id) if (not isinstance(grouped, dict) and r.id in grouped.groups) else pd.DataFrame()
            drop = False
            if not dfq.empty:
                cond = (dfq['pident'] / 100.0 >= identity_thr) & (dfq['qcov'] >= coverage_thr)
                if cond.any():
                    drop = True
            if not drop:
                keep_records.append(r)
                keep_rows.append({'slice_id': r.id, 'length': len(r.seq)})
        return keep_records, keep_rows

    kept_records, kept_rows = [], []
    progress = st.progress(0.0, text="Filtering slices against non-target DB…")
    completed = 0
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = {ex.submit(worker, i, ch): i for i, ch in enumerate(chunks)}
        for fut in as_completed(futs):
            recs, rows = fut.result()
            kept_records.extend(recs)
            kept_rows.extend(rows)
            completed += 1
            progress.progress(completed / len(chunks), text=f"Filtering {completed}/{len(chunks)} chunks…")

    sp_fa = workdir / 'specific_slices.fasta'
    sp_csv = workdir / 'specific_slices.csv'
    with open(sp_fa, 'w') as f:
        SeqIO.write(kept_records, f, 'fasta')
    pd.DataFrame(kept_rows).to_csv(sp_csv, index=False)
    return sp_fa, sp_csv


# -----------------------------
# Coverage computation
# -----------------------------
def compute_coverage(markers_fa: Path, target_fastas: List[Path], identity_thr: float, coverage_thr: float, threads: int, workdir: Path) -> Tuple[pd.DataFrame, Path]:
    tgt_concat = workdir / 'targets.concat.fa'
    concat_fastas(target_fastas, tgt_concat)
    db_prefix = workdir / 'targets_db'
    make_blastdb(tgt_concat, db_prefix)
    tsv = workdir / 'markers_vs_targets.tsv'
    blastn_query(markers_fa, db_prefix, tsv, evalue=1e-5, max_target_seqs=100000)
    hits = parse_hits(tsv)
    if hits.empty:
        return pd.DataFrame({'marker': [], 'genomes_hit': [], 'hit_fraction': []}), tsv

    hits['ok'] = (hits['pident'] / 100.0 >= identity_thr) & (hits['qcov'] >= coverage_thr)
    okhits = hits[hits['ok']]
    okhits['sseqid_clean'] = okhits['sseqid'].astype(str).str.split().str[0]
    grp = okhits.groupby('qseqid')['sseqid_clean'].nunique().reset_index(name='genomes_hit')
    total_genomes = len(target_fastas)
    grp['hit_fraction'] = grp['genomes_hit'] / total_genomes
    grp = grp.rename(columns={'qseqid': 'marker'})

    all_markers = [rec.id for rec in SeqIO.parse(str(markers_fa), 'fasta')]
    missing = set(all_markers) - set(grp['marker'])
    if missing:
        grp = pd.concat([grp, pd.DataFrame({'marker': list(missing), 'genomes_hit': 0, 'hit_fraction': 0.0})], ignore_index=True)

    grp = grp.sort_values(['genomes_hit', 'marker'], ascending=[False, True]).reset_index(drop=True)
    return grp, tsv


# -----------------------------
# Streamlit UI
# -----------------------------
st.title("🔬 Species-Specific Marker Finder & Coverage Analyzer")
st.caption("Upload target and non-target genomes (FASTA/FNA). The app will find species-specific slices and compute coverage across target strains.")

with st.sidebar:
    st.header("⚙️ Parameters")
    win = st.number_input("Slice length (bp)", min_value=50, max_value=5000, value=500, step=50)
    step = st.number_input("Step (bp)", min_value=10, max_value=2000, value=100, step=10)
    ident_thr = st.slider("Non-target filter: identity ≥", 0.1, 1.0, 0.3, 0.05)
    cov_thr = st.slider("Non-target filter: coverage ≥", 0.1, 1.0, 0.3, 0.05)
    threads = st.number_input("Threads", 1, 32, 8, 1)
    st.markdown("---")
    cov_ident_thr = st.slider("Coverage calc: identity ≥", 0.5, 1.0, 0.9, 0.05)
    cov_cov_thr = st.slider("Coverage calc: coverage ≥", 0.3, 1.0, 0.9, 0.05)

st.subheader("1) Upload Inputs")
col1, col2 = st.columns(2)
with col1:
    target_files = st.file_uploader("Target genomes (same species)", type=["fa", "fasta", "fna"], accept_multiple_files=True)
with col2:
    nontarget_files = st.file_uploader("Non-target genomes (other species)", type=["fa", "fasta", "fna"], accept_multiple_files=True)

work_root = Path(tempfile.gettempdir()) / "marker_finder_app"
work_root.mkdir(parents=True, exist_ok=True)

if st.button("🚀 Run Specificity Screening & Coverage", type="primary"):
    if not target_files or not nontarget_files:
        st.error("Please upload both target and non-target genomes.")
        st.stop()

    run_id = time.strftime("%Y%m%d_%H%M%S")
    rundir = work_root / run_id
    tgt_dir = rundir / 'targets'
    non_dir = rundir / 'non_targets'
    for d in (tgt_dir, non_dir):
        d.mkdir(parents=True, exist_ok=True)

    st.write("Saving uploaded files…")
    tgt_paths, non_paths = [], []
    for uf in target_files:
        p = tgt_dir / uf.name
        with open(p, 'wb') as f:
            f.write(uf.getbuffer())
        tgt_paths.append(p)
    for uf in nontarget_files:
        p = non_dir / uf.name
        with open(p, 'wb') as f:
            f.write(uf.getbuffer())
        non_paths.append(p)

    st.subheader("2) Build Non-target Database")
    non_concat = rundir / 'non_targets.concat.fa'
    st.write("Concatenating non-target genomes…")
    nrecs = concat_fastas(non_paths, non_concat)
    st.write(f"Non-target total records: **{nrecs}**")
    st.write("Creating BLAST database…")
    non_db = rundir / 'non_target_db'
    make_blastdb(non_concat, non_db)
    st.success("Non-target BLAST database ready.")

    st.subheader("3) Slice Target Genomes")
    tgt_concat = rundir / 'targets.concat.fa'
    nrecs_t = concat_fastas(tgt_paths, tgt_concat)
    st.write(f"Target total records: **{nrecs_t}**")
    slices_fa = rundir / 'target_slices.fa'
    slices_info = rundir / 'target_slices.csv'
    nslices = write_slices(tgt_concat, win, step, slices_fa, slices_info)
    st.success(f"Total slices: {nslices}")

    st.subheader("4) Specificity Filtering vs Non-target DB")
    sp_fa, sp_csv = filter_specific_slices(
        slices_fa=slices_fa,
        nontarget_db=non_db,
        identity_thr=ident_thr,
        coverage_thr=cov_thr,
        threads=int(threads),
        chunk_size=2000,
        workdir=rundir,
    )
    st.success(f"Specific markers retained: {sum(1 for _ in SeqIO.parse(str(sp_fa), 'fasta'))}")
    st.download_button("Download specific markers (FASTA)", data=open(sp_fa, 'rb').read(), file_name="specific_markers.fasta")
    st.download_button("Download marker table (CSV)", data=open(sp_csv, 'rb').read(), file_name="specific_markers.csv")

    st.subheader("5) Coverage Across Target Genomes")
    cov_df, hits_tsv = compute_coverage(
        markers_fa=sp_fa,
        target_fastas=tgt_paths,
        identity_thr=cov_ident_thr,
        coverage_thr=cov_cov_thr,
        threads=int(threads),
        workdir=rundir / 'coverage'
    )
    st.dataframe(cov_df.head(50))
    st.download_button("Download coverage table (CSV)", data=cov_df.to_csv(index=False).encode(), file_name="marker_coverage.csv")
    st.download_button("Download raw BLAST hits (TSV)", data=open(hits_tsv, 'rb').read(), file_name="markers_vs_targets.tsv")

st.markdown("""
---
### 📝 Notes
- Requires **BLAST+** (`makeblastdb`, `blastn`) in PATH.
- Specificity rule: slice is **discarded** if any hit in non-target DB has `pident ≥ threshold` and `coverage ≥ threshold`.
- Coverage = number of target genomes where each marker aligns with given thresholds.
- Temporary results are stored under `/tmp` (auto-created).
""")

