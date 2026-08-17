"""Stage 1 of the Decipher-marker surrogate pipeline.

Stream the Durham expression Excel and dump the full gene x sample matrix (all
887 assayed samples) to NPZ. This is the full Durham expression SOURCE from
which Stage 2 extracts the fitted marker rows. Stage 2 does not
full-transcriptome-normalize Durham; it subsets to the 555 clinically eligible
patients and rank-maps each fitted marker across that cohort to its stored
training reference.

Fails closed on: a Durham sample count other than 887, blank/duplicate sample
IDs, and duplicated / unparseable / nonfinite rows for markers required by
config/decipher_marker_surrogate_v1.csv. Blank, None, NA and 'nan' gene symbols
and sample IDs are treated as missing, not as valid values.
"""
import os, csv, time, numpy as np, openpyxl as oxl

EXPECTED_DURHAM_N = 887


def _find_met_pca_root():
    env = os.environ.get("MET_PCA_ROOT")
    if env and os.path.isdir(env):
        return env
    d = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    while d != os.path.dirname(d):
        if os.path.isdir(os.path.join(d, "code")) and (
            os.path.isdir(os.path.join(d, "outs")) or os.path.isdir(os.path.join(d, "output"))):
            return d
        d = os.path.dirname(d)
    raise FileNotFoundError("Could not locate Met_PCa project root; set MET_PCA_ROOT.")


ROOT = _find_met_pca_root()
OUT = os.path.join(ROOT, "outs", "Decipher", "_durham_full_expr.npz")
CFG = os.path.join(ROOT, "config", "decipher_marker_surrogate_v1.csv")
os.makedirs(os.path.dirname(OUT), exist_ok=True)


def _missing_symbol(s):
    if s is None:
        return True
    t = str(s).strip()
    return t == "" or t.lower() in ("na", "nan")


def _required_symbols(path):
    """Gene-level-mappable canonical genes and their aliases (the surrogate's
    required markers)."""
    req = set()
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if str(r.get("gene_level_mappable", "")).strip().upper() != "TRUE":
                continue
            g = (r.get("canonical_gene") or "").strip()
            if g:
                req.add(g)
            for a in (r.get("gene_aliases") or "").split(";"):
                a = a.strip()
                if a:
                    req.add(a)
    return req


REQUIRED = _required_symbols(CFG)

t0 = time.time()
wb = oxl.load_workbook(os.path.join(ROOT, "data", "Durham_cohort_and_GRID_cohort",
                                     "Durham_cohort_011526.xlsx"),
                       read_only=True, data_only=True)
ws = wb["eset_gene_filtered"]
rows = ws.iter_rows(values_only=True)
header = next(rows)
raw_samples = header[2:]
# validate raw header cells BEFORE str() so a None cell is not masked as "None"
if any(_missing_symbol(s) for s in raw_samples):
    raise ValueError("Durham header contains a blank/None/NA sample ID")
samples = [str(s).strip() for s in raw_samples]
if len(samples) != len(set(samples)):
    raise ValueError("Durham header contains duplicate sample IDs")
if len(samples) != EXPECTED_DURHAM_N:
    raise ValueError(f"Expected {EXPECTED_DURHAM_N} Durham samples, got {len(samples)}")
n_cols = len(samples)

genes, mat_rows = [], []
n_skip_sym = n_skip_parse = 0
req_counts = {}
for row in rows:
    sym = row[1]
    if _missing_symbol(sym):
        n_skip_sym += 1
        continue
    sym = str(sym).strip()
    is_req = sym in REQUIRED
    try:
        vals = np.asarray(row[2:], dtype=float)
    except Exception:
        if is_req:
            raise ValueError(f"Required marker {sym}: unparseable expression row")
        n_skip_parse += 1
        continue
    if is_req:
        req_counts[sym] = req_counts.get(sym, 0) + 1
        if not np.all(np.isfinite(vals)):
            raise ValueError(f"Required marker {sym}: nonfinite expression values")
    genes.append(sym)
    mat_rows.append(vals)

if not genes:
    raise ValueError("No gene rows parsed from eset_gene_filtered; check the sheet layout.")
dup_req = sorted(g for g, c in req_counts.items() if c > 1)
if dup_req:
    raise ValueError(f"Duplicated required-marker rows: {dup_req}")

mat = np.vstack(mat_rows).astype(np.float64)   # float64 so R/Python agree exactly
if mat.shape != (len(genes), n_cols):
    raise ValueError(f"matrix shape {mat.shape} != ({len(genes)}, {n_cols})")

print(f"Rows skipped: {n_skip_sym} (missing symbol), {n_skip_parse} (unparseable non-required)")
print(f"Required markers found in Durham: {sorted(req_counts)} ({len(req_counts)}/{len(REQUIRED)} config symbols)")
print(f"Loaded full Durham expression {mat.shape} in {time.time()-t0:.1f}s")
# atomic write via a same-directory temp file (large NPZ)
tmp = OUT + ".tmp"
with open(tmp, "wb") as fh:
    np.savez_compressed(fh, mat=mat, genes=np.array(genes, dtype=object),
                        samples=np.array(samples, dtype=object))
os.replace(tmp, OUT)
print(f"Saved -> {OUT}  ({os.path.getsize(OUT)/1e6:.1f} MB)")
