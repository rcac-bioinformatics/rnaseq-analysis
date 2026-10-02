#!/usr/bin/env python3
"""build_kit.py: generate the Negishi run kit from the lesson episodes.

Run from anywhere:  python3 negishi-run/build_kit.py
Writes negishi-run/generated/ and negishi-run/staged-scripts/ (both replaced on every run).

Code is never hand-copied: every learner command the kit runs is extracted from the
episodes (same fence parsing as local-test/harness/extract.py) and located by a content
anchor, so the build fails loudly if an episode changes shape. Output is deterministic:
no timestamps, sorted listings, fixed line endings. Rerunning on unchanged episodes gives
byte-identical files (checked with --check).
"""
import hashlib
import pathlib
import re
import shlex
import shutil
import sys

KIT = pathlib.Path(__file__).resolve().parent
ROOT = KIT.parent
GEN = KIT / "generated"
STAGED_OUT = KIT / "staged-scripts"
EPISODES = ROOT / "episodes"
SETUP = ROOT / "learners" / "setup.md"

sys.dont_write_bytecode = True  # keep the kit directory free of __pycache__
sys.path.insert(0, str(KIT / "py"))
import summarize  # noqa: E402  (metric handlers; also checked for coverage below)

# ---------------------------------------------------------------- parsing
OPEN_RE = re.compile(r'^(\s*)```\s*([A-Za-z]*)\s*$')
CLOSE_RE = re.compile(r'^\s*```\s*$')
DIV_OPEN = re.compile(r'^:{3,}\s*\{?\.?\s*([\w-]+)')
DIV_CLOSE = re.compile(r'^:{3,}\s*$')
FIG_RE = re.compile(r'include_graphics\("([^"]+)"\)')
MARK_RE = re.compile(r'<!-- NEGISHI:([a-z0-9-]+) -->')


class Block(object):
    def __init__(self, ep, nn, lang, start, end, context, text):
        self.ep, self.nn, self.lang = ep, nn, lang
        self.start, self.end, self.context, self.text = start, end, context, text

    @property
    def ref(self):
        return "%s/%02d" % (self.ep, self.nn)


def parse(path, ep):
    """Return (blocks, figures[(line, path)], markers[(line, id, own_line)])."""
    lines = path.read_text().splitlines()
    blocks, figs, marks, divs = [], [], [], []
    n, i = 0, 0
    while i < len(lines):
        L = lines[i]
        for m in MARK_RE.finditer(L):
            marks.append((i + 1, m.group(1), L.strip() == m.group(0)))
        if DIV_OPEN.match(L):
            divs.append(DIV_OPEN.match(L).group(1)); i += 1; continue
        if DIV_CLOSE.match(L):
            if divs:
                divs.pop()
            i += 1; continue
        if L.lstrip().startswith("```{"):  # knitr chunk: only include_graphics() here
            i += 1
            while not lines[i].strip().startswith("```"):
                f = FIG_RE.search(lines[i])
                if f:
                    figs.append((i + 1, f.group(1)))
                i += 1
            i += 1; continue
        m = OPEN_RE.match(L)
        if m:
            indent, lang = m.group(1), m.group(2).lower()
            start, body = i + 1, []
            i += 1
            while not CLOSE_RE.match(lines[i]):
                body.append(lines[i][len(indent):] if lines[i].startswith(indent) else lines[i])
                i += 1
            n += 1
            blocks.append(Block(ep, n, lang or "plain", start, i + 1, "/".join(divs) or "-",
                                "\n".join(body) + "\n"))
        i += 1
    return blocks, figs, marks


EPS = {}
for rmd in sorted(EPISODES.glob("*.Rmd")):
    key = rmd.stem.split("-")[0]  # 01 02 03 04a 04b 05 05b 06
    EPS[key] = (rmd,) + parse(rmd, key)
SETUP_BLOCKS, _, SETUP_MARKS = parse(SETUP, "setup")


def find(ep, anchor, lang=None):
    """The single code block of episode ep whose text contains anchor."""
    blocks = SETUP_BLOCKS if ep == "setup" else EPS[ep][1]
    hits = [b for b in blocks if anchor in b.text and (lang is None or b.lang == lang)]
    if len(hits) != 1:
        sys.exit("build_kit: anchor %r in %s matched %d blocks (need exactly 1)" % (anchor, ep, len(hits)))
    return hits[0]


# ---------------------------------------------------------------- the plan
# Steps run in this order; deps are afterok unless noted. kind:
#   kitjob   batch job written by the kit; runs learner blocks in a learner-like shell
#   learner  the learner's own SLURM script, submitted with sbatch from scripts/
#   ood      R episode in the OOD image (run_blocks.R), one fresh R session
# Resources for kitjob steps that replace an interactive session come from the
# episode's own sinteractive line. est = expected minutes on Negishi (local measurement
# x1.5, rounded up; used only for the preflight core-hour estimate).
STEPS = []


def step(name, kind, track, deps, cpus, time, est, desc, **kw):
    d = dict(name=name, kind=kind, track=track, deps=deps, cpus=cpus, time=time, est=est,
             desc=desc, mem=kw.pop("mem", ""), ntasks=kw.pop("ntasks", 1), array=kw.pop("array", 1))
    d.update(kw)
    STEPS.append(d)
    return d


def sint(block):
    """Parse the sinteractive line of a block -> (resources dict, block text without that line)."""
    out, res = [], None
    for line in block.text.splitlines():
        if line.strip().startswith("sinteractive "):
            a = shlex.split(line)
            res = {}
            k = 1
            while k < len(a):
                tok = a[k]
                if tok in ("-A", "-q", "-p", "-N", "-n"):
                    res[tok] = a[k + 1]; k += 2; continue
                if tok.startswith("--time="):
                    res["--time"] = tok.split("=", 1)[1]
                k += 1
            continue
        out.append(line)
    if res is None:
        sys.exit("build_kit: no sinteractive line in %s" % block.ref)
    return res, "\n".join(out).strip("\n") + "\n"


def learner_script(ep, anchor, fname):
    b = find(ep, anchor, "bash")
    if not b.text.startswith("#!/bin/bash"):
        sys.exit("build_kit: %s is not a SLURM script" % b.ref)
    return b, fname


LEARNER_SCRIPTS = [
    learner_script("04a", "--runMode genomeGenerate", "index_genome.sh"),
    learner_script("04a", "--job-name=read_mapping", "map_reads.sh"),
    learner_script("04a", "--job-name=featurecounts", "count_features.sh"),
    learner_script("04b", "--job-name=kallisto_index", "index_kallisto.sh"),
    learner_script("04b", "--job-name=kallisto_quant", "quant_kallisto.sh"),
]


def header_of(script_text):
    h = {}
    for line in script_text.splitlines():
        m = re.match(r'#SBATCH\s+--([\w-]+)(?:=(\S+))?', line)
        if m:
            h[m.group(1)] = m.group(2) or ""
    return h


def check_learner_headers():
    for b, fname in LEARNER_SCRIPTS:
        h = header_of(b.text)
        for k, v in (("account", "rcac-rnaseq"), ("qos", "standby"), ("partition", "cpu")):
            if h.get(k) != v:
                sys.exit("build_kit: %s (%s) has --%s=%s, expected %s" % (fname, b.ref, k, h.get(k), v))
        for k in ("time", "cpus-per-task"):
            if k not in h:
                sys.exit("build_kit: %s (%s) has no --%s" % (fname, b.ref, k))
        # Negishi allocates memory by cores; scripts must not request memory directly
        if "mem" in h or "mem-per-cpu" in h:
            sys.exit("build_kit: %s (%s) sets --mem; request more cores instead" % (fname, b.ref))


check_learner_headers()

# --- session pieces: ("block", Block, text) or ("extra", label, shell text)
S = {}

# 02: directory creation in the learner directory; reference download in a fresh
# directory; subsample spoiler on tiny inputs
b02_mkdir = find("02", "mkdir -p rnaseq-workshop/{data,scripts,results}", "bash")
b02_dl = find("02", 'GTFlink="https://ftp.ebi.ac.uk', "bash")
b02_clean = find("02", "transcripts-clean.fa", "bash")
b02_sub = find("02", "seqtk sample -s 42", "bash")
S["02"] = [
    ("block", b02_mkdir, b02_mkdir.text),
    ("extra", "fresh-dir reference download (02 download block; SCRATCH points to an empty test root)",
     'export SCRATCH="$RUN/ep02_fresh"; mkdir -p "$SCRATCH/rnaseq-workshop/data"'),
    ("block", b02_dl, b02_dl.text),
    ("block", b02_clean, b02_clean.text),
    ("extra", "compare the fresh downloads with the staged copies",
     'kit_compare_refs "$RUN/ep02_fresh/rnaseq-workshop/data" "$STAGED/data" > "$RUN/records/listings/02-ref-compare.tsv"'),
    ("extra", "subsample spoiler on 10,000-read inputs (20000000 -> 5000)",
     'export SCRATCH="$KIT_LEARNER_SCRATCH"; kit_subsample_inputs "$RUN/ep02_subsample"; cd "$RUN/ep02_subsample"'),
    ("block", b02_sub, b02_sub.text.replace(" 20000000 ", " 5000 ")),
    ("extra", "check the subsample result",
     'kit_check_subsample "$RUN/ep02_subsample" > "$RUN/records/listings/02-subsample-check.txt"'),
]
if b02_sub.text.count(" 20000000 ") != 2:
    sys.exit("build_kit: 02 subsample block no longer has two ' 20000000 ' arguments")

# 03: FastQC session, MultiQC
b03_fq = find("03", "fastqc data/*.fastq.gz", "bash")
b03_mq = find("03", "multiqc results/qc_fastq/", "bash")
b03_fastp = find("03", "fastp \\", "bash")
res03, txt03 = sint(b03_fq)
S["03"] = [
    ("block", b03_fq, txt03),
    ("extra", "listing for 03-out-fastqc-ls",
     'ls -1 "$W/results/qc_fastq" > "$RUN/records/listings/03-fastqc-ls.txt"'),
    ("block", b03_mq, b03_mq.text),
    ("extra", "MultiQC plot export for the 02_qc figures (kit output dir, not the learner dir)",
     'kit_multiqc_export "$W/results/qc_fastq" "$RUN/records/multiqc_export_qc"'),
    ("extra", "fastp spoiler on WT_Bcell_IR_rep1 (placeholder names substituted)",
     'kit_fastp_test'),
]
FASTP_TEXT = (b03_fastp.text.replace("SRRXXXXXXX_1.fastq.gz", "\"$W/data/WT_Bcell_IR_rep1_R1.fastq.gz\"")
              .replace("SRRXXXXXXX_2.fastq.gz", "\"$W/data/WT_Bcell_IR_rep1_R2.fastq.gz\"")
              .replace("SRRXXXXXXX_1.trimmed.fastq.gz", "\"$RUN/tmp/fastp/IR_rep1_R1.trimmed.fastq.gz\"")
              .replace("SRRXXXXXXX_2.trimmed.fastq.gz", "\"$RUN/tmp/fastp/IR_rep1_R2.trimmed.fastq.gz\"")
              .replace("fastp_report.html", "\"$RUN/records/fastp/fastp_report.html\"")
              .replace("fastp_report.json", "\"$RUN/records/fastp/fastp_report.json\""))
if "SRRXXXXXXX" in FASTP_TEXT:
    sys.exit("build_kit: 03 fastp block has placeholder names the kit does not substitute")

# 04a
b04a_salmon = find("04a", "salmon index --transcripts", "bash")
res04a_salmon, txt04a_salmon = sint(b04a_salmon)
b04a_samples = find("04a", "ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//'", "bash")
b04a_sub_index = find("04a", "sbatch index_genome.sh", "bash")
b04a_sub_map = find("04a", "sbatch map_reads.sh", "bash")
b04a_mqmap = find("04a", "multiqc results/mapping", "bash")
b04a_mkcounts = find("04a", "mkdir -p $SCRATCH/rnaseq-workshop/results/counts", "bash")
b04a_sub_count = find("04a", "sbatch count_features.sh", "bash")
b04a_clean = find("04a", "gene_counts_clean.txt", "bash")
b04a_mqcount = find("04a", "multiqc results/counts", "bash")
b04a_star_ex = find("04a", "SRR1234567", "bash")
S["04a-salmon"] = [("block", b04a_salmon, txt04a_salmon)]
S["04a-mapqc"] = [
    ("block", b04a_mqmap, b04a_mqmap.text),
    ("extra", "listing for 04a-out-star-index-ls",
     'ls -1 "$W/data/star_index" > "$RUN/records/listings/04a-star-index-ls.txt"'),
]
S["04a-post"] = [
    ("extra", "USER points at the test root so the literal /scratch/negishi/$USER path in the next block resolves to it",
     'export USER="$KIT_USER_OVERRIDE"'),
    ("block", b04a_clean, b04a_clean.text),
    ("extra", "restore USER", 'export USER="$KIT_REAL_USER"'),
    ("block", b04a_mqcount, b04a_mqcount.text),
]
# submit-time (login node) blocks, run by submit_all.sh just before the sbatch
SUBMIT_TIME = {
    "04a-map": [b04a_samples],
    "04a-count": [b04a_mkcounts],
}
SUBMIT_BLOCK = {"04a-index": b04a_sub_index, "04a-map": b04a_sub_map, "04a-count": b04a_sub_count}
# --prebuilt-index: the 04a callout's ln -s instead of building the STAR index
b04a_prebuilt = find("04a", "ln -s /scratch/negishi/aseethar/rnaseq-workshop_results/data/star_index", "bash")

# 04b
b04b_sub_index = find("04b", "sbatch index_kallisto.sh", "bash")
b04b_mkdir = find("04b", "mkdir -p $SCRATCH/rnaseq-workshop/results/kallisto_quant", "bash")
b04b_example = find("04b", "-o results/kallisto_quant/WT_Bcell_mock_rep1", "bash")  # not run: the array runs it
b04b_samples = find("04b", "ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//'", "bash")
b04b_sub_quant = find("04b", "sbatch quant_kallisto.sh", "bash")
b04b_tx2gene = find("04b", "> tx2gene.tsv", "bash")
b04b_rstart = find("04b", "module load r-rnaseq", "bash")
res04b_r, txt04b_r = sint(b04b_rstart)
b04b_rblocks = [b for b in EPS["04b"][1] if b.lang == "r"]
b04b_mq = find("04b", "multiqc results/kallisto_quant", "bash")
if not txt04b_r.rstrip().endswith("\nR"):
    sys.exit("build_kit: 04b R-session block no longer ends with a bare 'R' line")
S["04b-tximport"] = [
    ("block", b04b_tx2gene, b04b_tx2gene.text),
    # the interactive 'R' becomes R reading the episode's R blocks on stdin, as if typed
    ("block", b04b_rstart, txt04b_r.rstrip()[:-1] + 'R --no-save --no-restore < "$KIT_GEN/rsessions/04b-tximport.R"\n'),
]
S["04b-qc"] = [
    ("block", b04b_mq, b04b_mq.text),
    ("extra", "listings for 04b-out-kallisto-dir and 04b-out-abundance",
     'ls -1 "$W/results/kallisto_quant/WT_Bcell_mock_rep1" > "$RUN/records/listings/04b-kallisto-dir.txt"; '
     'head -6 "$W/results/kallisto_quant/WT_Bcell_mock_rep1/abundance.tsv" > "$RUN/records/listings/04b-abundance-head.txt"'),
]
SUBMIT_TIME["04b-quant"] = [b04b_mkdir, b04b_samples]
SUBMIT_BLOCK["04b-quant"] = b04b_sub_quant
SUBMIT_BLOCK["04b-index"] = b04b_sub_index

# R episodes: which blocks to skip in the main run, and why
R_SKIP = {
    "05": {find("05", "useMart(").nn: "biomaRt spoiler; run separately in step 05-biomart",
           find("05", "design = ~ batch + condition").nn: "fragment inside a callout, not runnable"},
    "05b": {}, "06": {},
}
# 05-biomart runs only the setup blocks plus the spoiler
B05_BIOMART = [find("05", "library(RColorBrewer)", "r").nn, find("05", 'read.delim(countsFile', "r").nn,
               find("05", '"data/mart.tsv"', "r").nn, find("05", "useMart(", "r").nn]
B05_BIOMART = sorted(set(B05_BIOMART))
# 06 transcript track: swap the two read_tsv lines, exactly as the episode tells 05b learners
SWAP_FROM = ('res <- read_tsv("results/deseq2/DESeq2_results_joined.tsv", show_col_types = FALSE)\n'
             '# Kallisto track (Episode 05b): use this line instead of the one above\n'
             '# res <- read_tsv("results/deseq2_kallisto/DESeq2_kallisto_results.tsv", show_col_types = FALSE)\n')
SWAP_TO = ('# res <- read_tsv("results/deseq2/DESeq2_results_joined.tsv", show_col_types = FALSE)\n'
           '# Kallisto track (Episode 05b): use this line instead of the one above\n'
           'res <- read_tsv("results/deseq2_kallisto/DESeq2_kallisto_results.tsv", show_col_types = FALSE)\n')
b06_swap = find("06", "# Kallisto track (Episode 05b): use this line instead", "r")
if SWAP_FROM not in b06_swap.text:
    sys.exit("build_kit: 06 track-switch lines changed; update SWAP_FROM")

# setup.md blocks for 01_learner_setup.sh
bS_echo = find("setup", "echo $RCAC_SCRATCH", "bash")
bS_rsync = find("setup", "rsync -avP /depot/DEPOT_PATH/rnaseq-workshop ${RCAC_SCRATCH}/", "bash")
bS_verify = find("setup", "ls ${RCAC_SCRATCH}/rnaseq-workshop/data/*.fastq.gz | wc -l", "bash")

# -------- step table
A, Q, P = "rcac-rnaseq", "standby", "cpu"
step("02", "kitjob", "common", [], 4, "1:00:00", 20, "02 directories, fresh reference download, subsample spoiler")
step("03", "kitjob", "common", [], int(res03["-n"]), res03["--time"], 45, "03 FastQC + MultiQC (+ fastp spoiler, MultiQC export)",
     ntasks=int(res03["-n"]))
step("04a-salmon", "kitjob", "genome", [], int(res04a_salmon["-n"]), res04a_salmon["--time"], 45,
     "04a Salmon index + strandedness check", ntasks=int(res04a_salmon["-n"]))
step("04a-index", "learner", "genome", [], 20, "", 65, "04a STAR index (index_genome.sh)", script="index_genome.sh")
step("04a-map", "learner", "genome", ["04a-index"], 20, "", 18, "04a STAR mapping array (map_reads.sh)",
     script="map_reads.sh", array=8)
step("04a-mapqc", "kitjob", "genome", ["04a-map"], 1, "0:30:00", 5, "04a MultiQC on mapping")
step("04a-count", "learner", "genome", ["04a-map"], 16, "", 6, "04a featureCounts (count_features.sh)",
     script="count_features.sh")
step("04a-post", "kitjob", "genome", ["04a-count"], 1, "0:30:00", 5, "04a clean count matrix + MultiQC on counts")
step("04b-index", "learner", "transcript", [], 48, "", 35, "04b kallisto index (index_kallisto.sh)",
     script="index_kallisto.sh")
step("04b-quant", "learner", "transcript", ["04b-index"], 16, "", 12, "04b kallisto quant array (quant_kallisto.sh)",
     script="quant_kallisto.sh", array=8)
step("04b-tximport", "kitjob", "transcript", ["04b-quant"], int(res04b_r["-n"]), res04b_r["--time"], 6,
     "04b tx2gene + tximport (r-rnaseq module)", ntasks=int(res04b_r["-n"]))
step("04b-qc", "kitjob", "transcript", ["04b-quant"], 1, "0:30:00", 5, "04b MultiQC on kallisto")
step("05", "ood", "genome", ["04a-post"], 4, "4:00:00", 10, "05 DESeq2 (featureCounts), OOD image", episode="05")
step("05-biomart", "ood", "genome", ["04a-post"], 4, "1:00:00", 10, "05 biomaRt spoiler alone, OOD image",
     episode="05")
step("05b", "ood", "transcript", ["04b-tximport"], 4, "4:00:00", 6, "05b DESeq2 (tximport), OOD image", episode="05b")
step("06-genome", "ood", "genome", ["05"], 4, "4:00:00", 15, "06 enrichment from the 05 table", episode="06")
step("06-transcript", "ood", "transcript", ["05b"], 4, "4:00:00", 15, "06 enrichment from the 05b table",
     episode="06")
step("metrics", "ood", "common", ["05", "05b", "05-biomart"], 4, "1:00:00", 8,
     "kit metrics (PCA, anchors, track comparison); not learner code", dep_type="afterany")

# learner-step headers come from the learner scripts themselves
for d in STEPS:
    if d["kind"] == "learner":
        b = [b for b, f in LEARNER_SCRIPTS if f == d["script"]][0]
        h = header_of(b.text)
        d["time"], d["cpus"] = h["time"], int(h["cpus-per-task"])

# ---------------------------------------------------------------- placeholders
# Sources for markers that are not "the console output of the code block just above".
FILE_SOURCES = {
    "02-out-data-tree": "records/listings/02-data-tree.txt",
    "03-out-fastqc-ls": "records/listings/03-fastqc-ls.txt",
    "04a-out-libformat": "learner:results/strand_check/lib_format_counts.json",
    "04a-out-star-index-ls": "records/listings/04a-star-index-ls.txt",
    "04b-out-kallisto-dir": "records/listings/04b-kallisto-dir.txt",
    "04b-out-abundance": "records/listings/04b-abundance-head.txt",
}
# Which step runs each episode's blocks (for console lookup)
BLOCK_STEP = {}
for name, pieces in S.items():
    for kind, b, _ in pieces:
        if kind == "block":
            BLOCK_STEP[b.ref] = name
for b in b04b_rblocks:
    BLOCK_STEP[b.ref] = "04b-tximport"
for ep, st in (("05", "05"), ("05b", "05b"), ("06", "06-genome")):
    for b in EPS[ep][1]:
        if b.lang == "r" and b.nn not in R_SKIP[ep]:
            BLOCK_STEP[b.ref] = st
BLOCK_STEP[find("05", "useMart(").ref] = "05-biomart"


def placeholder_rows():
    rows, seen = [], set()
    for key in sorted(EPS):
        rmd, blocks, figs, marks = EPS[key]
        for line, mid, own in marks:
            if mid in seen:
                sys.exit("build_kit: duplicate placeholder id %s" % mid)
            seen.add(mid)
            relf = rmd.relative_to(ROOT).as_posix()
            if mid in FILE_SOURCES:
                rows.append((mid, relf, line, "file", FILE_SOURCES[mid]))
            elif own and "-out-" in mid:
                prev = [b for b in blocks if b.end < line and b.lang in ("bash", "r")]
                if not prev:
                    sys.exit("build_kit: %s has no code block above it" % mid)
                b = prev[-1]
                if b.ref not in BLOCK_STEP:
                    sys.exit("build_kit: %s comes from %s, which the kit does not run" % (mid, b.ref))
                src = "console:%s:%s" % (BLOCK_STEP[b.ref], b.ref)
                if key == "06":
                    src += ";console:06-transcript:%s" % b.ref
                rows.append((mid, relf, line, "block", src))
            else:
                if mid not in summarize.METRICS:
                    sys.exit("build_kit: no metric handler for %s in py/summarize.py" % mid)
                rows.append((mid, relf, line, "metric", mid))
    for line, mid, own in SETUP_MARKS:
        if mid not in summarize.METRICS:
            sys.exit("build_kit: no metric handler for %s in py/summarize.py" % mid)
        rows.append((mid, SETUP.relative_to(ROOT).as_posix(), line, "metric", mid))
    # After the placeholders are filled, the metrics are still computed by summarize.py and
    # reported under "Measured values", so a rerun can be compared with the published text.
    return rows


# ---------------------------------------------------------------- writers
OUT = {}  # relative path -> text


def put(rel, text):
    OUT[rel] = text if text.endswith("\n") else text + "\n"


def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def render_kitjob(d):
    name = d["name"]
    lines = [
        "#!/bin/bash",
        "# GENERATED by build_kit.py from the episodes; do not edit. Step: %s" % name,
        "# %s" % d["desc"],
        "#SBATCH --job-name=kit-%s" % name,
        "#SBATCH --account=%s" % A,
        "#SBATCH --qos=%s" % Q,
        "#SBATCH --partition=%s" % P,
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=%d" % d["ntasks"],
        "#SBATCH --time=%s" % d["time"],
        "set -euo pipefail",
        'source "$KIT_DIR/lib.sh"',
        'kit_job_begin "%s"' % name,
        'kit_run_session "%s" "$KIT_GEN/sessions/%s.sh"' % (name, name),
        'kit_job_end "%s"' % name,
    ]
    return "\n".join(lines)


def render_oodjob(d):
    name = d["name"]
    if name == "metrics":
        run = 'kit_ood_run metrics "$KIT_GEN/rplans/05.tsv" "$KIT_USER_OVERRIDE" "$KIT_DIR/R/metrics.R"'
    elif name == "06-transcript":
        run = ('kit_transcript_dir\n'
               'kit_ood_run "06-transcript" "$KIT_GEN/rplans/06-transcript.tsv" "$KIT_USER_OVERRIDE_T"')
    else:
        run = 'kit_ood_run "%s" "$KIT_GEN/rplans/%s.tsv" "$KIT_USER_OVERRIDE"' % (name, name)
    lines = [
        "#!/bin/bash",
        "# GENERATED by build_kit.py; do not edit. Step: %s" % name,
        "# %s" % d["desc"],
        "# Resources match the Open OnDemand form in learners/setup.md (4 cores, standby).",
        "#SBATCH --job-name=kit-%s" % name,
        "#SBATCH --account=%s" % A,
        "#SBATCH --qos=%s" % Q,
        "#SBATCH --partition=%s" % P,
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=%d" % d["cpus"],
        "#SBATCH --time=%s" % d["time"],
        "set -euo pipefail",
        'source "$KIT_DIR/lib.sh"',
        'kit_job_begin "%s"' % name,
        run,
    ]
    return "\n".join(lines)


def render_session(name, pieces):
    """Learner-like shell: no -e/-u/pipefail; errors are trapped and logged, not fatal."""
    out = ["# GENERATED by build_kit.py; learner session for step %s." % name,
           "# Run by kit_run_session in a child bash without set -euo pipefail, like a learner's shell.",
           "# Lmod's module function, as in a learner's login shell (sbatch normally exports it)",
           'type module >/dev/null 2>&1 || source "${LMOD_PKG:-/opt/lmod/lmod}/init/bash" 2>/dev/null || source /etc/profile.d/lmod.sh',
           "set -E",
           "trap 'echo \"### KIT-ERR rc=$? line=$LINENO cmd=$BASH_COMMAND\" >&2' ERR",
           'cd "$W" 2>/dev/null || cd "$KIT_LEARNER_SCRATCH"']
    for kind, b, text in pieces:
        if kind == "block":
            out.append('echo "### KIT-BLOCK %s begin $(date +%%s) %s:%d"' % (b.ref, EPS[b.ep][0].name, b.start))
            out.append(text.rstrip("\n"))
            out.append('echo "### KIT-BLOCK %s end $(date +%%s)"' % b.ref)
        else:
            out.append("# KIT EXTRA: %s" % b)
            out.append("echo '### KIT-EXTRA begin' \"$(date +%%s)\" '%s'" % b.replace("'", ""))
            out.append(text)
            out.append('echo "### KIT-EXTRA end $(date +%s)"')
    return "\n".join(out)


def render_rsession_04b():
    out = []
    for b in b04b_rblocks:
        out.append('cat("### KIT-BLOCK %s begin", format(Sys.time(), "%%s"), "\\n")' % b.ref)
        out.append(b.text.rstrip("\n"))
        out.append('cat("### KIT-BLOCK %s end", format(Sys.time(), "%%s"), "\\n")' % b.ref)
    return "\n".join(out)


def rplan(ep, skip, only=None):
    rows = ["episode\tblock\tfile\taction\tline\tcontext"]
    for b in EPS[ep][1]:
        if b.lang != "r":
            continue
        if only is not None:
            act = "run" if b.nn in only else "skip:not part of this step"
        else:
            act = ("skip:" + skip[b.nn]) if b.nn in skip else "run"
        rows.append("%s\t%02d\tblocks/%s/%02d.R\t%s\t%d\t%s" % (ep, b.nn, ep, b.nn, act, b.start, b.context))
    return "\n".join(rows)


def build():
    # episode checksums
    sums = []
    for key in sorted(EPS):
        rmd = EPS[key][0]
        sums.append("%s  %s" % (sha256(rmd), rmd.relative_to(ROOT).as_posix()))
    sums.append("%s  %s" % (sha256(SETUP), SETUP.relative_to(ROOT).as_posix()))
    put("episodes.sha256", "\n".join(sums))

    # every block, verbatim
    idx = ["episode\tblock\tlang\tstart\tend\tcontext\tfile\tstep"]
    for key in sorted(EPS):
        for b in EPS[key][1]:
            ext = {"bash": "sh", "r": "R"}.get(b.lang, "txt")
            rel = "blocks/%s/%02d.%s" % (key, b.nn, ext)
            put(rel, b.text)
            idx.append("%s\t%02d\t%s\t%d\t%d\t%s\t%s\t%s" % (key, b.nn, b.lang, b.start, b.end, b.context, rel,
                                                             BLOCK_STEP.get(b.ref, "-")))
    put("index.tsv", "\n".join(idx))
    # 06 transcript-track copy of the read block
    for b in EPS["06"][1]:
        if b.lang == "r":
            t = b.text.replace(SWAP_FROM, SWAP_TO) if b is b06_swap else b.text
            put("blocks/06@transcript/%02d.R" % b.nn, t)

    # learner scripts (saved by the learner as the episode says) and staged copies
    for b, fname in LEARNER_SCRIPTS:
        put("learner-scripts/" + fname, b.text)

    # sessions and job scripts
    for d in STEPS:
        if d["kind"] == "kitjob":
            put("sessions/%s.sh" % d["name"], render_session(d["name"], S[d["name"]]))
            put("jobs/%s.sh" % d["name"], render_kitjob(d))
        elif d["kind"] == "ood":
            put("jobs/%s.sh" % d["name"], render_oodjob(d))
    put("rsessions/04b-tximport.R", render_rsession_04b())
    put("sessions/03-fastp.sh", "# GENERATED: 03 fastp spoiler with placeholder names substituted\n"
        'mkdir -p "$RUN/tmp/fastp" "$RUN/records/fastp"\n' + FASTP_TEXT)
    for d in STEPS:
        if d["name"] in SUBMIT_TIME:
            put("submit-time/%s.sh" % d["name"], render_session(d["name"] + " (submit time, login node)",
                                                               [("block", b, b.text) for b in SUBMIT_TIME[d["name"]]]))
    for name, b in sorted(SUBMIT_BLOCK.items()):
        put("submit-blocks/%s.sh" % name, b.text)
    put("submit-time/04a-index-prebuilt.sh", render_session("04a-index (prebuilt, login node)",
                                                            [("block", b04a_prebuilt, b04a_prebuilt.text)]))

    # R plans
    put("rplans/05.tsv", rplan("05", R_SKIP["05"]))
    put("rplans/05-biomart.tsv", rplan("05", {}, only=B05_BIOMART))
    put("rplans/05b.tsv", rplan("05b", R_SKIP["05b"]))
    put("rplans/06-genome.tsv", rplan("06", R_SKIP["06"]))
    put("rplans/06-transcript.tsv", rplan("06", R_SKIP["06"]).replace("blocks/06/", "blocks/06@transcript/"))

    # packages the OOD episodes need: every library() call, plus packages used without one
    pk = set()
    for ep in ("05", "05b", "06"):
        for b in EPS[ep][1]:
            if b.lang == "r" and b.nn not in R_SKIP[ep]:
                pk.update(re.findall(r'^\s*library\(([\w.]+)\)', b.text, re.M))
    pk.update(["apeglm", "hexbin", "matrixStats", "jsonlite", "biomaRt"])  # lfcShrink(apeglm), meanSdPlot(), kit metrics, 05 spoiler
    put("ood_packages.txt", "\n".join(sorted(pk)))
    pk4b = set(re.findall(r'^\s*library\(([\w.]+)\)', "".join(b.text for b in b04b_rblocks), re.M)) | {"rhdf5"}
    put("r-rnaseq_packages.txt", "\n".join(sorted(pk4b)))

    # setup.md blocks for 01_learner_setup.sh
    put("setup/01-echo.sh", bS_echo.text)
    put("setup/02-rsync.sh", bS_rsync.text)
    put("setup/03-verify.sh", bS_verify.text)

    # step table
    st = ["step\tkind\ttrack\tdeps\tdep_type\tcpus\tarray\tmem\ttime\test_min\tscript\tdescription"]
    for d in STEPS:
        st.append("\t".join([d["name"], d["kind"], d["track"], ",".join(d["deps"]) or "-",
                             d.get("dep_type", "afterok"), str(d["cpus"]), str(d["array"]), d["mem"] or "-",
                             d["time"], str(d["est"]), d.get("script", "-"), d["desc"]]))
    put("steps.tsv", "\n".join(st))

    # figure map: R blocks -> committed figure names (next include_graphics before the next code block)
    fm = ["step\tepisode\tblock\tk\tfigure"]
    for ep, steps in (("05", ["05"]), ("05b", ["05b"]), ("06", ["06-genome", "06-transcript"])):
        rmd, blocks, figs, marks = EPS[ep]
        code = [b for b in blocks if b.lang in ("r", "bash")]
        for i, b in enumerate(code):
            if b.lang != "r":
                continue
            nxt = code[i + 1].start if i + 1 < len(code) else 10 ** 9
            targets = [f for (ln, f) in figs if b.end < ln < nxt]
            for k, f in enumerate(targets, 1):
                for s in steps:
                    fm.append("%s\t%s\t%02d\t%d\t%s" % (s, ep, b.nn, k, f))
    put("figmap.tsv", "\n".join(fm))

    # placeholders
    ph = ["id\tfile\tline\tkind\tsource"]
    for r in placeholder_rows():
        ph.append("%s\t%s\t%d\t%s\t%s" % r)
    put("placeholders.tsv", "\n".join(ph))


def write_all(check=False):
    staged = {}
    for b, fname in LEARNER_SCRIPTS:
        staged[fname] = b.text
    samples_csv = find("05", "sample,condition", "plain").text
    staged["samples.csv"] = samples_csv
    # samples.txt as `ls *_R1.fastq.gz | sed ...` produces it (IR sorts before mock)
    names = [l.split(",")[0] for l in samples_csv.strip().splitlines()[1:]]
    staged["samples.txt"] = "\n".join(sorted(names, key=lambda s: s.lower())) + "\n"
    if check:
        diffs = []
        for rel, text in sorted(OUT.items()):
            p = GEN / rel
            if not p.exists() or p.read_text() != text:
                diffs.append("generated/" + rel)
        for f, text in sorted(staged.items()):
            p = STAGED_OUT / f
            if not p.exists() or p.read_text() != text:
                diffs.append("staged-scripts/" + f)
        have = {p.relative_to(GEN).as_posix() for p in GEN.rglob("*") if p.is_file()} if GEN.exists() else set()
        diffs += ["generated/%s (stale)" % r for r in sorted(have - set(OUT))]
        if diffs:
            print("build_kit --check: differs from a fresh build:\n  " + "\n  ".join(diffs))
            sys.exit(1)
        print("build_kit --check: generated/ and staged-scripts/ match a fresh build (%d files)" % (len(OUT) + len(staged)))
        return
    for d in (GEN, STAGED_OUT):
        if d.exists():
            shutil.rmtree(d)
    for rel, text in sorted(OUT.items()):
        p = GEN / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
        if rel.endswith(".sh"):
            p.chmod(0o755)
    STAGED_OUT.mkdir()
    for f, text in sorted(staged.items()):
        p = STAGED_OUT / f
        p.write_text(text)
        if f.endswith(".sh"):
            p.chmod(0o755)
    print("build_kit: wrote %d files to %s and %d to %s" % (len(OUT), GEN.relative_to(ROOT),
                                                         len(staged), STAGED_OUT.relative_to(ROOT)))


if __name__ == "__main__":
    build()
    write_all(check="--check" in sys.argv[1:])
