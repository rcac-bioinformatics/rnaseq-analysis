#!/usr/bin/env python3
"""summarize.py RUN_DIR KIT_GEN OUT_DIR

Reads a kit run directory and writes, into OUT_DIR:
  SUMMARY.md               per-step and per-episode status, anchors, runtimes, placeholder values
  placeholder_values.tsv   one row per <!-- NEGISHI:id --> marker with the value or source file
  outputs/<id>.txt         console output (or file content) for every output-block marker
Python 3.6 compatible, standard library only (Negishi login nodes).
METRICS maps every non-output placeholder id to a function; build_kit.py checks coverage.
"""
import csv
import glob
import io
import json
import os
import re
import statistics
import sys
import zipfile

SAMPLES = ["WT_Bcell_IR_rep1", "WT_Bcell_IR_rep2", "WT_Bcell_IR_rep3", "WT_Bcell_IR_rep4",
           "WT_Bcell_mock_rep1", "WT_Bcell_mock_rep2", "WT_Bcell_mock_rep3", "WT_Bcell_mock_rep4"]


class Ctx(object):
    def __init__(self, run, gen):
        self.run, self.gen = run, gen
        env = {}
        p = os.path.join(run, "kit.env")
        if os.path.exists(p):
            for line in open(p):
                if "=" in line:
                    k, v = line.rstrip("\n").split("=", 1)
                    env[k] = v.strip("'\"")
        self.W = env.get("W", os.path.join(run, "scratch", "rnaseq-workshop"))
        self.W_T = env.get("W_T", os.path.join(run, "scratch_t", "rnaseq-workshop"))
        self.rec = os.path.join(run, "records")
        kenv = os.path.join(run, "kit.env")
        self.start = os.path.getmtime(kenv) if os.path.exists(kenv) else 0
        self.used = []  # output files read by the current metric (freshness check)
        self.jobs = read_tsv(os.path.join(run, "jobs.tsv"))
        self.sacct = parse_sacct(os.path.join(self.rec, "sacct.txt"))
        self.steps = read_tsv(os.path.join(gen, "steps.tsv"))

    def path(self, *a):
        return os.path.join(self.W, *a)

    def use(self, p):
        """record an output file a metric reads; returns p"""
        self.used.append(p)
        return p

    def stale(self, paths):
        """output files older than this run: copied in with the data, not produced by the run"""
        return [p for p in paths if os.path.exists(p) and os.path.getmtime(p) < self.start]

    def jobids(self, step):
        return [r["jobid"] for r in self.jobs if r["step"] == step and r["jobid"] not in ("", "-")]

    def job_stats(self, step):
        """list of (task_id, elapsed_s, maxrss_gb, state) for a step's job(s), last submission only"""
        ids = self.jobids(step)
        if not ids:
            return []
        jid = ids[-1]
        return [(k, v["elapsed"], v["maxrss"], v["state"]) for k, v in sorted(self.sacct.items())
                if k == jid or k.startswith(jid + "_")]


def read_tsv(p):
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def hms(s):
    """SLURM elapsed [D-]HH:MM:SS -> seconds"""
    if not s or s in ("Unknown",):
        return None
    d = 0
    if "-" in s:
        d, s = s.split("-", 1)
        d = int(d)
    parts = [float(x) for x in s.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    return d * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def rss_gb(s):
    if not s:
        return None
    m = re.match(r'([\d.]+)([KMGT]?)', s)
    if not m:
        return None
    f = {"": 1.0 / 1024 ** 3, "K": 1.0 / 1024 ** 2, "M": 1.0 / 1024, "G": 1.0, "T": 1024.0}[m.group(2)]
    return float(m.group(1)) * f


def parse_sacct(p):
    """sacct -P -n --format=JobID,JobName,Elapsed,MaxRSS,ReqMem,State,Submit,Start,End,AllocCPUS,Timelimit,ExitCode"""
    out = {}
    if not os.path.exists(p):
        return out
    for line in open(p):
        f = line.rstrip("\n").split("|")
        if len(f) < 12:
            continue
        jid = f[0]
        root = jid.split(".")[0]
        e = out.setdefault(root, {"elapsed": None, "maxrss": None, "state": "", "submit": "", "start": "",
                                  "end": "", "cpus": "", "limit": "", "reqmem": "", "name": ""})
        if "." not in jid:
            e.update(elapsed=hms(f[2]), state=f[5], submit=f[6], start=f[7], end=f[8], cpus=f[9],
                     limit=f[10], reqmem=f[4], name=f[1])
        r = rss_gb(f[3])
        if r is not None and (e["maxrss"] is None or r > e["maxrss"]):
            e["maxrss"] = r
    return out


def mins(s):
    return "%.0f" % (s / 60.0) if s is not None else "?"


def rng(vals, fmt="%.1f"):
    vals = [v for v in vals if v is not None]
    if not vals:
        return "n/a"
    return (fmt + " to " + fmt) % (min(vals), max(vals))


def short(s):
    return s.replace("WT_Bcell_", "")


# ------------------------------------------------------------ parsers of tool outputs
def fastqc(ctx):
    res = {}
    for z in sorted(glob.glob(ctx.path("results", "qc_fastq", "*_fastqc.zip"))):
        ctx.use(z)
        name = os.path.basename(z)[:-len("_fastqc.zip")]
        with zipfile.ZipFile(z) as zf:
            base = [n for n in zf.namelist() if n.endswith("/summary.txt")][0].rsplit("/", 1)[0]
            summ = {}
            for line in io.TextIOWrapper(zf.open(base + "/summary.txt")):
                st, mod, _ = line.rstrip("\n").split("\t")
                summ[mod] = st
            dedup = None
            for line in io.TextIOWrapper(zf.open(base + "/fastqc_data.txt")):
                if line.startswith("#Total Deduplicated Percentage"):
                    dedup = float(line.split("\t")[1])
        res[name] = (summ, dedup)
    return res


def count_status(ctx, module, status):
    fq = fastqc(ctx)
    if not fq:
        return None
    return "%d of %d" % (sum(1 for s, _ in fq.values() if s.get(module) == status), len(fq))


def star_logs(ctx):
    res = {}
    for f in sorted(glob.glob(ctx.path("results", "mapping", "*Log.final.out"))):
        ctx.use(f)
        s = os.path.basename(f)[:-len("Log.final.out")]
        d = {}
        for line in open(f):
            if "|" in line:
                k, v = line.split("|", 1)
                d[k.strip()] = v.strip()
        res[s] = d
    return res


def pct(v):
    return float(v.rstrip("%"))


def fc_summary(ctx):
    p = ctx.use(ctx.path("results", "counts", "gene_counts.txt.summary"))
    if not os.path.exists(p):
        return {}
    rows = list(csv.reader(open(p), delimiter="\t"))
    names = [short(os.path.basename(x).replace("Aligned.sortedByCoord.out.bam", "")) for x in rows[0][1:]]
    res = {n: {} for n in names}
    for r in rows[1:]:
        for n, v in zip(names, r[1:]):
            res[n][r[0]] = int(v)
    return res


def kallisto(ctx):
    res = {}
    for s in SAMPLES:
        p = ctx.use(ctx.path("results", "kallisto_quant", s, "run_info.json"))
        if not os.path.exists(p):
            continue
        ri = json.load(open(p))
        fl = None
        lg = ctx.path("results", "kallisto_quant", s, s + ".log")
        if os.path.exists(lg):
            m = re.search(r'estimated average fragment length:\s*([\d.]+)', open(lg).read())
            fl = float(m.group(1)) if m else None
        res[s] = (ri, fl)
    return res


def metric_tsv(ctx, name):
    return read_tsv(os.path.join(ctx.rec, "metrics", name))




def timing(ctx, step):
    st = ctx.job_stats(step)
    if not st:
        return None
    el = [x[1] for x in st if x[1] is not None]
    rs = [x[2] for x in st if x[2] is not None]
    return el, rs


# ------------------------------------------------------------ metric handlers
def m_star_map_time(ctx):
    t = timing(ctx, "04a-map")
    if not t:
        return None
    el = [e for e in t[0]]
    tasks = [x for x in ctx.job_stats("04a-map") if "_" in x[0]]
    el = [x[1] for x in tasks if x[1] is not None] or el
    return "%s to %s min per sample (8 tasks); peak memory %.1f GB" % (mins(min(el)), mins(max(el)), max(t[1] or [0]))


def m_single_time(step, what):
    def f(ctx):
        t = timing(ctx, step)
        if not t or not t[0]:
            return None
        return "%s: %s min elapsed, peak memory %.1f GB" % (what, mins(max(t[0])), max(t[1] or [0]))
    return f


def m_isf_isr(ctx):
    p = ctx.use(ctx.path("results", "strand_check", "lib_format_counts.json"))
    if not os.path.exists(p):
        return None
    j = json.load(open(p))
    return "ISF %.1f million (%d), ISR %.1f million (%d), format %s" % (
        j["ISF"] / 1e6, j["ISF"], j["ISR"] / 1e6, j["ISR"], j["expected_format"])


def m_bias(ctx):
    p = ctx.use(ctx.path("results", "strand_check", "lib_format_counts.json"))
    return "%.3f" % json.load(open(p))["strand_mapping_bias"] if os.path.exists(p) else None


def m_unique(ctx):
    s = star_logs(ctx)
    v = {k: pct(d["Uniquely mapped reads %"]) for k, d in s.items()}
    return ("%s percent; " % rng(list(v.values()))) + ", ".join(
        "%s %.1f" % (short(k), x) for k, x in sorted(v.items())) if v else None


def m_unique_max(ctx):
    s = star_logs(ctx)
    if not s:
        return None
    k = max(s, key=lambda k: pct(s[k]["Uniquely mapped reads %"]))
    return "%s (%.1f percent)" % (short(k), pct(s[k]["Uniquely mapped reads %"]))


def m_unique_outliers(ctx):
    s = star_logs(ctx)
    if not s:
        return None
    v = {short(k): pct(d["Uniquely mapped reads %"]) for k, d in s.items()}
    mu, sd = statistics.mean(v.values()), statistics.pstdev(v.values())
    low = ["%s %.1f" % (k, x) for k, x in sorted(v.items()) if x < mu - sd]
    ir = statistics.mean([x for k, x in v.items() if k.startswith("IR")])
    mock = statistics.mean([x for k, x in v.items() if k.startswith("mock")])
    return "mean %.1f, sd %.1f; below mean-1sd: %s; IR mean %.1f vs mock mean %.1f" % (
        mu, sd, ", ".join(low) or "none", ir, mock)


def m_too_short(ctx):
    s = star_logs(ctx)
    v = [int(d.get("Number of reads unmapped: too short", "0")) / 1e6 for d in s.values()]
    return "%s million pairs" % rng(v) if v else None


def assigned_pct(fc):
    return {k: 100.0 * d["Assigned"] / sum(d.values()) for k, d in fc.items()}


def m_assigned_range(ctx):
    fc = fc_summary(ctx)
    if not fc:
        return None
    a = assigned_pct(fc)
    return "%s percent; " % rng(list(a.values())) + ", ".join("%s %.1f" % (k, x) for k, x in sorted(a.items()))


def m_assigned_min(ctx):
    fc = fc_summary(ctx)
    if not fc:
        return None
    a = assigned_pct(fc)
    k = min(a, key=a.get)
    return "%s (%.1f percent)" % (k, a[k])


def m_unassigned_largest(ctx):
    fc = fc_summary(ctx)
    if not fc:
        return None
    out = []
    for k, d in sorted(fc.items()):
        tot = float(sum(d.values()))
        un = {c: v for c, v in d.items() if c.startswith("Unassigned") and v > 0}
        top = sorted(un.items(), key=lambda x: -x[1])[:3]
        out.append("%s: %s" % (k, ", ".join("%s %.1f%%" % (c.replace("Unassigned_", ""), 100 * v / tot) for c, v in top)))
    return "; ".join(out)


def m_assigned_vs_depth(ctx):
    fc = fc_summary(ctx)
    if not fc:
        return None
    a = assigned_pct(fc)
    most = max(fc, key=lambda k: fc[k]["Assigned"])
    best = max(a, key=a.get)
    return "most assigned: %s (%d); highest percent: %s (%.1f); same sample: %s" % (
        most, fc[most]["Assigned"], best, a[best], "yes" if most == best else "no")


def m_pseudo_range(ctx):
    k = kallisto(ctx)
    if not k:
        return None
    v = {short(s): r["p_pseudoaligned"] for s, (r, _) in k.items()}
    return "%s percent; " % rng(list(v.values())) + ", ".join("%s %.1f" % (s, x) for s, x in sorted(v.items()))


def m_pseudo_consistency(ctx):
    k = kallisto(ctx)
    if not k:
        return None
    v = [r["p_pseudoaligned"] for r, _ in k.values()]
    n = [r["n_processed"] for r, _ in k.values()]
    return "sd %.1f points across samples; n_processed %s; kallisto %s, index %s, n_bootstraps %s" % (
        statistics.pstdev(v), rng(n, "%d"), list(k.values())[0][0].get("kallisto_version"),
        list(k.values())[0][0].get("index_version"), list(k.values())[0][0].get("n_bootstraps"))


def m_fraglen(ctx):
    k = kallisto(ctx)
    v = {short(s): fl for s, (_, fl) in k.items() if fl is not None}
    return ("%s bp; " % rng(list(v.values())) + ", ".join("%s %.1f" % kv for kv in sorted(v.items()))) if v else None


def m_kallisto_quant(ctx):
    tasks = [x for x in ctx.job_stats("04b-quant") if "_" in x[0] and x[1] is not None]
    if not tasks:
        return None
    return "%s to %s min per sample (16 cores, -b 0); peak memory %.1f GB" % (
        mins(min(x[1] for x in tasks)), mins(max(x[1] for x in tasks)), max(x[2] or 0 for x in tasks))


def anchors_text(ctx, gene=None):
    rows = metric_tsv(ctx, "anchors.tsv")
    if gene:
        rows = [r for r in rows if r["external_gene_name"] == gene]
    return "; ".join("%s %s LFC %.2f padj %s (%s)" % (r["track"], r["external_gene_name"], float(r["log2FoldChange"]),
                                                       r["padj"], r["call"]) for r in rows) or None


def m_top(direction):
    def f(ctx):
        rows = [r for r in metric_tsv(ctx, "top_genes.tsv") if r["track"] == "genome" and r["direction"] == direction]
        return ", ".join("%s (%.2f)" % (r["external_gene_name"], float(r["log2FoldChange"])) for r in rows) or None
    return f


def m_pca(track):
    def f(ctx):
        p = [r for r in metric_tsv(ctx, "pca.tsv") if r["track"] == track]
        lib = [r for r in metric_tsv(ctx, "libsize.tsv") if r["track"] == track]
        if not p:
            return None
        s = "PC1 %s%%, PC2 %s%%, PC1 separates IR from mock: %s" % (p[0]["PC1_pct"], p[0]["PC2_pct"], p[0]["PC1_separates_condition"])
        if lib:
            c = [float(r["counts"]) / 1e6 for r in lib]
            sf = [float(r["size_factor"]) for r in lib]
            s += "; library sizes %s million; size factors %s" % (rng(c), rng(sf, "%.2f"))
        return s
    return f


def m_pca_batch(ctx):
    p = [r for r in metric_tsv(ctx, "pca.tsv") if r["track"] == "genome"]
    return ("PC1 separates condition: %s (PC1 %s%%)" % (p[0]["PC1_separates_condition"], p[0]["PC1_pct"])) if p else None


def m_compare(ctx):
    r = metric_tsv(ctx, "track_comparison.tsv")
    return ", ".join("%s=%s" % kv for kv in sorted(r[0].items())) if r else None


def m_ora_gsea(ctx):
    d = os.path.join(ctx.W, "results", "enrichment")
    g = read_csv(os.path.join(d, "GSEA_GO_BP.csv"))
    up = {r["ID"] for r in read_csv(os.path.join(d, "GO_BP_up_regulated.csv"))}
    dn = {r["ID"] for r in read_csv(os.path.join(d, "GO_BP_down_regulated.csv"))}
    if not g:
        return None
    both = [r for r in g if r["ID"] in up or r["ID"] in dn]
    agree = sum(1 for r in both if (float(r["NES"]) > 0) == (r["ID"] in up))
    top = sorted(g, key=lambda r: float(r["p.adjust"]))[:5]
    return "%d GSEA terms; %d also in direction-specific ORA, sign agrees for %d; top: %s" % (
        len(g), len(both), agree, "; ".join("%s NES %.2f" % (r["Description"], float(r["NES"])) for r in top))


def read_csv(p):
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def m_ties(ctx):
    for r in read_tsv(os.path.join(ctx.rec, "events", "06-genome.tsv")):
        m = re.search(r'ties in the preranked stats \(([\d.]+)% of the list\)', r["message"])
        if m:
            return m.group(1)
    return None


def m_fastqc_dup(ctx):
    fq = fastqc(ctx)
    if not fq:
        return None
    per = {}
    for name, (_, dd) in fq.items():
        per.setdefault(short(name.rsplit("_R", 1)[0]), []).append(dd)
    return "deduplicated percentage (mean of R1, R2): " + ", ".join(
        "%s %.0f" % (k, statistics.mean(v)) for k, v in sorted(per.items()))


def m_tile_n(ctx):
    a = count_status(ctx, "Per tile sequence quality", "FAIL")
    b = count_status(ctx, "Per base N content", "FAIL")
    return "per tile FAIL %s; per base N content FAIL %s" % (a, b) if a else None


def m_tx2gene_coverage(ctx):
    """transcripts in the basic GTF vs in the transcript FASTA the Kallisto index is built from"""
    gtf = ctx.path("data", "gencode.vM38.primary_assembly.basic.annotation.gtf")
    fa = ctx.path("data", "gencode.vM38.transcripts.fa")
    if not (os.path.exists(gtf) and os.path.exists(fa)):
        return None
    n_gtf = sum(1 for line in open(gtf) if not line.startswith("#") and line.split("\t", 3)[2:3] == ["transcript"])
    n_fa = sum(1 for line in open(fa) if line.startswith(">"))
    return "%d of %d transcripts (%.0f%%)" % (n_gtf, n_fa, 100.0 * n_gtf / n_fa)


def m_pkg(ctx):
    p = os.path.join(ctx.rec, "preflight", "pkgcheck.txt")
    return open(p).read().strip().replace("\n", " | ") if os.path.exists(p) else None


def m_du(ctx):
    p = os.path.join(ctx.rec, "listings", "01-du.txt")
    return open(p).read().strip() if os.path.exists(p) else None


METRICS = {
    "02-time-star-map-per-sample": m_star_map_time,
    "03-sol-perbase-fail": lambda c: count_status(c, "Per base sequence quality", "FAIL"),
    "03-sol-tile-n-fail": m_tile_n,
    "03-sol-adapter-pass": lambda c: count_status(c, "Adapter Content", "PASS"),
    "03-sol-duplication": m_fastqc_dup,
    "03-sol-gc-fail": lambda c: count_status(c, "Per sequence GC content", "FAIL"),
    "04a-time-salmon": m_single_time("04a-salmon", "Salmon index + quant session, 4 cores"),
    "04a-txt-isf-isr": m_isf_isr,
    "04a-txt-strand-bias": m_bias,
    "04a-time-star-index": m_single_time("04a-index", "STAR genomeGenerate, 48 cores"),
    "04a-time-star-map": m_star_map_time,
    "04a-sol-unique-range": m_unique,
    "04a-sol-unique-max": m_unique_max,
    "04a-sol-unique-outliers": m_unique_outliers,
    "04a-sol-too-short": m_too_short,
    "04a-time-featurecounts": m_single_time("04a-count", "featureCounts, 16 cores"),
    "04a-sol-assigned-range": m_assigned_range,
    "04a-sol-assigned-min": m_assigned_min,
    "04a-sol-unassigned-largest": m_unassigned_largest,
    "04a-sol-assigned-vs-depth": m_assigned_vs_depth,
    "04b-time-quant-intro": m_kallisto_quant,
    "04b-time-index": m_single_time("04b-index", "kallisto index, 48 cores (single-threaded)"),
    "04b-time-quant": m_kallisto_quant,
    "04b-sol-pseudo-range": m_pseudo_range,
    "04b-txt-tx2gene-coverage": m_tx2gene_coverage,
    "04b-sol-pseudo-consistency": m_pseudo_consistency,
    "04b-sol-fraglen": m_fraglen,
    "05-sol-exploratory": m_pca("genome"),
    "05-txt-pca-batch": m_pca_batch,
    "05-sol-top-up": m_top("up"),
    "05-sol-top-down": m_top("down"),
    "05-sol-gadd45a": lambda c: anchors_text(c, "Gadd45a"),
    "05b-sol-exploratory": m_pca("transcript"),
    "05b-sol-compare-tracks": m_compare,
    "06-txt-gsea-ties": m_ties,
    "06-sol-ora-gsea": m_ora_gsea,
    "setup-pkg-test": m_pkg,
    "setup-data-size": m_du,
}


# ------------------------------------------------------------ outputs and status
def console_for(ctx, step, ref):
    ep, nn = ref.split("/")
    p = os.path.join(ctx.rec, "console", step, "%s_%s.txt" % (ep, nn))
    if os.path.exists(p):
        return open(p, errors="replace").read()
    p = os.path.join(ctx.run, "logs", "session-%s.log" % step)
    if os.path.exists(p):
        # "begin <epoch>" is printed output; R also echoes the cat() command, which has no digits there
        m = re.search(r'### KIT-BLOCK %s begin \d+[^\n]*\n(.*?)### KIT-BLOCK %s end \d+' % (re.escape(ref), re.escape(ref)),
                      open(p, errors="replace").read(), re.S)
        if m:
            lines = m.group(1).split("\n")
            while lines and (lines[-1].startswith('> cat("### KIT-BLOCK') or not lines[-1].strip()):
                lines.pop()
            return "\n".join(lines) + "\n"
    return None


def first_error(ctx, step, kind):
    if kind == "ood":
        for r in read_tsv(os.path.join(ctx.rec, "events", step + ".tsv")):
            if r["type"] == "ERROR":
                return "block %s (line %s, %s): %s" % (r["block"], r["line"], r["context"], r["message"][:200])
        return ""
    p = os.path.join(ctx.run, "logs", "session-%s.log" % step)
    if os.path.exists(p):
        for line in open(p, errors="replace"):
            if line.startswith("### KIT-ERR"):
                return line.strip()[:240]
    return ""


def status_of(ctx, step):
    p = os.path.join(ctx.run, "markers", step + ".done")
    st = ctx.job_stats(step)
    if os.path.exists(p):
        m = open(p).read().split()[0]
        if m == "RUNNING" and st and not {x[3].split()[0] for x in st} & {"RUNNING", "PENDING"}:
            return "FAIL(" + ",".join(sorted({x[3].split()[0] for x in st})) + ")"
        return m
    if st:
        states = {x[3].split()[0] for x in st}
        if states == {"COMPLETED"}:
            return "PASS"
        return "FAIL(" + ",".join(sorted(states)) + ")"
    return "NOT RUN"


def main():
    run, gen, out = sys.argv[1:4]
    ctx = Ctx(run, gen)
    os.makedirs(os.path.join(out, "outputs"), exist_ok=True)
    md = []
    env = open(os.path.join(run, "kit.env")).read() if os.path.exists(os.path.join(run, "kit.env")) else ""
    md.append("# Negishi run kit summary: %s\n" % os.path.basename(run))
    md.append("Kit generated from episodes with checksums in `generated/episodes.sha256` "
              "(prefix %s). See `checks/episode_checksums.txt` for the comparison done at collection time.\n"
              % (re.search(r"KIT_BUILT_FROM=(\S+)", env).group(1) if "KIT_BUILT_FROM" in env else "?"))

    # per-step table
    md.append("## Steps\n")
    md.append("| Step | Status | Job | Elapsed (min) | Limit | Peak mem (GB) | Requested mem | Queue wait (min) | First error |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    ep_status = {}
    for s in ctx.steps:
        st = status_of(ctx, s["step"])
        stats = ctx.job_stats(s["step"])
        el = [x[1] for x in stats if x[1] is not None]
        rs = [x[2] for x in stats if x[2] is not None]
        if s["kind"] == "ood":
            # sacct does not see processes inside the container reliably; use R's own VmHWM
            hw = [r["message"] for r in read_tsv(os.path.join(ctx.rec, "events", s["step"] + ".tsv"))
                  if r["block"] == "ALL" and r["type"] == "MAXRSS"]
            rs = [rss_gb(hw[0].replace(" kB", "K"))] if hw and rss_gb(hw[0].replace(" kB", "K")) else []
        ids = ctx.jobids(s["step"])
        e = ctx.sacct.get(ids[-1], {}) if ids else {}
        wait = ""
        try:
            import datetime
            fmt = "%Y-%m-%dT%H:%M:%S"
            wait = "%.0f" % ((datetime.datetime.strptime(e["start"], fmt) - datetime.datetime.strptime(e["submit"], fmt))
                             .total_seconds() / 60.0)
        except Exception:
            pass
        err = first_error(ctx, s["step"], s["kind"]) if not st.startswith("PASS") else ""
        md.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            s["step"], st, ids[-1] if ids else "-", rng([x / 60.0 for x in el], "%.0f") if el else "-",
            s["time"], "%.1f" % max(rs) if rs else "-", s["mem"], wait or "-", err.replace("|", "/")))
        ep = s["step"].split("-")[0] if s["step"] != "metrics" else "kit"
        ep_status.setdefault(ep, []).append(st)
    md.append("")
    md.append("## Episodes\n")
    md.append("| Episode | Result |")
    md.append("|---|---|")
    for ep, sts in sorted(ep_status.items()):
        res = "FAIL" if any(x.startswith("FAIL") for x in sts) else (
            "NOT RUN" if all(x == "NOT RUN" for x in sts) else (
                "PASS" if all(x == "PASS" for x in sts) else "PASS with warnings / partial"))
        md.append("| %s | %s (%s) |" % (ep, res, ", ".join(sts)))
    md.append("")

    # data checks written by kit extras
    md.append("## Data checks\n")
    for f, title in (("02-ref-compare.tsv", "Fresh GENCODE vM38 downloads vs the staged references (match = md5 equal)"),
                     ("02-subsample-check.txt", "Subsample spoiler on 10,000-read inputs (expect 5000 reads, names paired, originals archived)")):
        rows = read_tsv(os.path.join(ctx.rec, "listings", f))
        md.append("%s: %s" % (title, "not run" if not rows else ", ".join(
            "%s %s" % (r.get("file", r.get("sample")), r.get("match", "%s/%s %s %s" % (r.get("R1_reads"), r.get("R2_reads"),
                       r.get("paired_names_match"), r.get("originals_archived")))) for r in rows)))
        md.append("")

    # anchors
    md.append("## Anchors (IR vs mock; shrunken LFC from the saved tables)\n")
    rows = metric_tsv(ctx, "anchors.tsv")
    if rows:
        md.append("| Track | Gene | log2FC | padj | Call | Passes padj<=0.05 and abs(LFC)>=log2(1.5) |")
        md.append("|---|---|---|---|---|---|")
        for r in rows:
            md.append("| %s | %s | %.2f | %s | %s | %s |" % (r["track"], r["external_gene_name"], float(r["log2FoldChange"]),
                                                          r["padj"], r["call"], r["passes_episode_cutoffs"]))
    else:
        md.append("No anchor table (metrics step did not run or the DE tables are missing).")
    md.append("")

    # runtime vs claims
    md.append("## Runtime against the lesson\n")
    md.append("Front matter minutes are teaching + exercises per episode. Compute time is the sum of the "
              "episode's job elapsed times (array tasks run in parallel, so the longest task counts).\n")
    fm = {}
    for f in sorted(glob.glob(os.path.join(os.path.dirname(gen.rstrip("/")), "..", "episodes", "*.Rmd"))):
        t = open(f).read()
        a = re.search(r'^teaching:\s*(\d+)', t, re.M)
        b = re.search(r'^exercises:\s*(\d+)', t, re.M)
        if a and b:
            fm[os.path.basename(f).split("-")[0]] = int(a.group(1)) + int(b.group(1))
    md.append("| Episode | Front matter (min) | Compute (min) | Queue wait (min, max) |")
    md.append("|---|---|---|---|")
    for ep in sorted(fm):
        comp, waits = 0.0, []
        for s in ctx.steps:
            if s["step"].split("-")[0] != ep:
                continue
            el = [x[1] for x in ctx.job_stats(s["step"]) if x[1] is not None]
            comp += max(el) if el else 0
        md.append("| %s | %d | %.0f | %s |" % (ep, fm[ep], comp / 60.0, "-"))
    md.append("")

    # placeholders
    vals = []
    for r in read_tsv(os.path.join(gen, "placeholders.tsv")):
        v, src = None, r["source"]
        if r["kind"] == "metric":
            ctx.used = []
            try:
                v = METRICS[r["id"]](ctx)
            except Exception as e:  # report, never crash the summary
                v = "ERROR computing metric: %s" % e
            old = ctx.stale(ctx.used)
            if v is not None and old:
                v = "STALE (from %d file(s) older than this run, e.g. %s): %s" % (len(old), os.path.relpath(old[0], ctx.W), v)
        elif r["kind"] == "file":
            p = src.replace("learner:", ctx.W + "/") if src.startswith("learner:") else os.path.join(run, src)
            if os.path.exists(p):
                text = open(p, errors="replace").read()
                open(os.path.join(out, "outputs", r["id"] + ".txt"), "w").write(text)
                v = "outputs/%s.txt" % r["id"]
                if ctx.stale([p]):
                    v = "STALE (file older than this run): " + v
        else:
            parts = []
            for s in src.split(";"):
                _, step, ref = s.split(":")
                text = console_for(ctx, step, ref)
                if text is not None:
                    fn = "%s%s.txt" % (r["id"], "" if step != "06-transcript" else "@transcript")
                    open(os.path.join(out, "outputs", fn), "w").write(text)
                    parts.append("outputs/" + fn)
            v = ", ".join(parts) or None
        vals.append((r["id"], r["file"], r["line"], r["kind"], v if v is not None else "MISSING"))
    with open(os.path.join(out, "placeholder_values.tsv"), "w") as fh:
        fh.write("id\tfile\tline\tkind\tvalue\n")
        for x in vals:
            fh.write("\t".join(str(y).replace("\t", " ").replace("\n", " | ") for y in x) + "\n")
    miss = [x[0] for x in vals if x[4] == "MISSING"]
    md.append("## Placeholder values\n")
    md.append("%d placeholders, %d missing%s. Full list: `placeholder_values.tsv`; outputs in `outputs/`.\n" % (
        len(vals), len(miss), (": " + ", ".join(miss)) if miss else ""))
    md.append("| Placeholder | Value |")
    md.append("|---|---|")
    for x in vals:
        if x[3] == "metric":
            md.append("| %s | %s |" % (x[0], str(x[4]).replace("|", "/")))
    md.append("")

    # every metric, with or without a placeholder: compare these with the numbers in the episodes
    md.append("## Measured values\n")
    md.append("All metrics, including those whose placeholders are already filled; compare with the episode text.\n")
    md.append("| Metric | Value |")
    md.append("|---|---|")
    for mid in sorted(METRICS):
        ctx.used = []
        try:
            v = METRICS[mid](ctx)
        except Exception as e:
            v = "ERROR computing metric: %s" % e
        old = ctx.stale(ctx.used)
        if v is not None and old:
            v = "STALE: %s" % v
        md.append("| %s | %s |" % (mid, str(v).replace("|", "/")))
    md.append("")

    # warnings digest from R steps
    md.append("## R warnings (unique, per step)\n")
    for f in sorted(glob.glob(os.path.join(ctx.rec, "events", "*.tsv"))):
        ws = sorted({r["message"][:160] for r in read_tsv(f) if r["type"] == "WARNING"})
        if ws:
            md.append("- `%s`: %d unique; %s" % (os.path.basename(f)[:-4], len(ws), " / ".join(w.replace("|", "/") for w in ws[:8])))
    md.append("")
    open(os.path.join(out, "SUMMARY.md"), "w").write("\n".join(md) + "\n")
    print("summarize: wrote %s (%d placeholders, %d missing)" % (os.path.join(out, "SUMMARY.md"), len(vals), len(miss)))


if __name__ == "__main__":
    main()
