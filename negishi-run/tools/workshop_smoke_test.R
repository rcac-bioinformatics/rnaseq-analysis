# workshop_smoke_test.R: check an Open OnDemand RStudio (bioconductor) session before teaching.
# In the RStudio console:  source("~/rnaseq-analysis/negishi-run/tools/workshop_smoke_test.R")
# Prints the R and Bioconductor versions, the library paths, whether every package the
# R episodes load is available, scratch access, and internet access for KEGG, MSigDB,
# and Ensembl BioMart. Ends with "All packages load" when the package check passes.
local({
  t0 <- Sys.time()
  cat("R:", R.version.string, "\n")
  cat("Bioconductor:", tryCatch(as.character(packageVersion("BiocVersion")), error = function(e) "unknown"), "\n")
  cat(".libPaths():\n"); for (p in .libPaths()) cat("  ", p, "\n")
  cat("cores available:", parallel::detectCores(), " TMPDIR:", tempdir(), "\n\n")

  pkgs <- c("tidyverse", "RColorBrewer", "pheatmap", "ggrepel", "reshape2", "hexbin", "msigdbr",
            "DESeq2", "apeglm", "vsn", "biomaRt", "clusterProfiler", "enrichplot", "org.Mm.eg.db",
            "dplyr", "ggplot2", "readr")
  loads <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  ver <- vapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_), "")
  for (p in pkgs) cat(sprintf("  %-16s %s\n", p, if (loads[[p]]) ver[[p]] else "DOES NOT LOAD"))
  missing <- pkgs[!loads]

  sc <- file.path("/scratch/negishi", Sys.getenv("USER"))
  cat("\nscratch", sc, if (dir.exists(sc) && file.access(sc, 2) == 0) "is writable" else "is NOT writable", "\n")

  net <- function(label, url) {
    ok <- tryCatch({
      h <- curl::new_handle(timeout = 20, followlocation = TRUE)   # GET: BioMart refuses HEAD
      curl::curl_fetch_memory(url, handle = h)$status_code
    }, error = function(e) conditionMessage(e))
    cat(sprintf("  %-22s %s\n", label, ok))
  }
  cat("internet from this session:\n")
  net("KEGG REST", "https://rest.kegg.jp/info/kegg")
  net("Ensembl BioMart", "https://www.ensembl.org/biomart/martservice?type=registry")
  net("MSigDB (zenodo)", "https://zenodo.org/")
  cat(sprintf("\nchecked in %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (length(missing)) message("Missing: ", paste(missing, collapse = ", ")) else message("All packages load")
})
