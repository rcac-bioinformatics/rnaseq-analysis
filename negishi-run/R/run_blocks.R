# run_blocks.R <plan.tsv>
# Runs the R blocks of one episode, in order, in this one fresh R session, the way a
# learner pastes them into the RStudio console. Per block it records the console output
# (commands echoed, printed values, messages, warnings, errors) to
# $KIT_RECORDS/console/<step>/<episode>_<block>.txt, plots to
# $KIT_RECORDS/figures/<step>/<episode>_<block>_NN.png, and events to
# $KIT_RECORDS/events/<step>.tsv. An error stops the rest of that block (as source()
# does) and the run continues with the next block, so later failures are visible too.
# Exit status: 0 all blocks ran, 3 only blocks inside a spoiler failed, 1 otherwise.
args <- commandArgs(trailingOnly = TRUE)
plan <- read.delim(args[1], colClasses = "character")
rec <- Sys.getenv("KIT_RECORDS"); step <- Sys.getenv("KIT_STEP"); gen <- Sys.getenv("KIT_GEN")
condir <- file.path(rec, "console", step); figdir <- file.path(rec, "figures", step)
dir.create(condir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(rec, "events"), showWarnings = FALSE)
dir.create(file.path(rec, "sessioninfo"), showWarnings = FALSE)
evf <- file.path(rec, "events", paste0(step, ".tsv"))
cat("block\tline\tcontext\ttype\tmessage\n", file = evf)
ev <- function(b, type, msg) {
  cat(sprintf("%s\t%s\t%s\t%s\t%s\n", b$block, b$line, b$context, type,
              gsub("[\t\r\n]+", " | ", msg)), file = evf, append = TRUE)
}

cat("R:", R.version.string, "\n")
cat("Bioconductor:", tryCatch(as.character(packageVersion("BiocVersion")), error = function(e) "BiocVersion not installed"), "\n")
cat(".libPaths():\n"); print(.libPaths())
cat("USER =", Sys.getenv("USER"), " HOME =", Sys.getenv("HOME"), " TMPDIR =", Sys.getenv("TMPDIR"), "\n")

t0 <- Sys.time(); n_err_required <- 0L; n_err_optional <- 0L
for (i in seq_len(nrow(plan))) {
  b <- as.list(plan[i, ])
  tag <- sprintf("%s_%s", b$episode, b$block)
  if (startsWith(b$action, "skip")) { ev(b, "SKIPPED", sub("^skip:", "", b$action)); next }
  f <- file.path(gen, b$file)
  cat(sprintf("\n######## %s block %s (Rmd line %s, %s) ########\n", b$episode, b$block, b$line, b$context))
  con <- file(file.path(condir, paste0(tag, ".txt")), open = "wt")
  sink(con); sink(con, type = "message")
  grDevices::png(file.path(figdir, paste0(tag, "_%02d.png")), width = 1400, height = 1000, res = 150)
  tb <- Sys.time()
  failed <- FALSE
  withCallingHandlers(
    tryCatch(
      source(f, echo = TRUE, print.eval = TRUE, max.deparse.length = Inf, local = globalenv(),
             spaced = FALSE, prompt.echo = "> ", continue.echo = "+ "),
      error = function(e) {
        failed <<- TRUE
        cat("Error:", conditionMessage(e), "\n")
        ev(b, "ERROR", conditionMessage(e))
      }),
    warning = function(w) {
      cat("Warning message:\n", conditionMessage(w), "\n", sep = "")
      ev(b, "WARNING", conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) ev(b, "MESSAGE", conditionMessage(m)))
  while (grDevices::dev.cur() > 1) grDevices::dev.off()
  sink(type = "message"); sink(); close(con)
  if (failed) {
    if (grepl("spoiler", b$context)) n_err_optional <- n_err_optional + 1L else n_err_required <- n_err_required + 1L
  }
  ev(b, "TIME_S", sprintf("%.1f", as.numeric(difftime(Sys.time(), tb, units = "secs"))))
}
# devices opened for blocks that drew nothing leave no file or an empty one
for (p in list.files(figdir, pattern = "\\.png$", full.names = TRUE)) if (file.size(p) < 3000) file.remove(p)
tot <- list(block = "ALL", line = "", context = "")
ev(tot, "TIME_S", sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
hwm <- tryCatch(sub("^VmHWM:\\s*", "", grep("^VmHWM", readLines("/proc/self/status"), value = TRUE)),
                error = function(e) "NA")
ev(tot, "MAXRSS", hwm)
writeLines(capture.output(sessionInfo()), file.path(rec, "sessioninfo", paste0(step, ".txt")))
cat(sprintf("\nrun_blocks: %s done, %d required and %d optional block(s) failed\n", step, n_err_required, n_err_optional))
quit(status = if (n_err_required > 0) 1L else if (n_err_optional > 0) 3L else 0L, save = "no")
