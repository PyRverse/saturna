# Genuine satuRn reference for the pasilla real-data benchmark.
# Usage: Rscript --vanilla benchmarks/pasilla/satuRn.R <raw counts.txt> <raw metadata.txt>
# Run from the repository root. Library locations belong in R_LIBS, never in outputs.
suppressPackageStartupMessages(library(satuRn))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
bench <- "benchmarks/pasilla"
out_dir <- file.path(bench, "reference")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.delim(args[2], stringsAsFactors = FALSE)
samples <- metadata$SampleName
stopifnot(identical(colnames(raw), c("feature_id", "gene_id", samples)))
rownames(raw) <- raw$feature_id

# Independent re-derivation of the shared filter (see prepare.py) from the raw table.
values <- as.matrix(raw[, samples])
gene_totals <- rowsum(values, raw$gene_id)[raw$gene_id, , drop = FALSE]
gene_ok <- rowSums(gene_totals >= 10) == length(samples)
transcript_ok <- rowSums(values >= 10) >= 3
passing <- gene_ok & transcript_ok
multi <- passing & ave(as.integer(passing), raw$gene_id, FUN = sum) >= 2
kept <- raw[multi, ]
distinct <- tapply(seq_len(nrow(kept)), kept$gene_id, function(i) {
    nrow(unique(as.matrix(kept[i, samples, drop = FALSE])))
})
uninformative <- names(distinct)[distinct == 1]
expected_set <- ifelse(kept$gene_id %in% uninformative, "uninformative", "primary")

inputs <- read.csv(file.path(bench, "inputs/counts.csv"), check.names = FALSE,
                   stringsAsFactors = FALSE)
design <- read.csv(file.path(bench, "inputs/design.csv"), row.names = 1)
stopifnot(identical(inputs$transcript_id, kept$feature_id),
          identical(inputs$gene_id, kept$gene_id),
          identical(inputs$analysis_set, expected_set),
          all(as.matrix(inputs[, samples]) == as.matrix(kept[, samples])),
          identical(rownames(design), samples),
          identical(names(design), c("intercept", "knockdown", "paired_end")),
          all(design$intercept == 1),
          identical(design$knockdown, as.integer(metadata$condition == "KD")),
          identical(design$paired_end, as.integer(metadata$LibraryLayout == "PAIRED")))

fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))

# Repeats the estimation steps of satuRn:::p.adjust_empirical, which returns only
# p-values, so the fitted null itself can be reported and compared.
empirical_null <- function(pval, t) {
    z <- qnorm(pval / 2) * sign(t)
    zm <- z[abs(z) < 10]
    zm <- zm[!is.na(zm)]
    b <- 4.3 * exp(-0.26 * log(length(zm), 10))
    sc <- diff(quantile(zm)[c(2, 4)]) / (2 * qnorm(0.75))
    first <- satuRn:::locfdr_locmle(zm, xlim = c(median(zm), b * sc))
    breaks <- seq(min(zm), max(zm), length = 120)
    x <- (breaks[-1] + breaks[-length(breaks)]) / 2
    y <- hist(zm, breaks = breaks, plot = FALSE)$counts
    cov_in <- list(x = x, X = cbind(1, poly(x, df = 7)),
                   f = glm(y ~ poly(x, df = 7), poisson)$fit, sw = 0)
    second <- satuRn:::locfdr_locmle(zm, xlim = c(first[1], b * first[2]),
                                     d = first[1], s = first[2], Cov.in = cov_in)$mle
    list(location = unname(second[1]), scale = unname(second[2]),
         observations = length(zm), width_factor = b,
         first_stage = list(location = unname(first[1]), scale = unname(first[2])), z = z)
}
contrast <- matrix(c(0, 1, 0), ncol = 1,
                   dimnames = list(names(design), "knockdown"))
formula <- ~ 0 + intercept + knockdown + paired_end
analyses <- list(primary = inputs$analysis_set == "primary",
                 unrestricted = rep(TRUE, nrow(inputs)))
summary <- list()
for (name in names(analyses)) {
    rows <- inputs[analyses[[name]], ]
    counts <- as.matrix(rows[, samples])
    rownames(counts) <- rows$transcript_id
    object <- SummarizedExperiment(assays = list(counts = counts),
        rowData = DataFrame(isoform_id = rows$transcript_id, gene_id = rows$gene_id),
        colData = DataFrame(design))
    fit_warnings <- character()
    object <- withCallingHandlers(
        fitDTU(object, formula = formula, verbose = FALSE),
        warning = function(w) {
            fit_warnings <<- c(fit_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    plot_file <- tempfile(fileext = ".pdf")
    pdf(plot_file)
    tested <- testDTU(object, contrast, diagplot1 = FALSE, diagplot2 = FALSE)
    dev.off()
    unlink(plot_file)
    models <- rowData(object)[["fitDTUModels"]]
    params <- function(field) vapply(models, function(x) {
        value <- x@params[[field]]
        if (length(value) == 1 && is.na(value)) NA_real_ else as.numeric(value)
    }, numeric(1))
    coefs <- t(vapply(models, function(x) {
        value <- x@params$coefficients
        if (length(value) == 3) unname(value) else rep(NA_real_, 3)
    }, numeric(3)))
    result <- as.data.frame(rowData(tested)[["fitDTUResult_knockdown"]])
    moderation <- limma::squeezeVar(params("dispersion"), params("df.residual"), robust = FALSE)
    null <- empirical_null(result$pval, result$t)
    replayed <- 2 * pnorm(-abs((null$z - null$location) / null$scale))
    stopifnot(isTRUE(all.equal(replayed, result$empirical_pval, tolerance = 1e-14)))
    out <- data.frame(
        transcript_id = rows$transcript_id,
        gene_id = rows$gene_id,
        model_type = vapply(models, function(x) x@type, character(1)),
        coef_intercept = fmt(coefs[, 1]),
        coef_knockdown = fmt(coefs[, 2]),
        coef_paired_end = fmt(coefs[, 3]),
        dispersion = fmt(params("dispersion")),
        residual_df = fmt(params("df.residual")),
        variance_posterior = fmt(vapply(models, function(x) x@varPosterior, numeric(1))),
        df_posterior = fmt(vapply(models, function(x) x@dfPosterior, numeric(1))),
        estimates = fmt(result$estimates),
        se = fmt(result$se),
        df = fmt(result$df),
        t = fmt(result$t),
        pval = fmt(result$pval),
        regular_FDR = fmt(result$regular_FDR),
        empirical_pval = fmt(result$empirical_pval),
        empirical_FDR = fmt(result$empirical_FDR))
    write.csv(out, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE, quote = FALSE)
    summary[[name]] <- list(
        transcripts = nrow(rows),
        genes = length(unique(rows$gene_id)),
        model_types = as.list(table(out$model_type)),
        fit_warnings = as.list(table(fit_warnings)),
        prior_df = moderation$df.prior,
        prior_variance = moderation$var.prior,
        empirical_null = null[c("location", "scale", "observations", "width_factor",
                                "first_stage")])
}

versions <- vapply(c("satuRn", "limma", "locfdr", "SummarizedExperiment", "jsonlite"),
                   function(pkg) as.character(packageVersion(pkg)), character(1))
write_json(list(R = as.character(getRversion()), packages = as.list(versions),
                analyses = summary),
           file.path(out_dir, "runtime.json"), digits = NA, auto_unbox = TRUE, pretty = TRUE)
