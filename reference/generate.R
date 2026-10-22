# Genuine R reference generation; run from the repository root.
# Library locations belong in R_LIBS / R_LIBS_USER, never fixture metadata.
suppressPackageStartupMessages(library(satuRn))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(jsonlite))

dir.create("reference/outputs", showWarnings = FALSE, recursive = TRUE)
write_json_strict <- function(value, path) {
    write_json(value, path, digits = NA, auto_unbox = TRUE, pretty = TRUE,
               na = "string", dataframe = "columns", matrix = "rowmajor")
}
versions <- vapply(c("satuRn", "limma", "locfdr", "jsonlite", "SummarizedExperiment"),
                   function(pkg) as.character(packageVersion(pkg)), character(1))
write_json_strict(list(R = as.character(getRversion()), packages = as.list(versions)),
                  "reference/outputs/runtime.json")

for (scenario in c("balanced", "unequal_df", "empirical", "covariate")) {
    counts_input <- read.csv(file.path("datasets", scenario, "counts.csv"),
                             row.names = 1, check.names = FALSE)
    design <- read.csv(file.path("datasets", scenario, "design.csv"), row.names = 1)
    counts <- as.matrix(counts_input[, -1])
    object <- SummarizedExperiment(assays = list(counts = counts),
        rowData = DataFrame(isoform_id = rownames(counts), gene_id = counts_input$gene_id),
        colData = DataFrame(design))
    formula <- as.formula(paste("~ 0 +", paste(names(design), collapse = " + ")))
    object <- fitDTU(object, formula = formula, verbose = FALSE)
    contrast_values <- if (scenario == "covariate") c(0, 1, .5) else c(0, 1)
    contrast <- matrix(contrast_values, ncol = 1,
                       dimnames = list(names(design), "condition"))
    plot_file <- tempfile(fileext = ".pdf")
    pdf(plot_file)
    tested <- testDTU(object, contrast, diagplot1 = FALSE, diagplot2 = FALSE)
    dev.off()
    unlink(plot_file)
    models <- rowData(object)[["fitDTUModels"]]
    extract <- function(name) unname(lapply(models, function(x) unname(x@params[[name]])))
    output <- list(
        transcript_ids = rownames(counts),
        contrast = contrast_values,
        coefficients = extract("coefficients"),
        covariance_unscaled = extract("vcovUnsc"),
        dispersion = unname(vapply(models, function(x) x@params$dispersion, numeric(1))),
        residual_df = unname(vapply(models, function(x) x@params$df.residual, numeric(1))),
        variance_posterior = unname(vapply(models, function(x) x@varPosterior, numeric(1))),
        df_posterior = unname(vapply(models, function(x) x@dfPosterior, numeric(1))),
        result = as.data.frame(rowData(tested)[["fitDTUResult_condition"]]))
    write_json_strict(output, file.path("reference/outputs", paste0(scenario, ".json")))
}

moderation_cases <- list(
    equal = list(var = c(.02, .2, 1, 2, 4, 16, 20, 80), df = rep(6, 8)),
    infinite_prior = list(var = c(.8, .9, 1, 1.1, 1.2), df = rep(8, 5)),
    unequal = list(var = c(.02, .2, 1, 2, 4, 16, 20, 80), df = c(4, 6, 8, 4, 6, 8, 4, 6)),
    zeros = list(var = c(0, 0, 1, 2, 20), df = rep(6, 5)),
    two = list(var = c(1, 4), df = c(4, 6)))
moderation <- lapply(moderation_cases, function(input) {
    list(input = input, output = limma::squeezeVar(input$var, input$df, robust = FALSE))
})
write_json_strict(moderation, "reference/outputs/moderation.json")

# Exercise a shifted and widened central null independently of transcript fitting.
z <- c(.25 + 1.3 * qnorm(((1:900) - .5) / 900), seq(3.8, 5.2, length.out = 50))
p <- 2 * pnorm(-abs(z))
t <- -sign(z) * abs(z)
corrected <- satuRn:::p.adjust_empirical(p, t, "reference", FALSE, FALSE)
write_json_strict(list(pvalues = p, statistics = t, result = corrected),
                  "reference/outputs/empirical_null.json")
