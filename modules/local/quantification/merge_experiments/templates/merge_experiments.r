#!/usr/bin/env Rscript --vanilla

library(SummarizedExperiment)

paths <- c('${experiments.join("\', \'")}')
experiments <- lapply(paths, readRDS)

phenotype <- read.csv('${phenotype}', stringsAsFactors = FALSE)
annotation <- rtracklayer::import('${gtf}')
tpm <- read.table('${tpm}', header=TRUE, row.names=1)[, -1]

se_assays <- list()

for (se in experiments) {
    assays <- assays(se)
    # Iterate over named list of assays
    for (assay_name in names(assays)) {
        assay <- assays[[assay_name]]

        # Add assay to se_assays for its name
        if (is.null(se_assays[[assay_name]])) {
            se_assays[[assay_name]] <- assay
        } else {
            se_assays[[assay_name]] <- cbind(se_assays[[assay_name]], assay)
        }
    }
}

se_cbind <- do.call(SummarizedExperiment::cbind, experiments)
se <- SummarizedExperiment(assays = se_assays, colData = colData(se_cbind), rowData = rowData(se_cbind))

# Join phenotype data
colData(se) <- merge(colData(se), phenotype, by.x="names", by.y=colnames(phenotype)[1])

# Convert string columns to factors
for (col in colnames(colData(se))) {
    if (is.character(colData(se)[[col]]) && !(col == "names")) {
        colData(se)[[col]] <- as.factor(colData(se)[[col]])
    }
}

# Add transcript annotation
annotation <- annotation[match(rownames(se), annotation\$transcript_id),]
rowData(se) <- annotation

# Add TPM
colnames(tpm) <- gsub("^X(\\\\d+)", "\\\\1", colnames(tpm)) # Remove ^X prefixes in tpm colnames if followed by a digit
assay(se, "tpm", withDimnames = FALSE) <- tpm[rownames(se), colData(se)\$names]

saveRDS(se, '${meta.id}.merged.rds')

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-summarizedexperiment:', packageVersion('SummarizedExperiment'))
    ),
'versions.yml')
