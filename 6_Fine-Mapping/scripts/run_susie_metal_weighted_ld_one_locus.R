suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: Rscript run_susie_metal_weighted_ld_one_locus.R <sumstats_file> <out_prefix> <L> <weights_file> <analysis_id> <ANC:ld_zst:vars> ...")
}

sumstats_file <- args[1]
out_prefix <- args[2]
L <- as.integer(args[3])
weights_file <- args[4]
analysis_id_arg <- args[5]
ld_args <- args[6:length(args)]

message("Reading sumstats: ", sumstats_file)
gwas <- fread(sumstats_file)

required <- c("variant_id", "CHR", "BP", "REF", "ALT", "BETA_REF", "Z_REF", "P", "N")
missing <- setdiff(required, names(gwas))
if (length(missing) > 0) {
  stop("Missing required sumstats columns: ", paste(missing, collapse=", "))
}

gwas[, Z_REF := as.numeric(Z_REF)]
gwas[, N := as.numeric(N)]

if (any(!is.finite(gwas$Z_REF))) {
  stop("Non-finite Z_REF in sumstats.")
}

weights <- fread(weights_file)
weights[, analysis_id := as.character(analysis_id)]
weights[, ancestry := as.character(ancestry)]

weights <- weights[analysis_id == analysis_id_arg]

if (nrow(weights) == 0) {
  stop("No weights found for analysis_id=", analysis_id_arg)
}

if ("weight" %in% names(weights)) {
  weights[, weight := as.numeric(weight)]
} else if ("N_ancestry" %in% names(weights)) {
  weights[, N_ancestry := as.numeric(N_ancestry)]
  weights[, weight := N_ancestry / sum(N_ancestry, na.rm=TRUE)]
} else {
  stop("Weights file must contain either weight or N_ancestry column.")
}

weights <- weights[is.finite(weight) & weight > 0]

if (nrow(weights) == 0) {
  stop("No valid LD weights.")
}

weights[, weight := weight / sum(weight)]

message("LD weights:")
print(weights[, .(ancestry, weight)])

read_ld <- function(spec) {
  parts <- strsplit(spec, ":", fixed=TRUE)[[1]]
  if (length(parts) != 3) {
    stop("LD spec must be ANC:ld_zst:vars, got: ", spec)
  }

  anc <- parts[1]
  ld_zst <- parts[2]
  vars_file <- parts[3]

  if (!file.exists(ld_zst)) {
    stop("Missing LD zst for ", anc, ": ", ld_zst)
  }

  if (!file.exists(vars_file)) {
    stop("Missing LD vars for ", anc, ": ", vars_file)
  }

  vars_raw <- fread(vars_file, header=FALSE)

  if (ncol(vars_raw) == 1) {
    vars <- as.character(vars_raw[[1]])
  } else if ("ID" %in% names(vars_raw)) {
    vars <- as.character(vars_raw$ID)
  } else {
    vars <- as.character(vars_raw[[ncol(vars_raw)]])
  }

  message("Reading LD matrix for ", anc, ": ", ld_zst)
  R <- as.matrix(fread(cmd=paste("zstdcat", shQuote(ld_zst)), header=FALSE))
  storage.mode(R) <- "double"

  if (nrow(R) != length(vars) || ncol(R) != length(vars)) {
    stop("LD dimension mismatch for ", anc)
  }

  R <- (R + t(R)) / 2
  diag(R) <- 1

  list(ancestry=anc, vars=vars, R=R)
}

ld_list <- lapply(ld_args, read_ld)
names(ld_list) <- sapply(ld_list, function(x) x$ancestry)

valid_ancestries <- intersect(names(ld_list), weights$ancestry)

if (length(valid_ancestries) == 0) {
  stop("No LD matrices overlap with valid weights.")
}

# Common variants across all ancestry LD matrices and METAL sumstats
common <- gwas$variant_id

for (anc in valid_ancestries) {
  common <- intersect(common, ld_list[[anc]]$vars)
}

if (length(common) < 2) {
  stop("Too few common variants across METAL sumstats and ancestry LD matrices.")
}

message("Original METAL SNPs: ", nrow(gwas))
message("Common SNPs after ancestry LD intersection: ", length(common))
message("Removed due to LD intersection: ", nrow(gwas) - length(common))

# Keep METAL order restricted to common variants
gwas2 <- gwas[variant_id %in% common]
gwas2 <- gwas2[match(common, variant_id)]

z <- as.numeric(gwas2$Z_REF)

if (any(!is.finite(z))) {
  stop("Non-finite z after intersection.")
}

# Re-normalize weights to available ancestries
w <- weights[ancestry %in% valid_ancestries]
w[, weight := weight / sum(weight)]

R_weighted <- NULL

for (i in seq_len(nrow(w))) {
  anc <- w$ancestry[i]
  wt <- w$weight[i]

  obj <- ld_list[[anc]]
  idx <- match(common, obj$vars)

  if (any(is.na(idx))) {
    stop("Internal error: missing variants after intersection for ", anc)
  }

  Ranc <- obj$R[idx, idx, drop=FALSE]
  storage.mode(Ranc) <- "double"

  n_nonfinite_anc <- sum(!is.finite(Ranc))
  if (n_nonfinite_anc > 0) {
    message("Replacing non-finite LD entries with 0 for ", anc, ": ", n_nonfinite_anc)
    Ranc[!is.finite(Ranc)] <- 0
  }

  Ranc <- (Ranc + t(Ranc)) / 2
  diag(Ranc) <- 1

  if (is.null(R_weighted)) {
    R_weighted <- wt * Ranc
  } else {
    R_weighted <- R_weighted + wt * Ranc
  }
}

R_weighted <- (R_weighted + t(R_weighted)) / 2

n_nonfinite_weighted <- sum(!is.finite(R_weighted))
if (n_nonfinite_weighted > 0) {
  message("Replacing non-finite weighted LD entries with 0: ", n_nonfinite_weighted)
  R_weighted[!is.finite(R_weighted)] <- 0
}

R_weighted <- (R_weighted + t(R_weighted)) / 2
diag(R_weighted) <- 1

if (length(z) < 2) {
  stop("Too few SNPs remain after weighted LD filtering.")
}

n_eff <- round(median(as.numeric(gwas2$N), na.rm=TRUE))

if (!is.finite(n_eff) || n_eff <= 0) {
  stop("Invalid n_eff from METAL N.")
}

message("Running weighted-LD SuSiE-RSS with n=", n_eff, ", L=", L, ", SNPs=", length(z))

fit <- susie_rss(
  z = z,
  R = R_weighted,
  n = n_eff,
  L = L,
  estimate_residual_variance = FALSE
)

cs <- susie_get_cs(fit, Xcorr=R_weighted, coverage=0.95)

pip_table <- data.table(
  variant_id = gwas2$variant_id,
  CHR = gwas2$CHR,
  BP = gwas2$BP,
  REF = gwas2$REF,
  ALT = gwas2$ALT,
  BETA_REF = gwas2$BETA_REF,
  Z_REF = gwas2$Z_REF,
  P = gwas2$P,
  N = gwas2$N,
  PIP = fit$pip
)

pip_table[, CS := NA_character_]

cs_table <- data.table(
  cs_id=character(),
  variant_id=character(),
  CHR=character(),
  BP=integer(),
  REF=character(),
  ALT=character(),
  PIP=numeric(),
  Z_REF=numeric(),
  P=numeric()
)

if (!is.null(cs$cs) && length(cs$cs) > 0) {
  for (cs_name in names(cs$cs)) {
    idx <- cs$cs[[cs_name]]
    pip_table[idx, CS := cs_name]

    tmp <- data.table(
      cs_id = cs_name,
      variant_id = gwas2$variant_id[idx],
      CHR = gwas2$CHR[idx],
      BP = gwas2$BP[idx],
      REF = gwas2$REF[idx],
      ALT = gwas2$ALT[idx],
      PIP = fit$pip[idx],
      Z_REF = gwas2$Z_REF[idx],
      P = gwas2$P[idx]
    )

    cs_table <- rbind(cs_table, tmp)
  }
}

summary_table <- data.table(
  out_prefix = out_prefix,
  analysis_id = analysis_id_arg,
  n_snps = length(z),
  n = n_eff,
  L = L,
  n_credible_sets = ifelse(is.null(cs$cs), 0, length(cs$cs)),
  max_pip = max(fit$pip, na.rm=TRUE),
  top_variant = gwas2$variant_id[which.max(fit$pip)]
)

weight_out <- copy(w)
weight_out[, out_prefix := out_prefix]

fwrite(pip_table, paste0(out_prefix, ".pip.tsv"), sep="\t")
fwrite(cs_table, paste0(out_prefix, ".credible_sets.tsv"), sep="\t")
fwrite(summary_table, paste0(out_prefix, ".summary.tsv"), sep="\t")
fwrite(weight_out, paste0(out_prefix, ".ld_weights_used.tsv"), sep="\t")

if (!is.null(cs$purity)) {
  purity_table <- as.data.table(cs$purity)
  purity_table[, cs_id := rownames(cs$purity)]
  fwrite(purity_table, paste0(out_prefix, ".credible_set_purity.tsv"), sep="\t")
}

saveRDS(fit, paste0(out_prefix, ".susie_fit.rds"))
saveRDS(cs, paste0(out_prefix, ".susie_cs.rds"))

message("Weighted-LD SuSiE completed successfully.")
