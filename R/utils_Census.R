
# census mode function
#' @importFrom stats density
.dmode <- function(x, breaks="Sturges") {
  if (length(x) < 2) return (0);
  den <- stats::density(x, kernel=c("gaussian"), na.rm = T)
  ( den$x[den$y==max(den$y)] )
}

# census estimate_t function
.estimate_t <- function(relative_expr_matrix,
                        relative_expr_thresh) {
  #apply each column
  unlist(apply(relative_expr_matrix, 2, function(relative_expr)
    10^mean(.dmode(log10(relative_expr[relative_expr > relative_expr_thresh])))))
}


# census estimatesizefactors function
.estimateSizeFactorsForDenseMatrix <- function(CM,
                                               locfunc,
                                               round_exprs,
                                               method){

  if (round_exprs)
    CM <- round(CM)
  if (method == "weighted-median"){
    log_medians <- apply(CM, 1, function(cell_expr) {
      log(locfunc(cell_expr))
    })

    weights <- apply(CM, 1, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })

    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <-  weights * (log(cnts) -  log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      exp( mean(norm_cnts) )
    })
  }else if (method == "median-geometric-mean"){
    log_geo_means <- rowMeans(log(CM))

    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <- log(cnts) -  log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      exp( locfunc( norm_cnts ))
    })
  }else if(method == "median"){
    row_median <- apply(CM, 1, median)
    sfs <- apply(t(t(CM) - row_median), 2, median)
  }else if(method == 'mode'){
    sfs <- .estimate_t(CM)
  }else if(method == 'geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- log(cell_total) / mean(log(cell_total))
  }else if(method == 'mean-geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }

  sfs[is.na(sfs)] <- 1
  return(sfs)
}

#' @importFrom stats ecdf
.calibrate_per_cell_total_proposal <- function(relative_exprs_matrix,
                                               t_estimate,
                                               expected_capture_rate,
                                               method = c('num_genes', 'tpm_fraction') ){
  split_relative_exprs <- split(relative_exprs_matrix,
                                rep(1:ncol(relative_exprs_matrix),
                                    each = nrow(relative_exprs_matrix)))

  proposed_totals <- unlist(lapply(1:length(split_relative_exprs), function(ind) {
    x <- split_relative_exprs[[ind]];
    x <- x[x > 0.1];
    if(method == 'num_genes') {
      P <- ecdf(x);
      frac_x <- P(t_estimate[ind]);
    }
    else if(method == 'tpm_fraction') {
      frac_x <- sum(x[x < t_estimate[ind]]) / sum(x)
    }
    num_single_copy_genes <- sum(x <= t_estimate[ind]);
    num_single_copy_genes / frac_x / expected_capture_rate
  }))
  return(proposed_totals)
}

#' @importFrom MASS rlm
#' @importFrom stats predict
#' @importFrom parallel mcmapply
.relative2abs <- function(relative_expr_matrix,
                          t_estimate,
                          relative_spike, # relative exprs matrix of spike-ins
                          spikeInfo, # spike input
                          verbose,
                          cores) {

  if(is.null(relative_spike)) {
    names(t_estimate) <- colnames(relative_expr_matrix)
    expected_total_mRNAs <- .calibrate_per_cell_total_proposal(relative_exprs_matrix = relative_expr_matrix,
                                                               t_estimate = t_estimate,
                                                               expected_capture_rate = 0.25,
                                                               method = "num_genes")

    expr_probs <-  t(t(relative_expr_matrix)/ colSums(relative_expr_matrix))
    census_transcript_counts <- t(t(expr_probs) * expected_total_mRNAs)

    return(list(norm_cds = census_transcript_counts,
                t_estimate=t_estimate,
                expected_total_mRNAs=expected_total_mRNAs))

  }

  if(!is.null(relative_spike)) {
    FPKM <- NULL

    # relative expression of spike-ins
    ERCC_controls <- relative_spike

    # spike input information
    ERCC_annotation <- spikeInfo
    valid_ids <- which(ERCC_annotation[, "SpikeInput"] >= 0)

    # robust linear regression
    if (verbose) {message("Performing robust linear regression for each cell based
                        on the spike-in data provided.")}
    molModels <- apply(ERCC_controls, 2, function(cell_exprs, input.ERCC.annotation, valid_ids) {
      spike_df <- input.ERCC.annotation
      spike_df <- cbind(spike_df, cell_exprs[row.names(spike_df)])
      colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
      spike_df$numMolecules <- spike_df$SpikeInput
      spike_df$rounded_numMolecules <- round(spike_df$numMolecules)
      if (is.null(valid_ids))
        spike_df <- subset(spike_df, FPKM >= 1e-10)
      else {
        spike_df <- spike_df[valid_ids, ]
        spike_df <- subset(spike_df, FPKM >= 1e-10)
      }
      spike_df$log_fpkm <- log10(spike_df$FPKM)
      spike_df$log_numMolecules <- log10(spike_df$numMolecules)
      molModel <- tryCatch({
        molModel <- MASS::rlm(log_numMolecules ~ log_fpkm,
                              data = spike_df)
        molModel
      }, error = function(e) {
        print(e)
        NULL
      })
      molModel
    }, ERCC_annotation, valid_ids)

    if (verbose) {message("Apply the fitted robust linear regression model
                        to recover the absolute copy number for all transcripts in each cell")}
    norm_fpkms <- parallel::mcmapply(function(cell_exprs, molModel) {
      tryCatch({
        norm_df <- data.frame(log_fpkm = log10(cell_exprs))
        res <- 10^stats::predict(molModel, type = "response",
                                 newdata = norm_df)
      }, error = function(e) {
        rep(NA, length(cell_exprs))
      })
    }, split(as.matrix(relative_expr_matrix), rep(1:ncol(relative_expr_matrix),
                                                  each = nrow(relative_expr_matrix))), molModels, mc.cores = cores)
    k_b_solution <- data.frame(b = unlist(lapply(molModels,
                                                 FUN = function(x) {
                                                   intercept = x$coefficients[1]
                                                 })), k = unlist(lapply(molModels, FUN = function(x) {
                                                   slope = x$coefficients[2]
                                                 })))
    kb_model <- MASS::rlm(b ~ k, data = k_b_solution)
    kb_slope <- kb_model$coefficients[2]
    kb_intercept <- kb_model$coefficients[1]

    rownames(norm_fpkms) <- rownames(relative_expr_matrix)
    colnames(norm_fpkms) <- colnames(relative_expr_matrix)

    # return object
    return(list(norm_cds = norm_fpkms,
                kb_slope = kb_slope,
                kb_intercept = kb_intercept,
                k_b_solution = k_b_solution))
  }

}
