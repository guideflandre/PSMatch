##' @title Extract feature columns from a PSM object
##'
##' @description
##' Returns the subset of columns in a [PSM] object that are not
##' registered as PSM variables (i.e. spectrum, peptide, protein,
##' decoy, rank, score, fdr). All remaining columns are treated as
##' numerical input features for the SVM.
##'
##' @param psm A [PSM] object with PSM variables set via [PSM()] or
##'     [psmVariables()].
##'
##' @return A `data.frame` with one row per PSM and one column per
##'     feature. Throws an error if no feature columns remain after
##'     removing PSM variables.
##'
##' @keywords internal
##'
##' @noRd
extract_features <- function(psm) {
    fixedCols <- psmVariables(psm)
    fixedCols <- unname(fixedCols[!is.na(fixedCols)])
    featureCols <- setdiff(colnames(psm), fixedCols)

    if (length(featureCols) == 0L)
        stop(
            "No feature columns found after excluding PSM variables. ",
            "Ensure the PSM object contains additional numeric columns ",
            "beyond those registered in psmVariables()."
        )

    as.data.frame(psm[, featureCols, drop = FALSE])
}


##' @title Prepare SVM training data from a PSM object
##'
##' @description
##' Extracts feature columns from a [PSM] object via
##' [extract_features()], removes zero-variance features (which are
##' uninformative), scales the retained features to zero mean and unit
##' variance, and encodes target/decoy labels as `1` / `-1`.
##'
##' The decoy column is resolved through [psmVariables()] so no
##' hard-coded column name is required.
##'
##' @param psm A [PSM] object. The decoy column is resolved via
##'     [psmVariables()].
##'
##' @return A `list` with three elements:
##' \describe{
##'   \item{`features`}{Scaled numeric matrix (rows = PSMs,
##'       columns = retained features).}
##'   \item{`labels`}{`factor` with levels `c(-1L, 1L)`, where
##'       `-1` encodes decoy and `1` encodes target PSMs.}
##'   \item{`scaleAttrs`}{`list` containing `scaled:center` and
##'       `scaled:scale` vectors, used to rescale new data
##'       consistently during prediction.}
##' }
##'
##' @keywords internal
##'
##' @noRd
prepare_svm_data <- function(psm) {
    features <- extract_features(psm)

    ## Remove zero-variance columns (uninformative for the SVM)
    colVar <- vapply(features, var, numeric(1L), na.rm = TRUE)
    features <- features[, colVar > 0, drop = FALSE]

    ## Scale to zero mean and unit variance
    scaled <- scale(features)

    ## Encode labels using the decoy column from PSM variables
    decoyCol <- psmVariables(psm)["decoy"]
    labels <- factor(
        ifelse(psm[[decoyCol]], -1L, 1L),
        levels = c(-1L, 1L)
    )

    list(
        features   = scaled,
        labels     = labels,
        scaleAttrs = attributes(scaled)
    )
}


##' @title Train an SVM classifier on prepared PSM data
##'
##' @description
##' Fits a Support Vector Machine (via `e1071`) to distinguish target
##' from decoy PSMs. The model is trained with probability estimation
##' enabled so that [score_psms()] can extract posterior target
##' probabilities as rescored PSM scores.
##'
##' Features must already be scaled (as returned by
##' [prepare_svm_data()]) because `scale = FALSE` is passed
##' internally.
##'
##' @param prepared A `list` as returned by [prepare_svm_data()],
##'     containing `features` and `labels`.
##'
##' @param kernel `character(1)` SVM kernel type passed to
##'     [e1071::svm()]. One of `"radial"` (default), `"linear"`,
##'     `"polynomial"`, or `"sigmoid"`.
##'
##' @param cost `numeric(1)` regularisation parameter `C` (cost of
##'     constraint violation). Default is `1`. Larger values impose
##'     a harder margin and may overfit.
##'
##' @param gamma `numeric(1)` or `NULL`. Kernel coefficient for
##'     `"radial"`, `"polynomial"`, and `"sigmoid"` kernels. When
##'     `NULL` (default), `1 / ncol(features)` is used.
##'
##' @return An object of class `svm` (see [e1071::svm()]).
##'
##' @keywords internal
##'
##' @noRd
train_svm <- function(prepared,
                      kernel = "radial",
                      cost   = 1,
                      gamma  = NULL) {
    features <- prepared$features
    labels   <- prepared$labels

    if (is.null(gamma))
        gamma <- 1 / ncol(features)

    e1071::svm(
        x           = features,
        y           = labels,
        kernel      = kernel,
        cost        = cost,
        gamma       = gamma,
        type        = "C-classification",
        probability = TRUE,
        scale       = FALSE  # already scaled in prepare_svm_data()
    )
}


##' @title Partition PSM indices into cross-validation folds
##'
##' @description
##' Randomly assigns each of `n` PSMs to one of `nFolds` folds,
##' returning a `list` of integer index vectors. Each vector contains
##' the row indices of the PSMs assigned to that fold.
##'
##' The assignment ensures that fold sizes differ by at most one row
##' (balanced folds) via `rep(..., length.out = n)`.
##'
##' @param n `integer(1)` total number of PSMs to partition.
##'
##' @param nFolds `integer(1)` number of folds. Default is `10L`.
##'
##' @return A named `list` of length `nFolds`. Each element is an
##'     `integer` vector of row indices belonging to that fold.
##'
##' @keywords internal
##'
##' @noRd
create_cv_folds <- function(n, nFolds = 10L) {
    foldAssignment <- sample(rep(seq_len(nFolds), length.out = n))
    split(seq_len(n), foldAssignment)
}


##' @title Score PSMs using a trained SVM model
##'
##' @description
##' Applies a trained SVM model to a [PSM] object and returns the
##' posterior probability of each PSM being a target. Feature columns
##' are extracted with [extract_features()] and rescaled using the
##' centre/scale parameters stored from training so that held-out
##' data are treated identically to the training data.
##'
##' This function returns raw scores only; use [assign_svm_results()]
##' to attach scores and ranks back to the [PSM] object.
##'
##' @param psm A [PSM] object to score (may be a held-out fold or
##'     the full dataset).
##'
##' @param svmModel An `svm` object as returned by [train_svm()].
##'
##' @param scaleAttrs A `list` with `scaled:center` and `scaled:scale`
##'     vectors, as stored in the `scaleAttrs` element returned by
##'     [prepare_svm_data()].
##'
##' @return A `numeric` vector of length `nrow(psm)` containing the
##'     posterior probability of being a target PSM for each row.
##'
##' @keywords internal
##'
##' @noRd
score_psms <- function(psm, svmModel, scaleAttrs) {
    features <- extract_features(psm)

    ## Retain only the columns the model was trained on
    trainedCols <- colnames(svmModel$SV)
    features <- features[, trainedCols, drop = FALSE]

    ## Apply identical scaling as used during training
    featuresScaled <- scale(
        features,
        center = scaleAttrs$`scaled:center`,
        scale  = scaleAttrs$`scaled:scale`
    )

    preds <- predict(svmModel, newdata = featuresScaled,
                     probability = TRUE)
    probs <- attr(preds, "probabilities")
    probs[, "1"]
}


##' @title Attach SVM scores and within-peptide ranks to a PSM object
##'
##' @description
##' Adds a `svmScore` column (posterior target probability) and a
##' `svmRank` column (within-peptide rank, 1 = best) to a [PSM]
##' object. PSMs are reordered by peptide group and descending score
##' before ranks are assigned.
##'
##' The peptide grouping column is resolved via [psmVariables()].
##'
##' @param psm A [PSM] object whose row order matches `svmScores`.
##'
##' @param svmScores `numeric` vector of length `nrow(psm)` containing
##'     the posterior target probability for each PSM, as returned by
##'     [score_psms()].
##'
##' @return The input [PSM] object with two additional columns:
##' \describe{
##'   \item{`svmScore`}{`numeric` posterior probability of being a
##'       target PSM (higher = more confident target).}
##'   \item{`svmRank`}{`integer` rank within each peptide group
##'       (1 = highest `svmScore`).}
##' }
##'
##' @keywords internal
##'
##' @noRd
assign_svm_results <- function(psm, svmScores) {
    psm[["svmScore"]] <- svmScores

    ## Re-rank PSMs within each peptide group by descending SVM score
    peptideCol <- psmVariables(psm)["peptide"]
    peptideVec <- as.character(psm[[peptideCol]])
    ord <- order(peptideVec, -psm[["svmScore"]])
    psm <- psm[ord, ]

    psm[["svmRank"]] <- as.integer(ave(
        psm[["svmScore"]],
        as.character(psm[[peptideCol]]),
        FUN = seq_along
    ))

    psm
}


##' @title Compute a ROC curve from PSM scores and decoy labels
##'
##' @description
##' Sorts PSMs by descending score and computes the cumulative true
##' positive rate (TPR) and false positive rate (FPR) at each
##' threshold, anchored at (0, 0) and ending at (1, 1).
##'
##' Targets are treated as positives, decoys as negatives.
##'
##' @param scores `numeric` vector of PSM scores (higher = target).
##'
##' @param isDecoy `logical` vector of the same length as `scores`;
##'     `TRUE` for decoy PSMs.
##'
##' @return A `data.frame` with columns `fpr` (false positive rate)
##'     and `tpr` (true positive rate).
##'
##' @keywords internal
##'
##' @noRd
compute_roc <- function(scores, isDecoy) {
    ord     <- order(scores, decreasing = TRUE)
    isDecoy <- isDecoy[ord]
    nTarget <- sum(!isDecoy)
    nDecoy  <- sum(isDecoy)
    tpr     <- c(0, cumsum(!isDecoy) / nTarget)
    fpr     <- c(0, cumsum(isDecoy)  / nDecoy)
    data.frame(fpr = fpr, tpr = tpr)
}


##' @title Compute AUC from a ROC curve
##'
##' @description
##' Applies the trapezoidal rule to the (fpr, tpr) pairs produced
##' by [compute_roc()]. `abs(diff(fpr))` handles any non-monotone
##' segments safely.
##'
##' @param fpr `numeric` false positive rate vector.
##'
##' @param tpr `numeric` true positive rate vector (same length as
##'     `fpr`).
##'
##' @return `numeric(1)` area under the ROC curve in `[0, 1]`.
##'
##' @keywords internal
##'
##' @noRd
compute_auc <- function(fpr, tpr) {
    sum(abs(diff(fpr)) * (tpr[-length(tpr)] + tpr[-1L]) / 2)
}


##' @title Compute target-decoy FDR along a score threshold
##'
##' @description
##' Sorts PSMs by descending score and computes the running
##' target-decoy FDR at each position as
##' \eqn{FDR(t) = n_{\text{decoy}}(t) \,/\, n_{\text{target}}(t)}.
##'
##' @param scores `numeric` vector of PSM scores.
##'
##' @param isDecoy `logical` vector of the same length; `TRUE` for
##'     decoy PSMs.
##'
##' @return A `data.frame` with four columns: `score` (threshold),
##'     `nTarget` (cumulative targets), `nDecoy` (cumulative decoys),
##'     and `fdr`.
##'
##' @keywords internal
##'
##' @noRd
compute_td_fdr <- function(scores, isDecoy) {
    ord     <- order(scores, decreasing = TRUE)
    scores  <- scores[ord]
    isDecoy <- isDecoy[ord]
    nTarget <- cumsum(!isDecoy)
    nDecoy  <- cumsum(isDecoy)
    fdr     <- nDecoy / pmax(nTarget, 1L)
    data.frame(
        score   = scores,
        nTarget = nTarget,
        nDecoy  = nDecoy,
        fdr     = fdr
    )
}


##' @title Evaluate SVM rescoring performance
##'
##' @description
##' Computes performance metrics comparing the original PSM score
##' (from [psmVariables()]) with the cross-validated SVM score
##' produced by [run_svm_rescoring()]. Metrics include:
##'
##' - AUC of the target-decoy ROC curve (targets = positives,
##'   decoys = negatives).
##' - Number of target PSMs accepted at 1 % and 5 % FDR, estimated
##'   by the target-decoy approach:
##'   \eqn{FDR(t) = n_{\text{decoy}}(t) / n_{\text{target}}(t)}.
##' - Full ROC curve `data.frame` (columns `fpr`, `tpr`).
##' - Full FDR curve `data.frame` (columns `score`, `nTarget`,
##'   `nDecoy`, `fdr`).
##'
##' If the `score` PSM variable is `NA` (not set), only SVM metrics
##' are returned. If it is set, an `original` sub-list is added for
##' comparison.
##'
##' @param psm A rescored [PSM] object, i.e. the `rescored` element
##'     of the `list` returned by [run_svm_rescoring()]. Must contain
##'     an `svmScore` column.
##'
##' @return A `list` with one or two named sub-lists (`svm` and,
##'     optionally, `original`). Each contains:
##' \describe{
##'   \item{`auc`}{`numeric(1)` area under the ROC curve.}
##'   \item{`nTargetFdr1`}{`integer(1)` number of target PSMs
##'       accepted at 1 % FDR.}
##'   \item{`nTargetFdr5`}{`integer(1)` number of target PSMs
##'       accepted at 5 % FDR.}
##'   \item{`roc`}{`data.frame` with columns `fpr` and `tpr`.}
##'   \item{`fdrDf`}{`data.frame` with columns `score`, `nTarget`,
##'       `nDecoy`, and `fdr`.}
##' }
##'
##' @examples
##' set.seed(42)
##' n <- 120L
##' psmDf <- data.frame(
##'     spectrumID      = paste0("sp", seq_len(n)),
##'     sequence        = sample(
##'         c("PEPTIDEA", "PEPTIDEB", "PEPTIDEC"),
##'         n, replace = TRUE),
##'     protein         = sample(
##'         paste0("Prot", LETTERS[seq_len(5)]),
##'         n, replace = TRUE),
##'     isDecoy         = sample(
##'         c(FALSE, TRUE), n,
##'         replace = TRUE, prob = c(0.7, 0.3)),
##'     rank            = rep(1L, n),
##'     score           = runif(n),
##'     massErrorPpm    = rnorm(n, 0, 5),
##'     missedCleavages = sample(0:2, n, replace = TRUE)
##' )
##' psm <- PSM(psmDf,
##'            spectrum = "spectrumID", peptide  = "sequence",
##'            protein  = "protein",   decoy    = "isDecoy",
##'            rank     = "rank",      score    = "score")
##'
##' result  <- run_svm_rescoring(psm)
##' metrics <- evaluate_rescoring(result$rescored)
##'
##' ## AUC comparison (SVM vs original score)
##' metrics$svm$auc
##' metrics$original$auc
##'
##' ## Target PSMs accepted at 1 % FDR
##' metrics$svm$nTargetFdr1
##' metrics$original$nTargetFdr1
##'
##' ## Plot ROC curves for both scores
##' plot(metrics$original$roc$fpr, metrics$original$roc$tpr,
##'      type = "l", col = "grey50", lwd = 2,
##'      xlab = "FPR", ylab = "TPR",
##'      main = "ROC: original vs SVM score")
##' lines(metrics$svm$roc$fpr, metrics$svm$roc$tpr,
##'       col = "#2166AC", lwd = 2)
##' legend("bottomright",
##'        legend = c("Original", "SVM"),
##'        col    = c("grey50", "#2166AC"),
##'        lwd    = 2, bty = "n")
##'
##' @export
evaluate_rescoring <- function(psm) {
    decoyCol  <- psmVariables(psm)["decoy"]
    scoreCol  <- psmVariables(psm)["score"]
    isDecoy   <- as.logical(psm[[decoyCol]])
    svmScores <- psm[["svmScore"]]

    if (is.null(svmScores))
        stop(
            "No 'svmScore' column found. ",
            "Run run_svm_rescoring() before evaluate_rescoring()."
        )

    svmRoc  <- compute_roc(svmScores, isDecoy)
    svmFdr  <- compute_td_fdr(svmScores, isDecoy)
    svmAuc  <- compute_auc(svmRoc$fpr, svmRoc$tpr)
    idx1    <- svmFdr$fdr <= 0.01
    idx5    <- svmFdr$fdr <= 0.05

    result <- list(
        svm = list(
            auc         = svmAuc,
            nTargetFdr1 = if (any(idx1)) max(svmFdr$nTarget[idx1]) else 0L,
            nTargetFdr5 = if (any(idx5)) max(svmFdr$nTarget[idx5]) else 0L,
            roc         = svmRoc,
            fdrDf       = svmFdr
        )
    )

    if (!is.na(scoreCol) && scoreCol %in% colnames(psm)) {
        origScores <- as.numeric(psm[[scoreCol]])
        origRoc    <- compute_roc(origScores, isDecoy)
        origFdr    <- compute_td_fdr(origScores, isDecoy)
        origAuc    <- compute_auc(origRoc$fpr, origRoc$tpr)
        oidx1      <- origFdr$fdr <= 0.01
        oidx5      <- origFdr$fdr <= 0.05

        result$original <- list(
            auc         = origAuc,
            nTargetFdr1 = if (any(oidx1)) max(origFdr$nTarget[oidx1])
                          else 0L,
            nTargetFdr5 = if (any(oidx5)) max(origFdr$nTarget[oidx5])
                          else 0L,
            roc         = origRoc,
            fdrDf       = origFdr
        )
    }

    result
}


##' @title Visualise the SVM decision boundary via PCA projection
##'
##' @description
##' Projects the scaled feature space onto the first two principal
##' components (PCA) and renders the SVM decision boundary in that
##' 2D plane using a dense prediction grid. PSMs are overlaid as
##' coloured points (target = blue circles, decoy = red triangles).
##'
##' **Interpretation note.** Because the boundary lives in a
##' high-dimensional feature space, the 2D projection is an
##' approximation: variation along PC3 and beyond is collapsed to
##' zero when reconstructing grid points. The plot is most reliable
##' when PC1 and PC2 together explain a large proportion of variance.
##' The percentage of variance explained is shown on each axis label.
##'
##' Grid reconstruction follows these steps:
##' 1. Scale features with `scaleAttrs` (training parameters).
##' 2. Run PCA on the scaled features.
##' 3. For each grid point (pc1, pc2): set all other PCs to 0,
##'    back-project to the scaled feature space via the PCA rotation.
##' 4. Predict the SVM class on the reconstructed grid.
##'
##' @param psm A [PSM] object. Used to extract features and decoy
##'     labels via [psmVariables()].
##'
##' @param svmModel The final `svm` object returned in `result$model`
##'     by [run_svm_rescoring()].
##'
##' @param scaleAttrs The `scaleAttrs` element returned by
##'     [run_svm_rescoring()], containing the `scaled:center` and
##'     `scaled:scale` vectors from the full-data training run.
##'
##' @param nGrid `integer(1)` number of grid points along each axis.
##'     Higher values give a smoother boundary at the cost of more
##'     SVM predictions. Default is `80L`.
##'
##' @return Invisibly returns a `list` with:
##' \describe{
##'   \item{`pca`}{The `prcomp` object from the PCA step.}
##'   \item{`scores`}{`matrix` of PC1/PC2 scores for each PSM.}
##' }
##' The function is called primarily for its side effect (the plot).
##'
##' @examples
##' set.seed(42)
##' n <- 120L
##' psmDf <- data.frame(
##'     spectrumID      = paste0("sp", seq_len(n)),
##'     sequence        = sample(
##'         c("PEPTIDEA", "PEPTIDEB", "PEPTIDEC"),
##'         n, replace = TRUE),
##'     protein         = sample(
##'         paste0("Prot", LETTERS[seq_len(5)]),
##'         n, replace = TRUE),
##'     isDecoy         = sample(
##'         c(FALSE, TRUE), n,
##'         replace = TRUE, prob = c(0.7, 0.3)),
##'     rank            = rep(1L, n),
##'     score           = runif(n),
##'     massErrorPpm    = rnorm(n, 0, 5),
##'     missedCleavages = sample(0:2, n, replace = TRUE)
##' )
##' psm <- PSM(psmDf,
##'            spectrum = "spectrumID", peptide  = "sequence",
##'            protein  = "protein",   decoy    = "isDecoy",
##'            rank     = "rank",      score    = "score")
##'
##' result <- run_svm_rescoring(psm)
##'
##' ## Plot the decision boundary using the final full-data model
##' plot_decision_boundary(psm, result$model, result$scaleAttrs)
##'
##' @export
plot_decision_boundary <- function(psm,
                                   svmModel,
                                   scaleAttrs,
                                   nGrid = 80L) {
    features <- extract_features(psm)

    ## Retain only the columns the model was trained on
    trainedCols <- colnames(svmModel$SV)
    features <- features[, trainedCols, drop = FALSE]

    if (ncol(features) < 2L)
        stop(
            "At least 2 feature columns are required ",
            "to plot a 2D decision boundary."
        )

    ## Apply training scaling to all PSMs
    featScaled <- scale(
        features,
        center = scaleAttrs$`scaled:center`,
        scale  = scaleAttrs$`scaled:scale`
    )

    ## PCA on the scaled features (centre within PCA space)
    pca    <- prcomp(featScaled, center = TRUE, scale. = FALSE)
    scores <- pca$x[, 1:2, drop = FALSE]

    ## Build a dense grid in the PC1 / PC2 plane
    r1   <- range(scores[, 1L])
    r2   <- range(scores[, 2L])
    buf1 <- diff(r1) * 0.1
    buf2 <- diff(r2) * 0.1
    pc1Grid <- seq(r1[1L] - buf1, r1[2L] + buf1, length.out = nGrid)
    pc2Grid <- seq(r2[1L] - buf2, r2[2L] + buf2, length.out = nGrid)
    gridDf  <- expand.grid(PC1 = pc1Grid, PC2 = pc2Grid)

    ## Reconstruct feature vectors: other PCs fixed at 0
    nPc        <- ncol(pca$rotation)
    gridPcFull <- matrix(0, nrow = nrow(gridDf), ncol = nPc)
    gridPcFull[, 1L] <- gridDf$PC1
    gridPcFull[, 2L] <- gridDf$PC2

    ## Back-project to scaled feature space and restore PCA centre
    gridFeat <- sweep(
        gridPcFull %*% t(pca$rotation),
        MARGIN = 2L, STATS = pca$center, FUN = "+"
    )

    ## Predict SVM class on the grid (-1 = decoy, 1 = target)
    gridPreds <- predict(svmModel, newdata = gridFeat)
    zMat <- matrix(
        as.integer(as.character(gridPreds)),
        nrow = nGrid, ncol = nGrid
    )

    ## Variance explained for axis labels
    varExp <- pca$sdev^2 / sum(pca$sdev^2)
    xlab   <- sprintf("PC1 (%.1f %%)", varExp[1L] * 100)
    ylab   <- sprintf("PC2 (%.1f %%)", varExp[2L] * 100)

    ## Background: light coral = decoy region, light blue = target
    image(pc1Grid, pc2Grid, zMat,
          col  = c("#FDDBC7", "#D1E5F0"),
          xlab = xlab, ylab = ylab,
          main = "SVM decision boundary (PC1 / PC2 projection)")

    ## Decision boundary contour
    contour(pc1Grid, pc2Grid, zMat,
            levels     = 0,
            add        = TRUE,
            lwd        = 2,
            drawlabels = FALSE)

    ## Overlay actual PSMs coloured by target / decoy status
    decoyCol <- psmVariables(psm)["decoy"]
    isDecoy  <- as.logical(psm[[decoyCol]])

    points(scores[!isDecoy, 1L], scores[!isDecoy, 2L],
           pch = 16L, col = "#2166ACAA", cex = 0.8)
    points(scores[isDecoy, 1L],  scores[isDecoy, 2L],
           pch = 17L, col = "#D6604DAA", cex = 0.8)

    legend("topright",
           legend = c("Target", "Decoy"),
           pch    = c(16L, 17L),
           col    = c("#2166AC", "#D6604D"),
           bty    = "n")

    invisible(list(pca = pca, scores = scores))
}


##' @title Re-score peptide-spectrum matches using 10-fold CV SVM
##'
##' @description
##' Top-level wrapper that rescores all PSMs in a [PSM] object using
##' a Support Vector Machine trained via k-fold cross-validation.
##'
##' **Why cross-validation?** Training the SVM and scoring the same
##' PSMs inflates scores for PSMs seen during training (overfitting).
##' In k-fold CV, each PSM is scored by a model trained on the
##' *other* k-1 folds, so every PSM receives an unbiased score.
##' This mirrors the approach used by Percolator.
##'
##' The column names used at each step (decoy flag, peptide grouping)
##' are resolved from the PSM variables stored in the [PSM] object
##' via [psmVariables()], so no hard-coded column names are assumed.
##'
##' The workflow proceeds as follows:
##' 1. Partition PSMs into `nFolds` balanced folds.
##' 2. For each fold k: train on the remaining k-1 folds; score
##'    fold k with the held-out model.
##' 3. Collect all cross-validated scores across folds.
##' 4. Attach `svmScore` and `svmRank` columns to the [PSM] object.
##' 5. Train a final model on **all** data (returned as `model`) for
##'    use on future, unseen PSMs.
##'
##' @param psm A [PSM] object. PSM variables (spectrum, peptide,
##'     protein, decoy, rank) must be set via [PSM()]; all remaining
##'     numeric columns are used as SVM features.
##'
##' @param nFolds `integer(1)` number of cross-validation folds.
##'     Default is `10L`.
##'
##' @param kernel `character(1)` SVM kernel. One of `"radial"`
##'     (default), `"linear"`, `"polynomial"`, or `"sigmoid"`.
##'     Passed to [e1071::svm()].
##'
##' @param cost `numeric(1)` regularisation parameter `C`. Default
##'     is `1`. Larger values impose a harder margin and may overfit.
##'
##' @param gamma `numeric(1)` or `NULL`. Kernel coefficient for
##'     radial, polynomial and sigmoid kernels. When `NULL` (default),
##'     `1 / number_of_features` is used.
##'
##' @return A `list` with three elements:
##' \describe{
##'   \item{`model`}{Final `svm` object trained on all PSMs
##'       (see [e1071::svm()]). Use this to score new, unseen PSMs.}
##'   \item{`rescored`}{The input [PSM] object with `svmScore`
##'       (`numeric`, cross-validated posterior target probability)
##'       and `svmRank` (`integer`, within-peptide rank) appended.
##'       Scores are unbiased because each PSM was scored by a model
##'       that did not train on it.}
##'   \item{`scaleAttrs`}{`list` of `scaled:center` and `scaled:scale`
##'       vectors from the final full-data training run. Pass this
##'       to [plot_decision_boundary()] or [score_psms()] when
##'       applying the model to new data.}
##' }
##'
##' @examples
##' ## ----------------------------------------------------------
##' ## Build a small synthetic PSM object with two feature columns
##' ## ----------------------------------------------------------
##' set.seed(42)
##' n <- 120L
##' psmDf <- data.frame(
##'     spectrumID      = paste0("sp", seq_len(n)),
##'     sequence        = sample(
##'         c("PEPTIDEA", "PEPTIDEB", "PEPTIDEC"),
##'         n, replace = TRUE
##'     ),
##'     protein         = sample(
##'         paste0("Prot", LETTERS[seq_len(5)]),
##'         n, replace = TRUE
##'     ),
##'     isDecoy         = sample(
##'         c(FALSE, TRUE), n,
##'         replace = TRUE, prob = c(0.7, 0.3)
##'     ),
##'     rank            = rep(1L, n),
##'     score           = runif(n),
##'     ## --- feature columns used by the SVM ---
##'     massErrorPpm    = rnorm(n, mean = 0, sd = 5),
##'     missedCleavages = sample(0:2, n, replace = TRUE)
##' )
##'
##' ## Wrap in a PSM object and register PSM variables
##' psm <- PSM(psmDf,
##'            spectrum = "spectrumID",
##'            peptide  = "sequence",
##'            protein  = "protein",
##'            decoy    = "isDecoy",
##'            rank     = "rank",
##'            score    = "score")
##' psm
##' psmVariables(psm)
##'
##' ## Run 10-fold CV SVM rescoring (radial kernel, default settings)
##' result <- run_svm_rescoring(psm)
##'
##' ## Inspect the rescored PSM object
##' ## svmScore is cross-validated: each PSM was scored by a model
##' ## trained on the other 9 folds, avoiding score inflation.
##' rescored <- result$rescored
##' rescored[, c("sequence", "isDecoy", "score",
##'              "svmScore", "svmRank")]
##'
##' ## Retrieve the best-ranked target PSM per peptide
##' decoyCol <- psmVariables(rescored)["decoy"]
##' topHits  <- rescored[
##'     !rescored[[decoyCol]] & rescored[["svmRank"]] == 1L,
##' ]
##' topHits  <- topHits[order(-topHits[["svmScore"]]), ]
##' topHits[, c("sequence", "protein", "score", "svmScore")]
##'
##' ## Evaluate rescoring performance
##' metrics <- evaluate_rescoring(rescored)
##' metrics$svm$auc
##' metrics$original$auc
##'
##' ## Visualise the SVM decision boundary
##' plot_decision_boundary(psm, result$model, result$scaleAttrs)
##'
##' ## Tune with a linear kernel and higher cost, 5 folds
##' result2 <- run_svm_rescoring(psm, nFolds = 5L,
##'                              kernel = "linear", cost = 10)
##'
##' @importFrom e1071 svm
##'
##' @export
run_svm_rescoring <- function(psm,
                              nFolds = 10L,
                              kernel = "radial",
                              cost   = 1,
                              gamma  = NULL) {
    n <- nrow(psm)
    folds <- create_cv_folds(n, nFolds)

    ## Pre-allocate score vector; order matches original PSM row order
    svmScores <- numeric(n)

    message(sprintf(
        "Running %d-fold cross-validation on %d PSMs...",
        nFolds, n
    ))

    for (k in seq_len(nFolds)) {
        testIdx  <- folds[[k]]
        trainIdx <- unlist(folds[-k], use.names = FALSE)

        prepared <- prepare_svm_data(psm[trainIdx, ])
        foldModel <- train_svm(
            prepared,
            kernel = kernel,
            cost   = cost,
            gamma  = gamma
        )

        svmScores[testIdx] <- score_psms(
            psm[testIdx, ],
            foldModel,
            prepared$scaleAttrs
        )

        message(sprintf("  Fold %d / %d complete.", k, nFolds))
    }

    ## Train a final model on all data for scoring future PSMs
    message("Training final model on all data...")
    preparedAll <- prepare_svm_data(psm)
    model <- train_svm(
        preparedAll,
        kernel = kernel,
        cost   = cost,
        gamma  = gamma
    )

    message("Assigning cross-validated scores and ranks...")
    rescored <- assign_svm_results(psm, svmScores)

    message("Done.")
    list(
        model      = model,
        rescored   = rescored,
        scaleAttrs = preparedAll$scaleAttrs
    )
}
