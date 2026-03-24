#' @rdname validatePSM
#'
#' @title Validate a PSM
#'
#' @description
#' Validate a PSM by checking some of its spectral features. Calling
#' `validatePSM()` applies the different checks mentioned below. Failing a check
#' does not guarantee a bad identification per se, but it does provide
#' information to assess the confidence in dubious PSMs or confirm already good
#' identifications.
#'
#' @param x A `Spectra` object with identifications. The identification
#' sequences need to be stocked in a variable name called with the parameter
#' `sequenceVariable`. Refer to `Spectra::joinSpectraData` for more
#' informations.
#'
#' @param sequenceVariable A `character` of length 1L, representing the peptide
#' sequence for that identification.
#'
#' @param ... Arguments passed down to `checkParentIonIntensity()`
#'
#' @returns A `data.frame()` that checks the different validation metrics.
#'
#' @importFrom Spectra spectraData peaksData
#'
#' @importFrom PTMods convertAnnotation
#'
#' @export
#'
#' @examples
#'
#' library(Spectra)
#'
#' sp <- Spectra("/home/guillaumedeflandre/Documents/drive_UCL/PHD/rformassspectrometry/PSMatch-oriented/PSMatch-paper/2025-PSMatch-MiMB/data/001_exp1-PO4_PO4_D20.mzML")
#'
#' sage_results <- read.delim("/home/guillaumedeflandre/Documents/drive_UCL/PHD/rformassspectrometry/PSMatch-oriented/PSMatch-paper/2025-PSMatch-MiMB/data/PXD039419-sage-results.tsv")
#'
#' sage_results[, "label"] <- sage_results[, "label"] < 0
#'
#' psms <- PSM(
#' sage_results,
#' spectrum = "scannr",
#' peptide = "peptide",
#' protein = "proteins",
#' decoy = "label",
#' rank = "rank",
#' score = "sage_discriminant_score",
#' fdr = "spectrum_q")
#'
#' psms <- filterPsmRank(psms)
#'
#' head(sp$pkey <- paste0(basename(sp$dataOrigin), sub("ˆ.+scan=", "::", sp$spectrumId)))
#' head(psms$pkey <- paste0(basename(psms$filename), sub("ˆ.+scan=", "::", psms$scannr)))
#'
#' sp <- joinSpectraData(sp, psms, by.x = "pkey")
#'
#' validatePSM(x = sp[3492:3500], sequenceVariable = "peptide", fdr = "spectrum_q")
validatePSM <- function(x, sequenceVariable = "peptide",
                        fdr = "spectrum_q", ...) {

    stopifnot(requireNamespace("Spectra"))
    stopifnot(inherits(x, "Spectra"))

    v <- Spectra::peaksData(x) ## all spectra (MS1 & MS2) peaks data

    sequenceIds <- Spectra::spectraData(x, sequenceVariable)[, 1L]
    identifications <- which(!is.na(sequenceIds))

    x$sequence <- sequenceIds ## Create 'sequence' variable for labelFragments

    x_sub <- x[identifications] ## Fetch only PSMs
    sub_v <- v[identifications] ## Identified spectra peak data
    stripped_sequence <- PTMods::getCanonicalSequence(sequenceIds[identifications])

    stopifnot(length(sub_v) != 0)

    frags <- suppressWarnings(labelFragments(x_sub,
        type = c("a", "b", "c", "x", "y", "z")))
    frags_mz <- suppressWarnings(labelFragments(x_sub,
        type = c("b", "y"), what = "mz"))

    ab <- checkABpresence(x_sub, fragments = frags)

    xy <- checkXYpresence(x_sub, fragments = frags)

    overlap <- checkOverlap(x_sub, fragments = frags)

    xrea <- checkXrea(x_sub, peaks = sub_v)

    shifts <- checkShiftConsistency(x_sub, fragments = frags)

    parent <- checkParentIonIntensity(x_sub, peaks = sub_v, ...)

    purity_ints <- checkPrecursorPurity(x)[identifications]

    res <- data.frame(spectrumId = Spectra::spectraData(x_sub, "spectrumId")[, 1L],
                      scanIndex = Spectra::scanIndex(x_sub),
                      peptide = sequenceIds[identifications],
                      canonical_seq = stripped_sequence,
                      fdr = Spectra::spectraData(x_sub, fdr)[, 1L],
                      a2b2_presence = ab,
                      x2y2_presence = xy,
                      by_overlap = overlap,
                      xrea = xrea,
                      shift_consistency = shifts,
                      parent_ion_int = parent,
                      precursor_purity = purity_ints)

    rownames(res) <- NULL
    return(res)
}

#' @rdname validatePSM
#'
#' @param fragments The result of `labelFragments()` on `x`.
#'
#' @returns `checkABpresence()` : `TRUE` if the a2-b2 fragments are both present.
#'
#' @export
checkABpresence <- function(x, fragments = NULL) {

    if (!length(fragments)) {
        labels <- labelFragments(x, type = c("a", "b"))
    } else {
        labels <- fragments
    }
    .abCouple <- function(labs) "a2" %in% labs & "b2" %in% labs
    unlist(lapply(labels, .abCouple))
}

#' @rdname validatePSM
#'
#' @param fragments The result of `labelFragments()` on `x`.
#'
#' @returns `checkXYpresence()` : `TRUE` if the x2-y2 fragments are both present.
#'
#' @export
checkXYpresence <- function(x, fragments = NULL) {

    if (!length(fragments)) {
        labels <- labelFragments(x, type = c("x", "y"))
    } else {
        labels <- fragments
    }
    .xyCouple <- function(labs) "x2" %in% labs & "y2" %in% labs
    unlist(lapply(labels, .xyCouple))
}

#' @rdname validatePSM
#'
#' @param fragments The result of `labelFragments()` on `x`.
#'
#' @returns `checkOverlap()` : Is there an overlap between b- and y-ions ? If
#' there is and no modification was identified, there might be an unsearched
#' for modification present. If there is not, this further confirms the
#' identification.  If the overlap checks out it returns `TRUE`, if not `FALSE`.
#'
#' @export
checkOverlap <- function(x, sequenceVariable = "peptide",
                         fragments = NULL, strippedSeq = NULL) {

    if (length(strippedSeq)) {
        stripped_sequence <- strippedSeq
    } else {
        sequence_ids  <- Spectra::spectraData(x, sequenceVariable)[, 1L]
        stripped_sequence <- gsub("[^A-Za-z]", "", sequence_ids)
    }

    if (!length(fragments)) {
        labels <- suppressWarnings(labelFragments(x))
    } else {
        labels <- fragments
    }

    index <- lapply(labels, function(x) as.integer(gsub("\\D", "", x)))
    b_ions <- lapply(labels, function(x) which(grepl("b", x)))
    y_ions <- lapply(labels, function(x) which(grepl("y", x)))

    ans <- vector(length = length(x))

    for (i in seq_along(x)) {
        max_b <- max(index[[i]][b_ions[[i]]], 0)
        max_y <- max(index[[i]][y_ions[[i]]], 0)
        ## n - max(y) > max(b), if respected = no overlap
        ans[i] <- nchar(stripped_sequence[i]) - max_y > max_b
    }
    return(ans)
}

#' @rdname validatePSM
#'
#' @param peaks The spectrum peak data (result from `Spectra::peaksData()`).
#'
#' @returns `checkXrea()` : Spectrum quality metric defined in
#' Na, Seungjin, and Eunok Paek. 2006. The higher the returned value, the better
#' the quality of the spectrum.
#'
#' @export
checkXrea <- function(x, peaks = NULL) {

    if (length(peaks)) {
        ints_list <- lapply(peaks, function(x) x[, "intensity"])
    } else {
        ints_list <- intensity(x)
    }

    xrea_list <- mapply(function(ints) {
        n <- length(ints)
        tic <- sum(ints)
        rel_ints <- ints / tic

        # Cumulative intensities (CIi)
        cum_ints <- cumsum(rel_ints)

        # Ideal cumulative line: linear from (0,0) to (1,1)
        ideal_cum <- seq_along(ints) / length(ints)

        # Area between real and ideal curve using strip method (fixed bin width = 1/n)
        diff_area <- sum(abs(cum_ints - ideal_cum)) * (1 / n)

        # Alpha = difference between top two cumulative intensities
        top_cis <- sort(cum_ints, decreasing = TRUE)
        alpha <- abs(top_cis[1] - top_cis[2])

        # Triangle area is fixed at 0.5 (maximum possible area)
        triangle_area <- 0.5

        # Final Xrea formula
        xrea <- diff_area / (triangle_area + alpha)
        return(xrea)
    }, ints_list, SIMPLIFY = TRUE)

    return(xrea_list)
}

#' @rdname validatePSM
#'
#' @param fragments The result of `labelFragments()` on `x`.
#'
#' @returns `checkShiftConsistency()` : In case of modifications present:
#' returns the percentage of potential mass shifts actually matched. If not
#' applicapble because there are no modifications: `NA`.
#'
#' @export
checkShiftConsistency <- function(x, fragments = NULL, strippedSeq = NULL) {

    if (length(fragments)) {
        labels = fragments
    } else {
        labels <- suppressWarnings(labelFragments(x))
    }

    if (length(strippedSeq)) {
        stripped_sequence <- strippedSeq
    } else {
        stripped_sequence <- PTMods::getCanonicalSequence(x$sequence)
    }

    index <- lapply(labels, function(x) as.integer(gsub("\\D", "", x)))
    b_ions <- lapply(labels, function(x) which(grepl("b", x)))
    y_ions <- lapply(labels, function(x) which(grepl("y", x)))

    ans <- vector(length = length(x))

    for (i in seq_along(x)) {

        if (is.na(x$sequence[i])) {
            ans[i] <- NA
        } else if (x$sequence[i] == stripped_sequence[i]) {
            ans[i] <- NA
        } else {
            parsed_mods <- PTMods:::.parseModifiedSequence(x$sequence[i])
            pep_len <- nchar(stripped_sequence[i])
            mod_pos <- which(parsed_mods != 0)
            max_b <- mod_pos[1] - 1
            max_y <- pep_len - mod_pos[length(mod_pos)]

            b_matches <- sum(unique(index[[i]][b_ions[[i]]]) > max_b)
            y_matches <- sum(unique(index[[i]][y_ions[[i]]]) > max_y)

            if (sum(max_b, max_y) != 0) {
                ans[i] <- sum(b_matches, y_matches)/sum(max_b, max_y)
                } else ans[i] <- 1
        }
    }
    return(ans)
}

#' @rdname validatePSM
#'
#' @param peaks The spectrum peak data (result from `Spectra::peaksData()`).
#'
#' @param tolerance `Numeric(1L)` The tolerance to use when matching peaks.
#'
#' @param ppm `Numeric(1L)` The ppm value to use when matching peaks that is added to `tolerance`.
#'
#' @importFrom MsCoreUtils common
#'
#' @returns `checkParentIonIntensity()` : The relative intensity of the parent
#' ion.
#'
#' @export
checkParentIonIntensity <- function (x,
                                     peaks = NULL,
                                     tolerance = 0,
                                     ppm = 20) {

    stopifnot(requireNamespace("Spectra"))
    stopifnot(inherits(x, "Spectra"))

    k <- numeric()
    if (length(peaks)) {
        v <- peaks
    } else {
        v <- Spectra::peaksData(x)
    }
    precursors <- Spectra::precursorMz(x)

    for (i in seq_along(x)) {

        peak <- v[[i]]
        precursor_index <- which(MsCoreUtils::common(peak[, "mz"],
                                                     precursors[i],
                                                     tolerance = tolerance,
                                                     ppm = ppm))
        if (length(precursor_index)) {
            ints <- peak[precursor_index, "intensity"]
        } else {ints <- 0}
        max_ints <- max(peak[, "intensity"])
        k[i] <- ints/max_ints
    }
    return(k)
}

#' @rdname validatePSM
#'
#' @param x A `Spectra` object containing **both MS1 and MS2 spectra** from
#'   the same run(s). MS1 spectra are used as the source for isolation window
#'   peak data; each MS2 is paired with the nearest preceding MS1 scan (sorted
#'   by retention time within each `dataOrigin`).
#'
#' @param tolerance `Numeric(1)` Absolute m/z half-width (in Da) of the
#'   isolation window used when `useReportedIsolationWindow = FALSE` (default
#'   `0.05`).
#'
#' @param ppm `Numeric(1)` Additional m/z-proportional tolerance (in ppm)
#'   added to `tolerance` when defining the isolation window
#'   (`default 0`).
#'
#' @param useReportedIsolationWindow `logical(1)` If `TRUE`, use the
#'   `isolationWindowLowerMz` / `isolationWindowUpperMz` metadata stored in the
#'   spectra rather than computing the window from `tolerance` and `ppm`.
#'   Defaults to `FALSE`.
#'
#' @param BPPARAM A `BiocParallelParam` instance controlling parallel
#'   evaluation (one job per unique `dataOrigin`). Defaults to
#'   `BiocParallel::SerialParam()`.
#'
#' @importFrom Spectra precursorPurity
#'
#' @returns `checkPrecursorPurity()` returns a `numeric` vector of length
#'   `length(x)`. Each value is the ratio of the most-intense peak to the
#'   total intensity within the isolation window of the corresponding MS1 scan,
#'   as computed by `Spectra::precursorPurity()`. MS1 spectra and MS2 spectra
#'   with no preceding MS1 scan return `NA`.
#'
#' @export
checkPrecursorPurity <- function(x, tolerance = 0.05, ppm = 0,
                                 useReportedIsolationWindow = FALSE,
                                 BPPARAM = BiocParallel::SerialParam()) {

    stopifnot(inherits(x, "Spectra"))

    Spectra::precursorPurity(x,
                            tolerance = tolerance,
                            ppm = ppm,
                            useReportedIsolationWindow =
                                useReportedIsolationWindow,
                            BPPARAM = BPPARAM)
}