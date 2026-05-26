#' Add interval, YSF/LBY, and last fire type (LFT) columns to an output data frame
#'
#' @description
#' Wraps the provided workflow:
#' - Identifies SEAS* and FireType* columns in OutDF
#' - Computes interval columns INT01..INTxx from SEAS years (0 treated as missing)
#' - Optionally updates FireType matrix using fireTypeLowToHigh() when max_interval > 0
#' - Computes LBY (last burnt year), YSF (years since fire), and LFT (last fire type) for each year
#'   in a calendar-year TimeSpan and appends them to OutDF.
#'
#' Assumptions:
#' - SEAS values are calendar years, with 0 meaning "no event"
#' - LBY_f(M, y) returns NA when no previous fire exists for year y (as in your implementation)
#'
#' @param OutDF data.frame containing SEAS* and FireType* columns
#' @param v terra spatVector fire history containing SEASON field (calendar years)
#' @param max_interval integer; if > 0 calls fireTypeLowToHigh(); if 0 skips; if < 0 errors
#' @param start.SEASON optional integer start year; if NA uses 2nd smallest unique v$SEASON (or smallest if only one)
#' @param end.SEASON optional integer end year; if NA uses max(v$SEASON); otherwise max(end.SEASON, max(v$SEASON))
#' @param seas_pattern regex/pattern used to find season columns in OutDF (default "SEAS")
#' @param firetype_pattern regex/pattern used to find fire type columns in OutDF (default "FireType")
#' @param LBY_f function(M, y) returning last burnt year (calendar year) per row; NA if none
#' @param fireTypeLowToHigh optional function(max_interval, Interval_Matrix, Firetype_Matrix) returning updated fire types
#' @param verbose logical; print progress messages
#'
#' @return Updated OutDF with INT*, YSF*, LBY*, LFT* columns appended
#' @export
add_fire_lft_lby_ysf <- function(
    OutDF,
    v,
    max_interval = 0L,
    start.SEASON = NA_integer_,
    end.SEASON   = NA_integer_,
    seas_pattern = "SEAS",
    firetype_pattern = "FireType",
    #LBY_f,
    #fireTypeLowToHigh = NULL,
    verbose = TRUE
) {
  # ---- checks ----
  #if (!is.data.frame(OutDF)) stop("OutDF must be a data.frame")
  if (missing(v) || is.null(v) || is.null(v$SEASON)) stop("v must be provided and contain column 'SEASON'.")

  if (length(max_interval) != 1L || is.na(max_interval)) stop("max_interval must be a non-NA scalar.")
  if (!is.numeric(max_interval)) stop("max_interval must be numeric/integer.")
  if (max_interval < 0) stop("max interval cannot be less than 0")
  max_interval <- as.integer(max_interval)

  #if (missing(LBY_f) || !is.function(LBY_f)) stop("LBY_f must be supplied as a function.")
  if (max_interval > 0 && (is.null(fireTypeLowToHigh) || !is.function(fireTypeLowToHigh))) {
    stop("fireTypeLowToHigh must be supplied as a function when max_interval > 0.")
  }

  # ---- identify columns ----
  cn <- colnames(OutDF)
  SEASNames <- cn[grep(pattern = seas_pattern, x = cn)]
  FTNames   <- cn[grep(pattern = firetype_pattern, x = cn)]

  if (length(SEASNames) < 2L) stop("Need at least 2 season columns (SEAS*) to compute intervals.")
  if (length(FTNames)  < 1L)  stop("No fire type columns found (FireType*).")

  SEAS_Matrix <- as.matrix(OutDF[, SEASNames, drop = FALSE])
  FT_matrix   <- as.matrix(OutDF[, FTNames,   drop = FALSE])

  # enforce numeric for year comparisons in LBY_f (and interval math)
  suppressWarnings(storage.mode(SEAS_Matrix) <- "numeric")

  # ---- compute intervals ----
  Cols <- ncol(SEAS_Matrix)

  # for interval calc, treat 0 as missing
  SEAS_Matrix_int <- SEAS_Matrix
  SEAS_Matrix_int[SEAS_Matrix_int == 0] <- NA_real_

  Interval <- SEAS_Matrix_int[, 2:Cols, drop = FALSE] - SEAS_Matrix_int[, 1:(Cols - 1), drop = FALSE]
  IntNames <- paste0("INT", sprintf("%02d", 1:(Cols - 1)))
  colnames(Interval) <- IntNames

  OutDF <- cbind(OutDF, as.data.frame(Interval, check.names = FALSE))

  # ---- optional fire type adjustment ----
  if (max_interval > 0) {
    FT_matrix <- fireTypeLowToHigh(
      max_interval     = as.integer(max_interval),
      Interval_Matrix  = Interval,
      Firetype_Matrix  = FT_matrix
    )
    OutDF[, FTNames] <- FT_matrix
  }

  # ---- derive calendar-year TimeSpan ----
  seas_v <- v$SEASON
  seas_v <- seas_v[!is.na(seas_v)]
  if (!length(seas_v)) stop("v$SEASON has no non-NA values.")
  if (any(abs(seas_v - round(seas_v)) > 0)) stop("v$SEASON must be integer-like (calendar years).")

  seas_unique <- sort(unique(as.integer(seas_v)))
  min.SEASON <- if (length(seas_unique) >= 2L) seas_unique[2L] else seas_unique[1L]

  if (is.na(start.SEASON)) {
    start.SEASON <- min.SEASON
  } else {
    start.SEASON <- as.integer(start.SEASON)
    if (start.SEASON < min.SEASON) start.SEASON <- min.SEASON
  }

  if (is.na(end.SEASON)) {
    max.SEASON <- max(seas_unique)
  } else {
    end.SEASON <- as.integer(end.SEASON)
    max.SEASON <- max(end.SEASON, max(seas_unique))
  }

  if (start.SEASON > max.SEASON) stop("start.SEASON is greater than max.SEASON after adjustment.")

  TimeSpan <- start.SEASON:max.SEASON
  LTR <- length(TimeSpan)

  # ---- compute LBY using your LBY_f (expects 0 == no event, returns NA if none) ----
  LBY <- matrix(NA_integer_, nrow(SEAS_Matrix), LTR)
  for (i in seq_len(LTR)) {
    y <- TimeSpan[i]
    LBY[, i] <- tryCatch(
      as.integer(LBY_f(M = SEAS_Matrix, y = y)),
      error = function(e) stop(sprintf("LBY_f failed for y=%s: %s", y, e$message), call. = FALSE)
    )
  }
  colnames(LBY) <- LBYNames<-paste0("LBY", TimeSpan)

  # ---- compute YSF explicitly (NA when LBY is NA) ----
  YSF <- matrix(NA_integer_, nrow(SEAS_Matrix), LTR)
  for (i in seq_len(LTR)) {
    y <- TimeSpan[i]
    ok <- !is.na(LBY[, i])
    YSF[ok, i] <- y - LBY[ok, i]
  }
  colnames(YSF) <- YSFNames<-paste0("YSF", TimeSpan)

  # ---- build lookup matrix for last fire type by calendar year ----
  # For lookup, treat 0 as missing
  SEAS_Matrix_lookup <- SEAS_Matrix
  SEAS_Matrix_lookup[SEAS_Matrix_lookup == 0] <- NA_real_

  years_from_seas <- as.integer(stats::na.omit(as.vector(SEAS_Matrix_lookup)))
  years_from_lby  <- as.integer(stats::na.omit(as.vector(LBY)))
  years_all <- c(years_from_seas, years_from_lby, TimeSpan)
  years_all <- years_all[!is.na(years_all)]
  if (!length(years_all)) stop("No valid years available to build the lookup matrix (LUM).")

  year_axis <- min(years_all):max(years_all)

  LUM <- matrix(NA, nrow(SEAS_Matrix_lookup), length(year_axis))
  colnames(LUM) <- as.character(year_axis)

  for (i in seq_len(nrow(SEAS_Matrix_lookup))) {
    C_years <- as.integer(stats::na.omit(SEAS_Matrix_lookup[i, ]))
    if (!length(C_years)) next

    idx <- match(C_years, year_axis)
    if (anyNA(idx)) stop(sprintf("Row %d has event years outside lookup range; unexpected.", i), call. = FALSE)

    if (ncol(FT_matrix) < length(C_years)) {
      stop(sprintf(
        "Row %d: need at least %d FireType columns but only %d found.",
        i, length(C_years), ncol(FT_matrix)
      ), call. = FALSE)
    }

    V <- FT_matrix[i, seq_len(length(C_years))]
    LUM[i, idx] <- V
  }

  if (verbose) message("calculating last fire type")

  LFT <- matrix(NA, nrow(SEAS_Matrix_lookup), LTR)
  for (i in seq_len(nrow(SEAS_Matrix_lookup))) {
    idx_lby <- match(as.integer(LBY[i, ]), year_axis)
    sel <- !is.na(idx_lby)
    if (any(sel)) LFT[i, sel] <- LUM[i, idx_lby[sel]]
  }
  colnames(LFT) <- LFTNames<-paste0("LFT", TimeSpan)

  # ---- bind results back ----
  OutDF <- cbind(
    OutDF,
    as.data.frame(YSF, check.names = FALSE),
    as.data.frame(LBY, check.names = FALSE),
    as.data.frame(LFT, check.names = FALSE)
  )

  results <-
    list(
      OutDF = OutDF,
      TimeSpan = TimeSpan,
      YSFNames = YSFNames,
      LBYNames = LBYNames,
      LFTNames = LFTNames
    )
}

