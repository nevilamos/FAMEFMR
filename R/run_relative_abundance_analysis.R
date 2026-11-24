#' Run the Relative Abundance Analysis
#'
#' Runs the complete species relative abundance analysis using the reactive
#' values associated with the current Shiny session.
#'
#' The function is designed to be called directly from an
#' \code{observeEvent()} or another expression evaluated within the Shiny
#' server environment. It retrieves the \code{rv} reactive values object from
#' the calling environment, so individual elements of \code{rv} do not need
#' to be passed as function arguments.
#'
#' The analysis includes:
#'
#' \itemize{
#'   \item reading the selected taxon list;
#'   \item constructing paths to habitat distribution model rasters;
#'   \item adjusting HDM paths for the selected raster resolution;
#'   \item checking that all required HDM rasters exist;
#'   \item reading and preparing species abundance response data;
#'   \item constructing species abundance lookup arrays;
#'   \item calculating species relative abundance;
#'   \item producing annual and grouped species abundance summaries;
#'   \item calculating changes in abundance relative to baseline years;
#'   \item writing analysis outputs to CSV files;
#'   \item saving the current FAME settings; and
#'   \item automatically saving the complete analysis as a QS file.
#' }
#'
#' @param .env Environment containing the Shiny reactive values object
#'   \code{rv}. The default is the calling environment returned by
#'   \code{parent.frame()}. This argument normally does not need to be
#'   supplied when the function is called directly from an
#'   \code{observeEvent()} inside the Shiny server function.
#'
#' @return Invisibly returns \code{NULL}. The function modifies the
#'   \code{rv} reactive values object directly and writes output files to
#'   \code{rv$resultsDir}.
#'
#' @details
#' The calling environment, or one of its parent environments, must contain:
#'
#' \itemize{
#'   \item \code{rv}, a \code{shiny::reactiveValues()} object;
#'   \item \code{EFG_TSF_4GS}, the growth-stage to time-since-fire lookup
#'     table; and
#'   \item \code{outputNames}, the object used when saving the settings table.
#' }
#'
#' The function reads and modifies the following elements of \code{rv}:
#'
#' \itemize{
#'   \item \code{FHAnalysis};
#'   \item \code{startBaseline};
#'   \item \code{endBaseline}, when available;
#'   \item \code{spListChoice};
#'   \item \code{customSpList};
#'   \item \code{TaxonListPath};
#'   \item \code{TaxonList};
#'   \item \code{RasterRes};
#'   \item \code{spResponseChoice};
#'   \item \code{customResponseFile};
#'   \item \code{SpGSResponses};
#'   \item \code{abundByGS};
#'   \item \code{allCombs};
#'   \item \code{makeRArasters};
#'   \item \code{resultsDir};
#'   \item \code{yearsForRasters};
#'   \item \code{cropRasters};
#'   \item \code{SpYearSumm};
#'   \item \code{grpSpYearSummLong};
#'   \item \code{raDeltaAbundWide};
#'   \item \code{raDeltaAbundLong}; and
#'   \item \code{outputSettings}.
#' }
#'
#' If \code{rv$endBaseline} is absent or \code{NULL}, the value of
#' \code{rv$startBaseline} is used as the end of the baseline period. This
#' preserves compatibility with the original implementation, which used the
#' start baseline value for both ends of the baseline period.
#'
#' The function is intended to be called directly from the Shiny server
#' environment. If it is called through another wrapper function, supply the
#' relevant server environment explicitly through \code{.env}.
#'
#' @section Shiny usage:
#'
#' A typical call from an \code{observeEvent()} is:
#'
#' \preformatted{
#' observeEvent(input$runRA, {
#'   run_relative_abundance_analysis()
#' })
#' }
#'
#' Busy indicators can remain in the \code{observeEvent()}:
#'
#' \preformatted{
#' observeEvent(input$runRA, {
#'   withBusyIndicatorServer("runRA", {
#'     withBusyIndicatorServer("runRA_TFI", {
#'       run_relative_abundance_analysis()
#'     })
#'   })
#' })
#' }
#'
#' @section Side effects:
#'
#' The function modifies the \code{rv} reactive values object directly.
#'
#' It writes the following files:
#'
#' \itemize{
#'   \item \code{SpYearSummLong.csv};
#'   \item \code{SpYearSummWide.csv};
#'   \item \code{grpSpYearSummLong.csv};
#'   \item \code{SppSummChangeRelativetoBaseline.csv}; and
#'   \item \code{autoSavedAnalysis.qs}.
#' }
#'
#' Additional settings files may be written by
#' \code{saveSettingsTable()}.
#'
#' @seealso
#' \code{\link{calc_SpeciesRA}},
#' \code{\link{calcDeltaAbund}},
#' \code{\link{make_Spp_LU_list}}
#'
#' @export
run_relative_abundance_analysis <- function(.env = parent.frame()) {

  # ---------------------------------------------------------------------------
  # Retrieve objects from the calling Shiny server environment
  # ---------------------------------------------------------------------------

  if (!exists("rv", envir = .env, inherits = TRUE)) {
    stop(
      paste0(
        "The reactive values object `rv` could not be found in the ",
        "calling environment."
      ),
      call. = FALSE
    )
  }

  rv <- get(
    "rv",
    envir = .env,
    inherits = TRUE
  )

  if (!exists("EFG_TSF_4GS", envir = .env, inherits = TRUE)) {
    stop(
      paste0(
        "The object `EFG_TSF_4GS` could not be found in the ",
        "calling environment."
      ),
      call. = FALSE
    )
  }

  EFG_TSF_4GS <- get(
    "EFG_TSF_4GS",
    envir = .env,
    inherits = TRUE
  )

  if (!exists("outputNames", envir = .env, inherits = TRUE)) {
    stop(
      paste0(
        "The object `outputNames` could not be found in the ",
        "calling environment."
      ),
      call. = FALSE
    )
  }

  outputNames <- get(
    "outputNames",
    envir = .env,
    inherits = TRUE
  )

  # ---------------------------------------------------------------------------
  # Check that a fire-history analysis is available
  # ---------------------------------------------------------------------------

  shiny::validate(
    shiny::need(
      rv$FHAnalysis,
      "You need to run or load a FH analysis first"
    )
  )

  # ---------------------------------------------------------------------------
  # Define baseline years
  # ---------------------------------------------------------------------------

  startBaseline <- as.integer(rv$startBaseline)

  if (is.null(rv$endBaseline)) {
    endBaseline <- as.integer(rv$startBaseline)
  } else {
    endBaseline <- as.integer(rv$endBaseline)
  }

  Baseline <- seq.int(
    from = startBaseline,
    to = endBaseline
  )

  # ---------------------------------------------------------------------------
  # Read the taxon list
  # ---------------------------------------------------------------------------

  if (isFALSE(rv$spListChoice)) {

    rv$TaxonListPath <- fame_data_path(
      "ReferenceTables",
      "FAME_TAXON_LIST.csv"
    )

  } else {

    rv$TaxonListPath <- rv$customSpList
  }

  rv$TaxonList <- readr::read_csv(
    rv$TaxonListPath,
    show_col_types = FALSE
  )

  # ---------------------------------------------------------------------------
  # Construct full paths to the HDM rasters
  # ---------------------------------------------------------------------------

  rv$TaxonList$HDMPath <- fame_data_path(
    rv$TaxonList$HDMPath
  )
  rv$TaxonList$HDMPath[1]
  # ---------------------------------------------------------------------------
  # Change HDM paths to reflect 75 m or 225 m raster versions
  # ---------------------------------------------------------------------------

  if (rv$RasterRes == 75) {

    rv$TaxonList$HDMPath <- gsub(
      "/225m/",
      "/75m/",
      rv$TaxonList$HDMPath
    )

    rv$TaxonList$HDMPath <- gsub(
      "_225m_",
      "_75m_",
      rv$TaxonList$HDMPath
    )

  } else if (rv$RasterRes == 225) {

    rv$TaxonList$HDMPath <- gsub(
      "/75m/",
      "/225m/",
      rv$TaxonList$HDMPath
    )
  }

  # ---------------------------------------------------------------------------
  # Check that required HDM rasters exist
  # ---------------------------------------------------------------------------

  HDMPathsToCheck <- rv$TaxonList$HDMPath[
    rv$TaxonList$Include == "Yes"
  ]

  MissingHDMS <- HDMPathsToCheck[
    !vapply(
      HDMPathsToCheck,
      file.exists,
      logical(1)
    )
  ]
  if (length(MissingHDMS) > 0) {

    msg <- paste(
      "The following HDM rasters are missing:\n\n",
      paste(MissingHDMS, collapse = "\n")
    )

    if (inherits(rv, "reactivevalues")) {

      shiny::showModal(
        shiny::modalDialog(
          title = "Missing HDM Rasters",
          HTML(
            paste(
              "The following HDM rasters are missing:<br><br>",
              paste(MissingHDMS, collapse = "<br>")
            )
          ),
          easyClose = TRUE,
          footer = shiny::modalButton("OK")
        )
      )

      return(invisible(rv))

    } else {

      stop(msg, call. = FALSE)

    }
  }



  # ---------------------------------------------------------------------------
  # Identify included species and species rasters to write
  # ---------------------------------------------------------------------------

  HDMSpp_NO <- rv$TaxonList$TAXON_ID[
    rv$TaxonList$Include == "Yes"
  ]

  writeSp <- rv$TaxonList$TAXON_ID[
    rv$TaxonList$WriteSpeciesRaster == "Yes"
  ]

  writeSp <- writeSp[
    writeSp %in% HDMSpp_NO
  ]

  # ---------------------------------------------------------------------------
  # Select species growth-stage response data
  # ---------------------------------------------------------------------------

  if (isFALSE(rv$spResponseChoice)) {

    rv$SpGSResponses <- fame_data_path(
      "ReferenceTables",
      "OrdinalExpertLong.csv"
    )

  } else {

    rv$SpGSResponses <- file.path(
      rv$customResponseFile
    )
  }

  # ---------------------------------------------------------------------------
  # Read and prepare species abundance data
  # ---------------------------------------------------------------------------

  if (isTRUE(rv$abundByGS)) {

    AbundDataByGS <- readr::read_csv(
      rv$SpGSResponses,
      show_col_types = FALSE
    ) |>
      dplyr::select(
        EFG_NO,
        GS4_NO,
        FireType,
        Abund,
        TAXON_ID
      )

    # Expand growth-stage abundance data to full time-since-fire format.
    #
    # The inner join prevents NA values for YSF where no YSF values exist for
    # a GS4_NO that is included in AbundDataByGS.

    AbundDataLong <- AbundDataByGS |>
      dplyr::mutate(
        FireTypeNo = dplyr::if_else(
          FireType == "High",
          2,
          dplyr::if_else(
            FireType == "Low",
            1,
            0
          )
        )
      ) |>
      dplyr::inner_join(
        EFG_TSF_4GS,
        by = c(
          "EFG_NO",
          "GS4_NO"
        ),
        relationship = "many-to-many"
      ) |>
      dplyr::arrange(
        TAXON_ID
      )

  } else {

    # Read abundance data already supplied in full long format.

    AbundDataLong <- readr::read_csv(
      rv$SpGSResponses,
      show_col_types = FALSE
    ) |>
      dplyr::arrange(
        TAXON_ID
      )
  }

  # ---------------------------------------------------------------------------
  # Construct lookup arrays for cell-by-cell abundance calculations
  # ---------------------------------------------------------------------------

  message("Making species abundance lookup list")

  LU_List <- make_Spp_LU_list(
    myHDMSpp_NO = HDMSpp_NO,
    myAbundDataLong = AbundDataLong
  )

  message("Finished species abundance lookup list")

  # ---------------------------------------------------------------------------
  # Run the main species relative abundance calculation
  # ---------------------------------------------------------------------------

  message("Making species annual summaries")

  rv$SpYearSumm <- calc_SpeciesRA(
    myFHAnalysis = rv$FHAnalysis,
    myAllCombs = rv$allCombs,
    myHDMSpp_NO = HDMSpp_NO,
    myWriteSpRasters = rv$makeRArasters,
    myResultsDir = rv$resultsDir,
    myLU_List = LU_List,
    myTaxonList = rv$TaxonList,
    writeYears = rv$yearsForRasters,
    myWriteSp = writeSp,
    myCropRasters = rv$cropRasters
  )


  # ---------------------------------------------------------------------------
  # Save species abundance summary outputs
  # ---------------------------------------------------------------------------

  readr::write_csv(
    rv$SpYearSumm$SpYearSummLong,
    file.path(
      rv$resultsDir,
      "SpYearSummLong.csv"
    )
  )

  readr::write_csv(
    rv$SpYearSumm$SpYearSummWide,
    file.path(
      rv$resultsDir,
      "SpYearSummWide.csv"
    )
  )

  # ---------------------------------------------------------------------------
  # Create grouped species abundance summaries
  # ---------------------------------------------------------------------------

  rv$grpSpYearSummLong <- rv$SpYearSumm$grpSpYearSumm |>
    dplyr::rename(
      Index_AllCombs = `myAllCombs$Index_AllCombs`
    ) |>
    tidyr::pivot_longer(
      cols = -tidyr::one_of(
        "TAXON_ID",
        "Index_AllCombs"
      ),
      names_to = "SEASON",
      values_to = "sumRA"
    )

  # ---------------------------------------------------------------------------
  # Check whether planning units or burn units are used
  # ---------------------------------------------------------------------------

  if (
    match(
      "PU",
      names(rv$allCombs$U_AllCombs_TFI),
      nomatch = 0
    ) > 0
  ) {

    rv$grpSpYearSummLong <- rv$grpSpYearSummLong |>
      dplyr::mutate(
        PU = rv$allCombs$U_AllCombs_TFI$PU[
          Index_AllCombs
        ],
        EFG_NAME = rv$allCombs$U_AllCombs_TFI$EFG_NAME[
          Index_AllCombs
        ]
      ) |>
      dplyr::group_by(
        TAXON_ID,
        PU,
        EFG_NAME,
        SEASON
      )

  } else {

    rv$grpSpYearSummLong <- rv$grpSpYearSummLong |>
      dplyr::mutate(
        EFG_NAME = rv$allCombs$U_AllCombs_TFI$EFG_NAME[
          Index_AllCombs
        ]
      ) |>
      dplyr::group_by(
        TAXON_ID,
        EFG_NAME,
        SEASON
      )
  }

  rv$grpSpYearSummLong <- rv$grpSpYearSummLong |>
    dplyr::summarise(
      sumRA = sum(sumRA),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      TAXON_ID = as.integer(TAXON_ID)
    )

  readr::write_csv(
    rv$grpSpYearSummLong,
    file.path(
      rv$resultsDir,
      "grpSpYearSummLong.csv"
    )
  )

  gc()

  message("Finished species annual summaries")

  # ---------------------------------------------------------------------------
  # Calculate abundance changes relative to the baseline years
  # ---------------------------------------------------------------------------

  rv$raDeltaAbundWide <- calcDeltaAbund(
    SpYearSumm = rv$SpYearSumm$SpYearSummWide,
    myFHAnalysis = rv$FHAnalysis,
    myBaseline = Baseline
  )

  readr::write_csv(
    rv$raDeltaAbundWide,
    file.path(
      rv$resultsDir,
      "SppSummChangeRelativetoBaseline.csv"
    )
  )

  gc()

  # ---------------------------------------------------------------------------
  # Create long-format abundance change data for plots
  # ---------------------------------------------------------------------------

  rv$raDeltaAbundLong <- rv$raDeltaAbundWide |>
    dplyr::select(
      -tidyr::one_of(
        "TAXON_ID",
        "DIVNAME",
        "FFG_ACT_STATUS",
        "EPBC_ACT_STATUS",
        "CombThreshold",
        "NoLessthanThreshhold",
        "LastLessThanThreshold"
      )
    ) |>
    tidyr::pivot_longer(
      cols = -tidyr::one_of(
        "COMMON_NAME",
        "SCIENTIFIC_NAME",
        "Baseline"
      ),
      names_to = "SEASON",
      values_to = "DeltaRA"
    ) |>
    dplyr::mutate(
      SEASON = as.integer(SEASON)
    )

  message("Finished relative abundance change calculations")

  # # ---------------------------------------------------------------------------
  # # Save fire-history and species settings
  # # ---------------------------------------------------------------------------
  #
  # rv$outputSettings <- saveFHsettings(
  #   rv$outputSettings,
  #   rv
  # )
  #
  # rv$outputSettings <- saveSPsettings(
  #   rv$outputSettings,
  #   rv
  # )
  #
  # saveSettingsTable(
  #   rv$outputSettings,
  #   rv$resultsDir,
  #   outputNames
  # )

  # ---------------------------------------------------------------------------
  # Save the complete analysis so that it can be restored
  # ---------------------------------------------------------------------------
  autoSavePath <- fame_data_path(
      "FH_Outputs",
      "autoSavedAnalysis.qs"
      )


  if (shiny::isRunning()) {
  myRvList <- shiny::reactiveValuesToList(
    rv
  )
  qs::qsave(
    myRvList,
    autoSavePath
  )
  }else{qs::qsave(
    rv,
    autoSavePath
  )}



  message("Relative abundance analysis complete")

  invisible(rv)
}
