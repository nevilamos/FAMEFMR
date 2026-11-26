#' Set up FAME data directories and input files
#'
#' Creates the requested FAME data directories and downloads and extracts
#' the specified input data if it is not already available.
#'
#' Directory names may include subdirectories, for example
#' `"results/dir1"` or `"results/dir1/subdir"`.
#'
#' @param newdirs Character vector containing directory paths, relative to
#'   the FAME data directory. Paths may contain subdirectories.
#' @param datazip Character. Name of the ZIP file containing the input data.
#' @param url Character. URL from which the ZIP file should be downloaded.
#' @param timeout Numeric. Maximum download time in seconds. Defaults to
#'   600 seconds (10 minutes).
#'
#' @return Invisibly returns `TRUE`.
#'
#' @export
setup_fame_data <- function(
    newdirs = c(
      "AdHocPolygons",
      "CustomCSV",
      "FH_Outputs",
      "HDMS",
      "InputGeneralRasters",
      "PUPolygons",
      "rawFH",
      "ReferenceShapefiles",
      "ReferenceTables",
      "results"
    ),
    datazip = "FAME4_VG2020_inputfiles.zip",
    url = "https://ecological-risk-analysis.s3.ap-southeast-2.amazonaws.com/FAME_FMR/FAME4_VG2020_inputfiles.zip",
    timeout = 600
) {

  # Convert directory names, including nested directories,
  # to full FAME data paths
  newDirPaths <- vapply(
    newdirs,
    function(x) {
      parts <- strsplit(x, "[/\\\\]")[[1]]
      do.call(fame_data_path, as.list(parts))
    },
    FUN.VALUE = character(1)
  )

  # Create missing directories
  for (mypath in newDirPaths) {

    if (!dir.exists(mypath)) {

      dir.create(
        mypath,
        recursive = TRUE,
        showWarnings = FALSE
      )

      if (!dir.exists(mypath)) {
        stop("Could not create directory: ", mypath)
      }
    }
  }

  zip_path <- fame_data_path(datazip)

  # Temporary download filename
  tmp_zip <- paste0(zip_path, ".download")

  # Clean up any partial download left from a previous failed attempt
  if (file.exists(tmp_zip)) {
    unlink(tmp_zip)
  }

  # Download and extract if the ZIP is not already present
  if (!file.exists(zip_path)) {

    # Ensure partial downloads are removed if anything below fails
    on.exit({
      if (file.exists(tmp_zip)) {
        unlink(tmp_zip)
      }
    }, add = TRUE)

    # Temporarily increase R's download timeout
    old_timeout <- getOption("timeout")

    on.exit(
      options(timeout = old_timeout),
      add = TRUE
    )

    options(
      timeout = max(timeout, old_timeout)
    )

    message("Downloading FAME data...")

    tryCatch(

      download.file(
        url = url,
        destfile = tmp_zip,
        mode = "wb"
      ),

      error = function(e) {

        if (file.exists(tmp_zip)) {
          unlink(tmp_zip)
        }

        stop(
          "Failed to download FAME data: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    # Make sure a file was actually downloaded
    if (!file.exists(tmp_zip)) {
      stop(
        "Download completed without creating the expected file: ",
        tmp_zip,
        call. = FALSE
      )
    }

    # Rename only after successful download
    if (!file.rename(tmp_zip, zip_path)) {

      unlink(tmp_zip)

      stop(
        "Download completed but the file could not be renamed to: ",
        zip_path,
        call. = FALSE
      )
    }

    message("Download complete.")

    # Extract ZIP
    message("Extracting FAME data...")

    extraction_successful <- FALSE

    tryCatch({

      unzip(
        zipfile = zip_path,
        exdir = fame_data_path()
      )

      extracted_dir <- fame_data_path(
        tools::file_path_sans_ext(datazip)
      )

      if (dir.exists(extracted_dir)) {

        files <- list.files(
          extracted_dir,
          all.files = TRUE,
          full.names = TRUE,
          no.. = TRUE
        )

        file.copy(
          files,
          fame_data_path(),
          recursive = TRUE,
          overwrite = TRUE
        )

        unlink(
          extracted_dir,
          recursive = TRUE
        )
      }

      extraction_successful <- TRUE

    }, error = function(e) {

      stop(
        "Failed to extract FAME data: ",
        conditionMessage(e),
        call. = FALSE
      )
    })

    # Remove ZIP once it is no longer needed
    if (file.exists(zip_path)) {
      unlink(zip_path)
    }

    if (extraction_successful) {
      message("FAME data setup complete.")
    }
  }

  invisible(TRUE)
}
