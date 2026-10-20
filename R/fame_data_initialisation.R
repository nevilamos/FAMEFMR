# ============================================================
# FAME external-data initialisation
#
# This file:
#   * establishes FAME_DATA_DIR;
#   * creates writable FAME working directories;
#   * checks whether required reference data already exist;
#   * downloads missing reference data from S3 when necessary;
#   * avoids downloading again when valid local data are present;
#   * supports Windows and Linux/Shiny Server.
#
# Expected S3 structure:
#
#   s3://<bucket>/<prefix>/
#       HDMS/
#       InputGeneralRasters/
#       ReferenceShapefiles/
#       ReferenceTables/
#       rawFH/
#
# Environment variables:
#
#   FAME_DATA_DIR
#       Optional explicit local data directory.
#
#   FAME_S3_BUCKET
#       Required only when missing data must be downloaded.
#
#   FAME_S3_PREFIX
#       Optional prefix within the bucket.
#
#   FAME_AWS_PROFILE
#       Optional AWS CLI profile. Normally leave unset on EC2
#       and use the EC2 instance role.
#
# ============================================================


# ------------------------------------------------------------
# Internal constants
# ------------------------------------------------------------

.fame_required_data_dirs <- c(
  "HDMS",
  "InputGeneralRasters",
  "ReferenceShapefiles",
  "ReferenceTables",
  "rawFH"
)


.fame_working_dirs <- c(
  "AdHocPolygons",
  "CustomCSV",
  "FH_Outputs",
  "results"
)


.fame_completion_file <- ".fame_s3_sync_complete"


.fame_lock_dir <- ".fame_s3_sync_lock"


# ------------------------------------------------------------
# Return the configured FAME data root
# ------------------------------------------------------------

#' Return the FAME external-data directory
#'
#' The function first checks the `FAME.data_dir` R option and
#' then the `FAME_DATA_DIR` environment variable.
#'
#' @param must_work Logical. Must the directory already exist?
#'
#' @return Normalised path to the FAME data directory.
#' @export
fame_data_dir <- function(must_work = TRUE) {

  data_dir <- getOption(
    "FAME.data_dir",
    ""
  )

  if (
    !is.character(data_dir) ||
    length(data_dir) != 1L ||
    is.na(data_dir) ||
    !nzchar(data_dir)
  ) {

    data_dir <- Sys.getenv(
      "FAME_DATA_DIR",
      unset = ""
    )
  }


  if (!nzchar(data_dir)) {

    stop(
      paste0(
        "FAME_DATA_DIR has not been configured.\n",
        "Run configure_fame_data_dir() or initialise_fame_data() ",
        "before accessing FAME data."
      ),
      call. = FALSE
    )
  }


  data_dir <- path.expand(data_dir)


  if (must_work && !dir.exists(data_dir)) {

    stop(
      "FAME data directory does not exist: ",
      data_dir,
      call. = FALSE
    )
  }


  normalizePath(
    data_dir,
    winslash = "/",
    mustWork = must_work
  )
}


# ------------------------------------------------------------
# Construct a path below FAME_DATA_DIR
# ------------------------------------------------------------

#' Construct a path within the FAME external-data directory
#'
#' @param ... Directory and filename components.
#' @param must_work Logical. Must the FAME data root exist?
#'
#' @return A platform-independent file path.
#' @export
fame_data_path <- function(..., must_work = TRUE) {

  file.path(
    fame_data_dir(must_work = must_work),
    ...
  )
}


# ------------------------------------------------------------
# Test whether a directory is writable
# ------------------------------------------------------------

.fame_directory_writable <- function(path) {

  if (!dir.exists(path)) {
    return(FALSE)
  }


  test_file <- tempfile(
    pattern = ".fame_write_test_",
    tmpdir = path
  )


  created <- tryCatch(
    file.create(test_file),
    warning = function(w) FALSE,
    error = function(e) FALSE
  )


  if (isTRUE(created)) {
    unlink(test_file, force = TRUE)
    return(TRUE)
  }


  FALSE
}


# ------------------------------------------------------------
# Choose and configure the local root directory
# ------------------------------------------------------------

#' Configure the FAME external-data directory
#'
#' An existing `FAME_DATA_DIR` setting always takes precedence.
#' On Windows, the default is `C:/data/FAME_data`. On Linux,
#' the default is a `FAME_data` directory beneath `app_dir`.
#'
#' @param app_dir Location of the running Shiny application.
#'   Required on Linux when `FAME_DATA_DIR` is unset.
#' @param data_dir Optional explicit local data directory.
#' @param create Logical. Create the directory if absent.
#'
#' @return Invisibly returns the normalised data-directory path.
#' @export
configure_fame_data_dir <- function(
    app_dir = NULL,
    data_dir = NULL,
    create = TRUE
) {

  # ----------------------------------------------------------
  # Explicit function argument has highest precedence
  # ----------------------------------------------------------

  if (
    !is.null(data_dir) &&
    length(data_dir) == 1L &&
    !is.na(data_dir) &&
    nzchar(data_dir)
  ) {

    selected_dir <- data_dir

  } else {

    # --------------------------------------------------------
    # Existing environment setting comes next
    # --------------------------------------------------------

    selected_dir <- Sys.getenv(
      "FAME_DATA_DIR",
      unset = ""
    )
  }


  # ----------------------------------------------------------
  # Select a default only when nothing was explicitly set
  # ----------------------------------------------------------

  if (!nzchar(selected_dir)) {

    if (.Platform$OS.type == "windows") {

      selected_dir <- "C:/data/FAME_data"

    } else {

      if (
        is.null(app_dir) ||
        length(app_dir) != 1L ||
        is.na(app_dir) ||
        !nzchar(app_dir)
      ) {

        stop(
          paste0(
            "On Linux, app_dir must be supplied when ",
            "FAME_DATA_DIR is unset."
          ),
          call. = FALSE
        )
      }


      app_dir <- normalizePath(
        app_dir,
        winslash = "/",
        mustWork = TRUE
      )


      selected_dir <- file.path(
        app_dir,
        "FAME_data"
      )
    }
  }


  selected_dir <- path.expand(selected_dir)


  # ----------------------------------------------------------
  # Create root when necessary
  # ----------------------------------------------------------

  if (!dir.exists(selected_dir) && create) {

    created <- dir.create(
      selected_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )


    if (!created && !dir.exists(selected_dir)) {

      stop(
        "Unable to create FAME data directory: ",
        selected_dir,
        call. = FALSE
      )
    }
  }


  if (!dir.exists(selected_dir)) {

    stop(
      "FAME data directory does not exist: ",
      selected_dir,
      call. = FALSE
    )
  }


  selected_dir <- normalizePath(
    selected_dir,
    winslash = "/",
    mustWork = TRUE
  )


  # ----------------------------------------------------------
  # Confirm real write access
  # ----------------------------------------------------------

  if (!.fame_directory_writable(selected_dir)) {

    stop(
      paste0(
        "The FAME data directory is not writable by the ",
        "current process:\n",
        selected_dir,
        "\n\n",
        "On Shiny Server, the directory must be writable by ",
        "the account running the Shiny application."
      ),
      call. = FALSE
    )
  }


  # ----------------------------------------------------------
  # Configure the current R process
  # ----------------------------------------------------------

  Sys.setenv(
    FAME_DATA_DIR = selected_dir
  )


  options(
    FAME.data_dir = selected_dir
  )


  invisible(selected_dir)
}


# ------------------------------------------------------------
# Determine whether a directory contains data
# ------------------------------------------------------------

.fame_directory_has_data <- function(path) {

  if (!dir.exists(path)) {
    return(FALSE)
  }


  contents <- list.files(
    path,
    all.files = TRUE,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )


  # Ignore FAME control files if encountered
  contents <- contents[
    !basename(contents) %in% c(
      .fame_completion_file,
      .fame_lock_dir
    )
  ]


  length(contents) > 0L
}


# ------------------------------------------------------------
# Report current data status
# ------------------------------------------------------------

#' Report the status of FAME external-data directories
#'
#' @return A data frame containing the path, role, existence,
#'   contents and readiness of each standard directory.
#' @export
fame_data_status <- function() {

  root <- fame_data_dir()


  directories <- c(
    .fame_required_data_dirs,
    .fame_working_dirs
  )


  role <- c(
    rep(
      "required input",
      length(.fame_required_data_dirs)
    ),
    rep(
      "working/output",
      length(.fame_working_dirs)
    )
  )


  paths <- file.path(
    root,
    directories
  )


  exists <- dir.exists(paths)


  contains_data <- vapply(
    paths,
    .fame_directory_has_data,
    logical(1)
  )


  required <- role == "required input"


  ready <- ifelse(
    required,
    exists & contains_data,
    exists
  )


  data.frame(
    directory = directories,
    role = role,
    path = paths,
    exists = exists,
    contains_data = contains_data,
    ready = ready,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Locate the AWS CLI
# ------------------------------------------------------------

.fame_aws_cli <- function() {

  aws <- Sys.which("aws")


  if (!nzchar(aws)) {

    stop(
      paste0(
        "The AWS CLI could not be found.\n",
        "It is required only because FAME reference data are ",
        "missing and must be downloaded from S3."
      ),
      call. = FALSE
    )
  }


  unname(aws)
}


# ------------------------------------------------------------
# Construct an S3 URI
# ------------------------------------------------------------

.fame_s3_uri <- function(
    bucket,
    prefix = "",
    directory = NULL
) {

  bucket <- sub(
    "^s3://",
    "",
    bucket
  )


  bucket <- sub(
    "/+$",
    "",
    bucket
  )


  prefix <- gsub(
    "^/+|/+$",
    "",
    prefix
  )


  parts <- c(
    paste0("s3://", bucket),
    if (nzchar(prefix)) prefix else NULL,
    directory
  )


  paste0(
    paste(parts, collapse = "/"),
    "/"
  )
}


# ------------------------------------------------------------
# Run an AWS CLI command
# ------------------------------------------------------------

.fame_run_aws <- function(args) {

  aws <- .fame_aws_cli()


  profile <- Sys.getenv(
    "FAME_AWS_PROFILE",
    unset = ""
  )


  if (nzchar(profile)) {

    args <- c(
      args,
      "--profile",
      profile
    )
  }


  output <- system2(
    command = aws,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )


  status <- attr(
    output,
    "status"
  )


  if (is.null(status)) {
    status <- 0L
  }


  list(
    status = as.integer(status),
    output = output,
    command = paste(
      c("aws", args),
      collapse = " "
    )
  )
}


# ------------------------------------------------------------
# Confirm that the configured bucket can be accessed
# ------------------------------------------------------------

.fame_check_s3_access <- function(
    bucket,
    prefix = ""
) {

  uri <- .fame_s3_uri(
    bucket = bucket,
    prefix = prefix
  )


  result <- .fame_run_aws(
    c(
      "s3",
      "ls",
      uri
    )
  )


  if (result$status != 0L) {

    stop(
      paste0(
        "FAME could not access the configured S3 location:\n",
        uri,
        "\n\nAWS CLI output:\n",
        paste(result$output, collapse = "\n")
      ),
      call. = FALSE
    )
  }


  invisible(TRUE)
}


# ------------------------------------------------------------
# Acquire a simple first-run lock
# ------------------------------------------------------------

.fame_acquire_sync_lock <- function(
    data_dir,
    timeout_seconds = 600,
    poll_seconds = 2
) {

  lock_path <- file.path(
    data_dir,
    .fame_lock_dir
  )


  completion_path <- file.path(
    data_dir,
    .fame_completion_file
  )


  started <- Sys.time()


  repeat {

    # dir.create() is used as the lock acquisition operation.
    acquired <- dir.create(
      lock_path,
      showWarnings = FALSE
    )


    if (isTRUE(acquired)) {

      owner_file <- file.path(
        lock_path,
        "owner.txt"
      )


      writeLines(
        c(
          paste("pid:", Sys.getpid()),
          paste("host:", Sys.info()[["nodename"]]),
          paste("started:", format(Sys.time(), tz = "UTC"))
        ),
        owner_file
      )


      return(lock_path)
    }


    # Another process may have completed initialisation.
    if (file.exists(completion_path)) {
      return(NULL)
    }


    elapsed <- as.numeric(
      difftime(
        Sys.time(),
        started,
        units = "secs"
      )
    )


    if (elapsed >= timeout_seconds) {

      stop(
        paste0(
          "Timed out waiting for another FAME process to ",
          "finish the initial S3 data synchronisation.\n",
          "Lock location: ",
          lock_path
        ),
        call. = FALSE
      )
    }


    Sys.sleep(poll_seconds)
  }
}


# ------------------------------------------------------------
# Release initialisation lock
# ------------------------------------------------------------

.fame_release_sync_lock <- function(lock_path) {

  if (
    !is.null(lock_path) &&
    dir.exists(lock_path)
  ) {

    unlink(
      lock_path,
      recursive = TRUE,
      force = TRUE
    )
  }


  invisible(NULL)
}


# ------------------------------------------------------------
# Download selected directories from S3
# ------------------------------------------------------------

#' Download FAME data directories from S3
#'
#' The AWS CLI must be installed. On EC2, authentication should
#' normally be provided by the instance IAM role.
#'
#' @param directories Directory names to synchronise.
#' @param destination Local FAME data directory.
#' @param bucket S3 bucket name. Defaults to `FAME_S3_BUCKET`.
#' @param prefix Optional key prefix. Defaults to
#'   `FAME_S3_PREFIX`.
#' @param delete Logical. Pass `--delete` to `aws s3 sync`.
#'   The default is FALSE to avoid removing local files.
#'
#' @return Invisibly returns TRUE.
#' @export
download_fame_data_from_s3 <- function(
    directories,
    destination = fame_data_dir(),
    bucket = Sys.getenv("FAME_S3_BUCKET", unset = ""),
    prefix = Sys.getenv("FAME_S3_PREFIX", unset = ""),
    delete = FALSE
) {

  if (!length(directories)) {
    return(invisible(TRUE))
  }


  if (!nzchar(bucket)) {

    stop(
      paste0(
        "FAME reference data are missing, but FAME_S3_BUCKET ",
        "has not been configured."
      ),
      call. = FALSE
    )
  }


  .fame_check_s3_access(
    bucket = bucket,
    prefix = prefix
  )


  for (directory in directories) {

    local_dir <- file.path(
      destination,
      directory
    )


    dir.create(
      local_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )


    source_uri <- .fame_s3_uri(
      bucket = bucket,
      prefix = prefix,
      directory = directory
    )


    message(
      "Synchronising ",
      source_uri,
      " -> ",
      local_dir
    )


    args <- c(
      "s3",
      "sync",
      source_uri,
      local_dir,
      "--only-show-errors"
    )


    if (isTRUE(delete)) {

      args <- c(
        args,
        "--delete"
      )
    }


    result <- .fame_run_aws(args)


    if (result$status != 0L) {

      stop(
        paste0(
          "Failed to synchronise FAME data directory '",
          directory,
          "' from S3.\n\n",
          "AWS CLI output:\n",
          paste(result$output, collapse = "\n")
        ),
        call. = FALSE
      )
    }


    if (!.fame_directory_has_data(local_dir)) {

      stop(
        paste0(
          "S3 synchronisation completed without an AWS CLI ",
          "error, but no files were found in:\n",
          local_dir,
          "\n\nExpected S3 source:\n",
          source_uri
        ),
        call. = FALSE
      )
    }
  }


  invisible(TRUE)
}


# ------------------------------------------------------------
# Write successful-initialisation marker
# ------------------------------------------------------------

.fame_write_completion_marker <- function(
    data_dir,
    bucket,
    prefix
) {

  marker <- file.path(
    data_dir,
    .fame_completion_file
  )


  writeLines(
    c(
      paste("completed_utc:", format(Sys.time(), tz = "UTC")),
      paste("bucket:", bucket),
      paste("prefix:", prefix),
      paste(
        "required_directories:",
        paste(
          .fame_required_data_dirs,
          collapse = ","
        )
      )
    ),
    marker
  )


  invisible(marker)
}


# ------------------------------------------------------------
# Complete FAME initialisation
# ------------------------------------------------------------

#' Initialise FAME external data
#'
#' Existing local data are retained and used. Missing required
#' input directories are synchronised from S3. Writable/output
#' directories are created locally and are not downloaded.
#'
#' @param app_dir Location of the running Shiny application.
#'   Used for the Linux default when FAME_DATA_DIR is unset.
#' @param data_dir Optional explicit local data directory.
#' @param download_missing Logical. Download missing required
#'   data from S3.
#' @param force_sync Logical. Synchronise all required input
#'   directories even when local data are present.
#' @param bucket S3 bucket name. Defaults to FAME_S3_BUCKET.
#' @param prefix Optional S3 prefix. Defaults to FAME_S3_PREFIX.
#' @param lock_timeout Maximum seconds to wait for another Shiny
#'   process performing first-run synchronisation.
#'
#' @return Invisibly returns the data status table.
#' @export
initialise_fame_data <- function(
    app_dir = NULL,
    data_dir = NULL,
    download_missing = TRUE,
    force_sync = FALSE,
    bucket = Sys.getenv("FAME_S3_BUCKET", unset = ""),
    prefix = Sys.getenv("FAME_S3_PREFIX", unset = ""),
    lock_timeout = 600
) {

  # ----------------------------------------------------------
  # 1. Configure local root
  # ----------------------------------------------------------

  root <- configure_fame_data_dir(
    app_dir = app_dir,
    data_dir = data_dir,
    create = TRUE
  )


  message(
    "FAME data directory: ",
    root
  )


  # ----------------------------------------------------------
  # 2. Create working/output directories
  # ----------------------------------------------------------

  for (directory in .fame_working_dirs) {

    path <- file.path(
      root,
      directory
    )


    dir.create(
      path,
      recursive = TRUE,
      showWarnings = FALSE
    )


    if (!dir.exists(path)) {

      stop(
        "Unable to create FAME working directory: ",
        path,
        call. = FALSE
      )
    }
  }


  # ----------------------------------------------------------
  # 3. Determine missing required inputs
  # ----------------------------------------------------------

  status <- fame_data_status()


  required_status <- status[
    status$role == "required input",
    ,
    drop = FALSE
  ]


  if (isTRUE(force_sync)) {

    directories_to_sync <- required_status$directory

  } else {

    directories_to_sync <- required_status$directory[
      !required_status$ready
    ]
  }


  # ----------------------------------------------------------
  # 4. Use existing local data when complete
  # ----------------------------------------------------------

  if (!length(directories_to_sync)) {

    message(
      "Required FAME reference data are already available locally."
    )


    return(
      invisible(status)
    )
  }


  # ----------------------------------------------------------
  # 5. Missing data without download permission
  # ----------------------------------------------------------

  if (!isTRUE(download_missing)) {

    stop(
      paste0(
        "Required FAME data are missing or empty:\n  ",
        paste(
          directories_to_sync,
          collapse = "\n  "
        ),
        "\n\nSet download_missing = TRUE or populate these ",
        "directories manually."
      ),
      call. = FALSE
    )
  }


  if (!nzchar(bucket)) {

    stop(
      paste0(
        "Required FAME data are missing or empty:\n  ",
        paste(
          directories_to_sync,
          collapse = "\n  "
        ),
        "\n\nFAME_S3_BUCKET has not been configured."
      ),
      call. = FALSE
    )
  }


  # ----------------------------------------------------------
  # 6. Prevent simultaneous first-run syncs
  # ----------------------------------------------------------

  lock_path <- .fame_acquire_sync_lock(
    data_dir = root,
    timeout_seconds = lock_timeout
  )


  # NULL means another process completed while this process was
  # waiting. Recheck the data before proceeding.
  if (is.null(lock_path)) {

    status <- fame_data_status()


    if (!all(
      status$ready[
        status$role == "required input"
      ]
    )) {

      stop(
        paste0(
          "Another FAME process completed initialisation, but ",
          "some required data are still missing."
        ),
        call. = FALSE
      )
    }


    return(
      invisible(status)
    )
  }


  on.exit(
    .fame_release_sync_lock(lock_path),
    add = TRUE
  )


  # ----------------------------------------------------------
  # 7. Recheck after acquiring the lock
  # ----------------------------------------------------------

  status <- fame_data_status()


  required_status <- status[
    status$role == "required input",
    ,
    drop = FALSE
  ]


  if (isTRUE(force_sync)) {

    directories_to_sync <- required_status$directory

  } else {

    directories_to_sync <- required_status$directory[
      !required_status$ready
    ]
  }


  # ----------------------------------------------------------
  # 8. Download only what is missing
  # ----------------------------------------------------------

  if (length(directories_to_sync)) {

    message(
      "Downloading missing FAME data from S3: ",
      paste(
        directories_to_sync,
        collapse = ", "
      )
    )


    download_fame_data_from_s3(
      directories = directories_to_sync,
      destination = root,
      bucket = bucket,
      prefix = prefix,
      delete = FALSE
    )
  }


  # ----------------------------------------------------------
  # 9. Validate final state
  # ----------------------------------------------------------

  final_status <- fame_data_status()


  missing_after_sync <- final_status$directory[
    final_status$role == "required input" &
      !final_status$ready
  ]


  if (length(missing_after_sync)) {

    stop(
      paste0(
        "FAME data initialisation did not produce all required ",
        "data directories:\n  ",
        paste(
          missing_after_sync,
          collapse = "\n  "
        )
      ),
      call. = FALSE
    )
  }


  # ----------------------------------------------------------
  # 10. Record successful S3 initialisation
  # ----------------------------------------------------------

  .fame_write_completion_marker(
    data_dir = root,
    bucket = bucket,
    prefix = prefix
  )


  message(
    "FAME external-data initialisation completed successfully."
  )


  invisible(final_status)
}

