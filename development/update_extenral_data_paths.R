# ============================================================
# Update FAMEshiny repository-relative external-data paths
#
# Repository:
#   C:/data/FAMEFMR
#
# This script:
#   1. Scans scripts under inst/shinyapp.
#   2. Excludes files beneath directories named "results".
#   3. Converts known repository-relative data paths into calls
#      to fame_data_path().
#   4. Creates R/fame_data_paths.R.
#   5. Creates backups before modifying files.
#   6. Produces reports of automatic changes and remaining
#      references requiring manual review.
#
# Recommended workflow:
#   1. Commit or stash existing changes.
#   2. Run with dry_run <- TRUE.
#   3. Inspect the reports.
#   4. Change dry_run to FALSE.
#   5. Run the script again.
#   6. Inspect git diff before committing anything.
#
# ============================================================


# ------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------

repo_dir <- "C:/data/FAMEFMR"

app_dir <- file.path(
  repo_dir,
  "inst",
  "shinyapp"
)

development_dir <- file.path(
  repo_dir,
  "development"
)

backup_dir <- file.path(
  development_dir,
  "path_migration_backup"
)

change_report_file <- file.path(
  development_dir,
  "external_path_changes.csv"
)

review_report_file <- file.path(
  development_dir,
  "external_path_manual_review.csv"
)

excluded_report_file <- file.path(
  development_dir,
  "external_path_excluded_files.csv"
)

helper_file <- file.path(
  repo_dir,
  "R",
  "fame_data_paths.R"
)


# ------------------------------------------------------------
# Preview or apply changes
# ------------------------------------------------------------

# TRUE:
#   Report proposed changes but do not change repository files.
#
# FALSE:
#   Back up and modify repository files.

dry_run <-FALSE


# ------------------------------------------------------------
# Known external-data subdirectories
# ------------------------------------------------------------

external_dirs <- c(
  "AdHocPolygons",
  "CustomCSV",
  "FH_Outputs",
  "HDMS",
  "InputGeneralRasters",
  "ReferenceShapefiles",
  "ReferenceTables",
  "rawFH"
)


# ------------------------------------------------------------
# Directory names to exclude
# ------------------------------------------------------------

excluded_dirs <- c(
  "results"
)


# ------------------------------------------------------------
# File types to process
# ------------------------------------------------------------

file_pattern <- "\\.(R|r|Rmd|rmd)$"


# ------------------------------------------------------------
# 2. Check repository structure
# ------------------------------------------------------------

if (!dir.exists(repo_dir)) {
  stop(
    "Repository directory does not exist: ",
    repo_dir,
    call. = FALSE
  )
}

if (!dir.exists(app_dir)) {
  stop(
    "Embedded Shiny application does not exist: ",
    app_dir,
    call. = FALSE
  )
}

if (!dir.exists(file.path(repo_dir, ".git"))) {
  warning(
    "No .git directory was found at: ",
    repo_dir,
    "\nConfirm that this is the root of the FAMEFMR repository."
  )
}

dir.create(
  development_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Package path-helper code
# ------------------------------------------------------------

helper_code <- c(
  "#' Return the root FAME data directory",
  "#'",
  "#' The directory is obtained first from the `FAME.data_dir`",
  "#' option and then from the `FAME_DATA_DIR` environment",
  "#' variable.",
  "#'",
  "#' @param must_work Logical. Should the directory already exist?",
  "#'",
  "#' @return Normalised path to the external FAME data directory.",
  "#' @export",
  "fame_data_dir <- function(must_work = TRUE) {",
  "",
  "  data_dir <- getOption(",
  "    \"FAME.data_dir\",",
  "    Sys.getenv(\"FAME_DATA_DIR\", unset = \"\")",
  "  )",
  "",
  "  if (!is.character(data_dir) ||",
  "      length(data_dir) != 1L ||",
  "      is.na(data_dir) ||",
  "      !nzchar(data_dir)) {",
  "",
  "    stop(",
  "      paste0(",
  "        \"The FAME external data directory has not been configured.\\n\",",
  "        \"Set it with either:\\n\\n\",",
  "        \"  Sys.setenv(FAME_DATA_DIR = 'C:/path/to/FAME_data')\\n\\n\",",
  "        \"or:\\n\\n\",",
  "        \"  options(FAME.data_dir = 'C:/path/to/FAME_data')\"",
  "      ),",
  "      call. = FALSE",
  "    )",
  "  }",
  "",
  "  data_dir <- path.expand(data_dir)",
  "",
  "  if (must_work && !dir.exists(data_dir)) {",
  "    stop(",
  "      \"FAME_DATA_DIR does not exist: \",",
  "      data_dir,",
  "      call. = FALSE",
  "    )",
  "  }",
  "",
  "  normalizePath(",
  "    data_dir,",
  "    winslash = \"/\",",
  "    mustWork = must_work",
  "  )",
  "}",
  "",
  "",
  "#' Construct a path within the FAME data directory",
  "#'",
  "#' @param ... Directory names or filenames below FAME_DATA_DIR.",
  "#' @param must_work Logical. Should the root data directory exist?",
  "#'",
  "#' @return A platform-independent file path.",
  "#' @export",
  "fame_data_path <- function(..., must_work = TRUE) {",
  "",
  "  file.path(",
  "    fame_data_dir(must_work = must_work),",
  "    ...",
  "  )",
  "}"
)


# ------------------------------------------------------------
# 4. Utility functions
# ------------------------------------------------------------

escape_regex <- function(x) {

  gsub(
    "([][{}()+*^$|\\\\?.])",
    "\\\\\\1",
    x
  )
}


normalise_file_path <- function(x) {

  normalizePath(
    x,
    winslash = "/",
    mustWork = FALSE
  )
}


relative_to <- function(path, root) {

  path <- normalise_file_path(path)
  root <- normalise_file_path(root)

  sub(
    paste0(
      "^",
      escape_regex(root),
      "/?"
    ),
    "",
    path
  )
}


is_excluded_file <- function(path, excluded_names) {

  path <- normalise_file_path(path)

  excluded_pattern <- paste0(
    "(^|/)(",
    paste(
      escape_regex(excluded_names),
      collapse = "|"
    ),
    ")(/|$)"
  )

  grepl(
    excluded_pattern,
    path,
    ignore.case = TRUE,
    perl = TRUE
  )
}


add_change <- function(
    report,
    file,
    line_number,
    change_type,
    original,
    replacement
) {

  if (identical(original, replacement)) {
    return(report)
  }

  rbind(
    report,
    data.frame(
      file = file,
      line = as.integer(line_number),
      change_type = change_type,
      original = original,
      replacement = replacement,
      stringsAsFactors = FALSE
    )
  )
}


# ------------------------------------------------------------
# 5. Conversion functions
# ------------------------------------------------------------

convert_file_path_calls <- function(line, directory) {

  escaped_dir <- escape_regex(directory)


  # file.path("HDMS", ...)
  double_pattern <- paste0(
    "file\\.path\\(\\s*\"(?:\\./)?",
    escaped_dir,
    "\"\\s*,"
  )

  double_replacement <- paste0(
    "fame_data_path(\"",
    directory,
    "\","
  )

  line <- gsub(
    double_pattern,
    double_replacement,
    line,
    perl = TRUE,
    ignore.case = TRUE
  )


  # file.path('HDMS', ...)
  single_pattern <- paste0(
    "file\\.path\\(\\s*'(?:\\./)?",
    escaped_dir,
    "'\\s*,"
  )

  single_replacement <- paste0(
    "fame_data_path('",
    directory,
    "',"
  )

  line <- gsub(
    single_pattern,
    single_replacement,
    line,
    perl = TRUE,
    ignore.case = TRUE
  )


  # file.path("HDMS")
  only_double_pattern <- paste0(
    "file\\.path\\(\\s*\"(?:\\./)?",
    escaped_dir,
    "\"\\s*\\)"
  )

  line <- gsub(
    only_double_pattern,
    paste0(
      "fame_data_path(\"",
      directory,
      "\")"
    ),
    line,
    perl = TRUE,
    ignore.case = TRUE
  )


  # file.path('HDMS')
  only_single_pattern <- paste0(
    "file\\.path\\(\\s*'(?:\\./)?",
    escaped_dir,
    "'\\s*\\)"
  )

  line <- gsub(
    only_single_pattern,
    paste0(
      "fame_data_path('",
      directory,
      "')"
    ),
    line,
    perl = TRUE,
    ignore.case = TRUE
  )

  line
}


convert_quoted_paths <- function(line, directory) {

  escaped_dir <- escape_regex(directory)


  # ----------------------------------------------------------
  # Avoid reprocessing lines already converted for this
  # directory.
  # ----------------------------------------------------------

  already_converted_double <- paste0(
    "fame_data_path\\(\\s*\"",
    escaped_dir,
    "\""
  )

  already_converted_single <- paste0(
    "fame_data_path\\(\\s*'",
    escaped_dir,
    "'"
  )

  if (
    grepl(
      already_converted_double,
      line,
      ignore.case = TRUE,
      perl = TRUE
    ) ||
    grepl(
      already_converted_single,
      line,
      ignore.case = TRUE,
      perl = TRUE
    )
  ) {
    return(line)
  }


  # ----------------------------------------------------------
  # "HDMS/file.tif"
  # "./HDMS/file.tif"
  # "HDMS\\file.tif"
  #
  # becomes:
  #
  # fame_data_path("HDMS", "file.tif")
  # ----------------------------------------------------------

  double_path_pattern <- paste0(
    "\"(?:\\./)?",
    escaped_dir,
    "[/\\\\]+([^\"\\r\\n]+)\""
  )

  double_path_replacement <- paste0(
    "fame_data_path(\"",
    directory,
    "\", \"\\1\")"
  )

  line <- gsub(
    double_path_pattern,
    double_path_replacement,
    line,
    perl = TRUE,
    ignore.case = TRUE
  )


  # ----------------------------------------------------------
  # 'HDMS/file.tif'
  # './HDMS/file.tif'
  #
  # becomes:
  #
  # fame_data_path('HDMS', 'file.tif')
  # ----------------------------------------------------------

  single_path_pattern <- paste0(
    "'(?:\\./)?",
    escaped_dir,
    "[/\\\\]+([^'\\r\\n]+)'"
  )

  single_path_replacement <- paste0(
    "fame_data_path('",
    directory,
    "', '\\1')"
  )

  line <- gsub(
    single_path_pattern,
    single_path_replacement,
    line,
    perl = TRUE,
    ignore.case = TRUE
  )


  # ----------------------------------------------------------
  # Directory-only string:
  #
  # "HDMS"
  # "./HDMS"
  #
  # becomes:
  #
  # fame_data_path("HDMS")
  # ----------------------------------------------------------

  directory_double_pattern <- paste0(
    "\"(?:\\./)?",
    escaped_dir,
    "[/\\\\]?\""
  )

  directory_double_replacement <- paste0(
    "fame_data_path(\"",
    directory,
    "\")"
  )

  line <- gsub(
    directory_double_pattern,
    directory_double_replacement,
    line,
    perl = TRUE,
    ignore.case = TRUE
  )


  directory_single_pattern <- paste0(
    "'(?:\\./)?",
    escaped_dir,
    "[/\\\\]?'"
  )

  directory_single_replacement <- paste0(
    "fame_data_path('",
    directory,
    "')"
  )

  line <- gsub(
    directory_single_pattern,
    directory_single_replacement,
    line,
    perl = TRUE,
    ignore.case = TRUE
  )

  line
}


convert_line <- function(line) {

  updated <- line

  for (directory in external_dirs) {

    # Convert complete file.path() expressions first.
    updated <- convert_file_path_calls(
      updated,
      directory
    )

    # Then convert directly quoted paths.
    updated <- convert_quoted_paths(
      updated,
      directory
    )
  }

  updated
}


# ------------------------------------------------------------
# 6. Locate candidate scripts
# ------------------------------------------------------------

all_files <- list.files(
  app_dir,
  pattern = file_pattern,
  recursive = TRUE,
  full.names = TRUE
)


# ------------------------------------------------------------
# Exclude every file beneath a directory named "results"
# ------------------------------------------------------------

excluded <- is_excluded_file(
  all_files,
  excluded_dirs
)

excluded_files <- all_files[excluded]
files <- all_files[!excluded]


# ------------------------------------------------------------
# Record excluded files
# ------------------------------------------------------------

excluded_files_report <- data.frame(
  file = relative_to(
    excluded_files,
    app_dir
  ),
  reason = rep(
    "File is beneath an excluded results directory",
    length(excluded_files)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  excluded_files_report,
  excluded_report_file,
  row.names = FALSE
)


cat(
  "\nFAME external-path migration\n",
  "----------------------------\n",
  "Repository:                    ",
  repo_dir,
  "\n",
  "Application:                   ",
  app_dir,
  "\n",
  "Candidate files found:         ",
  length(all_files),
  "\n",
  "Files excluded under results:  ",
  length(excluded_files),
  "\n",
  "Files to process:              ",
  length(files),
  "\n",
  "Dry run:                       ",
  dry_run,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Prepare change report
# ------------------------------------------------------------

change_report <- data.frame(
  file = character(),
  line = integer(),
  change_type = character(),
  original = character(),
  replacement = character(),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 8. Process each application script
# ------------------------------------------------------------

for (file in files) {

  relative_file <- relative_to(
    file,
    app_dir
  )

  original_lines <- tryCatch(
    readLines(
      file,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    error = function(e) {

      warning(
        "Could not read ",
        file,
        ": ",
        conditionMessage(e)
      )

      character()
    }
  )

  if (!length(original_lines)) {
    next
  }

  updated_lines <- original_lines

  for (i in seq_along(original_lines)) {

    updated_lines[i] <- convert_line(
      original_lines[i]
    )

    if (!identical(
      original_lines[i],
      updated_lines[i]
    )) {

      change_report <- add_change(
        report = change_report,
        file = relative_file,
        line_number = i,
        change_type = "External path conversion",
        original = original_lines[i],
        replacement = updated_lines[i]
      )
    }
  }


  # ----------------------------------------------------------
  # Back up and write changed files only when dry_run is FALSE
  # ----------------------------------------------------------

  if (
    !dry_run &&
    !identical(original_lines, updated_lines)
  ) {

    backup_file <- file.path(
      backup_dir,
      relative_file
    )

    dir.create(
      dirname(backup_file),
      recursive = TRUE,
      showWarnings = FALSE
    )

    copied <- file.copy(
      from = file,
      to = backup_file,
      overwrite = TRUE,
      copy.date = TRUE
    )

    if (!copied) {
      stop(
        "Failed to create backup for: ",
        file,
        call. = FALSE
      )
    }

    writeLines(
      updated_lines,
      file,
      useBytes = TRUE
    )
  }
}


# ------------------------------------------------------------
# 9. Create or update package data-path helpers
# ------------------------------------------------------------

helper_change_required <- TRUE

if (file.exists(helper_file)) {

  existing_helper <- readLines(
    helper_file,
    warn = FALSE,
    encoding = "UTF-8"
  )

  helper_change_required <- !identical(
    existing_helper,
    helper_code
  )
}


if (helper_change_required) {

  change_report <- add_change(
    report = change_report,
    file = relative_to(
      helper_file,
      repo_dir
    ),
    line_number = 1L,
    change_type = "Create or update package helper",
    original = if (file.exists(helper_file)) {
      "Existing helper file"
    } else {
      "File does not exist"
    },
    replacement = paste(
      "Create or update",
      relative_to(
        helper_file,
        repo_dir
      )
    )
  )


  if (!dry_run) {

    dir.create(
      dirname(helper_file),
      recursive = TRUE,
      showWarnings = FALSE
    )


    if (file.exists(helper_file)) {

      helper_backup <- file.path(
        backup_dir,
        "R",
        basename(helper_file)
      )

      dir.create(
        dirname(helper_backup),
        recursive = TRUE,
        showWarnings = FALSE
      )

      copied <- file.copy(
        from = helper_file,
        to = helper_backup,
        overwrite = TRUE,
        copy.date = TRUE
      )

      if (!copied) {
        stop(
          "Failed to back up: ",
          helper_file,
          call. = FALSE
        )
      }
    }


    writeLines(
      helper_code,
      helper_file,
      useBytes = TRUE
    )
  }
}


# ------------------------------------------------------------
# 10. Save automatic-change report
# ------------------------------------------------------------

if (nrow(change_report)) {

  change_report <- change_report[
    order(
      change_report$file,
      change_report$line
    ),
  ]

  rownames(change_report) <- NULL
}


write.csv(
  change_report,
  change_report_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Scan proposed or updated text for manual-review cases
#
# Only the filtered "files" object is used here. Therefore,
# results directories remain excluded from manual review.
# ------------------------------------------------------------

review_results <- lapply(
  files,
  function(file) {

    original <- tryCatch(
      readLines(
        file,
        warn = FALSE,
        encoding = "UTF-8"
      ),
      error = function(e) character()
    )

    if (!length(original)) {
      return(NULL)
    }


    # In dry-run mode, scan the proposed converted text.
    # Otherwise, scan the actual updated text.

    txt <- if (dry_run) {

      vapply(
        original,
        convert_line,
        character(1)
      )

    } else {

      original
    }


    directory_match <- vapply(
      txt,
      function(line) {

        hits <- external_dirs[
          vapply(
            external_dirs,
            function(directory) {

              grepl(
                directory,
                line,
                fixed = TRUE,
                ignore.case = TRUE
              )
            },
            logical(1)
          )
        ]

        if (!length(hits)) {
          return(NA_character_)
        }

        paste(
          hits,
          collapse = "; "
        )
      },
      character(1)
    )


    # Retain directory references that were not converted to
    # fame_data_path(). These commonly identify:
    #
    #   * dynamically assembled paths;
    #   * paste() or paste0() paths;
    #   * paths stored in variables;
    #   * comments or labels;
    #   * cases requiring manual judgement.

    remaining <- !is.na(directory_match) &
      !grepl(
        "fame_data_path\\(",
        txt,
        perl = TRUE
      )


    if (!any(remaining)) {
      return(NULL)
    }


    data.frame(
      file = relative_to(
        file,
        app_dir
      ),
      line = which(remaining),
      directory = directory_match[remaining],
      code = trimws(txt[remaining]),
      stringsAsFactors = FALSE
    )
  }
)


review_results <- review_results[
  lengths(review_results) > 0
]


if (length(review_results)) {

  manual_review <- do.call(
    rbind,
    review_results
  )

  manual_review <- manual_review[
    order(
      manual_review$directory,
      manual_review$file,
      manual_review$line
    ),
  ]

  rownames(manual_review) <- NULL

} else {

  manual_review <- data.frame(
    file = character(),
    line = integer(),
    directory = character(),
    code = character(),
    stringsAsFactors = FALSE
  )
}


write.csv(
  manual_review,
  review_report_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Print summary
# ------------------------------------------------------------

cat(
  "\nMigration scan complete\n",
  "-----------------------\n",
  "Files processed:             ",
  length(files),
  "\n",
  "Files excluded:              ",
  length(excluded_files),
  "\n",
  "Proposed changes:            ",
  nrow(change_report),
  "\n",
  "Manual-review records:       ",
  nrow(manual_review),
  "\n\n",
  "Change report:\n  ",
  change_report_file,
  "\n\n",
  "Manual-review report:\n  ",
  review_report_file,
  "\n\n",
  "Excluded-files report:\n  ",
  excluded_report_file,
  "\n",
  sep = ""
)


if (dry_run) {

  cat(
    "\nNo repository files were changed because dry_run is TRUE.\n",
    "\nInspect the reports, then change:\n\n",
    "  dry_run <- TRUE\n\n",
    "to:\n\n",
    "  dry_run <- FALSE\n\n",
    "and run this script again.\n",
    sep = ""
  )

} else {

  cat(
    "\nRepository files were updated.\n",
    "\nOriginal versions of changed files are stored under:\n  ",
    backup_dir,
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 13. Display reports in RStudio
# ------------------------------------------------------------

if (interactive()) {

  View(change_report)

  if (nrow(manual_review)) {
    View(manual_review)
  }

  if (nrow(excluded_files_report)) {
    View(excluded_files_report)
  }
}
