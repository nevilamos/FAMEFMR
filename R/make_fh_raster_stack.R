make_fh_raster_stack <- function(
    r_template = template_r,
    vector_file = infile,
    fields = c("SEASON", "FIRETYPE_NO"),
    combField_multiplier = 10L,

    background = 0,
    overwrite = TRUE,
    datatype = "INT2U",
    compress = "LZW",
    memfrac = NULL,
    progress = TRUE,
    gc_every = 0,
    OtherAndUnknown = 2
) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required.")
  }
  stopifnot(length(fields) == 2)

  # --- template checks
  if (!inherits(r_template, "SpatRaster")) {
    stop("r_template must be a terra::SpatRaster.")
  }

  # --- terra memory option (optional)
  if (!is.null(memfrac)) {
    terra::terraOptions(memfrac = memfrac)
  }

  # --- output dirs
  out_dir<-file.path(tempdir(),"FIRE_RASTER_STACK")
  if (dir.exists(out_dir)) {
    unlink(out_dir, recursive = TRUE)
  }
  dir.create(out_dir,     recursive = TRUE, showWarnings = FALSE)

  # --- read vector
  v <- terra::vect(vector_file)
  nms <- names(v)

  # --- field checks
  if ("FIRETYPE" %in% nms) {
  #Recode FIRETYPE -> FIRETYPE_NO
  v$FIRETYPE_NO[v$FIRETYPE == "BURN"]     <- 1L
  v$FIRETYPE_NO[v$FIRETYPE == "BUSHFIRE"] <- 2L
  v$FIRETYPE_NO[v$FIRETYPE %in% c("OTHER", "UNKNOWN")] <- OtherAndUnknown
  nms <- names(v)
  }



  if (!all(fields %in% nms)) {
    missing <- fields[!fields %in% nms]
    stop(sprintf("Missing field(s) in vector attributes: %s", paste(missing, collapse = ", ")))
  }



  # --- compute combFields = SEASON*100 + FIRETYPE_NO
  mult <- as.integer(combField_multiplier)
  if (is.na(mult) || mult <= 0L) stop("combField_multiplier must be a positive integer.")

  v$combFields <- as.integer(v$SEASON * mult + v$FIRETYPE_NO)

  # --- split once (fast)
  vlist <- split(v[,"combFields"], v$combFields)
  group_names <- names(vlist)

  # --- progress
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = length(vlist), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  # --- write options
  wopt <- list(gdal = sprintf("COMPRESS=%s", compress), datatype = datatype)

  # --- main loop
  for (k in seq_along(group_names)) {
    nm <- group_names[[k]]

    # Rasterize
    rr<-terra::rasterize(
      x = vlist[[nm]],
      y = r_template,
      field = "combFields",
      background = background)
    plot(rr)
    writeRaster(x = rr,
      filename = file.path(out_dir, paste0(nm, ".tif")),
      overwrite = TRUE,
      wopt = wopt
    )
    rm(rr)



    if (progress) utils::setTxtProgressBar(pb, k)

    if (gc_every > 0 && (k %% gc_every) == 0) gc()
  }

  invisible(list(
    groups = group_names,
    n_groups = length(group_names),
    out_dir_season = out_dir_season,
    out_dir_ft = out_dir_ft,
    fields = fields,
    combField_multiplier = mult
  ))
}
