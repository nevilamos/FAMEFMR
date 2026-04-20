#' Create per-(SEASON × FIRETYPE_NO) fire-history raster layers (terra)
#'
#' @description
#' Builds a set of raster layers from a polygon/feature vector file by:
#' \enumerate{
#'   \item Reading the vector file with \pkg{terra}.
#'   \item Ensuring the required fields exist (default: \code{SEASON} and \code{FIRETYPE_NO}).
#'   \item Optionally recoding a text \code{FIRETYPE} field into \code{FIRETYPE_NO} if \code{FIRETYPE}
#'         is present (BURN = 1, BUSHFIRE = 2, OTHER/UNKNOWN = \code{OtherAndUnknown}).
#'   \item Creating a combined integer field \code{combFields = SEASON * combField_multiplier + FIRETYPE_NO}.
#'   \item Splitting features by \code{combFields} and rasterizing each group onto the template raster.
#'   \item Writing each raster layer as an individual GeoTIFF into a temporary output directory.
#' }
#'
#' The function uses a template \code{SpatRaster} to define the target grid (extent/resolution/CRS).
#'
#' @param r_template A \code{\link[terra:SpatRaster-class]{terra::SpatRaster}} template defining
#'   output extent, resolution and CRS.
#' @param vector_file Character. Path to the input vector dataset (e.g., \code{.shp}, \code{.gpkg})
#'   readable by \code{\link[terra:vect]{terra::vect}}.
#' @param fields Character vector of length 2 giving the names of the fields used to form the
#'   combined code. Defaults to \code{c("SEASON","FIRETYPE_NO")}.
#' @param combField_multiplier Integer scalar. Multiplier used when combining the two fields:
#'   \code{combFields = SEASON * combField_multiplier + FIRETYPE_NO}. Defaults to \code{10L}.
#'   Ensure this multiplier is large enough that \code{FIRETYPE_NO < combField_multiplier} to avoid
#'   code collisions.
#' @param background Numeric. Background value for cells not covered by any features. Defaults to 0.
#' @param overwrite Logical. If \code{TRUE}, overwrite existing GeoTIFFs in the output directory.
#'   Defaults to \code{TRUE}.
#' @param datatype Character. GDAL datatype for output rasters passed via \code{wopt}. Defaults to
#'   \code{"INT2U"}.
#' @param compress Character. GDAL compression option (e.g., \code{"LZW"}). Defaults to \code{"LZW"}.
#' @param memfrac Numeric in (0,1], optional. If supplied, passed to
#'   \code{\link[terra:terraOptions]{terra::terraOptions}} as \code{memfrac} to control memory use.
#' @param progress Logical. If \code{TRUE}, show a \code{\link[utils:txtProgressBar]{txtProgressBar}}.
#'   Defaults to \code{TRUE}.
#' @param OtherAndUnknown Integer. Value used when recoding \code{FIRETYPE} levels OTHER/UNKNOWN
#'   into \code{FIRETYPE_NO}. Defaults to 2.
#'
#' @details
#' \strong{Output directory:} The function writes GeoTIFFs into
#' \code{file.path(tempdir(), "FIRE_RASTER_STACK")}. If the directory already exists, it is deleted
#' and recreated.
#'
#' \strong{Recoding FIRETYPE:} If a field named \code{FIRETYPE} is present, the function writes
#' into \code{FIRETYPE_NO} (creating/overwriting values for matched rows).
#'
#' \strong{Performance note:} The line \code{plot(rr)} is retained from the supplied code and will
#' substantially slow large runs. Consider removing it or gating it behind a flag for production.
#'
#' @return Invisibly returns a list with:
#' \describe{
#'   \item{groups}{Character vector of \code{combFields} group labels (as returned by \code{split()}).}
#'   \item{n_groups}{Number of groups/layers written.}
#'   \item{out_dir}{Output directory containing the GeoTIFF files.}
#'   \item{fields}{The \code{fields} argument used.}
#'   \item{combField_multiplier}{The integer multiplier used to compute \code{combFields}.}
#' }
#'
#' @examples
#' \dontrun{
#' library(terra)
#' v <- vect("FH_State.shp")
#' r <- rast(ext = ext(v), res = 75, crs = crs(v))
#'
#' make_fh_raster_stack(
#'   r_template = r,
#'   vector_file = "FH_State.shp",
#'   fields = c("SEASON", "FIRETYPE_NO"),
#'   combField_multiplier = 100L,
#'   datatype = "INT2U"
#' )
#' }
#'
#' @importFrom terra vect rast rasterize writeRaster terraOptions
#' @export
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

  # --- read vector
  v <- terra::vect(vector_file)
  nms <- names(v)

  # --- field checks / optional recode
  if ("FIRETYPE" %in% nms) {
    # Recode FIRETYPE -> FIRETYPE_NO
    v$FIRETYPE_NO[v$FIRETYPE == "BURN"]     <- 1L
    v$FIRETYPE_NO[v$FIRETYPE == "BUSHFIRE"] <- 2L
    v$FIRETYPE_NO[v$FIRETYPE %in% c("OTHER", "UNKNOWN")] <- OtherAndUnknown
    nms <- names(v)
  }

  if (!all(fields %in% nms)) {
    missing <- fields[!fields %in% nms]
    stop(sprintf("Missing field(s) in vector attributes: %s", paste(missing, collapse = ", ")))
  }

  # --- compute combFields = SEASON*mult + FIRETYPE_NO
  mult <- as.integer(combField_multiplier)
  if (is.na(mult) || mult <= 0L) stop("combField_multiplier must be a positive integer.")

  v$combFields <- as.integer(v$SEASON * mult + v$FIRETYPE_NO)

  # --- split once (fast)
  vlist <- split(v[, "combFields"], v$combFields)
  group_names <- names(vlist)

  # --- progress
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = length(vlist), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  # --- write options
  wopt <- list(gdal = sprintf("COMPRESS=%s", compress), datatype = datatype)

  # --- main loop
  stackfiles<-character(0)
  for (k in seq_along(group_names)) {
    nm <- group_names[[k]]

    # Rasterize
    rr <- terra::rasterize(
      x = vlist[[nm]],
      y = r_template,
      field = "combFields",
      background = background
    )
    names(rr)<-nm
    f<-tempfile(pattern = paste0(nm, "_"), fileext = ".tif")
    stackfiles<-c(stackfiles,f)


    terra::writeRaster(
      x = rr,
      filename = f,
      overwrite = overwrite,
      wopt = wopt
    )

    rm(rr)

    if (progress) utils::setTxtProgressBar(pb, k)
    }

  # invisible(list(
  #   groups = group_names,
  #   n_groups = length(group_names),
  #   out_dir = out_dir,
  #   fields = fields,
  #   combField_multiplier = mult,
  #   stackfiles =stackfiles
  # ))

  s <- terra::rast(stackfiles)

  out_file <- tempfile(fileext = ".tif")
  out <- terra::writeRaster(s, out_file, overwrite = TRUE)

  # now safe to remove intermediates
  unlink(stackfiles)



out
}
