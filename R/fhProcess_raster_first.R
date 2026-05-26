#' Fire history raster workflow (wrapper): rasterize -> pack -> rebuild -> analyse
#'
#' @description
#' High-level wrapper that runs an end-to-end fire history raster workflow using
#' a raster template and a polygon vector file. The function:
#' \enumerate{
#'   \item reads the vector fire-history layer with [terra::vect()]
#'   \item generates a fire-history raster stack (via `make_fh_raster_stack()`)
#'   \item block-packs the raster stack to a binary representation
#'     (via `process_blocks_to_bin()`)
#'   \item reconstructs a packed raster from the binary files
#'     (via `write_raster_from_bin()`)
#'   \item assigns unique fire-history IDs, splits encoded integers into
#'     season and fire-type matrices (via `split_integer_first_last()`),
#'     then
#'   \item creates a list of outputs for downstream processing in FAME
#'     `add_fire_lft_lby_ysf()`.
#' }
#'
#' @details
#' This is a **workflow wrapper** around several lower-level
#' FAMEFMR utilities. It expects that the internal helper
#' functions `make_fh_raster_stack()`, `process_blocks_to_bin()`,
#' `write_raster_from_bin()`, `split_integer_first_last()`, and
#' `add_fire_lft_lby_ysf()` are available in the namespace.
#'
#'
#' **Encoding assumption**: the reconstructed raster values are assumed to be
#' integers that encode season and fire type digits; the call
#' `split_integer_first_last(as.integer(um), 5, 4, 1)` is currently hard-coded.
#'
#' @param r_template A `SpatRaster` template defining target extent/resolution/CRS.
#'   Defaults to `template_r`.
#' @param vector_file Path to a vector file readable by [terra::vect()] (e.g.
#'   shapefile, geopackage), or an object accepted by [terra::vect()].
#'   Defaults to `rawFH`.
#' @param fields Character vector of attribute names in `vector_file` used to
#'   derive the fire-history coding (e.g. `c("SEASON", "FIRETYPE_NO")`).
#' @param combField_multiplier Integer multiplier used when combining fields
#'   into a single code (e.g. `SEASON * combField_multiplier + FIRETYPE_NO`).
#'   (Currently passed through to downstream helpers; not used directly here.)
#' @param background Background value for rasterization (typically 0).
#' @param overwrite Logical; whether to overwrite existing outputs (passed to
#'   downstream helpers).
#' @param datatype Output raster datatype (e.g. `"INT2U"`). Passed to downstream
#'   helpers.
#' @param compress Compression method for written rasters (e.g. `"LZW"`).
#' @param memfrac Fraction of available memory to allow for processing.
#'   Passed to downstream helpers.
#' @param progress Logical; show progress messages/bars in downstream helpers.
#' @param OtherAndUnknown Integer code for "Other/Unknown" fire type handling in
#'   downstream logic.
#' @param start.SEASON Integer; first season/year for analyses in
#'   `add_fire_lft_lby_ysf()`.
#' @param end.SEASON Integer or `NA`; optional final season/year for analyses.
#' @param max_interval Integer; maximum interval argument forwarded to
#'   `add_fire_lft_lby_ysf()`.
#'
#' @return A list containing: \itemize{
#' \item OutDF sf polygons dataframe containing all the fire history attributes
#' each row represents a unique sequence of fire events
#' \item TimeSpan integer vector sequence of SEASONS to in the analysis output
#' \item YSFNames character vector of names of TSF years in output, needed by downstream functions
#' \item LBYNames character vector of names of LBY years in output, needed by downstream functions
#' \item LFTNames character vector of names of LBY years in output, needed by downstream functions
#' \item FH_ID integer vector giving the cell values of a raster of same extent resolution and crs
#'  as the input r_template, these values are the indeces and ID values of the rows in OutDF and
#'  used in downstream processes}
#'
#' @seealso
#' [terra::vect()], [terra::rast()], [terra::subst()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' # r_template <- rast("template.tif")
#' # vector_file <- "fire_history.gpkg"
#' res <- fhProcess_raster_first(
#'   r_template = r_template,
#'   vector_file = vector_file,
#'   fields = c("SEASON", "FIRETYPE_NO"),
#'   start.SEASON = 1980
#' )
#' }
fhProcess_raster_first <- function(r_template,
                                   v,
                                   fields = c("SEASON", "FIRETYPE_NO"),
                                   combField_multiplier = 10L,
                                   background = 0,
                                   overwrite = TRUE,
                                   datatype = "INT2U",
                                   compress = "LZW",
                                   memfrac = .8,
                                   progress = TRUE,
                                   OtherAndUnknown = 2,
                                   start.SEASON = 1980,
                                   end.SEASON = NA,
                                   max_interval = 0) {

  s <- make_fh_raster_stack(
    r_template = r_template,
    v = v,
    fields = fields,
    combField_multiplier = combField_multiplier,
    background = background,
    overwrite = overwrite,
    datatype = datatype,
    compress = compress,
    memfrac = memfrac,
    progress = progress,
    OtherAndUnknown = OtherAndUnknown
  )



  p1<-process_blocks_to_bin(s)

  r2<-write_raster_from_bin(r_template = r_template,bin_file = p1$bin,index_file = p1$index,max_ncol = p1$max_ncol)

  um<-as.matrix(terra::unique(r2,na.rm = FALSE))
  FH_ID<-1:nrow(um)
  FH_IDr<-terra::subst(x = r2,um,FH_ID)
  names(FH_IDr)<-"ID" #  need to rename this year so that fhAnalysis$FH_ID also has name ID

  ums<-split_integer_first_last(as.integer(um),5,4,1 )


  FTm<-SEASm<-um
  SEASm[]<-ums$first
  colnames(SEASm)<-paste0("SEAS",sprintf("%02d",1:ncol(um)))
  FTm[]<-ums$last
  colnames(FTm)<-paste0("FireType",sprintf("%02d",1:ncol(um)))

  OutDF<-cbind(FH_ID,SEASm,FTm)
  rm(FH_ID)
  fhAnalysis<-add_fire_lft_lby_ysf (OutDF = OutDF,max_interval = max_interval,start.SEASON = start.SEASON,end.SEASON = end.SEASON,v=v)
  fhAnalysis$FH_ID<-values(FH_IDr)

  return(fhAnalysis)
}
