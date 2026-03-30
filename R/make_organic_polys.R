#' Create organic "forest-cover" polygons and clip to a bounding box (terra-only)
#'
#' @description
#' Generates \code{n} organic-looking polygon "blobs" using only \pkg{terra}
#' (no \pkg{sf}). Each polygon is generated around a centroid that is sampled
#' uniformly inside a user-specified bounding box (\code{ext}). Polygons are
#' allowed to extend beyond the bounding box during generation, and are then
#' clipped to the bounding box using \code{terra::intersect()}.
#'
#' The organic boundary is created by sampling points on a noisy ellipse in
#' polar angle, smoothing the noise on a circle, and scaling the resulting
#' polygon to approximately match a target area.
#'
#' @param n Integer. Number of polygons to generate.
#'
#' @param ext A \code{terra::ext()} extent (or something coercible via
#'   \code{terra::ext()}) giving the bounding box to use for centroid placement
#'   and clipping. Default matches your earlier Vicgrid extent.
#'
#' @param crs Character. Coordinate reference system for the output, passed to
#'   \code{terra::vect(..., crs=)}. Default \code{"EPSG:3111"} (GDA94 / Vicgrid).
#'
#' @param area_prop Numeric length-2 vector. Target polygon area as a proportion
#'   of bounding-box area \emph{before clipping}. Defaults to \code{c(0.05, 0.50)}.
#'
#' @param aspect Numeric length-2 vector. Range for ellipse axis ratio \code{a/b}
#'   used to shape blobs. Defaults to \code{c(0.6, 1.8)}.
#'
#' @param n_vertices Integer length-2 vector. Range for number of boundary
#'   vertices per polygon. More vertices gives smoother boundaries but larger
#'   objects. Defaults to \code{c(40L, 120L)}.
#'
#' @param noise Numeric length-2 vector. Range for boundary roughness (fraction
#'   of radius). Higher values produce more ragged edges. Defaults to
#'   \code{c(0.12, 0.45)}.
#'
#' @param smooth_k Integer. Smoothing window size (odd is best) for the circular
#'   moving-average applied to boundary noise. Larger values yield smoother blobs.
#'   Defaults to \code{7L}.
#'
#' @param integer_coords Logical. If \code{TRUE}, rounds vertices to whole units
#'   \emph{before clipping}. Note: clipping will generally introduce new vertices
#'   at the bbox edge which may be non-integers, but the pre-clip vertices will
#'   be integers. Defaults to \code{FALSE}.
#'
#' @param recalc_metrics Logical. If \code{TRUE}, adds \code{area_m2} and
#'   \code{area_prop_postclip} after clipping. Adds \code{perim_m} if supported
#'   by your \pkg{terra} version. Defaults to \code{TRUE}.
#'
#' @param seed Optional integer. If supplied, sets \code{set.seed(seed)} inside
#'   the function for reproducible output.
#'
#' @return A \code{SpatVector} of clipped polygons with attribute columns:
#'   \describe{
#'     \item{poly_id}{1..n}
#'     \item{centroid_x, centroid_y}{sampled centroid coordinates (pre-clip)}
#'     \item{area_prop_target}{target area proportion (pre-clip)}
#'     \item{area_prop_preclip}{actual area proportion achieved (pre-clip)}
#'     \item{aspect_target}{target aspect ratio \code{a/b}}
#'     \item{vertices}{number of vertices used}
#'     \item{roughness}{boundary roughness used}
#'     \item{area_m2}{(optional) clipped polygon area in map units squared}
#'     \item{area_prop_postclip}{(optional) clipped area proportion of bbox area}
#'     \item{perim_m}{(optional) clipped polygon perimeter in map units}
#'   }
#'
#' @examples
#' \dontrun{
#' library(terra)
#' v10 <- make_organic_polys(seed = 1)
#' plot(v10, col = 1:nrow(v10), border = "darkgreen")
#' }
#'
#' @importFrom terra ext xmin xmax ymin ymax vect expanse as.polygons intersect crs
#'
#' @export
make_organic_polys <- function(
    n = 10,
    ext = terra::ext(2550385.6607, 2564188.7741, 2385877.2636, 2406129.8666),
    crs = "EPSG:3111",
    area_prop = c(0.05, 0.50),
    aspect = c(0.6, 1.8),
    n_vertices = c(40L, 120L),
    noise = c(0.12, 0.45),
    smooth_k = 7L,
    integer_coords = FALSE,
    recalc_metrics = TRUE,
    seed = NULL
) {
  stopifnot(requireNamespace("terra", quietly = TRUE))
  if (!is.null(seed)) set.seed(seed)

  # ---- robust extent extraction ----
  e <- terra::ext(ext)
  xminE <- terra::xmin(e); xmaxE <- terra::xmax(e)
  yminE <- terra::ymin(e); ymaxE <- terra::ymax(e)
  if (any(!is.finite(c(xminE, xmaxE, yminE, ymaxE))) || xmaxE <= xminE || ymaxE <= yminE) {
    stop("Invalid extent.")
  }

  W <- xmaxE - xminE
  H <- ymaxE - yminE
  bbox_area <- W * H

  stopifnot(
    length(area_prop) == 2, area_prop[1] > 0, area_prop[2] <= 1, area_prop[1] < area_prop[2],
    length(aspect) == 2, aspect[1] > 0, aspect[1] < aspect[2],
    length(n_vertices) == 2, n_vertices[1] >= 10L, n_vertices[1] <= n_vertices[2],
    length(noise) == 2, noise[1] >= 0, noise[1] < noise[2]
  )

  # ---- circular moving-average smoother ----
  smooth_circular <- function(z, k) {
    k <- as.integer(k)
    if (k <= 1L) return(z)
    if (k %% 2L == 0L) k <- k + 1L
    pad <- (k - 1L) %/% 2L
    zz <- c(tail(z, pad), z, head(z, pad))
    filt <- rep(1 / k, k)
    out <- as.numeric(stats::filter(zz, filt, sides = 2))
    out <- out[(pad + 1):(pad + length(z))]
    out[is.na(out)] <- z[is.na(out)]
    out
  }

  # ---- blob boundary around centroid using noisy ellipse radius ----
  blob_mat <- function(cx, cy, a, b, nv, rough, k_smooth) {
    theta <- seq(0, 2*pi, length.out = nv + 1L)[-(nv + 1L)]
    z <- smooth_circular(rnorm(nv), k_smooth)

    r <- 1 + rough * z / max(1e-9, max(abs(z)))
    r <- pmax(r, 0.15)

    x <- cx + (a * r) * cos(theta)
    y <- cy + (b * r) * sin(theta)

    m <- cbind(x, y)
    rbind(m, m[1, , drop = FALSE])  # close ring
  }

  # ---- scale polygon about centroid ----
  scale_about <- function(m, cx, cy, s) {
    mm <- m
    mm[, 1] <- cx + (m[, 1] - cx) * s
    mm[, 2] <- cy + (m[, 2] - cy) * s
    mm
  }

  # ---- axes from target area and aspect ratio ----
  sample_axes <- function(A_tgt, r_as) {
    a <- sqrt(A_tgt * r_as / pi)
    b <- sqrt(A_tgt / (r_as * pi))
    c(a = a, b = b)
  }

  mats <- vector("list", n)

  atts <- data.frame(
    poly_id = seq_len(n),
    centroid_x = NA_real_,
    centroid_y = NA_real_,
    area_prop_target = NA_real_,
    area_prop_preclip = NA_real_,
    aspect_target = NA_real_,
    vertices = NA_integer_,
    roughness = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n)) {
    cx <- runif(1, xminE, xmaxE)
    cy <- runif(1, yminE, ymaxE)

    a_prop <- runif(1, area_prop[1], area_prop[2])
    A_tgt  <- a_prop * bbox_area

    r_as <- runif(1, aspect[1], aspect[2])
    axes <- sample_axes(A_tgt, r_as)
    a0 <- as.numeric(axes["a"]); b0 <- as.numeric(axes["b"])

    nv <- sample.int(n_vertices[2] - n_vertices[1] + 1L, 1L) + n_vertices[1] - 1L
    rough <- runif(1, noise[1], noise[2])

    m <- blob_mat(cx, cy, a0, b0, nv, rough, smooth_k)

    # scale to match target area (roughness changes area)
    cand <- terra::vect(list(m), type = "polygons", crs = crs)
    A0 <- as.numeric(terra::expanse(cand))
    if (!is.finite(A0) || A0 <= 0) {
      # rare fallback
      m <- blob_mat(cx, cy, a0, b0, nv, rough = 0, k_smooth = 1L)
      cand <- terra::vect(list(m), type = "polygons", crs = crs)
      A0 <- as.numeric(terra::expanse(cand))
    }

    s <- sqrt(A_tgt / A0)
    m2 <- scale_about(m, cx, cy, s)

    if (isTRUE(integer_coords)) {
      m2 <- round(m2)
      if (!all(m2[1, ] == m2[nrow(m2), ])) m2 <- rbind(m2, m2[1, , drop = FALSE])
    }

    mats[[i]] <- m2

    cand2 <- terra::vect(list(m2), type = "polygons", crs = crs)
    A2 <- as.numeric(terra::expanse(cand2))

    atts$centroid_x[i] <- cx
    atts$centroid_y[i] <- cy
    atts$area_prop_target[i] <- a_prop
    atts$area_prop_preclip[i] <- A2 / bbox_area
    atts$aspect_target[i] <- r_as
    atts$vertices[i] <- nv
    atts$roughness[i] <- rough
  }

  # Build SpatVector, then attach attributes (works on older terra versions)
  v <- terra::vect(mats, type = "polygons", crs = crs)
  v <- cbind(v, atts)

  # ---- CLIP to bounding box ----
  bb_poly <- terra::as.polygons(e, crs = terra::crs(v))
  v_clip <- terra::intersect(v, bb_poly)

  # Recompute metrics after clipping (optional)
  if (isTRUE(recalc_metrics) && nrow(v_clip) > 0) {
    v_clip$area_m2 <- as.numeric(terra::expanse(v_clip))
    v_clip$area_prop_postclip <- v_clip$area_m2 / bbox_area

    # perim() is not present in very old terra; guard
    if ("perim" %in% getNamespaceExports("terra")) {
      v_clip$perim_m <- as.numeric(terra::perim(v_clip))
    }
  }

  v_clip
}
