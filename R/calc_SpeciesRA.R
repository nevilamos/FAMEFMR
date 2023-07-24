#' Species' relative abundance calculation and summary
#' @details Calculates the relative abundance of species for each raster cell in
#'   analysis and summaries these as summed abundance each season. Optionally it
#'   also write relative abundance rasters for species to disk
#' @param myFHAnalysis list containing all the fire history spatial attributes
#'   created by function fhProcess
#' @param myHDMSpp_NO vector of TAXON_IDs for species to be included in output
#' @param myWriteSpRasters logical: whether to also write species abundance
#'   rasters to disk
#' @param myLU_List list of species abundance lookup arrays created by function
#'   make_Spp_LU_list()
#' @param myResultsDir path of directory where raster results will be written usually
#'   generated  by FAME script
#' @param myTaxonList data.frame of species attributes (read from default or
#' user provided .csv) including field HDMPaths with the relative or absolute
#' path to the HDM raster used in the analysis.
#' @param writeYears  vector for SEASONS for which rasters are to be written
#'   otherwise if writeSpRasters == TRUE, if writeYears == NULL then all SEASONS
#'   are written out
#' @param myWriteSp vector of TAXON_IDs provided if only subset of species
#'   rasters are required as output.
#' @param myAllCombs all combinations of raster values object produced by
#'   function calc_U_AllCombs
#' @param myCropRasters  object containing cell indices and values of input rasters cropped to area of interest created by cropToOutput()
#'
#' @return list of two or three data frames:
#' \itemize{
#' \item SpYearSummWide summary of relative abundance of species by pivoted wide by SEASONS
#' \item SpYearSummLong Long Format summary of relative abundance of species by SEASONS
#' \item grpSpYearSumm summary of abundance of species by pivoted wide by SEASONS grouped by
#' myAllCombs
#' }

#' @export
calc_SpeciesRA <- function(myFHAnalysis,
                           myAllCombs = allCombs,
                           myHDMSpp_NO = HDMSpp_NO,
                           myWriteSpRasters = TRUE,
                           myLU_List = LU_List,
                           myResultsDir = ResultsDir,
                           myTaxonList = TaxonList,
                           writeYears = NULL,
                           myWriteSp = writeSp,
                           myCropRasters = cropRasters
                           ) {
  #change the directory for HDMs to the correct resolution
  if(myFHAnalysis$RasterRes == 225){}else{myTaxonList$HDMPath = gsub("/225m/","/75m/",myTaxonList$HDMPath)}
  myIDX=myCropRasters$IDX
  # set time range for analysis
  TimeRange <- as.integer(myFHAnalysis$TimeSpan)
  TimeNames <- as.character(myFHAnalysis$TimeSpan)
  if(length(grep("NoBurn",myFHAnalysis$YSFNames))>0){TimeNames=c(TimeNames,"NoBurn")}
  LTR <- length(TimeNames)

  # reads in raster from fhAnalysis as Template
  r <- eval(myCropRasters$rasterDef)

  # determine what years to write out (each year or TimeSpan)
  TimeSpan <- myFHAnalysis$TimeSpan
  myDF <- myFHAnalysis$OutDF
  if (is.null(writeYears))
  {
    writeYears = TimeSpan
  } else{
    writeYears <- writeYears[writeYears %in% TimeSpan]
  }

  # remove geometry myFHAnalysis DF to create a standard dataframe so columns
  # can be subset without sticky geometry of original spatial Features data
  # frame
  sf::st_geometry(myDF) <- NULL

  # create empty matrix with rownames containing HDM TAXON ID numbers and column
  # names of the years (SEASONS) used to house the Species abundance data for
  # each year
  SpYearSumm <-
    matrix(
      NA,
      nrow = (length(myHDMSpp_NO)),
      ncol = LTR,
      dimnames = list(as.character(myHDMSpp_NO), TimeNames)
    )
  grpSpYearSumm <- NULL
  # loop through calculation of per cell species abundance values
  # output raster values of flagged spp
  for (sp in myHDMSpp_NO) {
    cat("\r", paste("calculating abundances for", sp))
    mySpp <- as.character(sp)

    #get the lookup array of abundance values from the list of species value lookup
    LU = myLU_List[[mySpp]]

    # gets the HDMVals for the relevant species from the  HDMvalues by
    # species. values multiplied by 100 so that they can later be converted to
    # integer if necessary without losing small values

    HDM_Vector<-terra::extract(terra::rast(myTaxonList[myTaxonList$TAXON_ID== sp,][["HDMPath"]]),myIDX)[,1]*100

    # makes matrices of YSF, EFG, and LFT to use in lookup from LU array
    YSF_M <-
      as.matrix(myDF[myAllCombs$U_AllCombs_TFI$ID, myFHAnalysis$YSFNames]) + 1
    LFT_M <-
      as.matrix(myDF[myAllCombs$U_AllCombs_TFI$ID, myFHAnalysis$LFTNames])
    EFG_M <-
      matrix(myAllCombs$U_AllCombs_TFI$EFG, nrow(YSF_M), ncol(YSF_M))

    # looks up the cell-wise species values for for the species abundance values
    # by indices in array
    Spp_M <- array(LU[cbind(as.vector(YSF_M),
                            as.vector(EFG_M),
                            as.vector(LFT_M))],
                   dim = dim(YSF_M))

    # Multiplies these values by cell-wise HDM values
    #effectively masking out values where the species does not occur.
    Spp_Val_Cell_Year <-
      Spp_M[myAllCombs$Index_AllCombs,] * HDM_Vector

    colnames(Spp_Val_Cell_Year) <- TimeNames



    # get the sum of cell values for each year for the species
    # put them in the compilation data frame
    SpYearSumm[mySpp,] <- colSums(Spp_Val_Cell_Year, na.rm = TRUE)

    grpSpYearSumm <-rbind(grpSpYearSumm,
                          data.table::as.data.table(Spp_Val_Cell_Year)%>%
                            dplyr::group_by(myAllCombs$Index_AllCombs)%>%
                            dplyr::summarise(dplyr::across(tidyselect::everything(),
                                                           \ (x) sum(x,na.rm=T)))%>%
                            dplyr::mutate(TAXON_ID = mySpp)
    )


    #clean memory
    gc()

    # if writing species rasters is desired (function setting myWriteSpRasters == TRUE)
    # write out species rasters
    # this is by far the most time consuming part of the FAME processing

    if (myWriteSpRasters == TRUE) {
      for (myYear in as.character(writeYears)) {
        cat("\r", paste("writing species abund rasters for", myYear))
        #print(paste("writing species abund rasters for", myYear))
        if (sp %in% myWriteSp | is.null(myWriteSp)) {
          OutTif <-
            file.path(myResultsDir,
                      "RA_Rasters",
                      paste0("Spp", sp, "_", myYear, ".tif"))
          print(OutTif)
          emptySpraster <- r
          terra::values(emptySpraster) <-
            Spp_Val_Cell_Year[, myYear]
          terra::writeRaster(
            x = emptySpraster,
            filename = OutTif,
            gdal = c("COMPRESS=LZW", "TFW=YES"),
            datatype = 'INT2U',
            overwrite = TRUE
          )
        }
      }
    }
  } # end species raster cell abundance calcs loop

  # join species details from TaxonList  to the output tables
  SpYearSumm <-
    tibble::rownames_to_column(as.data.frame(SpYearSumm))
  names(SpYearSumm)[1] <- "TAXON_ID"
  TL <- myTaxonList %>%
    dplyr::mutate(TAXON_ID = as.character(TAXON_ID))
  SpYearSummWide <- dplyr::right_join(TL,
                                      SpYearSumm)
  SpYearSummLong <- dplyr::right_join(
    TL,
    SpYearSumm %>%
      tidyr::pivot_longer(-TAXON_ID,
                          names_to = "SEASON",
                          values_to = "SUM_RAx100"),
    multiple = "all") %>%
    dplyr::mutate(SEASON = ifelse(SEASON == "NoBurn","9999",SEASON))%>%
    dplyr::mutate(SEASON = as.integer(SEASON))



  # return calc_SpeciesRA function
  print("finished calc_SpeciesRA")
  return(list("SpYearSummWide" = SpYearSummWide, "SpYearSummLong" = SpYearSummLong,"grpSpYearSumm"=grpSpYearSumm))
}
