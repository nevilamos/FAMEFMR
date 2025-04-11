#' Main Fire History Fire Sequence analysis function
#'
#' @details The function takes sf polygon collection of fire history
#'   containing polygons with two fields: FIRETYPE and SEASON Where polygons of
#'   different FIRETYPE or SEASON overlap the function constructs unique
#'   non-overlapping polygon of their intersections ( and non intersecting areas
#'   ) and attributes each polygon with sequential fire SEASON (SEAS01, SEAS02
#'   ...) and corresponding FIRETYPE (TYPE01,TYPE02 ...)
#'
#' @param x input Fire History sf polygons with fields SEASON and FIRETYPE
#' @seealso [fhProcess1()] mainFHProcess() normally run from within it
#' @return output processed fire history with unique polygons each with sequence
#'  of fires and firetypes
#' @export
#'
#'
#' @examples
#' # randomFH<-generate_random_fire_history(20)
#' # plot(randomFH)
#' # #the function is run from within fhProcess1()
#' # outFH1<-fhProcess1(randomFH)
#' # plot(outFH1)
#' # #complete all fields in FH analysis
#' # outFH2<-fhProcess2(outFH1)
#' #
#' # plot(outFH2$OutDF,max.plot = 20)
#'
#'
mainFHProcess<-function (x = inFH)  {
  mySF<-x %>%
    dplyr::select(SEASON,FIRETYPE_NO) %>%
    sf::st_make_valid(
      oriented = FALSE,
      #s2_options = s2::s2_options(snap = s2::s2_snap_precision(1e+07), ...),
      geos_method = "valid_structure",
      geos_keep_collapsed = FALSE
    )  %>%
    dplyr::arrange(SEASON,FIRETYPE_NO)

  SEASONS<-mySF$SEASON
  FIRETYPE_NOS<-mySF$FIRETYPE_NO

  outSF<-sf::st_intersection(mySF) %>%
    sf::st_collection_extract()%>%
    sf::st_cast(.,"POLYGON") %>%
    tidyr::unnest_wider(origins,names_sep = "_")




  Vnames<-names(outSF)[grep("origins",names(outSF))]

  originIndices<-outSF[,Vnames] %>%
    sf::st_drop_geometry() %>%
    as.matrix()

  seasonSequence<-as.data.frame(matrix(SEASONS[as.matrix(originIndices)],dim(originIndices)))

  paddedSequence<-formatC(1:dim(originIndices)[2],1,flag="0")

  names(seasonSequence)<-paste0("SEAS",paddedSequence)

  typeSequence<-as.data.frame(matrix(FIRETYPE_NOS[as.matrix(originIndices)],dim(originIndices)))

  names(typeSequence)<-paste0("FireType",paddedSequence)
  outSF<-dplyr::bind_cols(outSF,seasonSequence,typeSequence) %>%
    sf::st_as_sf() %>%
    dplyr::select(-all_of(c(Vnames,"SEASON","FIRETYPE_NO","n.overlaps")))%>%
    dplyr::group_by(across(c(-geometry))) %>%
    dplyr::summarise() %>%
    dplyr::ungroup()



  sf::st_crs(outSF)<-sf::st_crs(mySF)

  return(outSF)

}

