#' PrepareFire History
#'
#'@description
#' A function that handles the subdivision of a Fire History Dataset  for
#' Parallel processing in the FAME fire history processing workflow.
#' It calculates the number of grid cells to be used in parallelizing function.
#' splitting the input fire history into separate sf datasets on a grid with
#' areas of approximately 100 square km
#'
#' @param inFH input fire history
#'
#' @return A list containing: \itemize{
#' \item hasData list of one or more sf Fire History sub units of the input Fire
#' History containing fire boundaries, defined by a grid over the bounding box
#' of inFH with a grid of squreside*squareside cells
#' \item hasIntersects a logical vector indicating which elements of "hasData"
#' have intersecting Fire boundaries within them
#' \item squareSide integer number of units of one side of the grid used to divide
#' the inFH to produce "hasData" for parrallel compute.
#' Calculated as the floor() of area of inFH bounding box/10^8
#'  (ie a 10km square expressed in square metres)}
#' @export
#'
prepFH<-function(inFH = myFH){

  #calculate dimensions of grid sides
  squareside<-floor(sqrt(as.numeric(sf::st_area(sf::st_as_sfc(sf::st_bbox(inFH))))/10^8))
  squareside<-ifelse(squareside<1,1,squareside)
  if (squareside ==1){
    hasData = list(inFH)
  } else {
    grid<-sf::st_make_grid(inFH ,n=c(squareside,squareside))
    grid<-sf::st_as_sf(grid)
    grid_ID<-1:nrow(grid)
    grid$grid_ID<-grid_ID
    gridded<-sf::st_intersection(inFH,grid) %>% sf::st_make_valid()
    griddedList<-list()
    for(i in grid_ID){griddedList[[i]]<-gridded%>%
      dplyr::filter(grid_ID == i)}

    hasData<-griddedList[lapply(griddedList,nrow)>0]




  }
  hasIntersects<-unlist(lapply(hasData,FUN = function(x){
    max(unlist(lapply(as.list(sf::st_intersects(x)),length)))>0}))
  return(list("hasData"=hasData,
              "squareside" =squareside,
              "hasIntersects"=hasIntersects))
}
