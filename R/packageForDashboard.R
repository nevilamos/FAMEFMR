#' packageForDashboard
#' takes a saved analysis file from FAME and packages it for use in FAME dashboard app
#' @param inPath path to saved analysis from fAME as .qs file
#'
#' @return None but saves the packaed file for use in dashboard to "./FH_Outputs/Dashboard_"
#' @export
#'
packageForDashboard <- function(inPath = "path to saved analysis from fAME as .qs file")
{

#  inPath <- "./FH_Outputs/savedAnalysis.qs"
inName <- basename(inPath)
rv <- qread(inPath)#qread("C:/Users/na03/Downloads/75mwith5sp.qs")

# to allow setting of whether or not to include PU index value in tables
PUcols <- ifelse("PU" %in% names(rv$allCombs$U_AllCombs_TFI), TRUE, FALSE)


if (is.null(rv$puPath)){
  puLUT<-foreign::read.dbf("./ReferenceShapefiles/LF_DISTRICT_with_PU_field.dbf")[,c("PU","DISTRICT_N")]
} else {
  puLUT<-foreign::read.dbf(gsub(".shp$",".dbf",rv$puPath))[,c("PU","DISTRICT_N")]
}


# check for whether any private land vegetation was included in calculations
if (
  (
    sum(
      !is.na(rv$allCombs$U_AllCombs_TFI$EFG)[is.na(rv$allCombs$U_AllCombs_TFI$PLM)]
    )> 0
  ) & rv$public == FALSE
) {

  #public and not public

} else if (  (
  sum(
    !is.na(rv$allCombs$U_AllCombs_TFI$EFG)[is.na(rv$allCombs$U_AllCombs_TFI$PLM)],na.rm=TRUE
  )== 0
) & rv$public == TRUE) {
  #publicOnly
} else{
  warning("PLM values and rv$public do not agree")
}



myColNames <- c("EFG", "DELWP", "FIRE_REG", "FIREFMZ", "PLM", "Index")
if (PUcols) myColNames <- c(myColNames, "PU")

allCombsAttr <- rv$allCombs$U_AllCombs_TFI[, ..myColNames]
allCombsAttr$PLM <- ifelse(is.na(allCombsAttr$PLM), 0, 1)
allCombsAttr <- data.table(sapply(allCombsAttr, as.integer))
x <- data.table(rv$SpYearSumm$grpSpYearSumm)
x <- data.table(sapply(x, as.integer))
names(x)[1] <- "Index"

gc()

seasNames <- names(x)[!names(x) %in% c("Index", "TAXON_ID")]

y <- melt(x,
          id.vars = c("TAXON_ID", "Index"),
          variable.name = "SEASON",
          value.name = "sumRA"
)

rm(x)
gc()

y <- y[allCombsAttr, on = "Index"]
y <- y[!is.na(EFG)]
y[, Index := NULL]
gc()
RA <- y[, .(sumRA = sum(sumRA)), by = c(names(y)[names(y) != "sumRA"])]

allZeros <- RA[, .(sumRA = sum(sumRA)), by = c(names(RA)[!names(RA) %in% c("sumRA", "SEASON")])][sumRA > 0][, sumRA := NULL]

RA <- dplyr::inner_join(RA, allZeros)
# need to add back category names
RA <- data.table(TFI_LUT[, c("EFG", "EFG_NAME")])[RA, on = "EFG"]
RA <- data.table(DELWP_LUT)[RA, on = "DELWP"]
RA <- data.table(REG_LUT)[RA, on = "FIRE_REG"]
RA <- data.table(FIREFMZ_LUT[, 1:2])[RA, on = "FIREFMZ"]
RA <- data.table(rv$TaxonList[, 1:2])[RA, on = "TAXON_ID"]
if(exists("puLUT") & ("PU" %in% names(RA))){
  RA <- data.table(puLUT)[RA, on = "PU"]
  rv$BBTFI$BBTFI_LONG<-data.table(puLUT)[data.table(rv$BBTFI$BBTFI_LONG), on = "PU"]
  rv$TFI<-data.table(puLUT)[data.table(rv$TFI), on = "PU"]
}
myList <- list("TFI" = rv$TFI, "BBTFI" = rv$BBTFI$BBTFI_LONG, "RA" = RA, "TaxonList" = rv$TaxonList)
savePath<-paste0("./FH_Outputs/Dashboard_",inName)
qsave(myList, paste0("./FH_Outputs/Dashboard_",inName))
print(paste("Saved dashboard input file to",savePath))
}

