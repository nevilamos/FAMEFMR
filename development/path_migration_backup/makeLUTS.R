## MAKE LOOK UP TABLES

# list of options for choice of region/ state or adhoc polygon for analysis
#ultimately should be possible to remove this and use  REG_LUT for the same purposes.
REG_NO <- c("USER DEFINED POLYGON" = 7,
  "WHOLE OF STATE" = 99,
  "BARWON SOUTH WEST" = 1,
  "GIPPSLAND" = 2 ,
  "GRAMPIANS" = 3,
  "HUME" = 4,
  "LODDON MALLEE" = 5,
  "PORT PHILLIP" = 6
)

REG_LUT <- data.frame(
  FIRE_REG = c(7,99, 1, 2, 3, 4, 5, 6),
  FIRE_REGION_NAME = c(
    "USER DEFINED POLYGON",
    "WHOLE OF STATE",
    "BARWON SOUTH WEST",
    "GIPPSLAND",
    "GRAMPIANS",
    "HUME",
    "LODDON MALLEE",
    "PORT PHILLIP"
    )
)

# Fire FMZ look up table
FIREFMZ_LUT <- data.frame(
  FIREFMZ = c(0, 1, 2, 3, 4, 5),
  FIRE_FMZ_NAME = c(
    "0 - Non FPA Land",
    "1 - Asset Protection Zone",
    "2 - Bushfire Moderation Zone",
    "3 - Landscape Management Zone",
    "4 - Planned Burn Exclusion Zone",
    "5 - Unknown. Contact Fire Management Officer"
  ),
  FIRE_FMZ_SHORT_NAME = c("Non FPA",
                          "APZ",
                          "BMZ",
                          "LMZ",
                          "PBEZ",
                          "UNK")
)

# DELWP region names look up table
DELWP_LUT <-  data.frame(
    DELWP = c(2, 4, 3, 5, 6, 7),
    #####-----LUT and then you could just slice that LUT.
    DELWP_REGION = c(
      "GIPPSLAND",
      "HUME",
      "PORT PHILLIP",
      "BARWON SOUTH WEST",
      "LODDON MALLEE",
      "GRAMPIANS"
    )
  )

# TFI_STATUS values lookup table and export to csv
TFI_STATUS_LUT <- structure(
  list(
    TFI_VAL = c(-99, 0, 1, 5, 6),
    TFI_STATUS = c(
      "NONE",
      "WITHIN_TFI",
      "BELOW_MIN_TFI",
      "ABOVE_MAX_TFI",
      "ABOVE_MAX_BELOW_MIN_HI_TFI"
    )
  ),
  row.names = c(NA, -5L),
  class = c("tbl_df", "tbl", "data.frame")
)


# Firetype look up table
FIRETYPE_LUT <- data.frame(TYPE = c(1, 2),
                           FIRETYPE = c("BURN", "BUSHFIRE"))

# Growth Stage look up table
GS_LUT <- data.frame(
  "GS" = c(0, 1, 2, 3, 4),
  "GROWTH_STAGE" = c("Unknown", "Juvenile", "Adolescent", "Mature", "Old")
)

# EFG to TFI attributes look up table
# read csv version of CGDL lookup table
TFI_LUT <-
  read_csv("./ReferenceTables/EFG_EVD_TFI.csv")[, c("EFG_NUM", "MIN_LO_TFI",
                                                    "MIN_HI_TFI", "MAX_TFI", "EFG_NAME")]
names(TFI_LUT)[1] <- "EFG"
