ui <- dashboardPage(  title = "FAME 4",
                      
  dashboardHeader(
    #formatting individual letters to different size and colour ----
    title = span(
      "F",
      style = "color: white; font-size: 40px",
      span("ire",
           style = "color: black; font-size: 30px"),
      span("A",
           style = "color: white; font-size: 40px"),
      span("nalysis",
           style = "color: black; font-size: 30px"),
      span("M",
           style = "color: white; font-size: 40px"),
      span("odule for",
           style = "color: Black; font-size: 30px"),
      span("E",
           style = "color: white; font-size: 40px"),
      span("cological values",
           style = "color: black; font-size: 30px")
    ),
    titleWidth = "100vw"
  ),
  #sidebar----
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      sidebarMenuOutput("tabsmenu"),
      # Button to open fame manual from pdf in www directory
      
      actionButton(inputId='ab1', label="Open manual in new tab", 
                   icon = icon("th"), 
                   onclick ="window.open('FAME_User_Manual.pdf', '_blank')"),
      
     
      
      
      # save and reload analysis  buttons ----
      
      box(
        width = 12,
        title = "Save or Load Analyses",
        withBusyIndicatorUI(
          shinyFilesButton(
            id = "loadAnalysis",
            label = "load existing analysis file",
            title = "select analysis to load",
            multiple = FALSE,
            buttonType = "default",
            class = NULL,
            icon = icon("folder-open"),
            style = NULL,
            viewtype = "icon"
          )
        ),
        withBusyIndicatorUI(
          shinySaveButton(
            id = "saveAnalysis",
            label = "save analysis to file",
            title = "save analysis file as...",
            filename = "savedAnalysis",
            filetype = list(qs = "qs"),
            icon = icon("file-export"),
            viewtype = "icon"
          )
        ),
        withBusyIndicatorUI(
          selectFileUI("fileForDashboard",
                       label = tags$div(HTML('<i class="fa-regular fa-box"></i> Select and package for dashboard'))
                        )
        ),

        uiOutput("openDashBoardBtn")

      ),
      box(
        width = 12,
        title = "WARNING BUTTON BELOW WILL SHUT DOWN SERVER",
        useShinyjs(),
        extendShinyjs(text = jscode, functions = c("closeWindow")),
        h5("Make sure you have downloaded all your data first"),
        actionButton(inputId="close", 
                     label =tags$div(HTML(
                       '<i class="fa-sharp fa-solid fa-circle-exclamation fa-2xl" style="color: #ff0000;"></i>  Shut Down Server')
                       )
                     )
      ),
      box(
        width = 12,
        h4("FAME version info:"),
        textOutput("FAMEGUIVersion"),
        textOutput("FAMEFMRVersion")
        
      )
      
      
    ),
    
    absolutePanel(
      conditionalPanel(
        # style = "background-color:#FFFFFF",
        condition = "($('html').hasClass('shiny-busy'))",
        img(
          src = "Fire-animation.gif",
          width = "300",
          height = "200"
        ),
        # Creative Commons Attribution-Share Alike 3.0 Unported license.https://commons.wikimedia.org/wiki/File:BurningFlame0.gif
        h6(
          "https://commons.wikimedia.org/wiki/File:BurningFlame0.gif"
        )
      )
    )
  ),
  #dashboard body ----
  dashboardBody(
    # this tags$head section changes the shinydashboard to DEECA corporate colours ----
    # DEECA Navy 100% =  hex #201547
    # DEECA Teal 100% =  hex #00B2A9
    # Vic Gov Blue = hex #004C97
    # Corporate Sky Blue = #88DBDF
    # ECA Lime = #CDDC29
    
    tags$head(tags$style(
      HTML(
        '
        .main-header .logo {
        font-weight: bold;
        font-size: 40px;
        font-family:Georgia;

      }
      .main-header .logo #mydiv b{
        font-weight: bold;
        font-size: 40px;
        font-family:Georgia;
      }


        /* button */
        .skin-blue .btn {
        background-color: #004C97;
        color: #ffffff;
        border-radius: 10px;
        }
        /* logo */
        .skin-blue .main-header .logo {
                              background-color: #00B2A9;
                              }

        /* logo when hovered */
        .skin-blue .main-header .logo:hover {
                              background-color: #00B2A9;
                              }

        /* navbar (rest of the header) */
        .skin-blue .main-header .navbar {
                              background-color: #00B2A9;
                              }

        /* main sidebar */
        .skin-blue .main-sidebar {
                              background-color: #00B2A9;
                              }

        /* active selected tab in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: #201547;
                              }

        /* other links in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                              background-color: #00B2A9;
                              color: #ffffff;
                              }

        /* other links in the sidebarmenu when hovered */
         .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: #004C97;
         }

        /* Box colours */
        .skin-blue .box {
                    color:#ffffff;
                    background:#201547;
                    border-bottom-color:#201547;
                    border-left-color:#201547;
                    border-right-color:#201547;
                    border-top-color:#201547;
                    background:#201547;
                    border-radius: 10px;
        }
       /* box title force to white */
        .skin-blue .box-header {
        color:#ffffff;

        }
      /* pagewell formatting */
      .skin-blue .well{
        color:#ffffff;
        background:#201547;
        border-bottom-color:#201547;
        border-left-color:#201547;
        border-right-color:#201547;
        border-top-color:#201547;
        background:#201547; 
        border-radius: 10px;
      }
      .skin-blue .sidebar-menu  {
        white-space: normal;
      }
      '
      )
    )),
    setBackgroundImage(src = "08732250_before_after_2014_fire.jpg", shinydashboard = TRUE),
    
    tabItems(
      # Upload and other utilities tab----
      tabItem(
        tabName = "util",
        fluidRow(column(
          6,
          box(
            width = 12,
            title = "Select and upload files",
            
            uploadFileUI(id="rawFH",label = "Select rawFH for upload geopackage (.gpkg) or shapefile (.dbf, .prj,.shp, and .shx)"),
            uploadFileUI(id="adHocPoly",label = "Select Area of Interest geopackage (.gpkg) or shapefile (.dbf, .prj,.shp, and .shx)"),
            uploadFileUI(id="puPoly",label = "Select PU/BU geopackage (.gpkg) or shapefile (.dbf, .prj,.shp, and .shx)"),
            
            uploadFileUI(id= "addCustomCSV",
              label = "Add Custom input CSV"),
            uploadFileUI(id=  "addCustomHDM225",
              label = "Add Custom 225m HDM files as .tif"
              #,
              #accept = c(".tif"),
              #multiple = T
            ),
            uploadFileUI(id ="addCustomHDM75",
                         label = "Add Custom 75m HDM files as .tif"
                         #,
                         #accept = c(".tif"),
                         #multiple = T
            ),
            
            fileInput("addCustomUnformattedHDM", NULL,accept=".tif", buttonLabel = "Load custom HDM .tif file for reformatting", multiple = FALSE),
            radioButtons("outRes",label = "Custom HDM output resolution",choices = list("225m"=225,"75m"=75),inline =TRUE ),
            textOutput(outputId = "formattedHDMPath"),
            
            uploadFileUI(id = "addSavedAnalysis",
                         label = "Upload saved analysis '.qs' file",
                         multiple = T)
          )
        ),
        column(
          6,
          box(
            width = 12,
            title = "Create draft species list or aspatial GSO inputs",
            selectInput("spREGION_NO", "Choose an area for species list",
                        choices = as.list(c(REG_NO))),
            conditionalPanel(
              condition = "input.spREGION_NO == '7'",
              selectInput(
                "spAdHocShape",
                "Select AdHoc Area shapefile",
                choice = c("", list.files("./AdHocPolygons/",
                                          pattern =
                                            ".shp$"))
              )
            ),
            radioButtons("sppublic", "Restrict to Public Land", c("Yes" = TRUE, "No" = FALSE)),
            
            
            # runscript button
            withBusyIndicatorUI(actionButton("runDSpList",
                                             label = "Run draft species list")),
            withBusyIndicatorUI(
              actionButton("runspEFGpList",
                           label = "Run Spp EFG LMU for list for GSO")
            )
          ),
          box(
            width = 12,
            title = "Download Files",
            fluidRow(column(4,
                            h5(
                              paste("Current Results Directory is:")
                            )),
                     column(8, h5(
                       textOutput("resultsDir")
                     ))),
            shinyFilesButton(
              "files",
              label = "Select file(s) to download",
              title = "",
              multiple = T
            ),
            tableOutput("filepaths"),
            downloadButton("downloadFiles", "Download Files")
          )
        )),
        fluidRow(column(
          4,
          box(
            width = 12,
            title = "Downloads",
            downloadButton("downloadTool", "Download FAME ArcGIS preproccessing tool")
            #downloadButton("downloadManual", "Download FAME manual"),
            # tags$h1("Revised, alternate 2"),
            
          )
        ),
        column(8,),)
      ),
      
      
      
      # fhAnalysis tab content----
      tabItem(
        tabName = "AnalysisSettings",
        fluidRow(
          column(
            6,
            box(
              title = "FH Analysis Settings",
              width = 12,
              #00B2A9
              splitLayout(
                selectFileUI(id = "rawFHPath",
                             label = "Select rawFH path"),
                textOutput("rawFHPath"),
                
                selectFileUI(id = "secondFH",
                             label = "Select secondFH path"),
                textOutput("secondFH")
              ),
              numericInput("baseFire",
                           "Include Base fire in FH analysis",
                           value =NA,width = "15%"),
              
              selectInput(
                "REGION_NO",
                "Choose a Region",
                choices = as.list(c(REG_NO)),
                width = "50%"
              ),
              conditionalPanel(condition = "input.REGION_NO == '7'",
                               splitLayout(
                                 selectFileUI(id = "AdHocPath",
                                              label = "Select AdHoc path"),
                                 textOutput("AdHocPath")
                                 
                               ),),
              fluidRow(
                column(
                  3,
                  radioButtons("RasterRes", "Select Raster Resolution", cellSizes),
                ),
                column(3,
                       radioButtons(
                         "public", "Restrict to Public Land", c("Yes" = TRUE, "No" = FALSE)
                       )),
                column(3,
                       radioButtons(
                         "max_interval", "Correct for Burn after Bushfire", c("Yes" = 5, "No" = 0)
                         )),
                
                column(
                  3,
                  radioButtons(
                    inputId = "otherUnknown",
                    label = "Value for other and unknown fires",
                    choices = c(
                      "Bushfire" = 2,
                      "Burn" = 1,
                      "NA" = NA
                    )
                  )
                )
              ),
              numericInput(
                "startTimespan",
                label = "First season for analysis output.<br />Value must be <= Last season for analysis output. ",
                value = 1980
              ),
              numericInput(
                "endTimespan",
                label = HTML("Last season for analysis output.<br />Value must be >= highest SEASON value in input fire mapping<br />or Run FH Analysis button  will return a 'subscript out of bounds' error"),
                value = get_fire_season(Sys.Date())
              ),
              
            ),
            box(
              title = "JFMP Settings",
              width = 12,
              #00B2A9
              checkboxInput(
                inputId = "usePUpolys",
                label = "Include burn unit/planning unit shapefile for JFMP analysis",
                value = FALSE
              ),
              
              # select planning unit shapefile----
              conditionalPanel(
                condition = "input.usePUpolys",
                splitLayout(
                  selectFileUI(id = "puPath",
                               label = "Select planning unit path"),
                  textOutput("puPath")
                ),
                numericInput("JFMPSeason0",
                             "JFMP SEASON 0",
                             as.integer(format(Sys.Date(
                             ), "%Y")))
              )
            )
          ),
          
          
          # Fauna abundance headings conditional on whether JFMP or "standard" analysis is been run ----
          column(
            6,
            box(
              width = 12,
              #00B2A9
              conditionalPanel(condition = "input.usePUpolys == 0",
                               h4("Fauna Abundance Calculation Choices")),
              conditionalPanel(condition = "input.usePUpolys == 1",
                               h4("JFMP Calculation Choices")),
              # inputs conditional on whether JFMP or "standard" analysis is been run ----
              textOutput("selected_FH_name"),
              splitLayout(
                numericInput(
                  "startBaseline",
                  "enter start season for abundance baseline",
                  "",
                  width = "40%"
                ),
                numericInput(
                  "endBaseline",
                  "enter end season for abundance baseline",
                  "",
                  width = "40%"
                )
              ),
              
              
              # use standard or choose custom species responses----
              
              fluidRow(
                column(
                  2,
                  checkboxInput(
                    inputId = "spListChoice",
                    label = "Use custom species list",
                    value = FALSE,
                    width = NULL
                  )
                ),
                conditionalPanel(
                  condition = "input.spListChoice",
                  selectFileUI(id = "customSpList",
                               label = "select custom species list"),
                  textOutput("customSpList")
                  
                ),
              ),
              splitLayout(
                checkboxInput(
                  inputId = "spResponseChoice",
                  label = "Use custom relative abundance table",
                  value = FALSE,
                  width = NULL
                ),
                conditionalPanel(
                  condition = "input.spResponseChoice",
                  checkboxInput(
                    inputId = "abundByGS",
                    label = "Relative abundance table by growth stage",
                    value = TRUE,
                    width = NULL
                  )
                )
              ),
              conditionalPanel(condition = "input.spResponseChoice",
                               fluidRow(
                                 selectFileUI(id = "customResponseFile",
                                              label = "Select user defined species response file"),
                                 textOutput("customResponseFile")
                                 
                               ),),
              # inputs conditional on whether JFMP or "standard" analysis is been run ----
              conditionalPanel(
                condition = "input.usePUpolys == 0",
                checkboxInput(
                  inputId = "makeRArasters",
                  label = "Make relative abundance rasters",
                  value = FALSE,
                  width = NULL
                ),
                conditionalPanel(
                  condition = 'input.makeRArasters',
                  radioButtons(
                    inputId = "allOrSomeYears",
                    label = "Years to make rasters",
                    inline = TRUE,
                    choices = c("all", "some")
                  ),
                  conditionalPanel(
                    condition = 'input.allOrSomeYears == "some"',
                    selectInput(
                      "yearsForRasters",
                      "Select one or more years for Rasters",
                      choices = "",
                      multiple = T
                    )
                  )
                ),
                conditionalPanel(condition = "input.usePUpolys == 0",
                                 withBusyIndicatorUI(
                                   actionButton("runRA",
                                                label = "Run fauna relative abundance calculations")
                                 )),
              ),
              conditionalPanel(
                condition = "input.usePUpolys == 1",
                splitLayout(
                  selectFileUI(id = "targetHaFilepath",
                               label = "select area target file for JFMP"),
                  textOutput("targetHaFilepath")
                ),
                splitLayout(
                  selectFileUI(id = "jfmpMetricWtFilePath",
                               label = "select file containing JFMP metric weights"),
                  textOutput("jfmpMetricWtFilePath")
                ),
                splitLayout(
                  selectFileUI(id = "zoneWtFile",
                               label = "select file containing JFMP zone weights"),
                  textOutput("zoneWtFile")
                ),
                withBusyIndicatorUI(actionButton("runJFMP1",
                                                 label = "Run JFMP calculations stage 1")),
              )
            ),
            conditionalPanel(
              condition = "input.usePUpolys == 1",
              box(
                width = 12,
                solidHeader = T,
                #00B2A9
                title = "Compare alternative JFMPs",
                splitLayout(
                  selectFileUI(id = "draftJFMPFile",
                               label = "select draft JFMP input"),
                  textOutput("draftJFMPFile")
                ),
                
                withBusyIndicatorUI(
                  actionButton(
                    "runCompareJFMP",
                    label = "Compare Draft JFMP"
                  )
                ),
              )
            ),
            conditionalPanel(
              condition = "input.usePUpolys == 0",
              box(
                title = "TFI and GS Calculations",
                width = 12,
                splitLayout(
                  checkboxInput(
                    inputId = "makeTFIrasters",
                    label = "Make TFIstatus maps for each year",
                    value = FALSE,
                    width = NULL
                  ),
                  checkboxInput(
                    inputId = "makeBBTFIrasters",
                    label = "Make BBTFIstatus maps for each year",
                    value = FALSE,
                    width = NULL
                  )
                ),
                splitLayout(
                  withBusyIndicatorUI(actionButton("runTFI",
                                                   label = "Run TFI calculations")),
                  withBusyIndicatorUI(actionButton("runGS",
                                                   label = "Run Growth Stage calculations"))
                )
              )
            )
          )
        ),
        fluidRow(
          column(6,
                 # runFH analysis action  button----
                 box(
                   width = 12, #00B2A9
                   withBusyIndicatorUI(actionButton("runFH",
                                                    label = "Run FH Analysis"))
                 )),
          conditionalPanel(condition = "input.usePUpolys == 0",
                           column(6,
                                  box(
                                    width = 6, #00B2A9
                                    withBusyIndicatorUI(actionButton("runRA_TFI",
                                                                     label = "Run all calculations")),
                                  )))
        )
      ),
      # fAbund_TFI tab content----
      tabItem(
        tabName = "fAbund_TFI",
        h2("Spatial TFI and Fauna Abundance Calculations "),
        fluidRow(wellPanel()),
        fluidRow(
          column(6,),
          # TFI choices inputs----
          conditionalPanel(condition = "input.usePUpolys == 0",
                           column(6,))
        ),
      ),
      # Tab  TFI charts ----
      
      tabItem(tabName = "TFIplots",
              fluidPage(
                wellPanel(h4("TFI Plots"),
                          fluidRow(
                            column(
                              6,
                              selectInput(
                                inputId = "EFGChoices",
                                label = "TFI EFG Choices",
                                choices = NULL,
                                selected = NULL,
                                multiple = F
                              )
                            ),
                            column(
                              6,
                              sliderInput(
                                "tfiSeasonChoices",
                                label = "Min and Max Season to plot",
                                min = 1980,
                                max = 2020,
                                value = c(1980, 2020),
                                sep = ""
                              )
                            )
                          ),),
                fluidRow(plotlyOutput("TFItrendPlot"),
                         plotlyOutput("BBTFIPlot"))
              )),
      # Tab  GS charts ----
      tabItem(tabName = "GSplots",
              fluidPage(
                wellPanel(h4("GS Plots"),
                          fluidRow(
                            column(
                              6,
                              selectInput(
                                "GSEFGChoices",
                                choices = NULL,
                                label = "Choose EFG Number  (to clear use backspace)",
                                selected = NULL,
                                multiple = F
                              )
                            ),
                            column(
                              6,
                              sliderInput(
                                "GSSeasonChoices",
                                label = "Min and Max Season to plot",
                                min = 1980,
                                max = 2020,
                                value = c(1980, 2020),
                                sep = ""
                              )
                            )
                          ),),
                fluidRow(plotlyOutput("GSPlot"),)
              )),
      # Tab Relative abundance plots----
      tabItem(tabName = "RAplots",
              fluidPage(
                wellPanel(
                  h4("Fauna Relative Abundance Plots"),
                  fluidRow(column(
                    6,
                    selectizeInput(
                      "raSpChoices",
                      choices = NULL,
                      label = "Choose 1-7 species  (to clear use backspace)",
                      selected = NULL,
                      multiple = T,
                      options = list(
                        maxItems = 7,
                        hideSelected = F,
                        placeholder = "Select a species"
                      )
                    )
                  ),
                  column(
                    6,
                    sliderInput(
                      "raSeasonChoices",
                      label = "Min and Max Season to plot",
                      min = 1980,
                      max = 2020,
                      value = c(1980, 2020),
                      sep = ""
                    )
                  )),
                ),
                fluidRow(plotlyOutput("RAtrendPlot"),
                         plotlyOutput("RADeltaPlot"))
              )),
      # tab GSO processing ----
      tabItem(tabName = "GSO",
              fluidRow(# h2("Aspatial GSO Input Selections"),
                column(
                  6,
                  box(
                    width = 12,
                    #00B2A9 title = "Select and upload GSO .csv files",
                    fileInput(
                      inputId = "addGSOCSV",
                      label = "gsofiles to upload",
                      accept = c("csv", ".csv")
                    )
                  ),
                  box(
                    width = 12,
                    #00B2A9 title = "Select GSO input tables",
                    selectInput(
                      "spEFGLMU",
                      "Select Spp_EFG_LMU.csv file",
                      choices = c(
                        "Spp_EFG_LMU.csv",
                        list.files("./GSOInputs/", pattern = "Spp_EFG_LMU.csv$")
                      )
                    ),
                    selectInput(
                      "lmuArea",
                      "LMU_Area.csv file",
                      choice = c(
                        "LMU_Area.csv",
                        list.files("./GSOInputs/", pattern = "LMU_Area.csv$")
                      )
                    ),
                    selectInput(
                      "lmuScenarios",
                      "LMU_Scenarios.csv file",
                      choice = c(
                        "LMU_Scenarios.csv",
                        list.files("./GSOInputs/", pattern = "LMU_Scenarios.csv$")
                      )
                    ),
                    selectInput(
                      "ObsData",
                      "ObsData.csv file",
                      choice = c(
                        "ObsData.csv",
                        list.files("./GSOInputs/", pattern = "ObsData.csv$")
                      )
                    )
                  )
                ),
                column(
                  6,
                  box(
                    width = 12,
                    #00B2A9 title = "Analysis options",
                    h4(
                      "Please note all inputs are case sensitive, do not include any spaces"
                    ),
                    selectInput(
                      "GSOFireType",
                      choices = c("High", "Low"),
                      label = "Low or High fire type."
                    ),
                    selectInput(
                      "GSOBaseLine",
                      label = "Baseline for comparisons.'Optimisation' or select from input scenarios",
                      choices = c("Optimisation")
                    ),
                    selectInput(
                      "GSOFaunaClasses",
                      label = "Which fauna classes to use",
                      choices = c("All", "Birds", "Mammals", "Reptiles", "Frogs"),
                      selected = "All",
                      # "Insects", "Fishes", "Crustaceans", "Annelids",Molluscs", "Nemerteans", "Flatworms", "Cnidarians","Echinoderms", "Zooplankton"
                      multiple = T
                    ),
                    selectInput(
                      "GSOrule",
                      "Select rule to use",
                      choices = c(
                        "Rule0",
                        "Rule1",
                        "Rule1a",
                        "Rule1b",
                        "Rule1c",
                        "Rule2",
                        "Rule2a",
                        "Rule2b",
                        "Rule2c",
                        "Rule3",
                        "Rule3a",
                        "Rule3b",
                        "Rule3c"
                      )
                    ),
                    numericInput("GSOdwt", "weight for option 2", 0.75),
                    numericInput("GSOnrep", "Number of iterations to run", 100),
                    numericInput("GSOnsim", "Number of simulations to generate 95% CI?", 5),
                    withBusyIndicatorUI(actionButton("runGSO",
                                                     "Run Aspatial GSO"))
                  )
                )))
    ),
  )
  
)
