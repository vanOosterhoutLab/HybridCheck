library(shinydashboard)

header <- dashboardHeader(title = "HybridCheck")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("DNA data", tabName = "dnaPage", icon = icon("dashboard")),
    menuItem("Four Taxon Tests", tabName = "fttPage", icon = icon("dashboard")),
    menuItem("Triplet Refinement", tabName = "tripletRefinement", icon = icon("dashboard")),
    menuItem("Triplet Analysis", tabName = "tripletAnalysis", icon = icon("dashboard")),
    menuItem("User Defined Blocks", tabName = "userBlocks", icon = icon("dashboard"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "dnaPage",
            fluidRow(
              column(width = 6,
                     box(title = "Sequence File", status = "warning", width = 12,
                         solidHeader = TRUE,
                         p("Upload your FASTA file containing sequences from your different
                    populations to get started."),
                         p("Note: Duplicate Sequences will be removed."),
                         fileInput('fastafile', 'Choose FASTA file to upload',
                                   accept = c(
                                     '.fas',
                                     '.fasta',
                                     '.FAS',
                                     '.FASTA'
                                   )
                         )
                     ),
                     box(title = "Define Populations", width = 12, status = "warning",
                         solidHeader = TRUE,
                         checkboxInput("oneSeqOnePop", "Define populations."),
                         conditionalPanel("input.oneSeqOnePop",
                                          numericInput("numPops", "How many populations?", value = 4, min = 1),
                                          uiOutput("populationsGen")
                         ),
                         actionButton(inputId = "setPops", label = "Set Populations")
                     )     
              ),
              column(width = 6,
                     fluidRow(
                       valueBox(uiOutput("seqNum"), "Sequences.", icon = icon("list"), 
                                width = 6),
                       valueBox(uiOutput("bpLen"), "Base pairs.", width = 6)
                     ),
                     fluidRow(
                       valueBox(uiOutput("infLen"), "Polymorphic sites.", 
                                color = "orange", width = 6),
                       valueBox(uiOutput("numPopulations"), "Populations.", 
                                color = "maroon", width = 6)
                     ),
                     fluidRow(
                       box(title = "Sequence Names", solidHeader = TRUE, status = "primary",
                           width = 12, uiOutput("seqNames")
                       ) 
                     ),
                     fluidRow(
                       box(title = "Populations", solidHeader = TRUE, status = "primary",
                           width = 12, uiOutput("popDetails")
                       )
                     )
                       
              )
            )       
    ),
    tabItem(tabName = "fttPage",
            fluidRow(
              box(title = "Four Taxon Tests", solidHeader = TRUE, status = "primary",
                  width = 6,
                  p("Four taxon tests accept four sequences or populations and calculates Patterson's D statistic.
                  The four taxon test is also known as the ABBA-BABA test. ABBA and BABA are the names given to
                  two distinctive patterns of mutation:"),
                  p("The ABBA SNP pattern means that P1 has the ancestral state at that site, whereas P2 and P3
                    have a derived state. The BABA pattern means that P1 and P3 share the derived state, and
                    P2 and the ancestral P4 are the same state."),
                  p("Notice in the images the test assumes an explicit evolutionary history: P1 and P2 coalesce first, and P4 is the ancestor / outgroup."),
                  p("Under a scenario of no introgression or no ancestral mutation, for all the two-state mutations in four taxa,
                               the two mutation patterns ABBA and BABA are expected to be about equal and Patterson's D is expected to be zero."),
                  p("Deviations from this exectation suggest introgression and hybridisation."),
                  fluidRow(
                    box(title = "The ABBA SNP pattern", solidHeader = TRUE,
                        width = 6, status = "primary", imageOutput("ABBAtree")),
                    box(title = "The BABA SNP pattern", solidHeader = TRUE,
                        width = 6, status = "primary", imageOutput("BABAtree"))
                  )   
              ),
              box(title = "Parameters", solidHeader = TRUE, status = "warning",
                  width = 6,
                  uiOutput("combosToAnalyze"),
                  hr(),
                  p("Define combinations of 4 taxa to be tested."),
                  actionButton("setCombinations", "Set Taxa Combinations"),
                  checkboxInput("fttAutoSets", "Automatically generate?"),
                  conditionalPanel("!input.fttAutoSets",
                                   numericInput("fttNumCombos", "Number of taxa combos:", 1),
                                   uiOutput("fttGen")
                  )
              )
            ),
          fluidRow(
            box(title = "Four Taxon Test Results", solidHeader = TRUE,
                status = "info", width = 12, fluidRow(uiOutput("combosToView")),
                downloadButton("saveFTT", "Save Test Results to CSV"),
                br(),br(),
                fluidRow(uiOutput("fttResults"))
            )
          )
    ),
    tabItem(tabName = "tripletRefinement",
            box(title = "Specify sequence combinations to scan", width = 6,
                "By default, every possible combination of 3 sequences will be analyzed.",
                "You can refine this by eliminating triplets containing sequences that are insufficiently diverged.",
                radioButtons("fttOrPops", "Triplet Generation Method:",
                             c("Generate triplets based on four taxon test results" = 'ftt',
                               "Generate triplets to scan for recombination between populations" = 'pop')),
                checkboxInput("distancebased", "Raw p-distance based"),
                conditionalPanel(condition = "input.fttOrPops == 'pop'",
                                 radioButtons("partitionStrictness", "How many sequences from the same group are allowed in a triplet?",
                                              c("One" = 1L,
                                                "Two" = 2L), selected=2L)),
                conditionalPanel("input.distancebased",
                                 h4("Distance based elimination options."),
                                 radioButtons("distancemethod", "Use automatic threshold detection?",
                                              c("Yes" = "yes",
                                                "No, specify a raw p-distance to use as a threshold manually." = "no")),
                                 conditionalPanel("input.distancemethod == 'no'",
                                                  sliderInput("mandistthreshold", "Raw p-distance threshold.", 0.01, 1, value=0.01, step=0.01))),
                actionButton("comboGen", "Generate Combinations to Scan")
            ),
            box(title = "Triplet combinations that will be scanned", width = 6,
                textOutput("NumCombos"),
                htmlOutput("GeneratedTriplets")
            )
    ),
    tabItem(tabName = "tripletAnalysis",
            box(title = "Parameters", status = "warning", solidHeader = TRUE,
             width = 12,
             fluidRow(
               box(title = "How to scan sequence similarity:", width = 4,
                   status = "danger",
                   numericInput("windowSize", "Size of sliding window (in bp):", min=1, value=100),
                   numericInput("stepSize", "Step size of sliding window (in bp):", min=1, value=1)
               ),
               box(title = "Block detection settings:", width = 4,
                   status = "danger",
                   radioButtons("detectionMethod", "Use a manual threshold?", c("Yes" = "yes", "No, automatically decide thresholds from data." = "no"), selected="no"),
                   checkboxInput("fallbackManual", "Fallback to a manual threshold?", value = TRUE),
                   conditionalPanel("input.detectionMethod == 'yes' || input.fallbackManual", 
                                    sliderInput("manBlockDetectDist", "Sequence Similarity Threshold:", 1, 100, value=95, step=1))
               ),
               box(title = "Block dating settings:", width = 4, status = "danger",
                   numericInput("mu", "Mutation Rate:", 10e-8),
                   numericInput("alpha", "Critical Value (alpha):", 0.05),
                   checkboxInput("bonf", "Bonferoni Correct Critical Value", TRUE),
                   checkboxInput("eliminateinsignificant", "Eliminate insignificant blocks", TRUE),
                   selectInput("correctionModel", "Mutation correction model:", c("JC69", "K80", "F81",
                                                                                  "K81", "F84", "BH87",
                                                                                  "T92", "TN93", "GG95"))
               )   
             ),
             fluidRow(
               box(title = "Start Analysis", width = 12,
                   status = "danger",
                   uiOutput("TripletToAnalyze"),
                   actionButton("analysisGO", "Run Analysis")
               )   
             ),
             fluidRow(
               box(title = "Plotting Settings", width = 12, status = "warning", collapsible = TRUE,
                   fluidRow(
                     box(width = 4,
                            h4("Plot title settings:"),
                            br(),
                            checkboxInput("plotTitle", "Include a title in the plot?", TRUE),
                            br(),
                            numericInput("plotTitleSize", "Font size of plot titles", 14),
                            br(),
                            textInput("plotTitleFace", "Font face of plot titles", "bold"),
                            br(),
                            textInput("plotTitleColour", "Colour of plot titles", "black"),
                            br(),
                            h4("Legend settings:"),
                            br(),
                            checkboxInput("plotLegends", "Add legends to the plots?", TRUE),
                            br(),
                            numericInput("plotLegendFontSize", "Font size for legend", 12)
                     ),
                     box(width = 4,
                            h4("X axis options:"),
                            br(),
                            checkboxInput("plotXTitle", "Include the title of the x-axis in the plots?", TRUE),
                            br(),
                            numericInput("plotXTitleFontSize", "Font size of x-axis title", 12),
                            br(),
                            textInput("plotXTitleColour", "Colour of x-axis title", "black"),
                            br(),
                            checkboxInput("plotXLabels", "Include the value labels of the x-axis in the plots?", TRUE),
                            br(),
                            numericInput("plotXLabelSize", "Font size of x-axis labels", 10),
                            br(),
                            textInput("plotXLabelColour", "Colour of x-axis labels", "black")
                     ),
                     box(width = 4,
                            h4("Y axis options:"),
                            br(),
                            checkboxInput("plotYTitle", "Include the title of the y-axis in the plots?", TRUE),
                            br(),
                            numericInput("plotYTitleFontSize", "Font size of y-axis title", 12),
                            br(),
                            textInput("plotYTitleColour", "Colour of y-axis title", "black"),
                            br(),
                            checkboxInput("plotYLabels", "Include the value labels of the y-axis in the plots?", TRUE),
                            br(),
                            numericInput("plotYLabelSize", "Font size of y-axis labels", 10),
                            br(),
                            textInput("plotYLabelColour", "Colour of y-axis labels", "black")
                     )
                   ),
                   fluidRow(
                     box(width = 12,
                         h4("Other Settings:"),
                         numericInput("plotMosaicScale", "Number of segments in RGB bars", 500),
                         numericInput("heightSave", "Height of plots when saved to file (inches)", 10),
                         numericInput("widthSave", "Width of plots when saved to file (inches)", 12),
                         numericInput("resSave", "Resolution of plots when saved (dpi)", 300),
                         radioButtons("whichPlotToSave", "Which plot to save?", list(Bars = "Bars", Lines = "Lines", Both = "Both"))
                     )
                   )
               )
             )
            ),
            box(title = "View Results", width = 12,
                uiOutput("TripletSelector"),
                downloadButton("saveBarPlots", "Save Plots"),
                downloadButton("saveTable", "Save Table"),
                hr(),
                fluidRow(
                box(title = "Heat Plot", width = 12, collapsable = TRUE,
                    plotOutput("barsPlot")
                )),
                fluidRow(
                box(title = "Lines Plot", width = 12, collapsable = TRUE,
                    plotOutput("linesPlot")
                )),
                dataTableOutput("blocksTable")
                
                  
                
            )
            
    ),
    tabItem(tabName = "userBlocks",
            fluidRow(
              box(title = "Add or clear user blocks", width = 6, status = "warning",
                  uiOutput("userBlocksSeqSelect"),
                  numericInput("startPosition", "Start of block in BP", value = NULL,
                               min = 1),
                  numericInput("endPosition", "End of block in BP", value = NULL,
                               min = 1),
                  actionButton("addUBButton", "Add user defined block between sequences"),
                  actionButton("clearUBButton", "Remove user defined blocks between sequences"),
                  actionButton("dateUBButton", "Test and date user defined blocks between sequences")
              ),
              box(title = "Block dating settings", width = 6, status = "warning",
                  numericInput("mu2", "Mutation Rate:", 10e-8),
                  numericInput("alpha2", "Critical Value (alpha):", 0.05),
                  checkboxInput("bonf2", "Bonferoni Correct Critical Value", TRUE),
                  checkboxInput("eliminateinsignificant2", "Eliminate insignificant blocks", TRUE),
                  selectInput("correctionModel2", "Mutation correction model:", c("JC69", "K80", "F81",
                                                                                  "K81", "F84", "BH87",
                                                                                  "T92", "TN93", "GG95"))
              )
            ),
            fluidRow(
              dataTableOutput("userBlocksTable")
            )
    )
  )
)

dashboardPage(header, sidebar, body, skin = "black")