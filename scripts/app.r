

library(shiny)
library(shinyjs)
library(shinyFiles)
library(tools)
library(data.table)
library(stringr)
library(igraph)


source("./scripts/utils.r")
source("./scripts/main.r")


server <- function(input, output, session) {


	shinyDirChoose(input,'cifdir',roots = c(home = '~'),filetypes = c('', 'cif'))

	dir <- reactive(input$cifdir)
	output$cifdir <- renderText({  # use renderText instead of renderPrint
		parseDirPath(c(home = '~'), dir())
	})

	observe(
    {
     if(is.null(input$annotationFile))
       {
        disable("script")
       }
     else
       {
         enable("script")
       }
     }
     )

	observeEvent(input$script, {

		home <- normalizePath("~")
        inDir <- file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))

		inFile <- input$annotationFile$datapath

	    # Function to toggle input elements. 
		# input_list: List of inputs, reactiveValuesToList(input)
		# enable_inputs: Enable or disable inputs?
		# Only buttons: Toggle all inputs, or only buttons?
		toggle_inputs <- function(input_list,enable_inputs=T,only_buttons=FALSE)
		{
		# Subset if only_buttons is TRUE.
		if(only_buttons){
		  buttons <- which(sapply(input_list,function(x) {any(grepl('Button',attr(x,"class")))}))
		  input_list = input_list[buttons]
		}

		# Toggle elements
		for(x in names(input_list))
		  if(enable_inputs){
		    shinyjs::enable(x)} else {
		      shinyjs::disable(x) }
		}

		updateProgress <- function(value = NULL, detail = NULL) 
		{
	      if (is.null(value)) {
	        value <- progress$getValue()
	        value <- value + (progress$getMax() - value) / 5
	      }
	      progress$set(value = value, detail = detail)
	    }

		input_list <- reactiveValuesToList(input)
		toggle_inputs(input_list,F,F)


	    # Create a Progress object
	    progress <- shiny::Progress$new()
	    progress$set(message = "Computing...", value = 0)
	    


		runmain(input$cutoff, input$aminoacid, input$task, inFile, inDir, updateProgress) 


		# Close the progress when this reactive exits (even if there's an error)
	    on.exit(progress$close())

		toggle_inputs(input_list,T,F)

	})


	session$onSessionEnded(function() {
        stopApp()
      })


}



ui <- fluidPage(


	titlePanel( h2("DynamicProtein", align = "center") ),

	sidebarLayout(
    
    sidebarPanel(

    	width = 3,

		numericInput("cutoff", ,label = h3("Distance cutoff"),value = 6, min = 4, max = 10),

		numericInput("aminoacid", label= h3("Number of amino acids"),value = 5, min = 5, max = 20),

		selectInput("task",label = h3("Select task"), choices = c("Mode 1","Mode 2","Mode 3")),

		shinyDirButton("cifdir", "Choose CIF folder", "Upload"),
		verbatimTextOutput("cifdir", placeholder = TRUE),

		fileInput("annotationFile", label= h3("Annotation file"),
		accept = c("text/tab-separated-values",".txt")),

		actionButton("script","Run program"),

		useShinyjs()

    ), 


    mainPanel(

		width = 7,
		tableOutput('upload'),

		strong("Distance cutoff:"),
		p("To make a PSN, amino acids are taken as nodes and any two amino acids that are close enough are connected by an edge. The 'distance cutoff' parameter defines the largest 3-dimensional (3D) distance (in Angstroms) between any two amino acids in a protein domain to be considered as close enough. The application currently has seven different values to choose from, i.e., 4,5,6,...,10. In the original study, the value of 6 was used."),
		br(),

		strong("Number of amino acids:"),
		p("A dynamic protein structure network (PSN) consists of two or more static PSN snapshots.
			The 'number of amino acids' parameter defines the number of amino acids that are added
		at each static protein structure network (PSN) snapshot, as the dynamic PSN is being created.
			In more detail, given a value 'X' for the 'number of amino acids' parameter, a dynamic PSN is created as follows. Starting from the begining of a given protein sequence (i.e., the first amino acid towards the N-terminal 
			of the protein), first X number of amino acids are converted into the first static PSN snapshot, then first 2X number of amino acids are converted into the second static PSN snapshot, then first 3X number of amino acids are converted into the thrid static PSN snapshot, and so on until all of the amino acids in the sequence are captured in the last snapshot of the corresponding dynamic PSN. In the original study, a value of 5 was used. In the application 16 different the values (5,6,7,...,20) can be used."),
		br(),

		strong("Select task:"),
		p("This parameter controls the extent of the application that the user wants to run.
			There are three options: Mode 1, Mode 2, and Mode 3.
			If a user wants to just create dynamic networks, then the option 'Mode 1' should be chosen.
			If a user wants to create dynamic networks and also extract dynamic graphlet-based features of the dynamic networks, then the option 'Mode 2' should be chosen.
			If a user wants to create dynamic networks, extract dynamic graphlet-based features of the dynamic networks, and classify the protein domains using the resulting dynamic features,
			 then the option 'Mode 3' should be chosen.
			"),
		br(),

		strong("Select CIF folder:"),
		p("This is a path to the directory where all of the Cyrstallographic Information File (CIF) of the protein 3-dimensional structures are stored in your local computing machine.
			"),
		br(),

		strong("Annotation file:"),
		p("This is a path to the file with a list of protein domains to be processed. The file should have a '.txt' extension.
			The contents of the file should be in a two column format, with the two columns being separated by a 'tab'.
			The first column should be the structual class of the corresponding protein domain in the second column. 
			This structural class information is only used in the protein structural classification part of the software. 
			If you do not have structural class information of the proteins in your dataset and if you do not wish to perform 
			protein structural classification, then please add a dummy structural class for each of the proteins in your list.
			In the second column, the name of the protein domain should be given as 'CIFID_chainID_startID1_endID1-startID2_endID2-...-startIDn_endIDn'. 
			Here, CIFID is the name of the CIF file. For example, for the CIF file '1fnn.cif' the CIFID is '1fnn'.
			chainID is the idenfication of the protein chain to which the corresponding protein domain belongs. startIDs and ednIDs represent the starting and ending sequence number for each of the contiguous amino acid sequence segments of the corresponding protein domain. For example, if a protein domain has two segments that forms that domain, then startID1 will represent the sequence number of the chainID where the corresponding first segment starts, and endID1 will represent the sequence number of the chainID where the corresponding first segment ends. Additionally, startID2 will represent the sequence number of the chainID where the corresponding second segment starts, and endID2 will represent the sequence number of the chainID where the corresponding second segment ends.
			For example, sequence number 1 to 17 in the chain A of the CIF file '1fnn.cif' belongs to the segment 1 and 
			sequence number 192 to 275 in the chain A of the CIF file '1fnn.cif' belongs to the segment 2 of the same protein domain protein that belongs to the structural class of 'alpha'. 
			In this case, the corresponsing line in the annotation file should be 'alpha	[tab separation] 1fnn_A_1_17-192_275'. As another example, sequence number 2 to 54 in the chain A of the CIF file '1aip.pdb' belongs to the segment 1, which is the only segment of contiugous amino acids of the protein domain, belongs to the structural class of 'alpha'. In this case, the corresponsing line in the annotation file should be 'alpha	[tab separation] 1aip_C_2_54'.
			If the protein domain you want to process spans the whole protein chain, then please write the first sequence number of the chain as startID and the last sequence number of the chain as the endID."
			
			),



		

    )

  )

)



shinyApp(ui = ui, server = server, options=list(port=6008))


