

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

		width = 5,
		tableOutput('upload'),

		strong("Distance cutoff:"),
		p("This parameter defines the largest 3-dimensional (3D) distance (in Angstroms) between any two amino acids
		 in a protein domain to be considered as close enough. To make a PSN, amino acids are taken as nodes and any 
		 two amino acids are are close enough are connected by an edge. The application currently has seven different 
		 values to choose from, i.e., 4,5,6,...,10. In the original study, the value of 6 was used."),
		br(),

		strong("Number of amino acids:"),
		p("A dynamic protein structure network (PSN) consists of two or more static PSN snapshots.
			The 'number of amino acids' parameter defines the number of amino acids that are added
		at each static protein structure network (PSN) snapshot, as the dynamic PSN is being created.
			In more detail, given a value 'X' for the 'number of amino acids' parameter, a dynamic PSN is created as follows. 
			Starting from the begining of a given protein sequence (i.e., the first amino acid towards the N-terminal 
			of the protein), first X number of amino acids are converted into the first static PSN snapshot, then first 2X 
			number of amino acids are converted into the second static PSN snapshot, then first 3X number of amino acids 
			are converted into the thrid static PSN snapshot, and so on until all of the amino acids are captured in the 
			last snapshot of the corresponding dynamic PSN. In the original study, a value of 5 was used. In the application
			 16 different the values (5,6,7,...,20) can be used."),
		br(),

		strong("Select task:"),
		p("This parameter controls the extent of the application that the user wants to run.
			There are three options: Mode 1, Mode 2, and Mode 3.
			If a user wants to just create the dynamic networks, then the option 'Mode 1' should be chosen.
			If a user wants to create the dynamic networks and also extract dynamic graphlet-based features of the dynamic networks,
			 then the option 'Mode 2' should be chosen.
			If a user wants to create the dynamic networks, extract dynamic graphlet-based features of the dynamic networks, and 
			classify the protein doamins using the resulting dynamic features,
			 then the option 'Mode 3' should be chosen.
			"),
		br(),

		strong("Select CIF folder:"),
		p("This is a path to the directory where all of the Cyrstallographic In formation File 
			(CIF) of the protein 3-dimansional structures are stored in your local computing machine.
			"),
		br(),

		strong("Annotation file:"),
		p("This is a path to the file with a list of protein domains to be processed. The file should have a '.txt' extension.
			The contents of the file should be in a two column format, with the two columns being separated by a 'tab'.
			The first column should be the structual class of the corresponding protein domain in the second column. 
			This structural class information is only used in the protein structural classification part of the software. 
			If you do not want the structural part of the proteins in your set and if you do not wish to perform 
			protein structural classification, then please add a dummy structural class for each of the proteins in your list.
			In the second column, the name of the protein domain should be given as 'PDBID_chainID_startID_endID'. 
			Here, PDBID is the name of the PDB file for example, for the PDB file '155c.pdb' the PDBID is '155c'.
			chainID is the idenfication of the protein chain to which the corresponding protein domain belongs.
			startID is the sequence number of the chainID where the corresponding protein domain starts.
			endID is the sequence number of the chainID where the corresponding protein domain ends. 
			For example, sequence number 1 to 134 in the chain A of the PDB file '155c.pdb' belongs to the protein structural class of 
			'alpha'. In this case, the corresponsing line in the annotation file would be 'alpha	[tab separation] 155c_A_1_134'.
			If the protein domain you want to process spans the whole protein chain, then please write the first sequence number 
			of the chain as startID and the last sequence number of the chain as the endID."
			
			),



		

    )

  )

)



shinyApp(ui = ui, server = server, options=list(port=6008))


