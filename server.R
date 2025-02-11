library(shiny)
library(shinybusy)
library(ggplot2)
library(ggpubr)  
library(ggtext)
library(eefAnalytics)
library(lme4)
source('Main_functions.R')

# Define the server
server <- function(input, output, session) {
  simulatedData <- eventReactive(input$simulate, {
    nt <- input$nt
    
    if (input$simType == "CRT") {
      np <- input$np
      ns <- input$ns
      n_schools_treated <- as.numeric(unlist(strsplit(as.character(input$n_schools_treated), ",")))
      sigma <- input$sigma
      ICC <- input$ICC
      sigmaPret <- input$sigmaPret
      B0 <- input$B0
      B1 <- input$B1
      es <- as.numeric(unlist(strsplit(as.character(input$es), ",")))
      seed <- input$seed
      attrition_rates <- as.numeric(unlist(strsplit(as.character(input$attrition_rates), ",")))
      
      crtdata_simulation(nt, n_schools_treated, np, ns, sigma, ICC, sigmaPret, B0, B1, es, seed, attrition_rates)
      
    } else if (input$simType == "MST") {
      np <- input$np
      ns <- input$ns
      tpi <- as.numeric(unlist(strsplit(as.character(input$tpi), ",")))
      sigma <- input$sigma
      sigmab0 <- input$sigmab0
      sigmab1 <- input$sigmab1
      sigmaPret <- input$sigmaPret
      B0 <- input$B0
      B1 <- input$B1
      es <- as.numeric(unlist(strsplit(as.character(input$es), ",")))
      seed <- input$seed
      attrition_rates <- as.numeric(unlist(strsplit(as.character(input$attrition_rates), ",")))
      
      mstdata_simulation(nt, tpi, np, ns, sigma, sigmab0, sigmab1, sigmaPret, B0, B1, es, seed, attrition_rates)
      
    } else if (input$simType == "SRT") {
      np <- input$np_srt
      tpi <- as.numeric(unlist(strsplit(as.character(input$tpi_srt), ",")))
      sigma <- input$sigma_srt
      sigmaPret <- input$sigmaPret_srt
      B0 <- input$B0_srt
      B1 <- input$B1_srt
      es <- as.numeric(unlist(strsplit(as.character(input$es_srt), ",")))
      seed <- input$seed_srt
      attrition_rates <- as.numeric(unlist(strsplit(as.character(input$attrition_rates_srt), ",")))
      
      srtdata_simulation(nt, tpi, np, sigma, sigmaPret, B0, B1, es, seed, attrition_rates)
    }
  })
  
  importedData <- reactive({
    req(input$dataset)
    read.csv(input$dataset$datapath)
  })
  
  observe({
    data <- importedData()
    updateSelectInput(session, "post_var", choices = names(data))
    updateSelectInput(session, "intervention_var", choices = names(data))
    updateSelectInput(session, "random_var", choices = names(data))
  })
  
  output$post_var_select <- renderUI({
    selectInput("post_var", "Select Post-test Outcome", choices = NULL, multiple = TRUE)
  })
  
  output$intervention_var_select <- renderUI({
    selectInput("intervention_var", "Select Intervention Variables", choices = NULL, multiple = TRUE)
  })
  
  output$random_var_select <- renderUI({
    selectInput("random_var", "Select Clustering Variable", choices = NULL)
  })
  
  output$covariate_select_fut <- renderUI({
    selectInput("covariates_fut", "Select Additional Covariates", choices = names(importedData()), multiple = TRUE)
  })
  
  futilityData <- eventReactive(input$analyzeFutility, {
    show_modal_spinner(spin = "circle", text = "Analyzing futility data...")
    
    tryCatch({
      data <- importedData()  # Load the data
      post_vars <- input$post_var  # Get post-test variable
      intervention_column <- input$intervention_var  # Get intervention variable
      covariates <- input$covariates_fut  # Get selected covariates (multiple)
      
      # Ensure that only one post-test variable and one intervention variable is selected
      if (length(post_vars) != 1 || length(intervention_column) != 1) {
        stop("Please select exactly one post-test and one intervention variable.")
      }
      
      # Call the appropriate futility function based on the simulation type
      result <- if (input$futSimType == "CRT") {
        crtfutility(
          data = data, 
          post_vars = post_vars,
          intervention_column = intervention_column, 
          Random = input$random_var, 
          Nsim = input$crtNsim, 
          Threshold = input$crtThreshold, 
          ProbThreshold = input$crtProbThreshold, 
          covariates = covariates  
        )
      } else if (input$futSimType == "MST") {
        mstfutility(
          data = data, 
          post_vars = post_vars,
          intervention_column = intervention_column, 
          Random = input$random_var, 
          Nsim = input$crtNsim, 
          Threshold = input$crtThreshold, 
          ProbThreshold = input$crtProbThreshold, 
          covariates = covariates  
        )
      } else {
        srtfutility(
          data = data, 
          post_vars = post_vars,
          intervention_column = intervention_column, 
          Nsim = input$crtNsim, 
          Threshold = input$crtThreshold, 
          ProbThreshold = input$crtProbThreshold, 
          covariates = covariates  
        )
      }
      
      remove_modal_spinner()
      return(result)
      
    }, error = function(e) {
      remove_modal_spinner()
      showNotification(paste("Error: ", e$message), type = "error")
      return(NULL)  
    })
  })
  
  
  # Render the futility table
  output$futilityTable <- renderTable({
    futility_result <- futilityData()  
    
    if (is.null(futility_result) || nrow(futility_result) == 0) {
      showNotification("Futility data is empty or invalid.", type = "error")
      return(NULL)
    }
    
    names(futility_result)[names(futility_result) == "Treatment"] <- "Intervention"
    names(futility_result)[names(futility_result) == "ProbES"] <- "P(Effect size > Threshold)"
    
    futility_result$Futility <- sapply(1:nrow(futility_result), function(i) {
      if (futility_result$Futility[i] == 1) {
        paste("Intervention", futility_result$Intervention[i], "is futile.")
      } else {
        paste("Intervention", futility_result$Intervention[i], "is not futile.")
      }
    })
    
    futility_result <- futility_result[, c("Intervention", "P(Effect size > Threshold)", "Futility")]
    
    return(futility_result)
  })
  
  # Add a download button for futility results
  output$downloadFutilityData <- downloadHandler(
    filename = function() {
      paste("futility_analysis_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      futility_result <- futilityData()
      
      if (is.null(futility_result) || nrow(futility_result) == 0) {
        showNotification("No futility analysis results available for download.", type = "error")
        return(NULL)
      }
      
      names(futility_result)[names(futility_result) == "Treatment"] <- "Intervention"
      names(futility_result)[names(futility_result) == "ProbES"] <- "P(Effect size > Threshold)"
      
      # Save the data as CSV
      write.csv(futility_result, file, row.names = FALSE)
    }
  )
  
  
  
  # Event to generate and preview the futility plot
  futilityPlotData <- reactiveVal(NULL)
  
  observeEvent(input$plotFutility, {
    show_modal_spinner(spin = "circle", text = "Plotting posterior probabilities...")
    
    data <- importedData()
    nt <- length(input$post_var)
    
    futilityPlotData(list(
      nt = nt,
      data = data,
      Random = input$random_var,
      Nsim = if (input$futSimType == "CRT") input$crtNsim else input$mstNsim,
      Threshold = input$crtThreshold,
      ProbThreshold = input$crtProbThreshold,
      covariates = input$covariates_fut, 
      threshold_range = input$threshold_range
    ))
  })
  
  output$futilityPlot <- renderPlot({
    plot_data <- futilityPlotData()
    req(plot_data)
    
    if (input$futSimType == "CRT") {
      plot_futility_decision_shiny(
        data = plot_data$data,
        post_vars = input$post_var,
        intervention_column = input$intervention_var,
        Random = plot_data$Random,
        Nsim = plot_data$Nsim,
        covariates = plot_data$covariates,  
        VerticalLine = input$crtVerticalLine,  
        ProbThreshold = input$crtProbThreshold,  
        threshold_range = plot_data$threshold_range
      )
    } else if (input$futSimType == "MST") {
      plot_futility_decision_shinymst(
        data = plot_data$data,
        post_vars = input$post_var,
        intervention_column = input$intervention_var,
        Random = plot_data$Random,
        Nsim = plot_data$Nsim,
        covariates = plot_data$covariates,  
        VerticalLine = input$crtVerticalLine,  
        ProbThreshold = input$crtProbThreshold,  
        threshold_range = plot_data$threshold_range
      )
    } else {
      plot_futility_decision_shinysrt(
        data = plot_data$data,
        post_vars = input$post_var,
        intervention_column = input$intervention_var,
        Random = plot_data$Random,
        Nsim = plot_data$Nsim,
        covariates = plot_data$covariates,  
        VerticalLine = input$crtVerticalLine,  
        ProbThreshold = input$crtProbThreshold, 
        threshold_range = plot_data$threshold_range
      )
    }
    
    remove_modal_spinner()  
  })
  
  output$downloadFutilityPlot <- downloadHandler(
    filename = function() {
      paste("futility_plot-", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      plot_data <- futilityPlotData()
      req(plot_data)
      
      tiff(file, width = 9, height = 7, units = "in", res = 300)
      
      # Call the appropriate plotting function based on the simulation type
      if (input$futSimType == "CRT") {
        plot_futility_decision_shiny(
          data = plot_data$data,
          post_vars = input$post_var,
          intervention_column = input$intervention_var,
          Random = plot_data$Random,
          Nsim = plot_data$Nsim,
          covariates = plot_data$covariates,  
          VerticalLine = plot_data$VerticalLine,  
          ProbThreshold = plot_data$ProbThreshold, 
          threshold_range = plot_data$threshold_range
        )
      } else if (input$futSimType == "MST") {
        plot_futility_decision_shinymst(
          data = plot_data$data,
          post_vars = input$post_var,
          intervention_column = input$intervention_var,
          Random = plot_data$Random,
          Nsim = plot_data$Nsim,
          covariates = plot_data$covariates,  
          VerticalLine = plot_data$VerticalLine,   
          ProbThreshold = plot_data$ProbThreshold,
          threshold_range = plot_data$threshold_range
        )
      } else if (input$futSimType == "SRT") {
        plot_futility_decision_shinysrt(
          data = plot_data$data,
          post_vars = input$post_var,
          intervention_column = input$intervention_var,
          Random = plot_data$Random,
          Nsim = plot_data$Nsim,
          covariates = plot_data$covariates, 
          VerticalLine = plot_data$VerticalLine,   
          ProbThreshold = plot_data$ProbThreshold, 
          threshold_range = plot_data$threshold_range
        )
      }
      
      dev.off()
    }
  )
  
  importedMultilevelData <- reactive({
    req(input$multilevelDataset)
    read.csv(input$multilevelDataset$datapath)
  })
  
  observe({
    data <- importedMultilevelData()
    updateSelectInput(session, "post_var_multilevel", choices = names(data))
    updateSelectInput(session, "intervention_var_multilevel", choices = names(data))
    updateSelectInput(session, "random_var_multilevel", choices = names(data))
  })
  
  output$post_var_select_multilevel <- renderUI({
    selectInput("post_var_multilevel", "Select Post-test Outcome", choices = NULL, multiple = TRUE)
  })
  
  output$intervention_var_select_multilevel <- renderUI({
    selectInput("intervention_var_multilevel", "Select Intervention Variables", choices = NULL, multiple = TRUE)
  })
  
  
  output$random_var_select_multilevel <- renderUI({
    conditionalPanel(
      condition = "input.method != 'srtBayes' && input.method != 'srtFREQ'",
      selectInput("random_var_multilevel", "Select Clustering Variable", choices = NULL)
    )
  })
  output$covariate_select_multilevel <- renderUI({
    selectInput("covariates_multilevel", "Select Additional Covariates", choices = names(importedMultilevelData()), multiple = TRUE)
  })
  
  importedSupData <- reactive({
    req(input$supDataset)
    read.csv(input$supDataset$datapath)
  })
  
  # Update the available columns for post_var, intervention_var, and covariates dynamically
  observe({
    data <- importedSupData()
    updateSelectInput(session, "post_var_sup", choices = names(data))  # Single select for post_var
    updateSelectInput(session, "intervention_var_sup", choices = names(data))  # Single select for intervention_var
    updateSelectInput(session, "random_var_sup", choices = names(data))  # Cluster variable
    updateSelectInput(session, "covariates_sup", choices = names(data))  # Multiple select for covariates
  })
  
  output$post_var_select_sup <- renderUI({
    selectInput("post_var_sup", "Select Post-test Outcome", choices = NULL)
  })
  
  output$intervention_var_select_sup <- renderUI({
    selectInput("intervention_var_sup", "Select Intervention Variable", choices = NULL)
  })
  
  output$random_var_select_sup <- renderUI({
    conditionalPanel(
      condition = "input.supSimType != 'SRT'",
      selectInput("random_var_sup", "Select Clustering Variable", choices = NULL)
    )
  })
  
  output$covariate_select_sup <- renderUI({
    selectInput("covariates_sup", "Select Additional Covariates", choices = NULL, multiple = TRUE)
  })
  
  superiorityData <- eventReactive(input$analyzeSuperiority, {
    show_modal_spinner(spin = "circle", text = "Analyzing superiority data...")
    
    tryCatch({
      data <- importedSupData()
      post_var <- input$post_var_sup  
      intervention_var <- input$intervention_var_sup 
      covariates <- input$covariates_sup  
      superiority_threshold <- input$crtSupSuperiorThreshold  
      
      if (is.null(post_var) || is.null(intervention_var)) {
        stop("Please select both a post-test variable and an intervention variable.")
      }
      
      # Superiority analysis based on the simulation type
      result <- if (input$supSimType == "CRT") {
        crtSuperiority(
          data = data, 
          post_var = post_var, 
          intervention_column = intervention_var, 
          Random = input$random_var_sup, 
          Nsim = input$crtSupNsim, 
          Threshold = input$crtSupThreshold, 
          reference_intervention = input$crtSupReference, 
          superiority_threshold = superiority_threshold,  
          covariates = covariates
        )
      } else if (input$supSimType == "MST") {
        mstSuperiority(
          data = data, 
          post_var = post_var, 
          intervention_column = intervention_var, 
          Random = input$random_var_sup, 
          Nsim = input$crtSupNsim, 
          Threshold = input$crtSupThreshold, 
          reference_intervention = input$crtSupReference, 
          superiority_threshold = superiority_threshold, 
          covariates = covariates
        )
      } else {
        srtSuperiority(
          data = data, 
          post_var = post_var, 
          intervention_column = intervention_var, 
          Nsim = input$crtSupNsim, 
          Threshold = input$crtSupThreshold, 
          reference_intervention = input$crtSupReference, 
          superiority_threshold = superiority_threshold,  
          covariates = covariates
        )
      }
      
      remove_modal_spinner() 
      return(result)
      
    }, error = function(e) {
      remove_modal_spinner()
      showNotification(paste("Error: ", e$message), type = "error")
      return(NULL)  # Return NULL if error occurs
    })
  })
  
  
  # Display the results of Superiority Analysis
  output$superiorityTable <- renderTable({
    superiority_result <- superiorityData()
    
    if (is.null(superiority_result) || nrow(superiority_result) == 0) {
      showNotification("No results available. Please check your input.", type = "error")
      return(NULL)
    }
    
    names(superiority_result)[names(superiority_result) == "Treatment"] <- "Intervention"
    
    # Rename the "ProbES" column to "P (Effect size > Threshold)"
    names(superiority_result)[names(superiority_result) == "ProbES"] <- "P(Effect size > Threshold)"
    
    # Modify the "Superiority" column:
    superiority_result$Superiority <- ifelse(
      superiority_result$Superiority == "Reference", 
      "Reference",  # Keep "Reference" unchanged
      ifelse(
        superiority_result$Superiority == "Superior",
        "Superior to the Reference Intervention",
        "Not Superior to the Reference Intervention"
      )
    )
    
    return(superiority_result)
  })
  
  # Download handler for CSV file
  output$downloadSuperiorityData <- downloadHandler(
    filename = function() {
      paste("superiority_analysis_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      superiority_result <- superiorityData()
      
      if (is.null(superiority_result) || nrow(superiority_result) == 0) {
        showNotification("No superiority analysis results available for download.", type = "error")
        return(NULL)
      }
      
      names(superiority_result)[names(superiority_result) == "Treatment"] <- "Intervention"
      names(superiority_result)[names(superiority_result) == "ProbES"] <- "P(Effect size > Threshold)"
      
      superiority_result$Superiority <- ifelse(
        superiority_result$Superiority == "Reference", 
        "Reference",  # Keep "Reference" unchanged
        ifelse(superiority_result$Superiority == "Superior", 1, 0)
      )
      
      write.csv(superiority_result, file, row.names = FALSE)
    }
  )
  
  
  
  
  output$dataTable <- renderTable({
    head(simulatedData(),n=10)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("simulated_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(simulatedData(), file, row.names = FALSE)
    }
  )
  
  analysisResult <- eventReactive(input$runAnalysis, {
    show_modal_spinner(spin = "circle", text = "Analyzing data...")
    
    tryCatch({
      data <- importedMultilevelData()
      post_vars <- input$post_var_multilevel
      intervention_column <- input$intervention_var_multilevel
      random_var <- input$random_var_multilevel
      covariates <- input$covariates_multilevel  # Handle multiple covariates
      
      if (length(post_vars) != 1 || length(intervention_column) != 1) {
        stop("Please select exactly one post-test and one intervention variable.")
      }
      
      if (is.null(covariates) || length(covariates) == 0) {
        covariates <- NULL  # No covariates provided
      }
      
      # Call the run_analysis function with user inputs
      result <- run_analysis(
        data = data, 
        post_vars = post_vars,
        intervention_column = intervention_column, 
        Random = random_var, 
        Nsim = input$multilevelNsim, 
        Threshold = input$multilevelThreshold, 
        method = input$method, 
        crtFREQoption = input$crtFREQoption, 
        nPerm = input$nPerm, 
        nBoot = input$nBoot, 
        bootType = input$bootType, 
        covariates = covariates 
      )
      
      remove_modal_spinner()
      return(result)
    }, error = function(e) {
      remove_modal_spinner()
      showNotification(paste("Error: ", e$message), type = "error")
      return(NULL)  # Return NULL in case of error
    })
  })
  
  
  output$multilevelPlot <- renderPlot({
    analysisResult()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("multilevel_analysis_plot-", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      tiff(file, width = 12, height = 9, units = 'in', res = 300)
      print(analysisResult())
      dev.off()
    }
  )
  
  plotTriggered <- reactiveVal(FALSE)
  
  # Plot Posterior Probabilities Tab
  importedPlotData <- reactive({
    req(input$plotDataset)
    read.csv(input$plotDataset$datapath)
  })
  
  observe({
    data <- importedPlotData()
    updateSelectInput(session, "post_var_plot", choices = names(data))
    updateSelectInput(session, "intervention_var_plot", choices = names(data))
    updateSelectInput(session, "random_var_plot", choices = names(data))
  })
  
  output$post_var_select_plot <- renderUI({
    selectInput("post_var_plot", "Select Post-test Outcome", choices = NULL)
  })
  
  output$intervention_var_select_plot <- renderUI({
    selectInput("intervention_var_plot", "Select Intervention Variable", choices = NULL)
  })
  
  output$random_var_select_plot <- renderUI({
    selectInput("random_var_plot", "Select Clustering Variable", choices = NULL)
  })
  
  output$covariate_select_plot <- renderUI({
    selectInput("covariates_plot", "Select Additional Covariates", choices = names(importedPlotData()), multiple = TRUE)
  })
  
  
  # Observe when the plotPosterior button is clicked
  observeEvent(input$plotPosterior, {
    plotTriggered(TRUE)
    show_modal_spinner(spin = "circle", text = "Plotting posterior probabilities...")
  })
  
  # Cache the plot data using a reactive expression
  plotPosteriorData <- reactive({
    req(plotTriggered())  # Ensure that this only runs after being triggered
    
    data <- importedPlotData()
    post_var <- input$post_var_plot
    intervention_var <- input$intervention_var_plot
    covariates <- input$covariates_plot  # List of covariates
    
    # Validate that the post-test and intervention variables are provided
    if (is.null(post_var) || is.null(intervention_var)) {
      stop("Please select both a post-test variable and an intervention variable.")
    }
    
    list(
      data = data,
      post_var = post_var,
      intervention_var = intervention_var,
      Random = if (input$plotSimType != "SRT") input$random_var_plot else NULL,
      Nsim = input$crtPlotNsim,
      VerticalLine = if (input$addVerticalLine) input$crtVerticalLine else NULL,  # Vertical line
      ProbThreshold = if (input$addHorizontalLine) input$crtPlotProbThreshold else NULL,  # Horizontal line
      covariates = covariates,  # Covariates
      threshold_range = input$plot_threshold_range  # Plot threshold range
    )
  })
  
  cachedPlot <- reactive({
    show_modal_spinner(spin = "circle", text = "Plotting posterior probabilities...")
    tryCatch({
      plot_data <- plotPosteriorData()
      req(plot_data)
      
      if (input$plotSimType == "CRT") {
        plot_futility_decision_shiny(
          data = plot_data$data,
          post_vars = plot_data$post_var,
          intervention_column = plot_data$intervention_var,
          Random = plot_data$Random,
          Nsim = plot_data$Nsim,
          covariates = plot_data$covariates,  # Covariates
          VerticalLine = plot_data$VerticalLine, 
          ProbThreshold = plot_data$ProbThreshold, 
          threshold_range = plot_data$threshold_range
        )
      } else if (input$plotSimType == "MST") {
        plot_futility_decision_shinymst(
          data = plot_data$data,
          post_vars = plot_data$post_var,
          intervention_column = plot_data$intervention_var,
          Random = plot_data$Random,
          Nsim = plot_data$Nsim,
          covariates = plot_data$covariates,  # Covariates
          VerticalLine = plot_data$VerticalLine, 
          ProbThreshold = plot_data$ProbThreshold, 
          threshold_range = plot_data$threshold_range
        )
      } else {
        plot_futility_decision_shinysrt(
          data = plot_data$data,
          post_vars = plot_data$post_var,
          intervention_column = plot_data$intervention_var,
          Nsim = plot_data$Nsim,
          covariates = plot_data$covariates,  # Covariates
          VerticalLine = plot_data$VerticalLine, 
          ProbThreshold = plot_data$ProbThreshold, 
          threshold_range = plot_data$threshold_range
        )
      }
      remove_modal_spinner()  # Remove spinner after the plot is generated
      recordPlot()  # Capture the plot  
    }, error = function(e) {
      remove_modal_spinner()  # Ensure the spinner is removed in case of an error
      showNotification(paste("Error: ", e$message), type = "error")  # Show error notification
      return(NULL)  # Return NULL in case of an error
    })
  })
  
  
  # Render the plot using the cached plot
  output$posteriorPlot <- renderPlot({
    plot <- cachedPlot()
    remove_modal_spinner()  # Remove the spinner after the plot is rendered
    plot  # Simply render the cached plot
  })
  
  output$downloadPosteriorPlot <- downloadHandler(
    filename = function() {
      paste("posterior_plot-", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      tiff(file, width = 7, height = 5, units = "in", res = 300)
      replayPlot(cachedPlot())  # Replay the cached plot for download
      dev.off()
    }
  )
  
  
  output$dataTable <- renderTable({
    req(simulatedData())
    head(simulatedData(), n = 10)
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("simulated_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(simulatedData(), file, row.names = FALSE)
    }
  )
  importedNewTreatmentData <- reactive({
    req(input$newTreatmentDataset)
    data <- read.csv(input$newTreatmentDataset$datapath)
    
    cat("Imported dataset columns:\n")
    print(names(data))
    
    return(data)
  })
  
  observe({
    data <- importedNewTreatmentData()
    
    updateSelectInput(session, "post_var_new", choices = names(data))
    updateSelectInput(session, "intervention_var_new", choices = names(data))
    updateSelectInput(session, "covariates_new", choices = names(data))
    
    # Pupils Variable (ID for SRT, Pupils for CRT/MST)
    updateSelectInput(session, "pupils_var_new", choices = names(data))
    
    # Only update school selection when the simulation type is CRT or MST
    if (input$newSimType %in% c("CRT", "MST")) {
      updateSelectInput(session, "schools_var_new", choices = names(data))
    }
  })
  
  output$post_var_select_new <- renderUI({
    selectInput("post_var_new", "Select Post-test Outcome", choices = NULL)
  })
  
  output$intervention_var_select_new <- renderUI({
    selectInput("intervention_var_new", "Select Intervention Variable", choices = NULL)
  })
  
  output$pupils_var_select_new <- renderUI({
    label_text <- if (input$newSimType == "SRT") "Select ID Variable" else "Select Pupils Variable"
    selectInput("pupils_var_new", label_text, choices = NULL)
  })
  
  output$covariate_select_new <- renderUI({
    selectInput("covariates_new", "Select Additional Covariates", choices = NULL, multiple = TRUE)
  })
  
  output$schools_var_select_new <- renderUI({
    if (input$newSimType %in% c("CRT", "MST")) {
      selectInput("schools_var_new", "Select Schools Variable", choices = NULL)
    }
  })
  
  validate_columns <- function(data, required_columns) {
    cat("Validating selected columns:\n")
    print(required_columns)
    
    missing_columns <- setdiff(required_columns, names(data))
    if (length(missing_columns) > 0) {
      stop(paste("The following columns are missing in the dataset:", paste(missing_columns, collapse = ", ")))
    }
  }
  
  # Event handler for adding a new treatment
  newTreatmentData <- eventReactive(input$addNewTreatment, {
    show_modal_spinner(spin = "circle", text = "Adding new treatment...")
    
    tryCatch({
      data <- importedNewTreatmentData()
      
      post_col <- input$post_var_new
      intervention_col <- input$intervention_var_new
      schools_col <- input$schools_var_new
      pupils_col <- input$pupils_var_new
      covariates <- input$covariates_new
      
      if (is.null(post_col) || is.null(intervention_col) || is.null(pupils_col)) {
        stop("Please make sure all column selections are made.")
      }
      
      required_columns <- c(post_col, intervention_col, pupils_col, covariates)
      validate_columns(data, required_columns)
      
      # Generate new intervention depending on the simulation type
      result <- switch(input$newSimType,
                       "CRT" = add_new_treatmentcrt(
                         existing_data = data,
                         new_schools = input$newSchools,
                         new_pupils_per_school = input$newPupilsPerSchool,
                         es = input$newEffectSize,
                         attrition_rate = input$newAttrition,
                         post_col = post_col,
                         intervention_col = intervention_col,
                         schools_col = schools_col,
                         pupils_col = pupils_col,
                         covariates = covariates
                       ),
                       "MST" = add_new_treatmentmst(
                         existing_data = data,
                         new_schools = input$newSchools,
                         new_pupils_per_school = input$newPupilsPerSchool,
                         es = input$newEffectSize,
                         attrition_rate = input$newAttrition,
                         treatment_percentage = input$treatmentPercentageMST,
                         post_col = post_col,
                         intervention_col = intervention_col,
                         schools_col = schools_col,
                         pupils_col = pupils_col,
                         covariates = covariates
                       ),
                       "SRT" = add_new_treatmentsrt(
                         existing_data = data,
                         new_pupils = input$newPupilsPerSchool * input$newSchools, # Assuming SRT adds all pupils at once
                         es = input$newEffectSize,
                         attrition_rate = input$newAttrition,
                         post_col = post_col,
                         intervention_col = intervention_col,
                         pupils_col = pupils_col,
                         covariates = covariates
                       )
      )
      
      if (nrow(result) > 1e6) {  # Example threshold of 1 million rows
        stop("The generated dataset is too large. Please adjust input parameters.")
      }
      
      remove_modal_spinner()
      return(result)
      
    }, error = function(e) {
      remove_modal_spinner()
      showNotification(paste("Error: ", e$message), type = "error")
      return(NULL)
    })
  })
  
  
  
  output$newTreatmentTable <- renderTable({
    req(newTreatmentData())
    head(newTreatmentData(), n = 10)
  })
  
  output$downloadNewTreatmentData <- downloadHandler(
    filename = function() {
      paste("new_treatment_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(newTreatmentData(), file, row.names = FALSE)
    }
  )
}