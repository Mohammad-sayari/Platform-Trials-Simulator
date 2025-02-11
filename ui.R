library(shiny)
library(shinybusy)
library(ggplot2)
library(ggpubr)  # For theme_pubclean()
library(ggtext)
library(eefAnalytics)
library(lme4)


# Define the UI
ui <- fluidPage(
  # Title Row
  fluidRow(
    column(
      width = 12,
      div(
        style = "text-align: left; margin-bottom: 10px;",  # Reduced margin-bottom for less spacing
        h1("Platform Trials Simulator", 
           style = "font-size: 38px; font-weight: bold; color: #333; margin: 0; padding: 0;")  # Inline style with reduced margins/padding
      )
    )
  ),
  
  # Tabs Row (Placed in a new row)
  fluidRow(
    column(
      width = 12,
      navbarPage(
        "",
        tabPanel(
          "Data Simulation",
          sidebarLayout(
            sidebarPanel(
              radioButtons("simType", "Simulation Type", choices = c("CRT", "MST", "SRT")),
              
              # Number of Interventions - This is common for all types
              numericInput("nt", "Number of Interventions (excluding control)", 2, min = 1),
              
              # Conditional inputs for SRT
              conditionalPanel(
                condition = "input.simType == 'SRT'",
                numericInput("np_srt", "Number of Participants", 100, min = 10),
                textInput("tpi_srt", "Percentage of Participants in Each Group (control, Intervention1, ...)", "50,30,20"),
                numericInput("sigma_srt", "Residual Standard Deviation", 1),
                numericInput("sigmaPret_srt", "Standard Deviation of Pre-test Scores", 1),
                numericInput("B0_srt", "Intercept", 0),
                numericInput("B1_srt", "Pre-test Coefficient", 0.5),
                textInput("es_srt", "Effect Sizes (comma-separated for each Intervention)", "0.2,0.3"),
                numericInput("seed_srt", "Random Seed", 1234),
                textInput("attrition_rates_srt", "Attrition Rates (comma-separated for control, Intervention1, ...)", "0.1,0.2,0.3")
              ),
              
              # Conditional inputs for CRT
              conditionalPanel(
                condition = "input.simType == 'CRT'",
                numericInput("np", "Number of Pupils Per School", 100, min = 10),
                numericInput("ns", "Number of Schools", 10, min = 1),
                textInput("n_schools_treated", "Number of schools treated (control, Intervention1, ...)", "5,3,2"),
                numericInput("sigma", "Residual Standard Deviation", 1),
                numericInput("ICC", "Intraclass Correlation Coefficient", 0.1),
                numericInput("sigmaPret", "Standard Deviation of Pre-test Scores", 1),
                numericInput("B0", "Intercept", 0),
                numericInput("B1", "Pre-test Coefficient", 0.5),
                textInput("es", "Effect Sizes (comma-separated for each Intervention)", "0.2,0.3"),
                numericInput("seed", "Random Seed", 1234),
                textInput("attrition_rates", "Attrition Rates (comma-separated for control, Intervention1, ...)", "0.1,0.1,0.1")
              ),
              
              # Conditional inputs for MST
              conditionalPanel(
                condition = "input.simType == 'MST'",
                numericInput("np", "Number of Pupils Per School", 100, min = 10),
                numericInput("ns", "Number of Schools", 10, min = 1),
                textInput("tpi", "Percentage of pupils in each group (control, Intervention1, ...)", "50,30,20"),
                numericInput("sigma", "Residual Standard Deviation", 1),
                numericInput("sigmab0", "Random Intercept Standard Deviation", 0.5),
                numericInput("sigmab1", "Random Slope Standard Deviation", 0.5),
                numericInput("sigmaPret", "Standard Deviation of Pre-test Scores", 1),
                numericInput("B0", "Intercept", 0),
                numericInput("B1", "Pre-test Coefficient", 0.5),
                textInput("es", "Effect Sizes (comma-separated for each Intervention)", "0.2,0.3"),
                numericInput("seed", "Random Seed", 1234),
                textInput("attrition_rates", "Attrition Rates (comma-separated for control, Intervention1, ...)", "0.1,0.1,0.1")
              ),
              
              # Simulate button
              actionButton("simulate", "Simulate Data")
            ),
            mainPanel(
              tabsetPanel(
                tabPanel(
                  "Data",
                  tableOutput("dataTable"),
                  downloadButton("downloadData", "Download Data")
                )
              )
            )
          )
        ),
        tabPanel(
          "Multilevel Analysis",
          sidebarLayout(
            sidebarPanel(
              radioButtons("method", "Analysis Method", choices = c("crtBayes", "crtFREQ", "mstBayes", "mstFREQ", "srtBayes", "srtFREQ")),
              fileInput("multilevelDataset", "Choose CSV File", accept = ".csv"),
              uiOutput("post_var_select_multilevel"),
              uiOutput("intervention_var_select_multilevel"),
              uiOutput("random_var_select_multilevel"),
              conditionalPanel(
                condition = "input.method == 'crtBayes' || input.method == 'mstBayes' || input.method == 'srtBayes'",
                numericInput("multilevelNsim", "Number of Simulations (>= 10000 is recommended)", 10000, min = 1),
                numericInput("multilevelThreshold", "Threshold", 0.05, min = 0, max = 1)
              ),
              conditionalPanel(
                condition = "input.method == 'crtFREQ' || input.method == 'mstFREQ' || input.method == 'srtFREQ'",
                selectInput("crtFREQoption", "Select Confidence Interval Calculation Method", choices = c("Analytic (Default)" = "Default", "Permutation", "Bootstrap")),
                conditionalPanel(
                  condition = "input.crtFREQoption == 'Permutation'",
                  numericInput("nPerm", "Number of Permutations", 1000, min = 1)
                ),
                conditionalPanel(
                  condition = "input.crtFREQoption == 'Bootstrap'",
                  numericInput("nBoot", "Number of Bootstraps", 1000, min = 1)
                ),
                
                # bootType appears only for crtFREQ and mstFREQ when Bootstrap is selected
                conditionalPanel(
                  condition = "input.crtFREQoption == 'Bootstrap' && (input.method == 'crtFREQ' || input.method == 'mstFREQ')",
                  selectInput("bootType", "Bootstrap Type", choices = c(
                    "Re-sampling at student level" = "case(1)",
                    "Re-sampling at school level" = "case(2)",
                    "Re-sampling at both levels" = "case(1,2)",
                    "Residual bootstrapping" = "residual"
                  ))
                )
              ),
              # Remove "Number of Covariates"
              # numericInput("num_covariates", "Number of Covariates", 1, min = 0), -- Removed
              
              # Allow selection of multiple covariates
              uiOutput("covariate_select_multilevel"),  # Updated to allow multiple selections
              
              actionButton("runAnalysis", "Run Analysis")
            ),
            mainPanel(
              plotOutput("multilevelPlot"),
              downloadButton("downloadPlot", "Download Plot")
            )
          )
        )  ,
        tabPanel(
          "Futility Analysis",
          sidebarLayout(
            sidebarPanel(
              radioButtons("futSimType", "Simulation Type", choices = c("CRT", "MST", "SRT")),
              fileInput("dataset", "Choose CSV File", accept = ".csv"),
              uiOutput("post_var_select"),
              uiOutput("intervention_var_select"),
              conditionalPanel(
                condition = "input.futSimType != 'SRT'",
                uiOutput("random_var_select")
              ),
              uiOutput("covariate_select_fut"),  # Updated to allow multiple selections
              numericInput("crtNsim", "Number of Simulations (>= 10000 is recommended)", 10000, min = 1),
              numericInput("crtThreshold", "Threshold", 0.05, min = 0, max = 1),
              numericInput("crtProbThreshold", "Futility Threshold", 0.8, min = 0, max = 1),
              actionButton("analyzeFutility", "Analyze Futility")
            ),
            mainPanel(
              tableOutput("futilityTable"),
              
              # Place the download button here
              downloadButton("downloadFutilityData", "Download Futility Results")
            )
          )
        ),
        tabPanel(
          "Superiority Analysis",
          sidebarLayout(
            sidebarPanel(
              radioButtons("supSimType", "Simulation Type", choices = c("CRT", "MST", "SRT")),
              fileInput("supDataset", "Choose CSV File", accept = ".csv"),
              
              uiOutput("post_var_select_sup"),      # Post-test variable (single select)
              uiOutput("intervention_var_select_sup"),  # Intervention variable (single select)
              
              conditionalPanel(
                condition = "input.supSimType != 'SRT'",
                uiOutput("random_var_select_sup")
              ),
              
              numericInput("crtSupReference", "Reference Intervention", 1, min = 1),
              numericInput("crtSupNsim", "Number of Simulations (>= 10000 is recommended)", 10000, min = 1),
              numericInput("crtSupThreshold", "Threshold", 0.05, min = 0, max = 1),
              
              # New Superiority Threshold input
              numericInput("crtSupSuperiorThreshold", "Superiority Threshold", 0.5, min = 0, max = 1),
              
              uiOutput("covariate_select_sup"),  # Covariates (multi-select)
              
              actionButton("analyzeSuperiority", "Analyze Superiority")
            ),
            mainPanel(
              tableOutput("superiorityTable"),
              
              # Place the download button here
              downloadButton("downloadSuperiorityData", "Download Superiority Results")
            )
          )
        ),
        tabPanel(
          "Add New Intervention",
          sidebarLayout(
            sidebarPanel(
              # Choose CRT, MST, or SRT option
              radioButtons("newSimType", "Simulation Type", choices = c("CRT", "MST", "SRT")),
              
              # File input to upload the existing dataset
              fileInput("newTreatmentDataset", "Choose Existing Dataset (CSV)", accept = ".csv"),
              
              # Show only for CRT and MST (Not for SRT)
              conditionalPanel(
                condition = "input.newSimType != 'SRT'",
                numericInput("newSchools", "Number of New Schools", 5, min = 1),
                numericInput("newPupilsPerSchool", "Pupils per New School", 100, min = 10)
              ),
              
              # Common Inputs for all types (CRT, MST, SRT)
              numericInput("newEffectSize", "Effect Size for New Intervention", 0.3),
              numericInput("newAttrition", "Attrition Rate for New Intervention", 0.1, min = 0, max = 1),
              
              # MST-specific input
              conditionalPanel(
                condition = "input.newSimType == 'MST'",
                sliderInput("treatmentPercentageMST", "Percentage of Pupils in New Intervention (MST)", 
                            min = 0, max = 1, value = 0.5, step = 0.1)
              ),
              
              # User selections for columns
              uiOutput("post_var_select_new"),
              uiOutput("intervention_var_select_new"),
              uiOutput("schools_var_select_new"),
              uiOutput("pupils_var_select_new"),
              
              # Covariate selection
              uiOutput("covariate_select_new"),
              
              # Button to initiate the addition of a new treatment
              actionButton("addNewTreatment", "Add New Intervention")
            ),
            
            mainPanel(
              # Display the modified dataset
              tableOutput("newTreatmentTable"),
              downloadButton("downloadNewTreatmentData", "Download New Dataset")
            )
          )
        ),
        
        tabPanel(
          "Plot Posterior Probabilities",
          sidebarLayout(
            sidebarPanel(
              radioButtons("plotSimType", "Simulation Type", choices = c("CRT", "MST", "SRT")),
              fileInput("plotDataset", "Choose CSV File", accept = ".csv"),
              
              uiOutput("post_var_select_plot"),      # Post-test variable
              uiOutput("intervention_var_select_plot"), # Intervention variable
              
              # Random variable selection (for CRT and MST, not for SRT)
              conditionalPanel(
                condition = "input.plotSimType != 'SRT'",
                uiOutput("random_var_select_plot")  # Random variable
              ),
              
              # Covariates selection (if any)
              uiOutput("covariate_select_plot"),
              
              numericInput("crtPlotNsim", "Number of Simulations (>= 10000 is recommended)", 10000, min = 1),
              
              # Checkbox to add vertical line
              checkboxInput("addVerticalLine", "Add a vertical line", value = FALSE),
              
              # Conditional numeric input for vertical line value
              conditionalPanel(
                condition = "input.addVerticalLine == true",
                numericInput("crtVerticalLine", "Value for Vertical Line", value = 0.05, min = 0, max = 1)
              ),
              
              # Checkbox to add horizontal line
              checkboxInput("addHorizontalLine", "Add a horizontal line", value = FALSE),
              
              # Conditional numeric input for horizontal line value
              conditionalPanel(
                condition = "input.addHorizontalLine == true",
                numericInput("crtPlotProbThreshold", "Value for Horizontal Line", value = 0.8, min = 0, max = 1)
              ),
              
              # Plot threshold range
              sliderInput("plot_threshold_range", "Plot Threshold Range", min = 0, max = 1, value = c(0, 1), step = 0.1),
              
              actionButton("plotPosterior", "Plot Posterior Probabilities")
            ),
            mainPanel(
              plotOutput("posteriorPlot"),
              downloadButton("downloadPosteriorPlot", "Download Plot")
            )
          )
        ),
        
        # Adding the "User Manual" tab
        tabPanel(
          "User Manual",
          fluidPage(
            tags$head(
              tags$style(HTML("
          /* Increase font size for the entire User Manual tab */
          #user-manual-content {
            font-size: 16px;  /* Increase the font size */
            line-height: 1.8;  /* Adjust line height for readability */
          }
          h1 {
            font-size: 32px;  /* Increase font size for main headers */
          }
          h2 {
            font-size: 26px;  /* Increase font size for section headers */
          }
          ul {
            margin-bottom: 17px;  /* Add some spacing between lists */
          }
        "))
            ),
            div(
              id = "user-manual-content",  # Apply custom styles to this div
              HTML("
          <h1>Platform Trials Simulator (PTS) Manual</h1>
          <h2>Summary</h2>
          <p>The PTS simulator is a web application developed using RShiny, a package within the R and RStudio statistical software environment. The software application includes functions for cluster-randomised, multisite, and simple randomised platform trial simulation, multilevel analysis, as well as futility and superiority analysis outputs. The software's rules are based on calculating Bayesian posterior probabilities of superiority and futility. To run the simulator, you need to input some data into the input bar manually. The software enables you to save and load simulation outputs. For all tabs, the input bar is located on the left side of the browser window, and the outputs are displayed on the right side. PTS consists of seven tabs: data simulation, multilevel analysis, futility analysis, superiority analysis, add new intervention, plot posterior probabilities and User Manual.</p>

          <h2>Data Simulation</h2>
          <p>The Data Simulation tab allows the user to simulate the CRT, MST and SRT.The simulation requires the user to input the following parameters: </p>
          <ul>
            <li><strong>Simulation Type (CRT, MST or SRT):</strong> Select the type of simulation: CRT (Cluster Randomized Trials), MST (Multi-site Trials) or SRT (Simple Randomized Trials).</li>
            <li><strong>Number of Interventions:</strong> The number of intervention arms in the trial, excluding the control arm. <br/>Example: 2 (1 control group and 2 intervention arms)</li>
            <li><strong>Number of Pupils Per School (for CRT/MST):</strong> The total number of pupils (or individuals) assigned to each school (or cluster).</li>
            <li><strong>Number of Schools (for CRT/MST):</strong> Total number of schools (or clusters) participating in the trial.</li>
            <li><strong>Number of Schools Treated (for CRT):</strong> A comma-separated list indicating the number of schools assigned to each group (control, Intervention 1, Intervention 2, etc.). This input applies only to CRT simulation.<br/>Example: '5, 3, 2' (5 control schools, 3 in Intervention 1, 2 in Intervention 2).</li>
            <li><strong>Percentage of Pupils in Each Group (for MST):</strong> A comma-separated list specifying the percentage of pupils in each group (control, Intervention 1, Intervention 2, etc.). This input applies only to MST simulation.<br/>Example: '50, 30, 20' (50% control, 30% in Intervention 1, 20% in Intervention 2).</li>
            <li><strong>Residual Standard Deviation:</strong> The standard deviation of residual errors.</li>
            <li><strong>Intraclass Correlation Coefficient (for CRT):</strong> The ICC measures the proportion of total variance attributed to differences between clusters (e.g., schools).</li>
            <li><strong>Random Intercept Standard Deviations (for MST):</strong> In MST simulations, it specifies the standard deviations for the random intercept, quantifying baseline heterogeneity between schools.</li>
            <li><strong>Random Slope Standard Deviations (for MST):</strong> In MST simulations, users must specify the standard deviations for the random slope, quantifying differential effects of the intervention across schools through school-by-intervention interactions.</li>
            <li><strong>Standard Deviation of Pre-test Scores:</strong> The standard deviation of pre-test scores for individuals.</li>
            <li><strong>Intercept (B0):</strong> The intercept (B0) represents the regression coefficient for intercept.</li>
            <li><strong>Pre-test Coefficient (B1):</strong> The coefficient (B1) represents the regression coefficient for the effect of pre-test.</li>
            <li><strong>Effect Sizes (es):</strong> A comma-separated list of effect sizes for each Intervention group. Example: '0.2, 0.3' (0.2 effect size for Intervention 1, 0.3 for Intervention 2).</li>
            <li><strong>Random Seed:</strong> A seed value used to ensure reproducibility of the simulation.</li>
            <li><strong>Attrition Rates:</strong> A comma-separated list specifying the attrition rate for the control group and each Intervention group. Attrition refers to participants dropping out of the trial.<br/>Example: '0.1, 0.2, 0.3' (10% attrition in control, 20% in Intervention 1, and 30% in Intervention 2).</li>
          </ul>
          
          <h2>Multilevel Analysis Tab</h2>
          <p>The Multilevel Analysis tab allows the user to analyze the CRT, MST, and SRT datasets using Frequentist and Bayesian multilevel models. To begin the analysis, users must upload a CSV file containing their trial data. The dataset should include post-test outcomes, Intervention group assignments, and pre-test scores or any covariates you want to include in the analysis.</p>
          <ul>
            <li><strong>Analysis Method:</strong> Select the type of analysis to be performed. Available options include:</li>
            <ul>
              <li>crtBayes: Bayesian analysis of cluster-randomized education trials using vague priors.</li>
              <li>crtFREQ: Frequentist analysis of cluster-randomized education trials using multilevel models.</li>
              <li>mstBayes: Bayesian analysis of multisite randomized education trials using vague priors.</li>
              <li>mstFREQ: Frequentist analysis of multisite randomized education trials using multilevel models.</li>
              <li>srtBayes: Bayesian analysis of simple randomized education trials using vague priors.</li>
              <li>srtFREQ: Frequentist analysis of simple randomized education trials using multilevel models.</li>
            </ul>
            <li><strong>Select Post-test Outcome:</strong> Select the variable from your dataset that represents the outcome or dependent variable.</li>
            <li><strong>Select Intervention Variables:</strong> Select the variable from your dataset that represents the intervention or Intervention groups.</li>
            <li><strong>Select Clustering Variable (for CRT/MST):</strong> Select the variable from your dataset that represents the clustering structure (e.g., schools, hospitals). This applies only for multilevel models (CRT and MST).</li>
            <li><strong>Number of Simulations (Bayesian Methods):</strong> Set the number of MCMC iterations per chain (>= 10000 is recommended). A minimum of 10,000 is recommended
to ensure convergence.</li>
            <li><strong>Threshold:</strong> A scalar pre-specified threshold for estimating Bayesian posterior probability that the observed effect size is greater than or equal to the threshold.</li>
            <li><strong>Frequentist Method Options (crtFREQ, mstFREQ, srtFREQ):</strong></li>
            <li><strong>Select Confidence Interval Calculation Method:</strong></li>
            <ul> 
              <li>Analytic (Default): Perform standard frequentist analysis.</li>
              <li>Permutation: Specify the number of permutations required to generate a permuted p-value.</li>
              <li>Bootstrap: Perform bootstrap analysis by specifying the number of bootstraps required to generate bootstrap confidence intervals and the bootstrap type.</li>
            </ul>
            <li><strong>Select Additional Covariates:</strong> Select one or more covariates from the dataset that should be included in the analysis model.</li>
            <li><strong>Run Analysis:</strong> After specifying all necessary parameters, click the 'Run Analysis' button to perform the chosen analysis method.</li>
            <li><strong>Results and Plot:</strong> Once the analysis is complete, a forest plot will be created to display the effect sizes using within and total variances. The plot can be downloaded as a TIFF image.</li>
          </ul>
          
          <h2>Futility Analysis Tab</h2>
          <p>The Futility Analysis tab allows users to assess futility in platform trials. Users must upload a CSV dataset and configure the simulation settings.</p>
          <ul>
            <li><strong>Simulation Type (CRT, MST, SRT):</strong> Select the type of trial simulation being analyzed for futility: CRT, MST, or SRT.</li>
            <li><strong>Threshold:</strong> The threshold for estimating the Bayesian Posterior Probability.</li>
            <li><strong>Futility Threshold:</strong> The threshold for the Bayesian Posterior Probability below which an intervention is considered futile.</li>
                        <li><strong>Select Additional Covariates:</strong> Select one or more covariates from the dataset that should be included in the analysis model.</li>
            <li><strong>Analyze Futility:</strong> Click the 'Analyze Futility' button to run the futility analysis based on the selected parameters.</li>
            <li><strong>Futility Decision:</strong> The futility analysis results will present the probability that an intervention's effect size exceeds a specified threshold, denoted as P(Effect Size > Threshold). Additionally, the analysis will determine and indicate whether each intervention is considered futile based on futility threshold.</li>
          </ul>

          <h2>Superiority Analysis Tab</h2>
          <p>The Superiority analysis tab is designed to compare the efficacy of a an Intervention against a reference Intervention. Users must upload a CSV dataset.</p>
          <ul>
           <li><strong>Threshold:</strong> The threshold for estimating the Bayesian Posterior Probability.</li>
            <li><strong>Superiority Threshold:</strong> The threshold for the Bayesian Posterior Probability below which an intervention is considered futile.</li>
            
            <li><strong>Reference Intervention:</strong> Select the reference Intervention against which all other Interventions will be compared.</li>
                        <li><strong>Select Additional Covariates:</strong> Select one or more covariates from the dataset that should be included in the analysis model.</li>
            <li><strong>Analyze Superiority:</strong> Click the 'Analyze Superiority' button to run the superiority analysis.</li>
            <li><strong>Superiority Decision:</strong> The superiority analysis results will present the probability that an intervention's effect size exceeds that of the reference intervention, denoted as P(Effect Size > Threshold). Additionally, the analysis will determine and indicate whether each intervention is considered superior based on superiority threshold.</li>
          </ul>
          
          <h2>Add New Intervention</h2>
        <p>The Add New Intervention tab allows users to introduce a new intervention into an existing dataset. Users must upload a CSV dataset and select the type of data (CRT, MST and SRT):</p>

        <li><strong>select required variables:</strong></li>
          <ul>
              <li><strong>Select Post-test Outcome:</strong> The outcome variable for analysis.</li>
              <li><strong>Select Intervention Variable:</strong> The column containing treatment assignments.</li>
              <li><strong>Select Pupils Variable (for CRT/MST):</strong> The pupil level identifier.</li>
             <li><strong>Select ID Variable (for SRT):</strong> The individual-level identifier.</li>

              <li><strong>Select Schools Variable (for CRT/MST):</strong> The cluster-level identifier..</li>
              <li><strong>Select Additional Covariates:</strong> Additional predictor variables (optional).</li>
            </ul>
          </li>
          <li><strong>Specify New Intervention Parameters:</strong></li>
          <ul>
            <li><strong>Number of New Schools (for CRT/MST):</strong> Specify how many new schools will be introduced.</li>
            <li><strong>Number of Pupils Per School (for CRT/MST):</strong> The number of pupils per new school.</li>
             <li><strong>Number of Participants (for SRT):</strong> The number of participants in new intervention.</li>

            <li><strong>Effect Size for New Intervention:</strong> The effect size of the new intervention.</li>
            <li><strong>Attrition Rate for New Intervention:</strong> The proportion of participants expected to drop out of the new intervention. Default: 0.1 (10% attrition in outcome)</li>
            <li><strong>Percentage of Pupils in New Treatment (for MST):</strong> The proportion of students in new intervention groups.</li>
          </ul>

          
          
          
          <h2>Plot Posterior Probabilities Tab</h2>
          <p>This tab allows users to visualize posterior probabilities across different thresholds for multiple interventions.</p>
          <ul>
            <li><strong>Add Vertical Line:</strong> Option to add a vertical reference line to the plot, typically representing a pre-specified threshold for estimating Bayesian posterior probability.</li>
            <li><strong>Add Horizontal Line:</strong> Option to add a horizontal reference line to the plot, typically representing a threshold of Bayesian posterior probability.</li>
            <li><strong>Plot Threshold Range:</strong> Adjust the range of threshold values for the posterior probability plot.</li>
            <li><strong>Plot Posterior Probabilities:</strong> Click the 'Plot Posterior Probabilities' button to generate the plot. The generated plot will display posterior probabilities across thresholds for each intervention group. You can also download the plot in TIFF format.</li>
          </ul>
          
          <h2>Associated Paper</h2>
          <p>The associated paper of this R Shiny application can be accessed via: <a href='#'>Link to Paper</a>.</p>
     ")
            )
          )
        )
      )
    )
  )
)