# Platform Trials Simulator (PTS)

**Software for Planning and Simulating Cluster-Randomised, Multisite, and Simple Randomised Platform Trials**

The PTS simulator is a web application developed using **RShiny**, a package within the **R** and RStudio statistical software environment. The software application includes functions for cluster-randomised, multisite, and simple randomised platform trial simulation, multilevel analysis, as well as futility and superiority analysis outputs. The software's rules are based on calculating Bayesian posterior probabilities of superiority and futility. To run the simulator, you need to input some data into the input bar manually. The software enables you to save and load simulation outputs. For all tabs, the input bar is located on the left side of the browser window, and the outputs are displayed on the right side. PTS consists of seven tabs: data simulation, multilevel analysis, futility analysis, superiority analysis, add new intervention, plot posterior probabilities, and User Manual.

### 1. **Data Simulation**
- The Data Simulation tab allows the user to simulate the CRT, MST, and SRT. Multilevel Analysis Tab
### 2. **Multilevel Analysis**
- The Multilevel Analysis tab allows the user to analyze the CRT, MST, and SRT datasets using Frequentist and Bayesian multilevel models. To begin the analysis, users must upload a CSV file containing their trial data. The dataset should include post-test outcomes, Intervention group assignments, and pre-test scores or any covariates you want to include in the analysis.

### 3. **Futility Analysis**
- The Futility Analysis tab allows users to assess futility in platform trials. Users must upload a CSV dataset and configure the simulation settings.

### 4. **Superiority Analysis**
- The Superiority analysis tab is designed to compare the efficacy of a an Intervention against a reference Intervention. Users must upload a CSV dataset.

### 5. **Add New Intervention**
- The Add New Intervention tab allows users to introduce a new intervention into an existing dataset. Users must upload a CSV dataset and select the type of data (CRT, MST and SRT)
### 6. **Plot Posterior Probabilities**
- This tab allows users to visualize posterior probabilities across different thresholds for multiple interventions.

### 7. **User Manual**
- The manual provides detailed descriptions of each input option and includes overviews of the outputs.

### Link to running app
https://pts-app.shinyapps.io/Platform-Trials-Simulator/

### Associated Paper
The associated paper of this **R Shiny** application can be accessed via: Link to Paper.
