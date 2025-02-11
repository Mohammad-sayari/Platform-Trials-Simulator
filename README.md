# Platform-Trials-Simulator
Platform Trials Simulator (PTS): Software for Planning and Simulating Cluster-Randomised, Multisite and Simple Randomised Platform Trials

The PTS simulator is a web application developed using RShiny, a package within the R and RStudio statistical software environment. The software application includes functions for cluster-randomised, multisite, and simple randomised platform trial simulation, multilevel analysis, as well as futility and superiority analysis outputs. The software's rules are based on calculating Bayesian posterior probabilities of superiority and futility. To run the simulator, you need to input some data into the input bar manually. The software enables you to save and load simulation outputs. For all tabs, the input bar is located on the left side of the browser window, and the outputs are displayed on the right side. PTS consists of seven tabs: data simulation, multilevel analysis, futility analysis, superiority analysis, add new intervention, plot posterior probabilities, and User Manual.

Data Simulation

The Data Simulation tab allows the user to simulate the CRT, MST, and SRT. Multilevel Analysis Tab

Futility Analysis Tab

The Futility Analysis tab allows users to assess futility in platform trials. Users must upload a CSV dataset and configure the simulation settings. The analysis will determine whether an intervention's effect size exceeds a specified threshold and indicate whether each intervention is considered futile.

Superiority Analysis Tab

The Superiority analysis tab is designed to compare the efficacy of an intervention against a reference intervention. Users must upload a CSV dataset. The results will present the probability that an intervention's effect size exceeds that of the reference intervention.

Add New Intervention

The Add New Intervention tab allows users to introduce a new intervention into an existing dataset. Users must upload a CSV dataset and select the type of data (CRT, MST, and SRT).

Plot Posterior Probabilities Tab

This tab allows users to visualize posterior probabilities across different thresholds for multiple interventions. Users can adjust the range of threshold values and add reference lines to the plot.

Associated Paper
The associated paper of this R Shiny application can be accessed via: Link to Paper.
