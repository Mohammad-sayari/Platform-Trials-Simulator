library(shiny)
library(shinybusy)
library(ggplot2)
library(ggpubr)  # For theme_pubclean()
library(ggtext)
library(eefAnalytics)
library(lme4)


ComparePlot5 <- function(eefAnalyticsList,group, Conditional=TRUE,ES_Total=TRUE,modelNames){
  if(!is(eefAnalyticsList,"list")){stop("eefAnalyticsList is not a list.")}
  
  if(!all(unlist(lapply(eefAnalyticsList,function(x) is(x,"eefAnalytics"))))){stop("Not all list objects are a eefAnalytics class object.")}
  
  if(missing(modelNames)){stop("modelNames must be specified.")}
  if(missing(group)){stop("group must be specified.")}
  
  plotObject(analyticObject=eefAnalyticsList,group=group, Conditional=Conditional,ES_Total=ES_Total, compare=TRUE,modelNames=modelNames)
  
}




##################################
#      Internal plot function    #
##################################


plotObject <- function(analyticObject,group, Conditional,ES_Total,slope, compare,modelNames,...){
  
  if(Conditional ==TRUE){analyticObject2=analyticObject; Condname="Conditional"}
  if(Conditional ==FALSE){analyticObject2=analyticObject$Unconditional; Condname="Unconditional"}
  if(ES_Total ==TRUE){ES_TW<-"Total"}
  if(ES_Total ==FALSE){ES_TW<-"Within"}
  
  if(compare==TRUE & !is.null(names(analyticObject))){stop("Specify the list of objects to compare")}
  if(!is.null(group)){
    trtname <-rownames(analyticObject2$ES)
    if(is.null(trtname)){trtname <-names(analyticObject2$ES)}
    trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
    trt <- trtname[trtpos]
  }
  #bootstrap, permutation and Pprobability plot for SRT model
  #---------------------------------------------------------
  if(sum(analyticObject$Method=="LM") ==1) {
    if(is.null(group)){stop("Group must be specified.")}
    if(!is.null(group)& sum(names(analyticObject)=="ProbES")== 0){
      if(sum(names(analyticObject)=="Bootstrap"|
             names(analyticObject)=="permES")==0){stop("Only relevant for bootstrapped or permutated values")}
      if(sum(names(analyticObject)=="Bootstrap")==1){
        ntp <- nrow(as.matrix(analyticObject$ES))
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
        
        obs.est <- analyticObject2$ES[trt,1]
        tmp2 <- as.numeric(analyticObject2$Bootstrap[,trt])
        xlabs=paste0(Condname," Bootstrap estimates")
        hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main="")
        abline(v=obs.est,col="red",lwd=2,lty=1)
        abline(v=0,col="grey48",lwd=2,lty=1)
        legend("topright",c("Observed Estimate","Zero-Line"),col=c("red","grey48"),bty="n",lty=1,lwd=2)
      }
      
      if(sum(names(analyticObject2)=="permES")==1){
        ntp <- nrow(as.matrix(analyticObject2$ES))
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
        
        
        Perm.names <-names(analyticObject2$permES)
        
        obs.est <- analyticObject2$ES[trt,1]
        tmp2 <- as.numeric(analyticObject2$permES[,trt])
        pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
        xlabs=paste0("Permutation values (PermES) based on ",Condname, " ES")
        
        hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main=paste("P(|PermES| > |ES|)=",pvalue,sep=""))
        abline(v=obs.est,col="red",lwd=2,lty=2)
        abline(v=-obs.est,col="red",lwd=2,lty=2)
        legend("topright",c("(-) Observed Estimate"),col=c("red"),bty="n",lty=2,lwd=2)
      }
    }
    if(sum(names(analyticObject)=="ProbES")> 0 ){
      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      tmp2 <- analyticObject2$ProbES[[which(trtpos==TRUE)]]
      
      thd0<- regmatches(rownames(tmp2),  gregexpr("[[:digit:]]+\\.*[[:digit:]]*",rownames(tmp2)))
      thd <- as.numeric(unlist(thd0))
      
      par_original <- par()[c("mar","xpd")]
      par_original0<- par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(thd, tmp2[,1], col="blue", type="b",ylim=c(0,max(tmp2)), xlab="Threshold", ylab="Posterior probability")
      on.exit(par(par_original0))
      on.exit(par(par_original))
    }
    
  }
  
  
  #bootstrap, permutation, and Pprobability plot for CRT and MST model
  #--------------------------------------------------------------------
  if(sum(analyticObject$Method=="MLM")==1){
    
    if(is.null(group)){
      
      tmp000  <- data.frame(analyticObject$SchEffects)#use analyticObject since both (un)condition has the same SchEffects object.
      if(slope==FALSE){
        tmp00 <- tmp000[,grep("Schools|Intercept|Estimate",names(tmp000))]
        mar11 <-  c(5, 4, 4, 2) + 0.1
      }
      if(slope==TRUE & dim(tmp000)[2] ==2){stop("x must be mstFREQ or mstBAyes object")}
      if(slope==TRUE & dim(tmp000)[2] >2){
        tmp00 <- tmp000[,!(names(tmp000) %in% "Intercept")]
        if(dim(tmp00)[2]==2){mar11 <- c(5, 4, 4, 2) + 0.1}
        if(dim(tmp00)[2]==3){mar11 <- c(5, 2, 4, 0) + 1.0}
        if(dim(tmp00)[2] >3){mar11 <- c(3, 2, 0, 0) + 1.0}}
      
      op <- par(mfrow = c(floor(dim(tmp00)[2]/2),round(dim(tmp00)[2]/2)),
                mar = mar11)
      
      
      for(i in 2:dim(tmp00)[2]){
        tmp <- data.frame(y=tmp00[,i],x=c(1:length(tmp00[,i])))
        tmp2 <- tmp[order(tmp$y),]
        ylabs=gsub("trt", "Intervention ",gsub("Estimate","Intercept", names(tmp00)[i]))
        barplot(tmp2$y,names.arg=tmp2$x,las=2,col="cornflowerblue",border="cornflowerblue")
        if(dim(tmp00)[2]<=2){mtext(ylabs, side = 2.5, line = 2)}
        if(dim(tmp00)[2] >2){mtext(ylabs, side = 2, line = 1.7, cex = 0.8)}
      }
      lines1=-2.5
      if(dim(tmp00)[2] >3){lines1=-1}
      title(xlab="School labels", outer = TRUE, line = lines1,cex.lab = 1.2)
      on.exit(par(op))
    }
    
    
    
    
    
    if( !is.null(group) & sum(names(analyticObject)=="Bootstrap")>0){
      ntp <- length(analyticObject2$ES)
      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      Boot.names <-names(analyticObject2$Bootstrap)
      obs.est <- analyticObject2$ES[[trt]][ES_TW,1]
      tmp2 <- as.numeric(analyticObject2$Bootstrap[,grep(ES_TW,grep(trt,Boot.names, ignore.case = T, value = T))])
      pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
      xlabs=paste0("Bootstrap estimates for ",Condname, " ES_",ES_TW)
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main="")
      abline(v=obs.est,col="red",lwd=2,lty=1)
      abline(v=0,col="grey48",lwd=2,lty=1)
      legend("topright",c("Observed Estimate","Zero-Line"),col=c("red","grey48"),bty="n",lty=1,lwd=2)
    }
    
    
    if( !is.null(group) & sum(names(analyticObject)=="permES")>0){
      ntp <- ifelse(is.list(analyticObject$ES),length(analyticObject$ES),1)
      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      Perm.names <-names(analyticObject2$permES)
      obs.est <- analyticObject2$ES[[trt]][ES_TW,1]
      tmp2 <- as.numeric(analyticObject2$permES[,grep(ES_TW,grep(trt,Perm.names, ignore.case = T, value = T))])
      pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
      xlabs=paste0("Permutation values(PermES) based on ",Condname, " ES_",ES_TW)
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main=paste("P(|PermES| > |ES|)=",pvalue,sep=""))
      abline(v=obs.est,col="red",lwd=2,lty=2)
      abline(v=-obs.est,col="red",lwd=2,lty=2)
      legend("topright",c("Observed Estimate"),col=c("red"),bty="n",lty=2,lwd=2)
    }
    
    
    if( !is.null(group) &sum(names(analyticObject)=="ProbES")> 0 ){
      
      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      tmp2 <- analyticObject2$ProbES[[which(trtpos==TRUE)]]
      tmp2.within<- tmp2[, grep("with",names(tmp2), ignore.case = TRUE)]
      tmp2.total <- tmp2[, grep("total",names(tmp2), ignore.case = TRUE)]
      thd0<- regmatches(rownames(tmp2),  gregexpr("[[:digit:]]+\\.*[[:digit:]]*",rownames(tmp2)))
      thd <- as.numeric(unlist(thd0))
      
      par_original <- par()[c("mar","xpd")]
      op<- par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(thd, tmp2.within, col="blue", type="b",ylim=c(0,max(tmp2)), xlab="Threshold", ylab="Posterior probability",...)
      lines(thd, tmp2.total, col="red", type="b", lty=2)
      legend("topright", legend=c("within", "total"), col=c("blue", "red"), lty=1:2, cex=0.8)
      on.exit(par(op))
      on.exit(par(par_original))
      
    }
    
  }
  
  # error bar for model comparing models
  #-------------------------------------
  if(is.null(names(analyticObject))){
    
    ltp <- names(analyticObject)
    if(!is.null(ltp)){stop("Specify list of eefAnalytics objects for comparison")}
    if(is.null(group)){stop("Group number must be defined")}
    ntp <- length(analyticObject)
    if(length(modelNames)!= ntp){stop("Names must be equal to the number of eefAnalytics objects")}
    
    es.mean <- es.lower <- es.upper <- p.name <- var.name <- NULL
    for(k in 1:ntp){
      tmp <- analyticObject[[k]]
      
      if(tmp$Method=="LM"){
        trtname <-rownames(tmp$ES)
        trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
        trt <- trtname[trtpos]
        
        if(Conditional==TRUE){tmp2 <- as.matrix(tmp$ES)}
        if(Conditional==FALSE){tmp2 <- as.matrix(tmp$Unconditional$ES)}
        trtname <-rownames(tmp$ES)
        if(is.null(trtname)){trtname <-names(tmp$ES)}
        trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
        trt <- trtname[trtpos]
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
        es.mean1 <-tmp2[trt,1]
        es.lower1 <-tmp2[trt,2]
        es.upper1 <-tmp2[trt,3]
        p.name1 <- rep(modelNames[k],length(es.mean1))
        var.name1 <- rep("Within",length(es.mean1))
        
      }
      
      if(tmp$Method=="MLM"){
        trtname <-names(tmp$ES)
        if(is.null(trtname)){trtname <-names(tmp$ES)}
        trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
        trt <- trtname[trtpos]
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
        if(Conditional==TRUE) {tmp2 <- tmp$ES[[trt]]}
        if(Conditional==FALSE){tmp2 <- tmp$Unconditional$ES[[trt]]}
        es.mean1 <- tmp2[,1]
        es.lower1 <-tmp2[,2]
        es.upper1 <-tmp2[,3]
        p.name1 <- rep(modelNames[k],length(es.mean1))
        var.name1 <- rownames(tmp2)
        
        
      }
      
      
      es.mean <- c(es.mean,es.mean1)
      es.lower <- c(es.lower,es.lower1)
      es.upper <- c(es.upper,es.upper1)
      p.name <- c(p.name,p.name1)
      var.name <- c(var.name,var.name1)
      
    }
    
    
    MyData1 <- data.frame(ES=es.mean,LB.95=es.lower,UB.95=es.upper,Variance=var.name,Name=p.name)
    MyData1<- MyData1[order(MyData1$Variance,decreasing = T),]
    MyData1$Anot <- paste0(MyData1$ES, " [", MyData1$LB.95,", ", MyData1$UB.95,"]")
    MyData1$Index <- 1:dim(MyData1)[1]
    MyData1$Xaxis <- (max(MyData1$UB.95))+0.05
    Mybreaks <-  round(c(min(MyData1$LB.95),(min(MyData1$LB.95)+(max(MyData1$UB.95)))/2,max(MyData1$UB.95)),2)
    xlimits <- c(min(min(MyData1$LB.95),0),(max(MyData1$UB.95))+((min(MyData1$LB.95)+(max(MyData1$UB.95)))/4)+0.6)
    
    
    Ann_text <- data.frame(Index = length(MyData1$Variance[MyData1$Variance=="Within"])+0.5,
                           ES = MyData1$Xaxis[1],LB.95=0, UB.95=0,lab = "Text",
                           Variance = factor("Total",levels = c("Within", "Total")))
    
    #ggplot
    p <- ggplot(data=MyData1, aes(x=ES, y=Name, xmin=LB.95, xmax=UB.95))
    p <- p + geom_point()
    p <- p + geom_errorbarh(height=.1)
    p <- p + scale_x_continuous(limits=xlimits ,breaks=Mybreaks, name=expression(paste("Hedges' ", italic("g"))))
    p <- p + geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)
    if(sum(unique(MyData1$Variance)%in% "Total")>0){p <- p + facet_grid(Variance~., scales= "free", space="free")}
    p <- p + ylab("Models")
    p <- p + theme_bw()
    p <- p + theme(axis.text.y =element_text(color="black"))
    p <- p + theme(text=element_text(size=16, color="black"))
    p <- p + theme(panel.spacing = unit(1, "lines"))
    p <- p + geom_text(aes(x = Xaxis,y = Name, label = Anot),hjust = 0, nudge_x = 0.04)
    p <- p + geom_text(data = Ann_text, y=Inf,label = "95% CI",hjust = -1.1,vjust = -0.5,size=4,fontface = "bold")
    p <- p + coord_cartesian(clip = "off")
    p <- p + theme(plot.margin = unit(c(30,5,5,5), "point"))
    p
    
  }
  
}


# Define the CRT and MST data simulation functions
crtdata_simulation <- function(nt, n_schools_treated, np, ns, sigma, ICC, sigmaPret, B0, B1, es, seed, attrition_rates) {
  
  # Error checking
  if (length(n_schools_treated) != nt + 1) {
    stop("Error: The length of n_schools_treated must be number of interventions + 1 (including control group).")
  }
  if (length(es) != nt) {
    stop("Error: The length of es must be equal to number of interventions.")
  }
  if (length(attrition_rates) != nt + 1) {
    stop("Error: The length of attrition_rates must be number of interventions + 1 (including control group).")
  }
  set.seed(seed)
  
  # Step 1: Create the base data structure with people and schools
  data <- expand.grid(pupils = 1:np, schools = 1:ns)
  data$pupils <- 1:nrow(data)
  
  # Step 2: Initialize treatment assignment and assign treatment groups to schools
  interventions <- "interventions"
  unique_schools <- unique(data$schools)
  
  # Initialize treatment column (0 = control)
  data[[interventions]] <- 0
  
  # Check if the total number of schools assigned matches the number of schools (ns)
  if (sum(n_schools_treated) != ns) {
    stop("Error: The total number of schools in n_schools_treated must be equal to the total number of schools (ns).")
  }
  
  # Assign schools to control group and treatment groups
  remaining_schools <- unique_schools
  
  # Step 2a: Assign schools to the control group
  control_schools <- sample(remaining_schools, n_schools_treated[1])
  remaining_schools <- setdiff(remaining_schools, control_schools)  # Remove control schools from remaining
  
  # Step 2b: Assign schools to each treatment group
  for (i in 1:nt) {
    treated_schools <- sample(remaining_schools, n_schools_treated[i + 1])
    data[[interventions]][data$schools %in% treated_schools] <- i  # Assign treatment group i
    remaining_schools <- setdiff(remaining_schools, treated_schools)  # Remove assigned schools
  }
  
  # Step 3: Generate pre-test scores for everyone
  prets <- "pretest"
  data[[prets]] <- rnorm(nrow(data), mean = 0, sd = sigmaPret)  # Pre-test scores for everyone
  
  # Initialize post-test scores
  posts <- "posttest"
  data[[posts]] <- NA  # Initially set all to NA
  
  # Step 4: Handle attrition for control and treatment groups
  non_attrition_idx <- 1:nrow(data)  # Start by assuming no one is attrited
  
  for (i in 0:nt) {
    group_idx <- which(data[[interventions]] == i)
    attrition_rate <- attrition_rates[i + 1]  # Use i+1 because first value is for control group
    attrition_size <- round(length(group_idx) * attrition_rate)
    
    if (attrition_size > 0) {
      attrition_idx <- sample(group_idx, attrition_size)
      data[attrition_idx, posts] <- NA  # Mark attrited participants
      non_attrition_idx <- setdiff(non_attrition_idx, attrition_idx)  # Update non-attrited participants
    }
  }
  
  # Step 5: Generate random individual errors and cluster-level random effects
  e_ij <- rnorm(length(non_attrition_idx), mean = 0, sd = sigma)  # Individual-level errors for non-attrited individuals
  sigmab <- sqrt(ICC / (1 - ICC) * sigma^2)  # School-level variance for all schools
  b_i_full <- rnorm(ns, mean = 0, sd = sigmab)  # Cluster-level random effects for all schools
  b_i <- b_i_full[data$schools[non_attrition_idx]]  # Get the random effects for relevant schools
  
  # Step 6: Create treatment effects for each group
  treatment_effects <- sapply(1:nt, function(i) as.numeric(data[non_attrition_idx, interventions] == i))
  
  treatment_effects_sizes <- sweep(treatment_effects, 2, es* sqrt(sigmab^2 + sigma^2), `*`)
  
  # Calculate the total treatment effect by combining individual treatment effects with corresponding effect sizes
  total_treatment_effect <- rowSums(treatment_effects_sizes)
  
  # Step 7: Compute post-test scores using the combined formula
  data[non_attrition_idx, posts] <- B0 + B1 * data[non_attrition_idx, prets] + total_treatment_effect + b_i + e_ij
  
  # Only keep relevant columns (pupils, schools, interventions, pretest, posttest)
  data <- data[, c("pupils", "schools", interventions, prets, posts)]
  
  return(data)
}

mstdata_simulation <- function(nt, tpi, np, ns, sigma, sigmab0, sigmab1, sigmaPret, B0, B1, es, seed, attrition_rates) {
  # Error checking: ensure tpi has length nt + 1 (control + treatment groups)
  if (length(tpi) != nt + 1) {
    stop("Error: 'tpi' must have length nt + 1 (first value for control group, remaining for treatment groups).")
  }
  
  # Error checking: ensure the sum of tpi is 100
  if (sum(tpi) != 100) {
    stop("Error: The sum of 'tpi' must be 100%.")
  }
  
  # Error checking: ensure es has length equal to nt
  if (length(es) != nt) {
    stop("Error: 'es' must have length nt (one for each treatment group).")
  }
  
  # Error checking: ensure attrition_rates has length equal to nt + 1
  if (length(attrition_rates) != nt + 1) {
    stop("Error: 'attrition_rates' must have length nt + 1 (one for control group and one for each treatment group).")
  }
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Step 1: Create the base data structure with people and schools
  data <- expand.grid(pupils = 1:np, schools = 1:ns)
  data$pupils <- 1:nrow(data)
  
  # Initialize treatment column (0 = control)
  interventions <- "interventions"
  data[[interventions]] <- 0
  
  # Total number of participants
  total_participants <- nrow(data)
  
  # Step 2: Assign control group based on tpi[1] (control percentage)
  control_size <- floor(total_participants * (tpi[1] / 100))
  control_pupils <- sample(1:total_participants, control_size)
  data[[interventions]][control_pupils] <- 0  # Assign control group (0)
  
  # Step 3: Assign treatment groups to the remaining pupils
  remaining_pupils <- setdiff(1:total_participants, control_pupils)
  
  # Calculate the number of participants for each treatment group except the last one
  treatment_sizes <- floor(total_participants * (tpi[2:(nt+1)] / 100))
  
  # Adjust the last treatment group to take the remaining participants
  treatment_sizes[nt] <- length(remaining_pupils) - sum(treatment_sizes[1:(nt-1)])
  
  for (i in 1:nt) {
    treated_pupils <- sample(remaining_pupils, treatment_sizes[i])
    data[[interventions]][treated_pupils] <- i  # Assign treatment group i to selected pupils
    remaining_pupils <- setdiff(remaining_pupils, treated_pupils)  # Remove assigned pupils
  }
  
  # Step 4: Generate pre-test scores for everyone
  prets <- "pret"
  data[[prets]] <- rnorm(nrow(data), mean = 0, sd = sigmaPret)
  
  # Initialize post-test scores
  posts <- "post"
  data[[posts]] <- NA  # Initially set all to NA
  
  # Step 5: Handle attrition for control and treatment groups
  non_attrition_idx <- 1:nrow(data)  # Start by assuming no one is attrited
  
  for (i in 0:nt) {
    group_idx <- which(data[[interventions]] == i)
    attrition_rate <- attrition_rates[i + 1]  # Use i+1 because first value is for control group
    attrition_size <- round(length(group_idx) * attrition_rate)
    
    if (attrition_size > 0) {
      attrition_idx <- sample(group_idx, attrition_size)
      data[attrition_idx, posts] <- NA  # Mark attrited participants
      non_attrition_idx <- setdiff(non_attrition_idx, attrition_idx)  # Update non-attrited participants
    }
  }
  
  # Step 6: Generate random individual errors and cluster-level random effects
  e_ij <- rnorm(length(non_attrition_idx), mean = 0, sd = sigma)  # Individual-level errors for non-attrited individuals
  
  # Cluster-level random effects
  b_i_full <- rnorm(ns, mean = 0, sd = sigmab0)  # Random intercepts for schools
  b1_i_full <- rnorm(ns, mean = 0, sd = sigmab1)  # Random slopes for treatment effect
  
  # Get the random effects for the relevant schools
  b_i <- b_i_full[data$schools[non_attrition_idx]]
  b1_i <- b1_i_full[data$schools[non_attrition_idx]]
  
  # Step 7: Use model.matrix to create treatment effect dummy variables matrix
  unique_interventions <- unique(data[non_attrition_idx, interventions])
  
  # Create the model matrix without intercept
  treatment_matrix <- model.matrix(~ factor(data[non_attrition_idx, interventions]) - 1)
  
  # Remove the first column
  treatment_matrix <- treatment_matrix[, -1]
  
  treatment_effects <-  sweep(treatment_matrix, 2, es* sqrt(sigmab0^2 + sigmab1^2 + sigma^2), `*`)
  
  # Calculate the total treatment effect by summing the contributions of each treatment group
  total_treatment_effect <- rowSums(treatment_effects) 
  
  random_slope <-  treatment_matrix*b1_i
  
  # Calculate the total treatment effect by summing the contributions of each treatment group
  total_random_slope_effect <- rowSums(random_slope) 
  
  # Step 8: Compute post-test scores for non-attrited participants
  data[non_attrition_idx, posts] <- B0 + 
    B1 * data[non_attrition_idx, prets] + 
    total_treatment_effect + 
    b_i + total_random_slope_effect + 
    e_ij
  
  # Only keep relevant columns (pupils, schools, interventions, pret, post)
  data <- data[, c("pupils", "schools", interventions, prets, posts)]
  
  return(data)
}

srtdata_simulation <- function(nt, tpi, np, sigma, sigmaPret, B0, B1, es, seed, attrition_rates) {
  # Error checking: ensure tpi has length nt + 1 (control + treatment groups)
  if (length(tpi) != nt + 1) {
    stop("Error: 'tpi' must have length nt + 1 (first value for control group, remaining for treatment groups).")
  }
  
  # Error checking: ensure the sum of tpi is 100
  if (sum(tpi) != 100) {
    stop("Error: The sum of 'tpi' must be 100%.")
  }
  
  # Error checking: ensure es has length equal to nt
  if (length(es) != nt) {
    stop("Error: 'es' must have length nt (one for each treatment group).")
  }
  
  # Error checking: ensure attrition_rates has length equal to nt + 1
  if (length(attrition_rates) != nt + 1) {
    stop("Error: 'attrition_rates' must have length nt + 1 (one for control group and one for each treatment group).")
  }
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Step 1: Create the base data structure with participants
  data <- data.frame(ID = 1:np)  # Assign unique IDs
  
  # Initialize treatment column (0 = control)
  interventions <- "interventions"
  data[[interventions]] <- 0
  
  # Step 2: Assign control group based on tpi[1] (control percentage)
  control_size <- floor(np * (tpi[1] / 100))
  control_pupils <- sample(1:np, control_size)
  data[[interventions]][control_pupils] <- 0  # Assign control group (0)
  
  # Step 3: Assign treatment groups to the remaining participants
  remaining_pupils <- setdiff(1:np, control_pupils)
  
  # Calculate the number of participants for each treatment group except the last one
  treatment_sizes <- floor(np * (tpi[2:(nt+1)] / 100))
  
  # Adjust the last treatment group to take the remaining participants
  treatment_sizes[nt] <- length(remaining_pupils) - sum(treatment_sizes[1:(nt-1)])
  
  for (i in 1:nt) {
    treated_pupils <- sample(remaining_pupils, treatment_sizes[i])
    data[[interventions]][treated_pupils] <- i  # Assign treatment group i to selected participants
    remaining_pupils <- setdiff(remaining_pupils, treated_pupils)  # Remove assigned participants
  }
  
  # Step 4: Generate pre-test scores for everyone
  prets <- "pret"
  data[[prets]] <- rnorm(nrow(data), mean = 0, sd = sigmaPret)
  
  # Initialize post-test scores
  posts <- "post"
  data[[posts]] <- NA  # Initially set all to NA
  
  # Step 5: Handle attrition for control and treatment groups
  non_attrition_idx <- 1:nrow(data)  # Start by assuming no one is attrited
  
  for (i in 0:nt) {
    group_idx <- which(data[[interventions]] == i)
    attrition_rate <- attrition_rates[i + 1]  # Use i+1 because first value is for control group
    attrition_size <- round(length(group_idx) * attrition_rate)
    
    if (attrition_size > 0) {
      attrition_idx <- sample(group_idx, attrition_size)
      data[attrition_idx, posts] <- NA  # Mark attrited participants
      non_attrition_idx <- setdiff(non_attrition_idx, attrition_idx)  # Update non-attrited participants
    }
  }
  
  # Step 6: Generate post-test scores for non-attrited participants
  # Create treatment matrix
  treatment_matrix <- model.matrix(~ factor(data[non_attrition_idx, interventions]) - 1)
  
  # Remove the first column (control group is baseline)
  if (ncol(treatment_matrix) > 1) {
    treatment_matrix <- treatment_matrix[, -1, drop = FALSE]
  }
  
  # Apply treatment effects using only individual-level variance (sigma)
  treatment_effects <- sweep(treatment_matrix, 2, es * sigma, `*`)
  
  # Calculate the total treatment effect by summing the contributions of each treatment group
  total_treatment_effect <- rowSums(treatment_effects)
  
  # Generate individual-level errors for non-attrited individuals
  e_ij <- rnorm(length(non_attrition_idx), mean = 0, sd = sigma)
  
  # Compute post-test scores
  data[non_attrition_idx, posts] <- B0 + 
    B1 * data[non_attrition_idx, prets] + 
    total_treatment_effect + 
    e_ij
  
  # Only keep relevant columns (ID, interventions, pret, post)
  data <- data[, c("ID", interventions, prets, posts)]
  
  return(data)
}


crtfutility <- function(data, post_vars = "post1", intervention_column = "treatments", Random = "schls", Nsim = 2000, 
                        Threshold = 0.05, ProbThreshold = 0.8, covariates = NULL) {
  
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  prob_es_values <- numeric(length(interventions))
  futility_decisions <- numeric(length(interventions))  # Initialize a single futility column
  
  output <- list()
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_vars, "~", intervention_column)
    } else {
      paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    output[[i]] <- crtBayes(
      as.formula(formula_str),
      random = Random,
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    prob_es_values[i] <- as.numeric(output[[i]]$ProbES[[1]]["Total1"])
    futility_decisions[i] <- if (prob_es_values[i] < ProbThreshold) 1 else 0
  }
  
  return(data.frame(
    Treatment = interventions,
    Futility = futility_decisions,  # Single futility column
    ProbES = prob_es_values
  ))
}


mstfutility <- function(data, post_vars = "post1", intervention_column = "treatments", Random = "schls", Nsim = 2000, 
                        Threshold = 0.05, ProbThreshold = 0.8, covariates = NULL) {
  
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  prob_es_values <- numeric(length(interventions))
  futility_decisions <- numeric(length(interventions))  # Initialize a single futility column
  
  output <- list()
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_vars, "~", intervention_column)
    } else {
      paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    output[[i]] <- mstBayes(
      as.formula(formula_str),
      random = Random,
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    prob_es_values[i] <- as.numeric(output[[i]]$ProbES[[1]]["Total1"])
    futility_decisions[i] <- if (prob_es_values[i] < ProbThreshold) 1 else 0
  }
  
  return(data.frame(
    Treatment = interventions,
    Futility = futility_decisions,  # Single futility column
    ProbES = prob_es_values
  ))
}


srtfutility <- function(data, post_vars = "post1", intervention_column = "treatments", Nsim = 2000, 
                        Threshold = 0.05, ProbThreshold = 0.8, covariates = NULL) {
  
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  prob_es_values <- numeric(length(interventions))
  futility_decisions <- numeric(length(interventions))  # Initialize a single futility column
  
  output <- list()
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_vars, "~", intervention_column)
    } else {
      paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    output[[i]] <- srtBayes(
      as.formula(formula_str),
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    prob_es_values[i] <- as.numeric(output[[i]]$ProbES[[1]]["Cond"])
    futility_decisions[i] <- if (prob_es_values[i] < ProbThreshold) 1 else 0
  }
  
  return(data.frame(
    Treatment = interventions,
    Futility = futility_decisions,  # Single futility column
    ProbES = prob_es_values
  ))
}



crtSuperiority <- function(data, post_var = "post1", intervention_column = "treatments", Random = "schls", 
                           Nsim = 2000, Threshold = 0.1, reference_intervention = 1, superiority_threshold = 0.5, covariates = NULL) {
  
  # Get unique interventions excluding the control group (0)
  interventions <- sort(unique(data[[intervention_column]]))
  interventions <- interventions[interventions != 0 & interventions != reference_intervention]  # Exclude control and reference intervention
  
  output <- list()
  prob_es_values <- numeric(length(interventions))
  sup_decisions <- numeric(length(interventions))
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    
    # Subset data for the current intervention and reference intervention
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == reference_intervention)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    # Build the formula
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_var, "~", intervention_column)
    } else {
      paste(post_var, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    # Run the Bayesian analysis
    output[[as.character(intervention)]] <- crtBayes(
      as.formula(formula_str),
      random = Random,
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    # Extract posterior probability of superiority
    prob_es_values[i] <- as.numeric(output[[as.character(intervention)]]$ProbES[[1]]["Total1"])
    sup_decisions[i] <- if (prob_es_values[i] > superiority_threshold) 1 else 0
  }
  
  # Add the reference treatment row
  reference_row <- data.frame(
    Treatment = reference_intervention,
    ProbES = NA,  # No comparison for reference intervention
    Superiority = "Reference"
  )
  
  # Create the result table
  result <- data.frame(
    Treatment = interventions,
    ProbES = prob_es_values,
    Superiority = ifelse(sup_decisions == 1, "Superior", "Not Superior")
  )
  
  # Combine the reference row with the results
  result <- rbind(reference_row, result)
  
  # Sort the table by Treatment
  result <- result[order(result$Treatment), ]
  
  return(result)
}



mstSuperiority <- function(data, post_var = "post1", intervention_column = "treatments", Random = "schls", 
                           Nsim = 2000, Threshold = 0.1, reference_intervention = 1, superiority_threshold = 0.5, covariates = NULL) {
  
  # Get unique interventions excluding the control group (0)
  interventions <- sort(unique(data[[intervention_column]]))
  interventions <- interventions[interventions != 0 & interventions != reference_intervention]  # Exclude control and reference intervention
  
  output <- list()
  prob_es_values <- numeric(length(interventions))
  sup_decisions <- numeric(length(interventions))
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    
    # Subset data for the current intervention and reference intervention
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == reference_intervention)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    # Build the formula
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_var, "~", intervention_column)
    } else {
      paste(post_var, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    # Run the Bayesian analysis
    output[[as.character(intervention)]] <- mstBayes(
      as.formula(formula_str),
      random = Random,
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    # Extract posterior probability of superiority
    prob_es_values[i] <- as.numeric(output[[as.character(intervention)]]$ProbES[[1]]["Total1"])
    sup_decisions[i] <- if (prob_es_values[i] > superiority_threshold) 1 else 0
  }
  
  # Add the reference treatment row
  reference_row <- data.frame(
    Treatment = reference_intervention,
    ProbES = NA,  # No comparison for reference intervention
    Superiority = "Reference"
  )
  
  # Create the result table
  result <- data.frame(
    Treatment = interventions,
    ProbES = prob_es_values,
    Superiority = ifelse(sup_decisions == 1, "Superior", "Not Superior")
  )
  
  # Combine the reference row with the results
  result <- rbind(reference_row, result)
  
  # Sort the table by Treatment
  result <- result[order(result$Treatment), ]
  
  return(result)
}


srtSuperiority <- function(data, post_var = "post1", intervention_column = "treatments", 
                           Nsim = 2000, Threshold = 0.1, reference_intervention = 1, superiority_threshold = 0.5, covariates = NULL) {
  
  # Get unique interventions excluding the control group (0)
  interventions <- sort(unique(data[[intervention_column]]))
  interventions <- interventions[interventions != 0 & interventions != reference_intervention]  # Exclude control and reference intervention
  
  output <- list()
  prob_es_values <- numeric(length(interventions))
  sup_decisions <- numeric(length(interventions))
  
  for (i in seq_along(interventions)) {
    intervention <- interventions[i]
    
    # Subset data for the current intervention and reference intervention
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == reference_intervention)
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    # Build the formula
    formula_str <- if (is.null(covariates) || length(covariates) == 0) {
      paste(post_var, "~", intervention_column)
    } else {
      paste(post_var, "~", intervention_column, "+", paste(covariates, collapse = " + "))
    }
    
    # Run the Bayesian analysis
    output[[as.character(intervention)]] <- srtBayes(
      as.formula(formula_str),
      intervention = intervention_column,
      nsim = Nsim,
      data = intervention_data,
      threshold = Threshold
    )
    
    # Extract posterior probability of superiority
    prob_es_values[i] <- as.numeric(output[[as.character(intervention)]]$ProbES[[1]]["Cond"])
    sup_decisions[i] <- if (prob_es_values[i] > superiority_threshold) 1 else 0
  }
  
  # Add the reference treatment row
  reference_row <- data.frame(
    Treatment = reference_intervention,
    ProbES = NA,  # No comparison for reference intervention
    Superiority = "Reference"
  )
  
  # Create the result table
  result <- data.frame(
    Treatment = interventions,
    ProbES = prob_es_values,
    Superiority = ifelse(sup_decisions == 1, "Superior", "Not Superior")
  )
  
  # Combine the reference row with the results
  result <- rbind(reference_row, result)
  
  # Sort the table by Treatment
  result <- result[order(result$Treatment), ]
  
  return(result)
}


# Function to add a new treatment to CRT data
add_new_treatmentcrt <- function(existing_data, new_schools, new_pupils_per_school, 
                                 es, attrition_rate, 
                                 post_col, intervention_col, 
                                 schools_col, pupils_col, covariates) {
  
  # Step 1: Fit a Linear Mixed-Effects Model using lmer
  formula <- as.formula(paste(post_col, "~", paste(covariates, collapse = " + "), "+ (1 |", schools_col, ")"))
  lmer_model <- lmer(formula, data = existing_data, REML = TRUE)
  
  # Extract fixed effect estimates
  fixed_effects <- fixef(lmer_model)
  
  # Extract variance components
  sigma <- sigma(lmer_model)  # Residual standard deviation
  school_sd <- as.numeric(VarCorr(lmer_model)[[schools_col]])^0.5  # School-level variance
  
  # Assign the new treatment number as the next available number
  new_treatment_num <- max(existing_data[[intervention_col]], na.rm = TRUE) + 1
  
  # Ensure school IDs continue sequentially
  max_existing_school_id <- max(existing_data[[schools_col]], na.rm = TRUE)
  new_school_ids <- seq(from = max_existing_school_id + 1, length.out = new_schools)
  
  # Generate the new dataset structure
  new_data <- expand.grid(
    schls = new_school_ids,
    ppls = 1:new_pupils_per_school
  )
  
  # Ensure pupil IDs start sequentially after the highest existing ID
  max_existing_pupil_id <- max(existing_data[[pupils_col]], na.rm = TRUE)
  new_data[[pupils_col]] <- seq(from = max_existing_pupil_id + 1, length.out = nrow(new_data))
  
  # Assign schools column properly
  new_data[[schools_col]] <- rep(new_school_ids, each = new_pupils_per_school)
  
  # Assign the new intervention to all new pupils
  new_data[[intervention_col]] <- new_treatment_num
  
  # Convert categorical covariates to factors in both datasets
  for (covariate in covariates) {
    if (is.character(existing_data[[covariate]]) || is.factor(existing_data[[covariate]])) {
      existing_data[[covariate]] <- as.factor(existing_data[[covariate]])
      new_data[[covariate]] <- as.factor(sample(levels(existing_data[[covariate]]), nrow(new_data), replace = TRUE))
    }
  }
  
  # Generate numerical covariates (including pretest)
  for (covariate in covariates) {
    if (is.numeric(existing_data[[covariate]])) {
      covariate_mean <- mean(existing_data[[covariate]], na.rm = TRUE)
      covariate_sd <- sd(existing_data[[covariate]], na.rm = TRUE)
      new_data[[covariate]] <- rnorm(nrow(new_data), mean = covariate_mean, sd = covariate_sd)
    }
  }
  
  # Generate school-level random effects
  new_school_effects <- rnorm(new_schools, mean = 0, sd = school_sd)
  
  # Assign school effects to students
  new_data$school_effects <- rep(new_school_effects, each = new_pupils_per_school)
  
  # Generate individual-level residuals
  new_data$individual_residuals <- rnorm(nrow(new_data), mean = 0, sd = sigma)
  
  # Compute baseline post-test score (excluding new treatment effect)
  new_data[[post_col]] <- fixed_effects[1]  # Intercept
  
  # Convert categorical variables to dummy variables using model.matrix
  cat_vars <- covariates[sapply(existing_data[covariates], is.factor) | sapply(existing_data[covariates], is.character)]
  if (length(cat_vars) > 0) {
    cat_matrix <- model.matrix(~ . - 1, data = new_data[, cat_vars, drop = FALSE])
    col_names <- colnames(cat_matrix)
    matched_cols <- intersect(col_names, names(fixed_effects))
    
    for (col in matched_cols) {
      new_data[[post_col]] <- new_data[[post_col]] + fixed_effects[col] * cat_matrix[, col]
    }
  }
  
  # Add continuous covariate effects dynamically
  num_vars <- setdiff(covariates, cat_vars)
  for (covariate in num_vars) {
    if (covariate %in% names(fixed_effects)) {
      new_data[[post_col]] <- new_data[[post_col]] + fixed_effects[covariate] * new_data[[covariate]]
    }
  }
  
  # Apply the effect size **only to the new treatment group**
  effect_size_adjustment <- (new_data[[intervention_col]] == new_treatment_num) * es * sqrt(sigma^2 + school_sd^2)
  
  # Final post-test scores
  new_data[[post_col]] <- new_data[[post_col]] + effect_size_adjustment + new_data$school_effects + new_data$individual_residuals
  
  # Remove auxiliary columns
  new_data <- new_data[, !names(new_data) %in% c("school_effects", "individual_residuals")]
  
  # Apply attrition randomly by setting some post-test scores to NA
  attrition_size <- round(nrow(new_data) * attrition_rate)
  attrition_idx <- sample(1:nrow(new_data), attrition_size)
  new_data[[post_col]][attrition_idx] <- NA
  
  # Ensure all necessary columns exist in both datasets before combining
  missing_in_new <- setdiff(names(existing_data), names(new_data))
  missing_in_existing <- setdiff(names(new_data), names(existing_data))
  
  # Add missing columns as NA to new data
  for (col in missing_in_new) {
    new_data[[col]] <- NA
  }
  
  # Add missing columns as NA to existing data
  for (col in missing_in_existing) {
    existing_data[[col]] <- NA
  }
  
  # Reorder new_data columns to match existing_data
  new_data <- new_data[names(existing_data)]
  
  # Combine the new treatment group with the existing data
  final_combined_data <- rbind(existing_data, new_data)
  
  # Keep only relevant columns
  final_combined_data <- final_combined_data[, c(pupils_col, schools_col, post_col, intervention_col, covariates)]
  
  return(final_combined_data)
}



# Function to add a new treatment to MST data
add_new_treatmentmst <- function(existing_data, new_schools, new_pupils_per_school, 
                                 es, attrition_rate, treatment_percentage,
                                 post_col, intervention_col, 
                                 schools_col, pupils_col, covariates) {
  
  # Fit hierarchical model with random slopes and intercepts
  formula <- as.formula(paste(post_col, "~", paste(covariates, collapse = " + "), "+ (1 +", intervention_col, "|", schools_col, ")"))
  lmer_model <- lmer(formula, data = existing_data, REML = TRUE)
  
  # Extract fixed effects and variance components
  fixed_effects <- fixef(lmer_model)
  random_effects <- VarCorr(lmer_model)
  
  sigma <- sigma(lmer_model)  # Residual standard deviation
  sigmab0 <- attr(random_effects[[schools_col]], "stddev")[1]  # Random intercept SD
  
  # Extract random slope SD (handle missing cases)
  sigmab1 <- ifelse(length(attr(random_effects[[schools_col]], "stddev")) > 1, 
                    attr(random_effects[[schools_col]], "stddev")[2], 
                    0)  # If missing, set to 0
  
  # Assign new treatment number
  new_treatment_num <- ifelse(is.na(max(existing_data[[intervention_col]], na.rm = TRUE)), 
                              1, 
                              max(existing_data[[intervention_col]], na.rm = TRUE) + 1)
  
  # Assign new school IDs
  last_school_id <- max(existing_data[[schools_col]], na.rm = TRUE)
  last_school_id <- ifelse(is.na(last_school_id), 0, last_school_id)
  new_school_ids <- seq(from = last_school_id + 1, length.out = new_schools)
  
  # Create new data structure
  new_data <- expand.grid(
    schls = new_school_ids,
    ppls = 1:new_pupils_per_school
  )
  
  # Rename the school column
  names(new_data)[names(new_data) == "schls"] <- schools_col
  
  # Assign unique pupil IDs
  max_existing_pupil_id <- max(existing_data[[pupils_col]], na.rm = TRUE)
  max_existing_pupil_id <- ifelse(is.na(max_existing_pupil_id), 0, max_existing_pupil_id)
  new_data[[pupils_col]] <- seq(from = max_existing_pupil_id + 1, length.out = nrow(new_data))
  
  # Assign treatment within schools
  new_data[[intervention_col]] <- unlist(lapply(new_school_ids, function(x) {
    treatment_group_size <- round(new_pupils_per_school * treatment_percentage)
    group_assignment <- c(rep(new_treatment_num, treatment_group_size), 
                          rep(0, new_pupils_per_school - treatment_group_size))
    sample(group_assignment)
  }))
  
  # Generate school-level random effects
  new_school_intercepts <- rnorm(new_schools, mean = 0, sd = sigmab0)
  new_treatment_slopes <- rnorm(new_schools, mean = 0, sd = sigmab1)
  
  # Handle covariates properly
  for (covariate in covariates) {
    if (covariate %in% names(existing_data)) {
      unique_values <- unique(na.omit(existing_data[[covariate]]))
      
      if (length(unique_values) == 2 && all(unique_values %in% c(0, 1))) {
        prob_1 <- mean(existing_data[[covariate]] == 1, na.rm = TRUE)
        new_data[[covariate]] <- rbinom(nrow(new_data), size = 1, prob = prob_1)
        
      } else if (is.numeric(existing_data[[covariate]]) && length(unique_values) > 10) {
        covariate_mean <- mean(existing_data[[covariate]], na.rm = TRUE)
        covariate_sd <- sd(existing_data[[covariate]], na.rm = TRUE)
        new_data[[covariate]] <- rnorm(nrow(new_data), mean = covariate_mean, sd = covariate_sd)
        
      } else {
        value_counts <- table(existing_data[[covariate]])
        category_probs <- value_counts / sum(value_counts)
        new_data[[covariate]] <- sample(names(category_probs), size = nrow(new_data), replace = TRUE, prob = category_probs)
      }
    } else {
      new_data[[covariate]] <- NA
    }
  }
  
  # Assign school effects to students
  new_data$school_effects <- rep(new_school_intercepts, each = new_pupils_per_school)
  new_data$treatment_slopes <- rep(new_treatment_slopes, each = new_pupils_per_school)
  new_data$individual_residuals <- rnorm(nrow(new_data), mean = 0, sd = sigma)
  
  # Compute post-test scores
  new_data[[post_col]] <- fixed_effects[1] +  # Intercept
    new_data$school_effects + 
    new_data$treatment_slopes * new_data[[intervention_col]] + 
    new_data$individual_residuals
  
  # Apply effect size **only to the new treatment group**
  effect_size_adjustment <- (new_data[[intervention_col]] == new_treatment_num) * 
    (es * sqrt(sigma^2 + sigmab0^2 + sigmab1^2))
  
  new_data[[post_col]] <- new_data[[post_col]] + effect_size_adjustment
  
  # Apply attrition
  attrition_size <- round(nrow(new_data) * attrition_rate)
  attrition_idx <- sample(1:nrow(new_data), attrition_size)
  new_data[[post_col]][attrition_idx] <- NA
  
  # Merge data
  new_data <- new_data[, names(existing_data)]
  final_combined_data <- rbind(existing_data, new_data)
  
  return(final_combined_data)
}

add_new_treatmentsrt <- function(existing_data, new_pupils, es, attrition_rate, 
                                 post_col, intervention_col, pupils_col, covariates) {
  
  # Extract standard deviation from existing data
  sigma <- sd(existing_data[[post_col]], na.rm = TRUE)
  
  # Determine new treatment number
  new_treatment_num <- ifelse(is.na(max(existing_data[[intervention_col]], na.rm = TRUE)), 
                              1, 
                              max(existing_data[[intervention_col]], na.rm = TRUE) + 1)
  
  # Select existing ID column
  last_pupil_id <- max(existing_data[[pupils_col]], na.rm = TRUE)
  last_pupil_id <- ifelse(is.na(last_pupil_id), 0, last_pupil_id)
  new_pupil_ids <- seq(from = last_pupil_id + 1, length.out = new_pupils)
  
  # Create new dataset, keeping only selected variables
  new_data <- existing_data[1:new_pupils, ]  # Duplicate the structure of existing data
  new_data[[pupils_col]] <- new_pupil_ids  # Assign new IDs from existing column
  new_data[[intervention_col]] <- new_treatment_num  # Assign new intervention
  
  # Handle covariates (keep same distribution as existing data)
  for (covariate in covariates) {
    if (covariate %in% names(existing_data)) {
      unique_values <- unique(na.omit(existing_data[[covariate]]))
      
      if (length(unique_values) == 2 && all(unique_values %in% c(0, 1))) {
        # Binary Variable: Sample using the same proportion from existing data
        prob_1 <- mean(existing_data[[covariate]] == 1, na.rm = TRUE)
        new_data[[covariate]] <- rbinom(nrow(new_data), size = 1, prob = prob_1)
        
      } else if (is.numeric(existing_data[[covariate]]) && length(unique_values) > 10) {
        # Continuous Variable: Generate using normal distribution
        covariate_mean <- mean(existing_data[[covariate]], na.rm = TRUE)
        covariate_sd <- sd(existing_data[[covariate]], na.rm = TRUE)
        new_data[[covariate]] <- rnorm(nrow(new_data), mean = covariate_mean, sd = covariate_sd)
        
      } else {
        # Categorical Variable: Sample based on observed proportions
        value_counts <- table(existing_data[[covariate]])
        category_probs <- value_counts / sum(value_counts)
        new_data[[covariate]] <- sample(names(category_probs), size = nrow(new_data), replace = TRUE, prob = category_probs)
      }
    } else {
      new_data[[covariate]] <- NA  # If covariate not in existing data, set NA
    }
  }
  
  # Generate post-test scores based on effect size and existing distribution
  effect_size_adjustment <- es * sigma
  new_data[[post_col]] <- mean(existing_data[[post_col]], na.rm = TRUE) + effect_size_adjustment + 
    rnorm(nrow(new_data), mean = 0, sd = sigma)
  
  # Apply attrition randomly by setting some post-test scores to NA
  attrition_size <- round(nrow(new_data) * attrition_rate)
  attrition_idx <- sample(1:nrow(new_data), attrition_size)
  new_data[[post_col]][attrition_idx] <- NA
  
  # Ensure consistency with existing data structure (remove unwanted columns)
  new_data <- new_data[names(existing_data)]
  
  # Combine with existing dataset
  final_combined_data <- rbind(existing_data, new_data)
  
  return(final_combined_data)
}




# Function to plot futility decision for CRT directly within Shiny
plot_futility_decision_shiny <- function(data, post_vars = "post1", intervention_column = "treatments", Random = "schls", Nsim = 2000, covariates = NULL, VerticalLine = NULL, ProbThreshold = NULL, threshold_range = c(0, 1.0)) {
  
  
  # Store the thresholds based on the provided range
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  
  thresholds <- seq(threshold_range[1], threshold_range[2], by = 0.1)
  
  # Initialize a list to store the crtBayes output for each threshold
  output <- list()
  
  # Prepare a matrix to hold the probabilities for each intervention and threshold
  probabilities_matrix <- matrix(NA, nrow = length(thresholds), ncol = length(interventions))
  
  for (i in seq_along(interventions)) {
    # Loop over each threshold value
    for (j in 1:length(thresholds)) {
      # Select intervention group and corresponding control group (treatments == 0)
      intervention <- interventions[i]
      intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
      
      # Change intervention column for the selected intervention group to 1
      intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
      
      # Construct formula
      if (is.null(covariates) || length(covariates) == 0) {
        formula_str <- paste(post_vars, "~", intervention_column)
      } else {
        formula_str <- paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
      }
      
      # Call the respective function based on the method
      output[[j]] <- crtBayes(
        as.formula(formula_str),
        random = Random,
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = thresholds[j]
      )
      
      # Extract the Total1 probability for the current threshold
      probabilities_matrix[j, i] <- as.numeric(output[[j]]$ProbES[[1]]["Total1"])
    }
  }
  
  # Convert to data frame for ggplot
  df <- data.frame(
    Threshold = rep(thresholds, length(interventions)),
    Probability = as.vector(probabilities_matrix),
    Intervention = factor(rep(paste("Intervention", 1:length(interventions)), each = length(thresholds)))
  )
  
  # Custom colors
  dim_blue <- rgb(0, 0, 1, alpha = 0.5)
  dim_red <- rgb(1, 0, 0, alpha = 0.5)
  
  # Default breaks for x and y axes
  default_x_breaks <- pretty(thresholds)
  default_y_breaks <- pretty(c(0, 1))
  
  # Ensure VerticalLine and ProbThreshold are included in the axis breaks
  if (!is.null(VerticalLine)) {
    default_x_breaks <- unique(c(VerticalLine, default_x_breaks))
  }
  if (!is.null(ProbThreshold)) {
    default_y_breaks <- unique(c(ProbThreshold, default_y_breaks))
  }
  
  # Define custom labels with color
  x_labels <- sapply(default_x_breaks, function(tick) {
    if (!is.null(VerticalLine) && tick == VerticalLine) {
      return(sprintf("<span style='color:%s;'>%s</span>", dim_blue, tick))
    } else {
      return(as.character(tick))
    }
  })
  
  y_labels <- sapply(default_y_breaks, function(tick) {
    if (!is.null(ProbThreshold) && tick == ProbThreshold) {
      return(sprintf("<span style='color:%s;'>%s</span>", dim_red, tick))
    } else {
      return(as.character(tick))
    }
  })
  
  # Create ggplot
  p <- ggplot(df, aes(x = Threshold, y = Probability, color = Intervention)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    labs(title = "Posterior Probabilities Across Thresholds",
         x = "Threshold",
         y = "Posterior Probability") +
    scale_y_continuous(limits = c(0, 1), 
                       breaks = default_y_breaks,
                       labels = y_labels) +
    scale_x_continuous(breaks = default_x_breaks,
                       labels = x_labels) +
    theme_pubclean() +
    theme(
      legend.position = "top",  # Move legend to the top
      legend.title = element_blank(),  # Remove legend title
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black border
      axis.title.x = element_text(margin = margin(t = 10)),  # Space for the x-axis title
      axis.title.y = element_text(margin = margin(r = 10)),  # Space for the y-axis title
      axis.text.x = element_markdown(),  # Enable HTML rendering for x-axis
      axis.text.y = element_markdown()   # Enable HTML rendering for y-axis
    ) +
    coord_cartesian(clip = "off")  # Disable clipping to allow text outside the panel
  
  # Add vertical line if specified
  if (!is.null(VerticalLine)) {
    p <- p + geom_vline(xintercept = VerticalLine, linetype = "dashed", color = dim_blue, size = 1)
  }
  
  # Add horizontal line if specified
  if (!is.null(ProbThreshold)) {
    p <- p + geom_hline(yintercept = ProbThreshold, linetype = "dashed", color = dim_red, size = 1)
  }
  
  # Print the plot
  print(p)
}


# Function to plot futility decision for CRT directly within Shiny
plot_futility_decision_shinymst <- function(data, post_vars = "post1", intervention_column = "treatments", Random = "schls", Nsim = 2000, covariates = NULL, VerticalLine = NULL, ProbThreshold = NULL, threshold_range = c(0, 1.0)) {
  
  
  # Store the thresholds based on the provided range
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  
  thresholds <- seq(threshold_range[1], threshold_range[2], by = 0.1)
  
  # Initialize a list to store the mstBayes output for each threshold
  output <- list()
  
  # Prepare a matrix to hold the probabilities for each intervention and threshold
  probabilities_matrix <- matrix(NA, nrow = length(thresholds), ncol = length(interventions))
  
  for (i in seq_along(interventions)) {
    # Loop over each threshold value
    for (j in 1:length(thresholds)) {
      # Select intervention group and corresponding control group (treatments == 0)
      intervention <- interventions[i]
      intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
      
      # Change intervention column for the selected intervention group to 1
      intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
      
      # Construct formula
      if (is.null(covariates) || length(covariates) == 0) {
        formula_str <- paste(post_vars, "~", intervention_column)
      } else {
        formula_str <- paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
      }
      
      # Call the respective function based on the method
      output[[j]] <- mstBayes(
        as.formula(formula_str),
        random = Random,
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = thresholds[j]
      )
      
      # Extract the Total1 probability for the current threshold
      probabilities_matrix[j, i] <- as.numeric(output[[j]]$ProbES[[1]]["Total1"])
    }
  }
  
  # Convert to data frame for ggplot
  df <- data.frame(
    Threshold = rep(thresholds, length(interventions)),
    Probability = as.vector(probabilities_matrix),
    Intervention = factor(rep(paste("Intervention", 1:length(interventions)), each = length(thresholds)))
  )
  
  # Custom colors
  dim_blue <- rgb(0, 0, 1, alpha = 0.5)
  dim_red <- rgb(1, 0, 0, alpha = 0.5)
  
  # Default breaks for x and y axes
  default_x_breaks <- pretty(thresholds)
  default_y_breaks <- pretty(c(0, 1))
  
  # Ensure VerticalLine and ProbThreshold are included in the axis breaks
  if (!is.null(VerticalLine)) {
    default_x_breaks <- unique(c(VerticalLine, default_x_breaks))
  }
  if (!is.null(ProbThreshold)) {
    default_y_breaks <- unique(c(ProbThreshold, default_y_breaks))
  }
  
  # Define custom labels with color
  x_labels <- sapply(default_x_breaks, function(tick) {
    if (!is.null(VerticalLine) && tick == VerticalLine) {
      return(sprintf("<span style='color:%s;'>%s</span>", dim_blue, tick))
    } else {
      return(as.character(tick))
    }
  })
  
  y_labels <- sapply(default_y_breaks, function(tick) {
    if (!is.null(ProbThreshold) && tick == ProbThreshold) {
      return(sprintf("<span style='color:%s;'>%s</span>", dim_red, tick))
    } else {
      return(as.character(tick))
    }
  })
  
  # Create ggplot
  p <- ggplot(df, aes(x = Threshold, y = Probability, color = Intervention)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    labs(title = "Posterior Probabilities Across Thresholds",
         x = "Threshold",
         y = "Posterior Probability") +
    scale_y_continuous(limits = c(0, 1), 
                       breaks = default_y_breaks,
                       labels = y_labels) +
    scale_x_continuous(breaks = default_x_breaks,
                       labels = x_labels) +
    theme_pubclean() +
    theme(
      legend.position = "top",  # Move legend to the top
      legend.title = element_blank(),  # Remove legend title
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black border
      axis.title.x = element_text(margin = margin(t = 10)),  # Space for the x-axis title
      axis.title.y = element_text(margin = margin(r = 10)),  # Space for the y-axis title
      axis.text.x = element_markdown(),  # Enable HTML rendering for x-axis
      axis.text.y = element_markdown()   # Enable HTML rendering for y-axis
    ) +
    coord_cartesian(clip = "off")  # Disable clipping to allow text outside the panel
  
  # Add vertical line if specified
  if (!is.null(VerticalLine)) {
    p <- p + geom_vline(xintercept = VerticalLine, linetype = "dashed", color = dim_blue, size = 1)
  }
  
  # Add horizontal line if specified
  if (!is.null(ProbThreshold)) {
    p <- p + geom_hline(yintercept = ProbThreshold, linetype = "dashed", color = dim_red, size = 1)
  }
  
  # Print the plot
  print(p)
}


# Function to plot futility decision for CRT directly within Shiny
plot_futility_decision_shinysrt <- function(data, post_vars = "post1", intervention_column = "treatments", Nsim = 2000, covariates = NULL, VerticalLine = NULL, ProbThreshold = NULL, threshold_range = c(0, 1.0)) {
  
  # Check if no covariates are selected
  if (is.null(covariates) || length(covariates) == 0) {
    # If no covariates are provided, add a column filled with zeros
    data$zero_covariate <- 0
    covariates <- "zero_covariate"  # Use this as the covariate in the model
  }
  
  # Store the thresholds based on the provided range
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  
  thresholds <- seq(threshold_range[1], threshold_range[2], by = 0.1)
  
  # Initialize a list to store the srtBayes output for each threshold
  output <- list()
  
  # Prepare a matrix to hold the probabilities for each intervention and threshold
  probabilities_matrix <- matrix(NA, nrow = length(thresholds), ncol = length(interventions))
  
  for (i in seq_along(interventions)) {
    # Loop over each threshold value
    for (j in 1:length(thresholds)) {
      # Select intervention group and corresponding control group (treatments == 0)
      intervention <- interventions[i]
      intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
      
      # Change intervention column for the selected intervention group to 1
      intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
      
      # Construct formula
      formula_str <- paste(post_vars, "~", intervention_column)
      if (!is.null(covariates) && length(covariates) > 0) {
        formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
      }
      interventioncall <- paste0(intervention_column, 1)
      
      # Call the respective function based on the method
      output[[j]] <- srtBayes(
        as.formula(formula_str),
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = thresholds[j]
      )
      
      # Extract the Total1 probability for the current threshold
      probabilities_matrix[j, i] <- as.numeric(output[[j]]$ProbES[[interventioncall]]["Cond"])
    }
  }
  
  # Convert to data frame for ggplot
  df <- data.frame(
    Threshold = rep(thresholds, length(interventions)),
    Probability = as.vector(probabilities_matrix),
    Intervention = factor(rep(paste("Intervention", 1:length(interventions)), each = length(thresholds)))
  )
  
  # Custom colors
  dim_blue <- rgb(0, 0, 1, alpha = 0.5)
  dim_red <- rgb(1, 0, 0, alpha = 0.5)
  
  # Default breaks for x and y axes
  default_x_breaks <- pretty(thresholds)
  default_y_breaks <- pretty(c(0, 1))
  
  # Ensure VerticalLine and ProbThreshold are included in the axis breaks
  if (!is.null(VerticalLine)) {
    default_x_breaks <- unique(c(VerticalLine, default_x_breaks))
  }
  if (!is.null(ProbThreshold)) {
    default_y_breaks <- unique(c(ProbThreshold, default_y_breaks))
  }
  
  # Define custom labels with color
  x_labels <- sapply(default_x_breaks, function(tick) {
    if (!is.null(VerticalLine) && tick == VerticalLine) {
      return(sprintf("<span style='color:%s;'>%s</span>", dim_blue, tick))
    } else {
      return(as.character(tick))
    }
  })
  
  y_labels <- sapply(default_y_breaks, function(tick) {
    if (!is.null(ProbThreshold) && tick == ProbThreshold) {
      return(sprintf("<span style='color:%s;'>%s</span>", dim_red, tick))
    } else {
      return(as.character(tick))
    }
  })
  
  # Create ggplot
  p <- ggplot(df, aes(x = Threshold, y = Probability, color = Intervention)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    labs(title = "Posterior Probabilities Across Thresholds",
         x = "Threshold",
         y = "Posterior Probability") +
    scale_y_continuous(limits = c(0, 1), 
                       breaks = default_y_breaks,
                       labels = y_labels) +
    scale_x_continuous(breaks = default_x_breaks,
                       labels = x_labels) +
    theme_pubclean() +
    theme(
      legend.position = "top",  # Move legend to the top
      legend.title = element_blank(),  # Remove legend title
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black border
      axis.title.x = element_text(margin = margin(t = 10)),  # Space for the x-axis title
      axis.title.y = element_text(margin = margin(r = 10)),  # Space for the y-axis title
      axis.text.x = element_markdown(),  # Enable HTML rendering for x-axis
      axis.text.y = element_markdown()   # Enable HTML rendering for y-axis
    ) +
    coord_cartesian(clip = "off")  # Disable clipping to allow text outside the panel
  
  # Add vertical line if specified
  if (!is.null(VerticalLine)) {
    p <- p + geom_vline(xintercept = VerticalLine, linetype = "dashed", color = dim_blue, size = 1)
  }
  
  # Add horizontal line if specified
  if (!is.null(ProbThreshold)) {
    p <- p + geom_hline(yintercept = ProbThreshold, linetype = "dashed", color = dim_red, size = 1)
  }
  
  # Print the plot
  print(p)
}



run_analysis <- function(data, post_vars = "post1", intervention_column = "treatments", Random = "schls", Nsim = 2000, 
                         Threshold = 0.05, ProbThreshold = 0.8, method = "crtBayes", crtFREQoption = "Default", 
                         nPerm = NULL, nBoot = NULL, bootType = NULL, covariates = NULL) {
  
  if (is.null(covariates) || length(covariates) == 0) {
    # No covariates provided; construct formula without the zero_covariate
    formula_str <- paste(post_vars, "~", intervention_column)
  } else {
    # Covariates provided; include them in the formula
    formula_str <- paste(post_vars, "~", intervention_column, "+", paste(covariates, collapse = " + "))
  }
  
  
  # Get unique interventions
  interventions <- sort(unique(data[[intervention_column]][data[[intervention_column]] != 0]))
  
  output <- list()
  
  for (i in seq_along(interventions)) {
    # Select intervention group and corresponding control group (treatments == 0)
    intervention <- interventions[i]
    intervention_data <- subset(data, data[[intervention_column]] == intervention | data[[intervention_column]] == 0)
    
    # Change intervention column for the selected intervention group to 1
    intervention_data[[intervention_column]] <- ifelse(intervention_data[[intervention_column]] == intervention, 1, 0)
    
    # Construct formula
    formula_str <- paste(post_vars, "~", intervention_column)
    if (!is.null(covariates) && length(covariates) > 0) {
      formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
    }
    
    # Call the respective function based on the method
    if (method == "crtBayes") {
      output[[i]] <- crtBayes(
        as.formula(formula_str),
        random = Random,
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = Threshold
      )
    } else if (method == "crtFREQ") {
      if (crtFREQoption == "Default") {
        output[[i]] <- crtFREQ(
          as.formula(formula_str),
          random = Random,
          intervention = intervention_column,
          data = intervention_data
        )
      } else if (crtFREQoption == "Permutation") {
        output[[i]] <- crtFREQ(
          as.formula(formula_str),
          random = Random,
          intervention = intervention_column,
          data = intervention_data,
          nPerm = nPerm
        )
      } else if (crtFREQoption == "Bootstrap") {
        output[[i]] <- crtFREQ(
          as.formula(formula_str),
          random = Random,
          intervention = intervention_column,
          data = intervention_data,
          nBoot = nBoot,
          type = bootType
        )
      }
    } else if (method == "mstBayes") {
      output[[i]] <- mstBayes(
        as.formula(formula_str),
        random = Random,
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = Threshold
      )
    } else if (method == "mstFREQ") {
      if (crtFREQoption == "Default") {
        output[[i]] <- mstFREQ(
          as.formula(formula_str),
          random = Random,
          intervention = intervention_column,
          data = intervention_data
        )
      } else if (crtFREQoption == "Permutation") {
        output[[i]] <- mstFREQ(
          as.formula(formula_str),
          random = Random,
          intervention = intervention_column,
          data = intervention_data,
          nPerm = nPerm
        )
      } else if (crtFREQoption == "Bootstrap") {
        output[[i]] <- mstFREQ(
          as.formula(formula_str),
          random = Random,
          intervention = intervention_column,
          data = intervention_data,
          nBoot = nBoot,
          type = bootType
        )
      }
    } else if (method == "srtBayes") {
      output[[i]] <- srtBayes(
        as.formula(formula_str),
        intervention = intervention_column,
        nsim = Nsim,
        data = intervention_data,
        threshold = Threshold
      )
    } else if (method == "srtFREQ") {
      if (crtFREQoption == "Default") {
        output[[i]] <- srtFREQ(
          as.formula(formula_str),
          intervention = intervention_column,
          data = intervention_data
        )
      } else if (crtFREQoption == "Permutation") {
        output[[i]] <- srtFREQ(
          as.formula(formula_str),
          intervention = intervention_column,
          data = intervention_data,
          nPerm = nPerm
        )
      } else if (crtFREQoption == "Bootstrap") {
        output[[i]] <- srtFREQ(
          as.formula(formula_str),
          intervention = intervention_column,
          data = intervention_data,
          nBoot = nBoot
        )
      }
    }
  }
  
  # Generate comparison plot using ComparePlot
  result <- ComparePlot5(output, modelNames = paste("Intervention", 1:length(interventions)), group = 1)
  
  # Return the result which includes the plot
  return(result)
}






