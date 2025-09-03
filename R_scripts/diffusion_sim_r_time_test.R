# Test Diffusion Simulation Runtime in R
# Tyler Barrett, August 27, 2025

################
#  Test Setup  #
################

# Clear Out Console Script
  cat("\014")
  rm(list = ls(all.names = TRUE))

    fp <- "C:/Users/tmbar/OneDrive/Documents/GitHub/NetEPI"
  }else{
    fp <- "/Users/tylerbarrett/Documents/GitHub/NetEPI/madagascar_network_sim"
  }

# Options
  options(stringsAsFactors = FALSE)
  options(mc.cores = parallel::detectCores())

# Load Packages
  library(expm)
  library(dplyr)
  library(ggplot2)
  library(lhs)

#################
#   Functions   #
#################

# Calculating the Mode
  arithmetic_mode <- function(x, na.rm = FALSE) {
  if (na.rm) { # if na.rm is TRUE, remove NA values from input x
    x <- x[!is.na(x)]
  }
  val <- unique(x)
  return(val[which.max(tabulate(match(x, val)))])
}

# Color Functions
  color_assignment <- function(prob_infection){
  # Creating Color Palette
  hist_data <- graphics::hist(prob_infection, plot=FALSE)
  cutcol <- cut(probact$prop_inf, breaks=hist_data$breaks, right=FALSE)
  base::levels(cutcol)
  colorpal <-   colorpal <- c("slateblue","royalblue1","lightblue2","#FFFFFF", "yellow1", "navajowhite", "lightsalmon", "orange", "darkorange1","coral", "firebrick1", "brown" )
  
  # Assign the data to the 12 color groups
  color_vector <- vector('character', length(cutcol))
  for (i in seq_along(cutcol)){
    cutcol_value <- as.character(cutcol[[i]])
    cutcol_value <- strsplit(cutcol_value, ',')[[1]]
    cutcol_value <- sub(")", "", cutcol_value, fixed=TRUE)
    cutcol_value <- gsub("\\[|\\]", "", cutcol_value)
    cutcol_value <- as.numeric(cutcol_value)
    
    if(cutcol_value[[1]] >= 0 & cutcol_value[[2]] <= 0.05){
      color_vector[[i]] <-  colorpal[[1]]
    }else if(cutcol_value[[1]] >= 0.05 & cutcol_value[[2]] <= 0.10){
      color_vector[[i]] <-  colorpal[[2]]
    }else if(cutcol_value[[1]] >= 0.10 & cutcol_value[[2]] <= 0.15){
      color_vector[[i]] <-  colorpal[[3]]
    }else if(cutcol_value[[1]] >= 0.15 & cutcol_value[[2]] <= 0.2){
      color_vector[[i]] <-  colorpal[[4]]
    }else if(cutcol_value[[1]] >= 0.2 & cutcol_value[[2]] <= 0.25){
      color_vector[[i]] <-  colorpal[[5]]
    }else if(cutcol_value[[1]] >= 0.25 & cutcol_value[[2]] <= 0.3){
      color_vector[[i]] <-  colorpal[[6]]
    }else if(cutcol_value[[1]] >= 0.3 & cutcol_value[[2]] <= 0.35){
      color_vector[[i]] <-  colorpal[[7]]
    }else if(cutcol_value[[1]] >= 0.35 & cutcol_value[[2]] <= 0.4){
      color_vector[[i]] <-  colorpal[[8]]
    }else if(cutcol_value[[1]] >= 0.4 & cutcol_value[[2]] <= 0.45){
      color_vector[[i]] <-  colorpal[[9]]
    }else if(cutcol_value[[1]] >= 0.45 & cutcol_value[[2]] <= 0.5){
      color_vector[[i]] <-  colorpal[[10]]
    }else if(cutcol_value[[1]] >= 0.5 & cutcol_value[[2]] <= 0.55){
      color_vector[[i]] <-  colorpal[[11]]
    }else{
      color_vector[[i]] <-  colorpal[[12]]
    }
  }
  
  # Return Color Vector
  return(color_vector)
}

# Function to Examine the Parameter Space in Preparation for the Simulation
  simparmspace <- function(b_int, b_cxn_symp, b_cxn_global, b_cxn_peers, 
                         b_close, b_cls_x_smp, edgemax, maxdeg, n) {
  
  # Defining parameters internally
  b_int <- b_int  # baseline activity
  b_cxn_symp <- b_cxn_symp  # activity lower if sick
  b_cxn_total <- b_cxn_global  # activity drop global sick
  b_cxn_peers <- b_cxn_peers  # activity drop local sick
  b_close <- b_close  # activity boost by edge weight
  b_cls_x_smp <- b_cls_x_smp
  n <- n
  
  # Create an empty dataframe for storing results
  probact <- data.frame(prob_act = numeric(), issymptomatic = integer(), ninf = integer(), edgwgt = integer(), peersinf = integer(), b_int = numeric(), b_cxn_symp = numeric(), b_cxn_total = numeric(), b_cxn_peers = numeric())
  
  # Looping Over Parameters
  for (issymptomatic in seq(0,1,1)) {
    for (ninf in seq(1, 300, 50)) {
      for (edgwgt in seq(1, edgemax, 1)) {
        for (peersinf in seq(0, maxdeg, 1)) {
          lwact <- b_int+((issymptomatic)*(b_cxn_symp))+(b_close*edgwgt)+(b_cxn_peers*peersinf) + b_cxn_total*(ninf/(n/3))+ b_cls_x_smp*(edgwgt)*issymptomatic
          
          prob_act <- exp(lwact) / (1 + exp(lwact))
          
          probact <- rbind(probact, data.frame(prob_act = prob_act, issymptomatic = issymptomatic, ninf = ninf, edgwgt = edgwgt, peersinf = peersinf, b_int = b_int, b_cxn_symp = b_cxn_symp, b_cxn_total = b_cxn_total, b_cxn_peers = b_cxn_peers))
        }
      }
    }
  }
  
  probact$prop_inf <- probact$ninf / 500
  
  # Creating Plotting Function
  function_plot <- function(){
    # Creating Layout
    viz_matrix <- matrix(c(7,7,7,7,7,7,7,7,7,7,
                           1,1,1,2,2,2,3,3,3,8,
                           1,1,1,2,2,2,3,3,3,8,
                           1,1,1,2,2,2,3,3,3,8,
                           4,4,4,5,5,5,6,6,6,8,
                           4,4,4,5,5,5,6,6,6,8,
                           4,4,4,5,5,5,6,6,6,8,
                           9,9,9,9,9,9,9,9,9,9), 
                         ncol  = 10, byrow = TRUE)
    layout(viz_matrix)
    
    
    # Looping Over Visualization Matrix
    for (issymptomatic in seq(0,1,1)) {
      r_colors <- color_assignment(probact$prop_inf)
      
      r_sample_data <- cbind(probact, r_colors)
      
      r_panel_data <- r_sample_data[(r_sample_data$issymptomatic == issymptomatic), ]
      for (edgwgt in seq(1, edgemax, 1)) {
        # Isolating Cell Data
        r_cell_data <- r_panel_data[(r_panel_data$edgwgt == edgwgt),]
        
        # Creating Base Plot
        par(mar=c(3, 4, 2, 2))
        plot(NA,type="n", main= paste('edgwgt =', unique(r_cell_data$edgwgt)), xlab = "", 
             ylab = "Probability of Interaction", family = 'serif',
             xlim= c(0,5), ylim=  c(0,1), bty="n", las=1)
        
        # Adding Reference Lines
        ref_lines <- seq(0,1,0.2)
        for (ref_line in seq_along(ref_lines)){
          abline(h = ref_lines[[ref_line]] , col = "gray60", lty=3)
        }
        
        # Plotting X-Axis Label
        mtext(side = 1, text = 'Number of Ill Peers', col = "black", line = 2, cex = 0.75, family='serif')
        
        # Plotting Lines
        for (ninf in seq(1, 300, 50)){
          r_line_data <- r_cell_data[(r_cell_data$ninf == ninf),]
          lines(r_line_data$peersinf, r_line_data$prob_act, lty=1, col=r_line_data$r_colors)
        }
        
        # Adding Symptomatic
        if(edgwgt == 3){
          mtext(side = 4, text = paste("issymptomatic =", issymptomatic), col = "black", line = 1, cex = 0.75, family='serif', font=2)
        }else{
          edgwgt = edgwgt
        }
      }
    }
    
    # Adding Plot Title
    par(mar = c(0,0,0,0), bty='n')
    plot(0, type='n', xlab=' ', ylab=' ', cex.axis=1.3,   xaxt='n', yaxt = 'n', family='serif', las=1, main=' ')
    text(x= 0.95, y= -0.1, labels = c("SIR Social Feedback Function Plot"), cex=2.5, family='serif', font=2)
    
    # Adding Color Gradient
    legend_image <- as.raster(matrix(rev(  unique(color_assignment(probact$prop_inf)))), ncol=1)
    
    par(mar = c(0,0,0,0), bty='n')
    plot(0, type='n', xlab=' ', ylab=' ', xaxt='n', yaxt = 'n', cex.axis=1.3, family='HersheySerif', las=1, main=' ')
    
    rasterImage(legend_image, 0.90, -0.75, 1.01, 0.75)
    text(x=1.15, y = seq(-0.7,0.7,l=12), labels = seq(0,0.55,0.05) , cex=1, family='serif')
    text(x= 1, y=0.8, labels = c("Prob. Inf."), cex=1.35, family='serif')
    
    # Adding Annotation
    par(mar = c(0,0,0,0), bty='n')
    plot(0, type='n', xlab=' ', ylab=' ', cex.axis=1.3,   xaxt='n', yaxt = 'n', family='HersheySerif', las=1, main=' ')
    text(x= 0.95, y= -0.1, labels = paste(c(paste("b_int =",b_int), paste("b_cxn_symp =", b_cxn_symp), 
                                            paste("b_cxn_global =",b_cxn_global), paste("b_cxn_peers =",b_cxn_peers),
                                            paste("b_close =",b_close)), collapse=" "), cex=1.3, family='serif', font=1)
    
  }
  
  g <- cowplot::as_grob(function_plot)
  p_1 <- cowplot::ggdraw(g)
  p_1 
  
  return_list <- vector('list', length=2)
  return_list[[1]] <- probact
  return_list[[2]] <- p_1
  return(return_list)
}

# Jim's Function to Convert Valued Edgelist to Adjacency Matrix
  el2adjval <- function(valid, el1, el2, elwgt) {
  nodeset <- unique(valid)
  n <- length(nodeset)
  adjmat <- matrix(0, nrow = n, ncol = n)
  
  for (i in seq_along(el1)) {
    iv <- el1[i]
    jv <- el2[i]
    wv <- elwgt[i]
    
    iloc <- which(nodeset == iv)
    jloc <- which(nodeset == jv)
    
    if (length(iloc) > 0 && length(jloc) > 0) {
      adjmat[iloc, jloc] <- wv
    }
  }
  
  rownames(adjmat) <- c(nodeset)
  colnames(adjmat) <- rownames(adjmat)
  
  return(adjmat)
}

# Function to Prepare Network Data for Simulation
  sim_prep <- function(network){
  
  # Extract Nodelist, Edgelist, and Weights from Graph Object
  nl <- as.integer(igraph::V(network))
  el <- igraph::as_edgelist(network)
  elwgt <- igraph::E(network)$weight
  
  # Convert Edgelist to a Binary Adjacency Matrix
  adj_mat <- el2adjval(nl, el[,1], el[,2], rep(1, length(elwgt)))
  
  # Loop through the Matrix and Replace 1s with Column Names
  for (i in 1:nrow(adj_mat)) {
    row_vector <- adj_mat[i,]
    row_vector <- row_vector * as.integer(colnames(adj_mat))
    adj_mat[i,] <- row_vector
  }
  
  # Convert Adjacency Matrix to Adjacency List
  # Change 0 Values to NAs
  adj_mat[adj_mat == 0] <- NA
  
  # Shift Non-NA Values to the Left
  alst <- as.data.frame(t(apply(adj_mat, 1, function(x) {return(c(x[!is.na(x)],x[is.na(x)]))})))
  colnames(alst) = colnames(adj_mat)
  
  # Make vlst by Shifting Non-NA Values to the Left
  adj_mat <- el2adjval(nl, el[,1], el[,2], elwgt)
  adj_mat[adj_mat == 0] <- NA
  
  vlst <- as.data.frame(t(apply(adj_mat, 1, function(x) {return(c(x[!is.na(x)],x[is.na(x)]))})))
  colnames(vlst) = colnames(adj_mat)
  
  # Convert Matrices into Lists
  adj_list <- vector('list', nrow(adj_mat))
  names(adj_list) <- row.names(adj_mat)
  
  edge_list <- vector('list', nrow(adj_mat))
  names(edge_list) <- row.names(adj_mat)
  for(i in seq_along(adj_mat[,1])){
    # Isolate row vector & Populate Lists
    adj_list[[i]] <- as.integer(alst[i, ])
    edge_list[[i]] <- as.numeric(vlst[i,])
  }
  
  # Return alst and vlst for Simulation Function
  return(list(alst = adj_list, vlst = edge_list))
}

# Social Feedback Simulation Function

# Arguments for Social Feedback Simulation Function
# alst = adjacency list
# vlst = edge weights for adjacency list
# infectedp = seed nodes (i.e., set of initial infected nodes)
# inf_r = infection rate
# rec_t = recovery time (days)
# maxtime = total number of days
# p_symp = probability of being symptomatic
# b_int = intercept / baseline activity
# b_close = effect of closeness /edge weights
# b_cxn_peers = change in activity if peers are infected
# b_cxn_total = change in activity based on global prevalence
# b_cxn_symp = change in activity if symptomatic
# b_cls_x_smp = interaction term for closeness and symptomatic

  transmission_prob <- function(inf_r, edgwgt){
  trans <- rbinom(1, 1, inf_r) # transmission
  trans <- as.logical(trans)
  tranprob <- 1 - ((1 - inf_r) ^ edgwgt)
  trans <- rbinom(1, 1, tranprob)
  
  return(trans)
}

  sirdif <- function(alst, vlst, infectedp, inf_r, rec_t, maxtime, p_symp,
                     b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, b_cls_x_smp) {
    # Getting Start Time
    start_time <- Sys.time()
    
    # Creating Function Input Object for Updating
    infection_update_list <- vector('list', 3)
    names(infection_update_list) <- c('infected', 'inftime', 'nbrsinf')
    
    # Creating Function to Replace in Place
    replace_in_place <- function(x, idx, alters) {
      s_alst[[idx]][!is.na(x) & x %in% alters] <<- NA
      return(invisible(NULL))
    }
    
    # Set Number of Nodes in Simulation
    n <- length(alst)
    
    # Pre-allocate with NAs or zeros
    infection_update_list$infected <- rep(NA, n)
    infection_update_list$inftime <- rep(0, n)
    
    
    # Assign Starting Values: Set Seed Nodes
    initial_infected_length <- length(infectedp)
    infection_update_list$infected[1:initial_infected_length] <- infectedp
    
    # Make a Copy of the Adjacency List and Edge Values List to Modify
    s_alst <- alst
    s_vlst <- vlst
    
    # Create an Empty Vector to Count Number of Each Node's Neighbors Who Have Been Infected
    infection_update_list$nbrsinf <- rep(0, length(alst))
    
    # Create a New List Limited to Infected Nodes
    inf_net <- s_alst[infection_update_list$infected]
    
    # Count the Number of Each Node's Neighbors Who Have Been Infected
    for (i in 1:length(inf_net)) {
      n_i <- sum(!is.na(inf_net[[i]]))
      infection_update_list$nbrsinf[n_i] <- infection_update_list$nbrsinf[n_i] + 1
    }
    
    # Schedule to Remove People, Initial Cases at Starting Value
    infection_update_list$inftime <- as.integer(matrix(data = rec_t, nrow = 1, ncol = length(infection_update_list$infected))) 
    infection_update_list$inftime[is.na(infection_update_list$infected)] <- NA
    
    # Compute Proportion Who are Infected
    pinf <- length(infection_update_list$infected[!is.na(infection_update_list$infected)])/n
    
    # Set Initial Time Record
    timesum <- data.frame(Time = 0, Number_Infected = length(infection_update_list$infected[!is.na(infection_update_list$infected)]), 
                          ProportionEverInfected = pinf, ProportionCurrentlyInfected = 0, 
                          Number_Recovered = 0, ProportionRecovered = 0)
    
    # Start Main Simulation Loop
    for (t in seq_along(seq(1,maxtime,1))) { # index time
      if (class(infection_update_list$infected) == "numeric" & length(infection_update_list$infected) != 0) {
        # Identify the number and time infected
        nctinf <- length(infection_update_list$inftime[!is.na(infection_update_list$infected)])
        ninf <- length(infection_update_list$infected[!is.na(infection_update_list$infected)])
        
        # Looping through alters
        for (j in seq_along(infection_update_list$infected[!is.na(infection_update_list$infected)])) { # search from infected
          ego <- infection_update_list$infected[[j]]
          isSymptomatic <- rbinom(1, 1, p_symp) # is ego symptomatic (random for now)
          
          jnbr <- which(s_alst[[ego]] != 0) # susceptible neighbors
          if(length(jnbr) > 0) {
            # Getting ID
            alters <- as.integer(s_alst[[ego]][!is.na(s_alst[[ego]])])
            
            # Getting Weights 
            edgwgts <- s_vlst[[ego]][alst[[ego]] %in% alters]
            
            # Checking for Infected Peers
            peersinf <- infection_update_list$nbrsinf[alters] # count of infected peers
            
            # Calculating score
            lwact <- b_int + (isSymptomatic * b_cxn_symp) +
              (b_close * edgwgts) +
              (b_cxn_peers * peersinf) +
              (b_cxn_total * (ninf/(n/3))) +
              (b_cls_x_smp * peersinf * isSymptomatic)
            
            # Adding Probability of Edge Activation
            prob_act <- exp(lwact) / (1 + exp(lwact))
            active_link <- rbinom(length(lwact), 1, prob_act) # is the edge active
            
            # Determining transmission
            alters <- alters[(active_link == 1)]
            edgwgts <-   edgwgts[(active_link == 1)]
            
            trans <- as.integer(lapply(edgwgts, function(x)  transmission_prob(inf_r, x)))
            alters <- alters[(trans == 1)]
            edgwgts <-   edgwgts[(trans == 1)]
            
            # Updating Infection & Infection Time Lists
            if (length(alters) > 0) {
              next_idx <- which(is.na(infection_update_list$infected))[1]
              infection_update_list$infected[next_idx:(next_idx + length(alters) - 1)] <- alters
              infection_update_list$inftime[next_idx:(next_idx + length(alters) - 1)] <- rep(rec_t, length(alters))
            }
            
            # Updating Susceptible Alters in Place (Infected People Can't be Infected While Infected)
            lapply(seq_along(s_alst), function(idx) replace_in_place(s_alst[[idx]], idx, alters))
            
            # Updating Neighbor List ???
            inf_nbrs <- alst[[ego]][c(2:length(alst[[ego]]))][is.na(alst[[ego]][c(2:length(alst[[ego]]))]) == FALSE]
            infection_update_list$nbrsinf[inf_nbrs] <- infection_update_list$nbrsinf[inf_nbrs] + 1 # update number of infected neighbors
          }else{
            jnbr <- jnbr
          }
        }
        
        # Recovery
        infection_update_list$inftime[!is.na(infection_update_list$inftime)] <- infection_update_list$inftime[!is.na(infection_update_list$inftime)] - 1
        
        # Remove Individuals Who Transitioned to Recovered
        any_recovered <- infection_update_list$infected[which(infection_update_list$inftime == 0)]
        infection_update_list$infected[infection_update_list$infected %in% any_recovered] <- NA
        infection_update_list$inftime[infection_update_list$inftime == 0] <- NA
        
        # Sorting List so the order of !NA values is preserved, but NA values are pushed to the back
        infection_update_list$infected <- infection_update_list$infected[order(is.na(infection_update_list$infected))]
        infection_update_list$inftime <- infection_update_list$inftime[order(is.na(infection_update_list$inftime))]
        
        # Calculating Number Infected & Number Recovered
        NI <- length(infection_update_list$infected[!is.na(infection_update_list$infected)])
        NR <- length(any_recovered)
        
        # Calculating Proportions
        propeverinf <- (NI + NR)/n
        propcurinf <- NI/n
        proprec <- (NI > 0) * (NR / (NI + NR))
        
        # Adding to Infection History in Place
        timesum <- rbind(timesum, c(t, NI, propeverinf, propcurinf, NR, proprec))
      }else { # no infected left, end simulation
        t <- maxtime
      }
    }
    rm(any_recovered, t, NI, propeverinf, propcurinf, NR, proprec)
    
    # Getting Stop Time
    stop_time <- Sys.time()
    total_time <- stop_time - start_time
    #print(total_time)
    
    # Return Infection Log
    output_list <- vector('list',2)
    names(output_list) <- c('infection_log', 'total_time')
    output_list[[1]] <- timesum
    output_list[[2]] <- total_time
    return(output_list)
}

# Function to Plot Simulation Results
  plot_sim_results <- function(infection_data) {
  # Combine rows using rbind
  infection_history <- do.call('rbind', infection_data$infection_data)
  
  # Extract iteration IDs
  iteration_ids <- gsub('Run_', '', rownames(infection_history))
  
  # Add iteration IDs as a new column
  infection_history <- cbind(as.integer(floor(as.numeric(iteration_ids))), infection_history)
  
  # Rename the first column to 'simulation'
  colnames(infection_history)[[1]] <- 'simulation'
  
  # Create the plot
  ggplot(infection_history, aes(x = Time, y = Number_Infected, group = simulation, color = as.factor(simulation))) +
    geom_line(alpha = 0.5) +
    theme_minimal() +
    labs(title = "Diffusion Simulation Results",
         x = "Time",
         y = "Value",
         color = "Simulation") +
    scale_color_viridis_d() +
    theme(legend.position = "none")
}

# Function to Extract Data for Analysis
  extract_data <- function(infection_data, df_name) {
  
  # Combine rows using rbind
  infection_history <- do.call('rbind', infection_data$infection_data)
  
  # Extract iteration IDs
  iteration_ids <- gsub('Run_', '', rownames(infection_history))
  
  # Add iteration IDs as a new column
  infection_history <- cbind(as.integer(floor(as.numeric(iteration_ids))), infection_history)
  
  # Rename the first column to 'simulation'
  colnames(infection_history)[[1]] <- 'simulation'
  
  # Add new first column with the network name
  infection_history <- cbind(network = df_name, infection_history)
}

# Function to Compute Outcomes
  compute_outcomes <- function(data) {
  
  # Compute total infected, peak prevalence, and area under the curve
  results_summary_1 <- data %>%
    group_by(network, simulation) %>%
    summarize(total_infected = sum(Number_Recovered),
              peak_prevalence = max(Number_Infected),
              prop_infected = sum(Number_Recovered) / pop_size)
  #auc = integrate(Vectorize(function(t) approx(Time, Number_Infected, xout = t)$y),
  #               min(Time), max(Time))$value)
  
  # Compute duration of epidemic
  results_summary_2 <- data %>%
    group_by(network, simulation) %>%
    filter(Number_Infected == 0) %>%
    summarize(duration = min(Time))
  
  # Compute timing of epidemic peak
  results_summary_3 <- data %>%
    group_by(network, simulation) %>%
    filter(Number_Infected == max(Number_Infected)) %>%
    summarize(peak_time = min(Time))
  
  # Combine outcomes into one data frame
  results_summary <- results_summary_1 %>%
    left_join(results_summary_2) %>%
    left_join(results_summary_3)
  
  # Return the outcome summary
  return(results_summary)
}

# Function to Rescale Latin Hypercube Sample Results
  rescale <- function(samples, min_val, max_val) {
  return(min_val + (max_val - min_val) * samples)
}
  
#################
#   Data Prep   #
#################
  
# Load Synthetic Networks
  synthetic_nets_village_a <- readRDS("/workspace/data/sim_test_data/synthetic_nets_village_a.rds")
  synthetic_nets_village_m <- readRDS("/workspace/data/sim_test_data/synthetic_nets_village_m.rds")
  synthetic_nets_village_s <- readRDS("/workspace/data/sim_test_data/synthetic_nets_village_s.rds")

# Select Five Networks for Each Village
  set.seed(123)
  test_nets <- list(
    test_nets_a = synthetic_nets_village_a[sample(1:length(synthetic_nets_village_a), 5)],
    test_nets_m = synthetic_nets_village_m[sample(1:length(synthetic_nets_village_m), 5)],
    test_nets_s = synthetic_nets_village_s[sample(1:length(synthetic_nets_village_s), 5)]
  )
  
# Prep Data for Simulation
  sim_data <- list(
    test_nets_a = lapply(test_nets$test_nets_a, sim_prep),
    test_nets_m = lapply(test_nets$test_nets_m, sim_prep),
    test_nets_s = lapply(test_nets$test_nets_s, sim_prep)
  )
  
# Latin Hypercube Sampling to Set Parameter Values
  set.seed(123)
  param_values <- as.data.frame(randomLHS(10, 4))
  param_values <- param_values %>%
    rename(inf_r = V1,
           p_symp = V2,
           b_cxn_symp = V3,
           b_cxn_peers = V4) %>%
    mutate(inf_r = rescale(inf_r, 0.001, 0.02)) %>%
    mutate(p_symp = rescale(p_symp, 0.1, 1)) %>%
    mutate(b_cxn_symp = rescale(b_cxn_symp, -2, 0)) %>%
    mutate(b_cxn_peers = rescale(b_cxn_peers, -1, 0)) %>%
    mutate(maxtime = 200) %>% # try 1000
    mutate(rec_t = 14) %>% # try 1000
    mutate(b_int = 1) %>%
    mutate(b_cxn_total = 0) %>%
    mutate(b_close = 1) %>%
    mutate(b_cls_x_smp = 0)
  
# Set Number of Iterations
  numitter = 10
  
#####################
#   Run Simulation  #
#####################
  
# Testing the runtime for 150 simulations: 10 runs for 15 networks (5 per village).
# The network size ranges from 900 to 2,500 nodes.
  
# Initialize list to store results for each village
  nets <- c(1:5)
  villages <- c("test_nets_a", "test_nets_m", "test_nets_s")
  infection_data_all_villages <- vector('list', length(villages))
  names(infection_data_all_villages) <- villages
  
# Run Simulation
  start_time <- Sys.time()
  for (village in villages){
    
    # Set Population Size for Each Village
    pop <- as.numeric(length(sim_data[[village]][[1]]$alst))
    
    # Compute 1% of Population for Index Cases
    index <- 0.01*length(sim_data[[village]][[1]]$alst)
    
    for (i in seq_along(nets)) {
      
      # Initialize infection_data and infection_times for each network
        infection_data <- vector('list', numitter)
        names(infection_data) <- paste0('Run_', seq(1, numitter, 1))
        
        infection_times <- vector('numeric', numitter)
        names(infection_times) <- paste0('Run_', seq(1, numitter, 1))
      
      # Calculate total simulations for progress tracking
        total_sims <- length(villages) * length(nets) * numitter
        current_village_index <- match(village, villages)
      
      for (j in seq_along(infection_data)) {
        
        # Calculate current simulation number
          current_sim <- ((current_village_index - 1) * length(nets) * numitter) + 
            ((i - 1) * numitter) + j
        
        # Print progress
          cat("Running simulation", current_sim, "of", total_sims, 
              "(Village:", village, "| Network:", i, "| Run:", j, ")\n")
        
        # Run Simulation
          infection_log <- sirdif(sim_data[[village]][[i]]$alst, sim_data[[village]][[i]]$vlst,
                                  infectedp = as.numeric(sample(1:pop, index)),
                                  inf_r = param_values$inf_r[j],
                                  rec_t = param_values$rec_t[j],
                                  maxtime = param_values$maxtime[j],
                                  p_symp = param_values$p_symp[j],
                                  b_int = param_values$b_int[j],
                                  b_close = param_values$b_close[j],
                                  b_cxn_peers = param_values$b_cxn_peers[j],
                                  b_cxn_total = param_values$b_cxn_total[j],
                                  b_cxn_symp = param_values$b_cxn_symp[j],
                                  b_cls_x_smp = param_values$b_cls_x_smp[j])
        
        # Populate Lists
          infection_data[[j]] <- infection_log[[1]]
          infection_times[[j]] <- infection_log[[2]]
      }
      
      # Store results for this network within the village
        infection_data_all_villages[[village]][[i]] <- list(
          infection_data = infection_data,
          infection_times = infection_times)
      
      # Clean up
        rm(infection_log)
    }
  }
  end_time <- Sys.time()
  total_runtime <- end_time - start_time
  print(paste("Total runtime:", total_runtime))

# First test run took approximately 45 minutes for 150 simulations.
  
# Rename Results for Analysis
  infection_data <- infection_data_all_villages
  
###############################
#   Check Simulation Results  #
###############################
  
# Initialize List to Store Results Plots
  results_plots <- vector('list', length(villages))
  results_plots <- names(villages)
  
# Plot Simulation Results for Each Network
  for (village in villages) {
    
    for (i in seq_along(nets)) {
      
      # Plot results
      p <- plot_sim_results(infection_data[[village]][[i]])
      
      # Store plots in list
      results_plots[[village]][[i]] <- p
      rm(p)
    }
  }
  
# Look at a Few Plots
  results_plots$test_nets_a[[1]]
  results_plots$test_nets_m[[1]]
  results_plots$test_nets_s[[1]]
  
# Create Empty Data Frame to Store Results Data
  analysis_data <- vector('list', length(villages))
  analysis_data <- names(villages)
  
# Extract Data for Analysis
  for (village in villages) {
    
    # Initialize an empty dataframe to store results
    village_df <- data.frame()
    
    for (i in seq_along(nets)) {
      
      # Extract data from results
      df <- extract_data(infection_data[[village]][[i]], paste0('Net_', nets[[i]]))
      
      # Add data to village dataframe
      village_df <- rbind(village_df, df)
    }
    # Add village dataframe to analysis_data list
    analysis_data[[village]] <- village_df
  }
  
  analysis_data$test_nets_a <- analysis_data$test_nets_a %>%
    mutate(pop_size = 1900)
  analysis_data$test_nets_m <- analysis_data$test_nets_m %>%
    mutate(pop_size = 2700)
  analysis_data$test_nets_s <- analysis_data$test_nets_s %>%
    mutate(pop_size = 900)
  
# Compute Outcomes for Each Network and village
  results_summary <- vector('list', length(villages))
  results_summary <- names(villages)
  
  for (village in villages) {
    
    df <- compute_outcomes(analysis_data[[village]])
    
    # Create a new variable with the village name
    df$village_name <- village
    
    results_summary[[village]] <- df
  }
  
# Combine All Village Dataframes
  combined_results <- do.call(rbind, results_summary)
  
# Move Village Column
  combined_results <- combined_results %>% relocate(village_name)
  
  combined_results$village_name <- recode(combined_results$village_name, test_nets_m = 1, test_nets_s = 2, test_nets_a = 3)
  combined_results <- combined_results %>%  
    rename(village = village_name) %>%
    mutate(village = as.factor(village))
  combined_results <- distinct(combined_results)
  
# Plot Outcomes by Village
  ggplot(combined_results, aes(x = village, y = total_infected, fill = village)) +
    geom_boxplot()
  
  ggplot(combined_results, aes(x = village, y = peak_prevalence, fill = village)) +
    geom_boxplot()
  
  ggplot(combined_results, aes(x = village, y = duration, fill = village)) +
    geom_boxplot()
  
  ggplot(combined_results, aes(x = village, y = peak_time, fill = village)) +
    geom_boxplot()
  
  ggplot(combined_results, aes(x = village, y = prop_infected, fill = village)) +
    geom_boxplot()
  