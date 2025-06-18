#Jim's New Diffusion Set-Up
#Jonathan H. Morgan & Tyler Barrett
#19 October 2023

# Clear Out Console Script
  cat("\014")
  rm(list = ls(all.names = TRUE))

  # Set Working Directory
  if (paste(Sys.info()[[1]], collapse=" ") == "Darwin"){
    setwd("/Users/tylerbarrett/Dropbox (Duke Bio_Ea)/Village_Simulations/villageData/village_networks_2023")
    getwd()
  }else{
    setwd("/mnt/d/Dropbox/DNAC/IDEANet/Village_Simulations/villageData/village_networks_2023")
    getwd()
  }

# Options
  options(stringsAsFactors = FALSE)
  options(mc.cores = parallel::detectCores())
  
################
#   PACKAGES   #
################
  
library(expm)

#################
#   FUNCTIONS   #
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
  # b_int = intercept
  # b_close = effect of closeness
  # b_cxn_peers = effect of infected peers
  # b_cxn_total = effect of total total infected nodes
  # b_cxn_symp = effect symptomatic peers
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
  
###############################
#   TEST DIFFUSION FUNCTION   #
###############################
  
# Read In Test Network
  test_net <- igraph::read_graph("WeakCore1_3.2_2.net", format = c("pajek"))
  
# Prep Data for Simulation
  sim_data <- sim_prep(test_net)

# Set Starting Values
  infectedp <- c(5, 10, 15, 20, 25)
  inf_r <- 0.33
  rec_t <- 14
  maxtime <- 200
  p_symp <- 0.5
  b_int <- -0.1
  b_cxn_symp <- -1.5
  b_cxn_total <- -3.5
  b_cxn_peers <- -0.5
  b_close <- 1.0
  b_cls_x_smp <- -0.1
  numitter <- 200
  
# Analyze the Simulation Function
  Rprof("Simulation_Profile.out")
    infection_log <- sirdif(sim_data$alst, sim_data$vlst, infectedp, inf_r, rec_t, maxtime, 
                            p_symp, b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, 
                            b_cls_x_smp)
  Rprof(NULL)
  summaryRprof("Simulation_Profile.out")

# Run Simulation
  infection_data <- vector('list', 100)
  names(infection_data) <- paste0('Run_', seq(1,100,1))
  
  infection_times <- vector('numeric', 100)
  names(infection_times) <- paste0('Run_', seq(1,100,1))
  for (i in seq_along(infection_data)){
    # Run Simulation
      infection_log <- sirdif(sim_data$alst, sim_data$vlst, infectedp, inf_r, rec_t, maxtime, 
                              p_symp, b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, 
                              b_cls_x_smp)
    
    # Populate Lists
      infection_data[[i]] <- infection_log[[1]]
      infection_times[[i]] <- infection_log[[2]]
  }
  rm(infection_log)
  
  infection_history <- do.call('rbind', infection_data)
  iteration_ids <- gsub('Run_','',rownames(infection_history))
  infection_history <- cbind(as.integer(floor(as.numeric(iteration_ids))), infection_history)
  colnames(infection_history)[[1]] <- c('simulation') 
  readr::write_csv(infection_history, 'Revised_Simulation_200Steps_100Runs_Check.csv')
  
# Comparing Simulation Performance
  r_sim_times <- infection_times
  
  o_x <- density(o_sim_times)$x
  r_x <- density(r_sim_times)$x
  x_coord <- unique(c(o_x, r_x))
  
  o_y <- density(o_sim_times)$y
  r_y <- density(r_sim_times)$y
  y_coord <- unique(c(o_y, r_y))
  
  plot(0, type='n', main='Simulation Time Comparisons', ylab=c('Density'),
       xlab=c('Elapsed Time'), xlim=c(min(x_coord), max(x_coord)), 
       ylim=c(min(y_coord), max(y_coord)), col='blue',las=1, bty='n')
  
  reference_lines <- c(10, 15, 20, 25)
  for (i in seq_along(reference_lines)){
    abline(v = reference_lines[[i]], col='grey', lwd = 2, lty=2)
  }
  
  lines(density(r_sim_times), col='blue')
  abline(v = mean(r_sim_times), col='blue')
  median(r_sim_times)
  arithmetic_mode(r_sim_times)
  min(r_sim_times)
  max(r_sim_times)
  
  lines(density(o_sim_times), col='brown')
  abline(v = mean(o_sim_times), col='brown')
  median(o_sim_times)
  arithmetic_mode(o_sim_times)
  min(o_sim_times)
  max(o_sim_times)
  
  legend("topright", legend=c("Simulation", "Revised Simulation"), col=c('brown', 'blue'), lty=1:1, cex=0.75, bty='n')
  
###########################################################
#   SAS VS. R COMPARATIVE ASSESSMENT: WEAK CORE NETWORK   #
###########################################################
  
# Load SAS
  sas_asim_data <-  readr::read_csv(file="SAS_Aggregate_Simulation_Results.csv")
  sas_full_data <- readr::read_csv(file="SAS_Simulation_Results.csv")
  sas_combinations <- sas_full_data[(is.na(sas_full_data$B_INT) != TRUE ),]
  
# Read In Test Network
  test_net <- igraph::read_graph("WeakCore1_3.2_2.net", format = c("pajek"))
  
# Prep Data for Simulation
  sim_data <- sim_prep(test_net)
  
# Specifying Comparison Function
  directory <-  "/mnt/d/Dropbox/DNAC/IDEANet/Village_Simulations/villageData/village_networks_2023/test_sim_results"
  sas_compare <- function(sas_combinations, sim_data, max_iteration=nrow(sas_combinations),
                          directory){
    # Getting Start Time
      start_time <- Sys.time()
      print(start_time)
    
    # Looping Over Combinations
      simulation_list <- vector('list', nrow(sas_combinations))
      names(simulation_list) <- seq(1, nrow(sas_combinations), 1)
      for (i in seq_along(simulation_list[c(1:max_iteration)])){
        # Isolate Combination Vector
          parameters <- sas_combinations[i,]
          
        # Populating Parameters
          infectedp <- parameters$SEEDNODE
          inf_r <- 0.33
          rec_t <- 14
          maxtime <- 200
          p_symp <- 0.5
          b_int <- parameters$B_INT
          b_cxn_symp <- parameters$B_SELF
          b_cxn_total <- parameters$B_CXN_GLOBAL
          b_cxn_peers <- -parameters$B_CXN_PEER
          b_close <- parameters$B_CLOSE
          b_cls_x_smp <- parameters$B_CLS_X_SMP
        
          iteration_data <- sirdif(sim_data$alst, sim_data$vlst, infectedp, inf_r, rec_t, maxtime, p_symp,
                                      b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, b_cls_x_smp)
          
          iteration_data$keep <- iteration_data[[2]] + iteration_data[[3]] + iteration_data[[4]]
          iteration_data <- rbind(iteration_data[(iteration_data$keep > 0),], iteration_data[201,])
          iteration_data <- iteration_data[c(1:4)]
          iteration_data$case <- rep(i, nrow(iteration_data))
          iteration_data <- iteration_data[c(5,1:4)]
          
        # Populate the File
        # simulation_list[[i]] <- iteration_data
          
        # Writing-Out File
          readr::write_csv(iteration_data, file= paste0(directory,"/iteration_",i,".csv"))
      }
  
    # Collapsing List into Data Fame
    # infection_log <- do.call('rbind', simulation_list)
      
    # Getting Stop Time
      stop_time <- Sys.time()
      total_time <- stop_time - start_time
      print(total_time)
        
    # Return Infection Log
    # return(infection_log)
  }
  
# Looping Over Parameters to Simulate Values: Time Check
  #for (i in seq_along(seq(1, 12, 1))){
  #  r_infection_log <- sas_compare(sas_combinations, sim_data, i, directory)
  #}
  
# Visualizing Time Plot
  r_sim_time_log <- readr::read_csv("R_Simulation_TimeLog.csv")
  plot(x=r_sim_time_log$`Number of Iterations`, y=r_sim_time_log[[5]], bty='n', las=1, type='l',
       ylab="Elapsed Time (s)", xlab="Number of Simulations", col='black')
  grid(lwd = 2)
  points(x=r_sim_time_log[[1]], y=r_sim_time_log[[5]], pch=16, col='blue')
  
# Running SAS Compare Where We Write-Out the Files
  sas_compare(sas_combinations, sim_data, 12, directory)

###################################
#   FUNCTION PLOT DEMONSTRATION   #
###################################
  
# Inputs & Parameters
  b_int <- -0.1           # Baseline activity
  b_cxn_symp <- -1.5      # Activity lower if sick
  b_cxn_global <- -3.5    # Activity drop global sick
  b_cxn_peers <- -0.5     # Activity drop local sick
  b_close <- 1.0          # Activity boost by edge weight
  b_cls_x_smp <- -0.1;
  edgemax = 3
  maxdeg = 5
  n <- 500
  
# Looking at Function Performance
  function_list <- simparmspace(b_int, b_cxn_symp, b_cxn_global, b_cxn_peers, 
                                b_close, b_cls_x_smp, edgemax, maxdeg, n)
  probact <- function_list[[1]]
  function_plot <- function_list[[2]]
  
  function_plot
  ggplot2::ggsave("R_SocialFeedback_Simulation_FunctionPlot.png", dpi=600, width = 10.6806, height = 7.30556)

#############
#   TESTS   #
#############
  
# Importing SAS Results
  sas_function_results <- readr::read_csv('/mnt/d/Dropbox/DNAC/IDEANet/Village_Simulations/SAS_Simulation_Outputs.csv')
  
# Checking Mean & SDs of Simulation Outputs Compared to SAS
  r_sim_means <- vector('numeric', ncol(probact))
  for (i in seq_along(r_sim_means)){
    r_sim_means[[i]] <- mean(probact[[i]])
    
  }
  
  sas_sim_means <- vector('numeric', ncol(sas_function_results))
  for (i in seq_along(sas_sim_means)){
    sas_sim_means[[i]] <- mean(sas_function_results[[i]])
  }
  
  r_sim_sd <- vector('numeric', ncol(probact))
  for (i in seq_along(probact)){
    r_sim_sd[[i]] <- sd(probact[[i]])
  }
  
  sas_sim_sd <- vector('numeric', ncol(sas_function_results))
  for (i in seq_along(sas_sim_sd)){
    sas_sim_sd[[i]] <- sd(sas_function_results[[i]])
  }
  
  comparative_results <- data.frame(r_mean = r_sim_means, r_sd = r_sim_sd, sas_mean=sas_sim_means, sas_sd=sas_sim_sd)
  prob_act_scores <- data.frame(sas = sas_function_results$prob_act, R = probact$prob_act)
  
# Visualizing Function Distributions
  comparative_function_plotter <- function(b_int, b_cxn_symp, b_cxn_global, b_cxn_peers, b_close, 
                                           probact, sas_function_results){
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
        sas_colors <- color_assignment(sas_function_results$prop_inf)
        
        r_sample_data <- cbind(probact, r_colors)
        sas_sample_data <- cbind(sas_function_results, sas_colors)
        
        r_panel_data <- r_sample_data[(r_sample_data$issymptomatic == issymptomatic), ]
        sas_panel_data <- sas_sample_data[(sas_sample_data$issymptomatic == issymptomatic), ]
        for (edgwgt in seq(1, edgemax, 1)) {
          # Isolating Cell Data
            r_cell_data <- r_panel_data[(r_panel_data$edgwgt == edgwgt),]
            sas_cell_data <- sas_panel_data[(sas_panel_data$edgwgt == edgwgt),]
          
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
              sas_line_data <- sas_cell_data[(sas_cell_data$ninf == ninf),]
            
              lines(r_line_data$peersinf, r_line_data$prob_act, lty=2, col=r_line_data$r_colors)
              lines(sas_line_data$peersinf, sas_line_data$prob_act, lty=1, col=sas_line_data$sas_colors)
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
  
  x11(width=10.6806, height=7.30556)
  comparative_function_plotter(b_int, b_cxn_symp, b_cxn_global, b_cxn_peers, b_close, probact, sas_function_results)
  