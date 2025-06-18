# Social Weighting Function
# Jonathan H. Morgan, Ph.D., Dana Pasquale, Ph.D., and Tyler Barrett
# 25 August 2023

# Clear Out Console Script
  cat("\014")

# Options
  options(stringsAsFactors = FALSE)
  options(scipen = 999)

# Set Working Directory
  if (paste(Sys.info()[[1]], collapse=" ") == "Darwin"){
    setwd("/mnt/d/Dropbox/DNAC/IDEANet/Village_Simulations/villageData/village_networks_2023")
    getwd()
  }else{
    setwd("C:/Users/tmbar/Duke Bio_Ea Dropbox/Tyler Barrett/Village_Simulations/villageData/village_networks_2023/")
    getwd()
  }

#################
#   FUNCTIONS   #
#################

`%notin%` <- Negate(`%in%`)

agg_function <- function(v, na.rm = TRUE) {
  c(mean = mean(v, na.rm = na.rm), sd = sd(v, na.rm = na.rm), median = median(v, na.rm = na.rm))
}

village_import <- function(directory){
  # Determining File Locations
    file_list <- base::list.files(path=directory,pattern =".csv", full.names = TRUE)
  
  # Generating Output-List
    output_list <- vector('list', length(file_list))
  
  # Looping Over 
    for(i in seq_along(output_list)){
      # Reading-In Data as Text
        data <- readLines(file_list[[i]])
        num_lines <- length(readLines(file_list[[i]])) - 1
    
      # Importing into a List of DataFrames
        data_lines <- vector('list', length=num_lines)
        for(j in 1:num_lines){
          data_line <- as.data.frame(  t(strsplit(data[[j+1]], ',')[[1]]) )
          data_lines[[j]] <- data_line
        }
    
      # Stacking DataFrames
        file_df <- do.call('rbind', data_lines)
        colnames(file_df) <- strsplit(data[[1]], ',')[[1]]
    
      # Populating Output List
        output_list[[i]] <- file_df
        rm(data, num_lines, file_df)
    }
  
  # Return List
    return(output_list)
}

village_cluster <- function(summary_network_list){
  # Transform Edgelist into Matrix for the Purposes Generating a Graph
    edges <- as.matrix(summary_network_list$edgelist[c(3,5,6)])
  
  # Creating igraph object
    colnames(summary_network_list$nodelist)[[2]] <- c('attr')
    g <- igraph::graph_from_data_frame(d = edges[,c(1,2)], directed = FALSE, vertices = summary_network_list$nodelist[c(1:2)]) 
  
  # compute weighted out degree
    igraph::degree(g, mode='out', loops=FALSE)
  
  # Adding edge weights
    igraph::edge.attributes(g)$weight <- edges[,3]
  
  # Perform Community Detection
    r <- 1
    ldc <- igraph::cluster_leiden(g, resolution_parameter=r)
  
  # Transform into a DataFrame of Values
    community_list <- data.frame(node_id = summary_network_list$nodelist$id, node_label = summary_network_list$nodelist$attr , 
                               community = ldc$membership)
  
  # Return community_list
    return(community_list)
}

village_cluster_unweighted <- function(network){
  
  # Convert Network Object to Igraph Object
  g <- intergraph::asIgraph(network)
  
  # Perform Community Detection
  r <- 0.3
  ldc <- igraph::cluster_leiden(igraph::as.undirected(g), resolution_parameter=r)
  
  # Transform into a DataFrame of Values
  community_list <- data.frame(node_id = igraph::V(g)$vertex.names,
                               community = ldc$membership)
  
  # Return community_list
  return(community_list)
}

################
#   PACKAGES   #
################

library(corrgram)
library(psych)
library(statnet)
library(network)
library(tidyverse)

######################
#   IMPORTING DATA   #
######################

directory <- getwd()
load(paste0(directory,'/','village_data_20July2023.Rda'))

#######################
#   FORMATTING DATA   #
#######################

test_cases <- c(83, 43, 58)
all_cases <- sort(unique(as.numeric(village_data[[2]]$villagecode)))
# test_cases <- all_cases
village_subestting <- function(community_list, village_data){
  # Isolating Nodelist from Data List
    node_list <- village_data[[2]]
  
  # Subset List by community_list
    node_list <- node_list[(node_list[[2]] %in% community_list),]
  
  # Isolating Raw Edgelist
    edge_list <- village_data[[1]]
  
  # Subset List by community_list
    edge_list <- edge_list[(edge_list[[6]] %in% community_list),]
  
  # Outputting a Network List
    network_list <- vector('list', 2)
    names(network_list) <- c('nodelist', 'edgelist')
    network_list[[1]] <- node_list
    network_list[[2]] <- edge_list
  
  # Returning Network List
    return(network_list)
}

network_data <- village_subestting(test_cases, village_data)

###########################
#   LOOKING AT THE DATA   #
###########################

# Summing Network Types
  tie_summation <- function(village_id){
    # Isolating Edges 
      edges <- network_data$edgelist
      edges <- edges[(edges$villagecode == village_id),]
  
    # Simplifying Edgelist 
      edges <- edges[c("ind_id", "tietype", "hh_id", "alter_hh_id")]
  
    # Loop Over Individual IDs
      individual_ids <- sort(unique(edges$ind_id))
      individual_ties <- vector('list', length(individual_ids))
      names(individual_ties) <- individual_ids
      for(i in seq_along(individual_ids)){
        # Isolate Individuals
          ego <- edges[(edges$ind_id == individual_ids[[i]]),]
    
        # Isolating Household ID
          household_id <- unique(ego$hh_id)
    
        # Isolate Alters
          alter_ids <- sort(unique(ego$alter_hh_id))
          tie_sums <- vector('numeric', length(alter_ids))
          tie_types <- vector('character', length(alter_ids))
          for(j in seq_along(alter_ids)){
            tie_types[[j]] <- paste(  t(unique(ego[(ego$alter_hh_id == alter_ids[[j]]), 2]))  , collapse=',')
            tie_sums[[j]] <- nrow(unique(ego[(ego$alter_hh_id == alter_ids[[j]]),]))
          }
    
        # Create Summary DataFrame
          ego_data <- data.frame(ego= rep(individual_ids[[i]], length(alter_ids)), hh_id= rep(household_id, length(alter_ids)),
                               alter_hh_id=alter_ids, tie_type=tie_types, weight=tie_sums)      
    
        # Populate List
          individual_ties[[i]] <- ego_data
      }
  
    # Collapse List into DataFrame
      tie_summary <- do.call('rbind', individual_ties)
  
    # Return Data
      return(tie_summary)
  }

# Specifying Function
  network_extractor <- function(village_id, tie_type, directed=TRUE){
    # Isolating Edges 
      edges <- network_data$edgelist
      edges <- edges[(edges$villagecode == village_id),]
  
    # Isolating tie specific edgelist at the household level
      edges <- edges[(edges$tietype == tie_type), ]
  
    # Simplifying Edgelist 
      edges <- edges[c("hh_id", "alter_hh_id")]
  
    # Getting Nodelist
      nodelist <- network_data$nodelist
      nodelist <- nodelist[(nodelist$villagecode == village_id),]
  
    # Creating ::Network (StatNet) Network Object
      g <- network::network.initialize(nrow(nodelist), directed = as.logical(directed))
  
    # Creating Sequential Edges
      nodes <- nodelist[c('id')]
      colnames(nodes) <- c('label')
      nodes <- cbind(seq(1, nrow(nodes), 1), as.data.frame(nodes))
      colnames(nodes)[[1]] <- c('id')
  
    # Mapping Sequential Node IDs to Edgelist
      edges <- cbind(seq(1, nrow(edges), 1), as.data.frame(edges))
      colnames(edges)[[1]] <- c('Obs_ID')
  
      senders <- edges[c(1,2)]
      colnames(senders)[[2]] <- c('label')
      senders <- dplyr::left_join(senders, nodes, by=c('label'))
      colnames(senders)[[3]] <- c('sender_id')
  
      targets <- edges[c(1,3)]
      colnames(targets)[[2]] <- c('label')
      targets <- dplyr::left_join(targets, nodes, by=c('label'))
      colnames(targets)[[3]] <- c('target_id')
  
      edgelist <- dplyr::left_join(senders, targets, by=c('Obs_ID'))
      rm(senders, targets)
  
    # Generating network data object
      g <- network::network.initialize(nrow(nodes), directed = as.logical(directed))
      el <- edgelist[,c('sender_id','target_id')]
      el[,1] <- as.character(el[,1])
      el[,2] <- as.character(el[,2])
      el <- na.omit(el)
      g <- network::add.edges(g, el[,1], el[,2])
  
    # Return
      return(g)
  }
  
# Generating a Dyad Dataset for the Purpose of Plotting Correlations
  edgelist <-network_data$edgelist
  tie_labels <- c('borrowing', 'lending', 'getting health advice', 'giving health advice', 'social/free time')
  tie_comparer <- function(edgelist, tie_labels){
    # Create Unique Edgelist ID with which Merge
      combined_edgelist <- edgelist[,c("hh_id", "alter_hh_id")]
  
    # Isolating Unique Set of Dyads
      dyad_id <- unique(combined_edgelist)
      dyad_id <- cbind(seq(1,nrow(dyad_id),1) ,as.data.frame(dyad_id))
      colnames(dyad_id)[[1]] <- c('dyad_id')
      dyad_id <- dplyr::left_join(dyad_id, edgelist[,c('villagecode','hh_id', 'alter_hh_id')], by=c('hh_id', 'alter_hh_id'))
  
    # Isolating Each Network and Appending from Left to Right
      tie_types <- unique(edgelist$tietype)
      for (i in seq_along(tie_types)){
        # Isolate Tie Data
          tie_data <- edgelist[(edgelist$tietype == tie_types[[i]]),c('hh_id', 'alter_hh_id')]
    
        # Adding Weight
          tie_data$weight <- rep(1,nrow(tie_data))
    
        # Creating Unique Name
          colnames(tie_data)[[3]] <- tie_labels[[i]]
    
        # Join with dyad_id to Append 
          dyad_id <- dplyr::left_join(dyad_id, tie_data, by=c('hh_id', 'alter_hh_id'))
    
        # Replacing NAs with 0s
          dyad_id[[ncol(dyad_id)]][is.na(dyad_id[[ncol(dyad_id)]])] <- 0 
    }
  
    # Return Dyad Data
    return(dyad_id)
  }
  
# Generate Dyad Dataset Comparing Tie Weights per Tie Type
  selected_villages_dyads <- tie_comparer(edgelist, tie_labels)

# Looking at Correlations
  contingency_table_constructor <- function(tie_labels, selected_villages_dyads){
    # Create Base Table to Populate
      tie_table <- matrix(ncol=length(unique(selected_villages_dyads$villagecode)), nrow = length(tie_labels))
      colnames(tie_table) <- unique(selected_villages_dyads$villagecode)
      rownames(tie_table) <- tie_labels
  
    # Populating Combinations
      for (i in seq_along(colnames(tie_table))){
        # Isolate Row Data
          village_data <- selected_villages_dyads[(selected_villages_dyads$villagecode == as.numeric(colnames(tie_table)[[i]])),]
    
        # Looping Through Tie Types to Populate Sums
          for (j in seq_along(tie_labels)){
            # Isolate Vector
              vector_sum <- sum(village_data[,c(tie_labels[[j]])])/nrow(village_data)
      
            # Populate Table
              tie_table[rownames(tie_table)[[j]],colnames(tie_table)[[i]]] <- vector_sum
          }
      }
  
    # Constructing Edge Tables for Village
  
    # Isolate Village Data
      village_data <- selected_villages_dyads[(selected_villages_dyads$villagecode == as.numeric(colnames(tie_table)[[i]])),]
  
    # Creating Base Table
      edge_table <- matrix(ncol=length(tie_labels), nrow = length(tie_labels))
      rownames(edge_table) <- tie_labels
      colnames(edge_table) <- tie_labels
  
    # Creating All Pairs
      edge_pairs <- as.data.frame(t(combn(tie_labels, 2)))
      colnames(edge_pairs) <- c('i', 'j')
  
    # Loop Through Pairs
      for (j in seq_along(edge_pairs$i)){
        # Isolate Data
          pair_data <- village_data[,c(edge_pairs[j,1], edge_pairs[j,2])]
    
        # Calculate Matching Pairs
          matched_pairs <- pair_data[(pair_data[[1]] == pair_data[[2]]),]
    
        # Get the Row Count of Matched Pairs
          edge_table[edge_pairs[j,1], edge_pairs[j,2]] <- nrow(matched_pairs)
      }
  
  }
  
# Pareto Plot of Node Attributes to Determine Inclusion in ERGM
  pareto_plot <- function(df, attribute, village){
    # Compute Proportion for Each Attribute Category
      df_summary <- df %>%
        select(villagecode, attribute) %>%
        group_by(villagecode) %>%
        mutate(total = n()) %>%
        group_by(across(c(villagecode, attribute, total))) %>%
        summarise(count = n()) %>%
        mutate(proportion = count / total) %>%
        mutate(proportion = round(proportion, digits = 2)) %>%
        pivot_wider(id_cols = villagecode,
                    names_from = attribute,
                    values_from = proportion) %>%
        mutate(across(everything(), ~replace_na(., 0)))
      
      # Convert to Long Format Dataframe
      df_summary <- pivot_longer(df_summary, cols = c(2:ncol(df_summary)),
                                      names_to = attribute,
                                      values_to = "proportion")
      
      # Filter to Village for Plotting and Compute Cumulative Percentage
      df_summary <- df_summary %>%
        filter(villagecode == village) %>%
        arrange(desc(proportion)) %>%
        mutate(cumulative_percent = cumsum(proportion))
      
      # Make Pareto Plot
      if (max(df_summary$proportion) > 0.5) {
        plot <- ggplot(df_summary, aes(x = reorder(df_summary[[attribute]], -proportion), y = proportion)) +
          geom_bar(stat = "identity", fill = "skyblue") +
          geom_line(aes(y = cumulative_percent / max(cumulative_percent), group = 1),
                    color = "red", size = 1) +
          geom_point(aes(y = cumulative_percent / max(cumulative_percent)), color = "red", size = 2) +
          scale_y_continuous(name = "Proportion", limits = c(0, 1),
                             sec.axis = sec_axis(trans = ~ ., 
                                                 name = "Cumulative Percent", 
                                                 labels = scales::percent)) +
          labs(title = paste0("Pareto Plot for Village ", village), x = attribute, y = "Proportion of Individuals in Category") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      } else {
        plot <- ggplot(df_summary, aes(x = reorder(df_summary[[attribute]], -proportion), y = proportion)) +
        geom_bar(stat = "identity", fill = "skyblue") +
        geom_line(aes(y = cumulative_percent / max(cumulative_percent) * 0.5, group = 1),
                  color = "red", size = 1) +
        geom_point(aes(y = cumulative_percent / max(cumulative_percent) * 0.5), color = "red", size = 2) +
        scale_y_continuous(name = "Proportion", limits = c(0, 0.5),
                           sec.axis = sec_axis(trans = ~ . * 2, 
                                               name = "Cumulative Percent", 
                                               labels = scales::percent)) +
        labs(title = paste0("Pareto Plot for Village ", village), x = attribute, y = "Proportion of Individuals in Category") +
        theme_bw()
      }
      
      # Print Plot
      print(plot)
      
      rm(df, df_summary)
  }

#####################################
#   SPECIFYING WEIGHTING FUNCTION   #
#####################################

# Notes
# Simulations are specified at the village level
# Nodes represent individuals, but individuals can have household characteristics.

# Identify Relevant Edges to apply edge weights
# 1) borrowing
# 2) lending
# 3) getting health advice
# 4) giving health advice
# 5) social/free time

# Relation Events: Counts (Navigating uncertainty in networks of social exchange, pg. 11)
# Exchange Type Repetition (ij)
# Partner repetition (ij)
# Reciprocity (ij)
# Transitive referral (I hear X is a good person to get advice from) (ij)
# Exchange Activity (Everyone Borrows Money from X) (j)

# Map Node-Level Characteristics: Intervention Phase
# Female Headed Household
# Subsistence Farming
# Has Children

# Punch List
# Look at tie type correlations
# We will use to specify weights

##############################
# ERGM ESTIMATION FUNCTION   #
##############################

networks_list <- test_cases
iteration <- 3

# Specifying a Tie Estimation Function that Can Be Looped Across the Entire Set of Networks
  tie_estimator <- function(networks_list, iteration){
    # Generate Tie-Specific Networks
       free_time <- network_extractor(networks_list[[iteration]], 5)
       lending <- network_extractor(networks_list[[iteration]], 2)
       borrowing <- network_extractor(networks_list[[iteration]], 1)
       health_receiving <- network_extractor(networks_list[[iteration]], 3)
       health_giving <- network_extractor(networks_list[[iteration]], 4)
      
    # Assign Attributes
      health_giving %v% "indegree" <- degree(health_giving,cmode="indegree")
      health_giving %v% "outdegree" <- degree(health_giving,cmode="outdegree")
      attributes <- network_data$nodelist[(network_data$nodelist$villagecode == networks_list[[iteration]]), ]
      attributes <- attributes[!duplicated(attributes$id), ]
      attributes$house_edu <- as.numeric(attributes$house_edu)
      attributes$house_edu[is.na(attributes$house_edu)] <- 0
      attributes$hhid_dist <- as.numeric(attributes$id)
      attributes$num_kid[is.na(attributes$num_kid)] <- 0
      health_giving %v% "male_hh" <- attributes$oldest_male
      health_giving %v% "kids" <- attributes$num_kid
      health_giving %v% "house_caste" <- attributes$house_caste
      health_giving %v% "house_edu" <- attributes$house_edu
      health_giving %v% "hhid_dist" <- attributes$hhid_dist
      
    # Assess Whether to Include Caste in ERGM
      
     # if (any(prop.table(table(free_time %v% "house_caste")) > 0.6) == TRUE) {
      
        # If a Network has More than 60% of Individuals in One Caste, Exclude Caste
        md = max(degree(free_time, cmode = "outdegree"))
        m1 <- ergm(free_time ~ edges +
                     absdiff("house_edu") +
                     gwidegree(.05,fixed=T) +
                     gwnsp(.10,fixed=T),
                   constraints = ~bd(maxout = md),
                   control = c(control.ergm(main.method = "MCMLE",
                                            MCMC.burnin = 800000, MCMC.interval = 1050,
                                            MCMC.samplesize = 80000, MCMLE.maxit = 30)),
                   verbose = FALSE)
        
        # Extract Feature Set
        m1_feature_set <- ergmMPLE(free_time ~ edges +
                                     absdiff("house_edu") +
                                     gwidegree(.05,fixed=T) +
                                     gwnsp(.10,fixed=T),
                                   constraints = ~bd(maxout = md),
                                   control = c(control.ergm(main.method = "MCMLE",
                                                            MCMC.burnin = 800000, MCMC.interval = 1050,
                                                            MCMC.samplesize = 80000, MCMLE.maxit = 30)),
                                   verbose = FALSE)
      } else {
        
        # If a Network has Less than 60% of Individuals in One Caste, Exclude Caste
        md = max(degree(free_time, cmode = "outdegree"))
        m1 <- ergm(free_time ~ edges +
                     absdiff("house_edu") +
                     nodematch("house_caste", diff = T) +
                     gwidegree(.05,fixed=T) +
                     gwnsp(.10,fixed=T),
                   constraints = ~bd(maxout = md),
                   control = c(control.ergm(main.method = "MCMLE",
                                            MCMC.burnin = 800000, MCMC.interval = 1050,
                                            MCMC.samplesize = 80000, MCMLE.maxit = 30)),
                   verbose = FALSE)
        
        # Extract Feature Set
        m1_feature_set <- ergmMPLE(free_time ~ edges +
                                     absdiff("house_edu") +
                                     nodematch("house_caste", diff = T) +
                                     gwidegree(.05,fixed=T) +
                                     gwnsp(.10,fixed=T),
                                   constraints = ~bd(maxout = md),
                                   control = c(control.ergm(main.method = "MCMLE",
                                                            MCMC.burnin = 800000, MCMC.interval = 1050,
                                                            MCMC.samplesize = 80000, MCMLE.maxit = 30)),
                                   verbose = FALSE)
      }
  
      md = max(degree(health_giving, cmode = "outdegree"))
      m1 <- ergm(health_giving ~ edges +
                   mutual + 
                   absdiff("house_edu") +
                   absdiff("hhid_dist") +
                   gwesp(.25,fixed=T),
                 constraints = ~bd(maxout = md),
                 control = c(control.ergm(main.method = "MCMLE",
                                          MCMC.burnin = 800000, MCMC.interval = 1050,
                                          MCMC.samplesize = 80000, MCMLE.maxit = 30)),
                 verbose = TRUE)
      
      
      # Extract Feature Set
      m1_feature_set <- ergmMPLE(free_time ~ edges +
                                   absdiff("house_edu") +
                                   gwidegree(.05,fixed=T) +
                                   gwnsp(.10,fixed=T),
                                 constraints = ~bd(maxout = md),
                                 control = c(control.ergm(main.method = "MCMLE",
                                                          MCMC.burnin = 800000, MCMC.interval = 1050,
                                                          MCMC.samplesize = 80000, MCMLE.maxit = 30)),
                                 verbose = FALSE)
        
    # Isolating Model Table
      response_vector <- m1_feature_set[[1]]
      covariates <- as.data.frame(m1_feature_set[[2]])
      model_table <- cbind(response_vector, covariates)
      colnames(model_table)[[1]] <- c('response')
      rm(response_vector, covariates)
          
    # Isolating Model Components
      model_coefficients <- m1$coefficients
      coefficient_index <- data.frame(variable = names(model_coefficients), coefficient=model_coefficients)
      
    # Calculating Values
      predicted_values <- vector('numeric', nrow(model_table))
      for (j in seq_along(predicted_values)) {
        # Isolate Model Table Row
          row_values <- model_table[j,]
          row_index <- data.frame(variable=names(row_values), value=as.numeric(row_values))
          
        # Map Coefficient Values to Row Values & Multiply Values
          row_index = dplyr::left_join(coefficient_index, row_index, by=c('variable'))
          row_index$model_term = row_index$coefficient * row_index$value
          
        # Calculate Predicted Value & Populating Output Vector
          predicted_value <- sum(row_index$model_term)
          predicted_values[[j]] <- predicted_value
      }
      
    # Adding Predicted Values to Model Table
      model_table <- cbind(predicted_values, model_table)
      colnames(model_table)[[1]] <- c('model_prediction')
    
    # Return Output List
      output_list <- vector('list', length=2)
      names(output_list) <- c('ergm', 'model_table')
      output_list[[1]] <- m1
      output_list[[2]] <- model_table
      
    # Return List
      return(output_list)
  
# Looking at Model Outputs
  model_1_list <- tie_estimator(test_cases, 1)
  model_1_table <- model_1_list$model_table

# base::dir.create("/Users/tylerbarrett/Dropbox (Duke Bio_Ea)/Village_Simulations/villageData/village_networks_2023/ergm_models")

# Looping Over Village Networks
  villages_estimation <- function(test_cases){
    # Creating directory objects to preserve working directory
      current_directory <- getwd()
      model_directory <- paste0(getwd(),"/ergm_models")

    # Setting Working Directory
      setwd(model_directory)

    # Getting Start Time
      start_time <- Sys.time()

    # Looping Over Test Cases
      for (iteration in seq_along(test_cases)){
        message("Network ", seq(1, length(test_cases), 1)[[iteration]], "/", length(test_cases))
        
        model_list <- tie_estimator(test_cases, iteration)
  
        write.csv(model_list$model_table, file = paste0("village_model_", test_cases[[iteration]], ".csv"),
                  row.names = FALSE)
  
        rm(model_list)
        gc(reset = TRUE)
        base::Sys.sleep(3)
      }

    # Getting Stop Time
      stop_time <- Sys.time()
      total_time <- stop_time - start_time
      print(total_time)

    # Restoring Working Directory
      setwd(current_directory)
  }

  villages_estimation(test_cases)

# Plotting Examples
  x11(width=10.6806, height=7.30556)
  dev.off()
  
  pdf("Village83_Predictions.pdf", width=8.5, height=5,onefile=T)
    layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
    layout(mat = layout.matrix) 
  
    plot(density(model_1_table$model_prediction), family='HersheySerif', las=1, main='Network 83: Total Values', bty='n', col='brown', axes=TRUE)
    grid(lwd = 2)
  
    plot(density(model_1_table$model_prediction[model_1_table$model_prediction > 0]), family='HersheySerif', las=1, main='Network 83: Positive Values', bty='n', col='brown', axes=TRUE)
    grid(lwd = 2)
  dev.off()
  
# Writing Model Table File Out for Evaluation
  readr::write_csv(model_1_table, file="Village_83_ModelTable.csv")
  
#############
#   TESTS   #
#############

# Plotting Component Networks
  par(mar=c(0,0,0,0))
  plot(free_time)
  plot(lending)
  plot(borrowing)
  plot(health_receiving)
  plot(health_giving)


