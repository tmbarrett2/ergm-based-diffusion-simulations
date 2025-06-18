#   Functions for the ERGM Diagnostic Pipeline
#   Tyler Barret & Jonathan H. Morgan, Ph.D.
#   23 April 2025

# Clear Out Console Script
  cat("\014")

# Options
  options(stringsAsFactors = FALSE)
  options(mc.cores = parallel::detectCores())

#################
#   Functions   #
#################
  
# EDGE CONCORDANCE CHECK
  
  # Function to Perform Edge Concordance Check
    ergm_concordance_check <- function(network_a, network_b) {
      
      # Convert Network Object to Matrix
        mat_a <- network::as.matrix.network.adjacency(network_a)
        mat_b <- network::as.matrix.network.adjacency(network_b)
      
      # Identify non-zero elements in both matrices 
        non_zero_a <- mat_a != 0 
        non_zero_b <- mat_b != 0 
      
      # Element-wise concordance where both mat_a and mat_b are non-zero 
        concordant_non_zero <- (mat_a == mat_b) & non_zero_a & non_zero_b 
      
      # Total number of non-zero comparisons 
        total_non_zero_comparisons <- sum(non_zero_a, non_zero_b) 
      
      # Number of concordant non-zero comparisons 
        num_concordant_non_zero <- sum(concordant_non_zero, na.rm = TRUE)
      
      # Calculate concordance rate for non-zero elements 
        concordance_rate_non_zero <- num_concordant_non_zero / total_non_zero_comparisons
  }
    
  # Function to Check Concordance on Multiple Networks
  # Meta Network List is a List of Lists, e.g., a List of Villages with Networks for Different Tie Types
    ergm_concordance_check_iterative <- function(meta_network_list) {
      
      # Village Output List
        network_list <- vector("list", length(meta_network_list))
        names(network_list) <- paste0("meta_network_", seq_along(meta_network_list))
      
        for(i in seq_along(network_list)){
        
          # Create Index of All Network Pairs
            comparison_index <- t(combn(names(meta_network_list[[i]]), m = 2))
            comparison_index <- as.data.frame(comparison_index)
          
          # Create a List of Comparison Networks Excluding Focal Network
            concordance_checks <- vector("numeric", nrow(comparison_index))
          
          # Loop Over Network and Alter List to Compute Concordance for Each Pair of Networks
            for(j in seq_along(comparison_index[,1])){
              
              i_element <- comparison_index[j,1]
              j_element <- comparison_index[j,2]
              
              concordance_checks[[j]] <- ergm_concordance_check(meta_network_list[[i]][i_element][[1]], meta_network_list[[i]][j_element][[1]])
            }
          
          # Return Results Data Frame
            comparison_index$concordance_value <- concordance_checks
            colnames(comparison_index)[c(1,2)] <- c("network_1", "network_2")
            
            network_list[[i]] <- comparison_index  
        }
      
      return(network_list)
    }

# DEGREE CORRELATION CHECKS

  # Function to Create a Dataframe with In and Outdegree for Each Network
    degree_comparer <- function(network_list){
      
        # Check if Network is Directed
          directed_check <- sapply(network_list, network::is.directed)
          
        # Create Empty Lists to Store Indegree and Outdegree 
          indegree_outputs <- vector("list", length(network_list))
          names(indegree_outputs) <- names(network_list)
            
          outdegree_outputs <- vector("list", length(network_list))
          names(outdegree_outputs) <- names(network_list)
            
          totaldegree_outputs <- vector("list", length(network_list))
          names(totaldegree_outputs) <- names(network_list)
            
          for (i in seq_along(directed_check)) {
            if (directed_check[[i]] == TRUE){
              # Compute Indegree and Outdegree for Each Network
                indegree_outputs[[i]] <- sna::degree(network_list[[i]], cmode = "indegree")
                  
                outdegree_outputs[[i]] <- sna::degree(network_list[[i]], cmode = "outdegree")
            }else{
              # Compute Toatal Degree
                totaldegree_outputs[[i]] <- sna::degree(network_list[[i]], cmode = "freeman")
            }
          }
      
        # Create Output Dataframe
          indegree_df <- data.frame(indegree_outputs)
          outdegree_df <- data.frame(outdegree_outputs)
          totaldegree_df <- data.frame(totaldegree_outputs)
          
          empty_check <- c(nrow(indegree_df), nrow(outdegree_df), nrow(totaldegree_df))
          if(empty_check[[1]] > 0){
            colnames(indegree_df) <- paste0(colnames(indegree_df), "_indegree")
          } 
          if(empty_check[[2]] > 0){
            colnames(outdegree_df) <- paste0(colnames(outdegree_df), "_outdegree")
          } 
          if(empty_check[[3]] > 0){
            colnames(totaldegree_df) <- paste0(colnames(totaldegree_df), "_totaldegree")
          }
          
          combined_df <- as.data.frame(matrix(ncol = 1, nrow = length(sna::degree(network_list[[1]], cmode = "freeman"))))
          df_list <- list(indegree_df = indegree_df,
                          outdegree_df = outdegree_df, 
                          totaldegree_df = totaldegree_df)
          df_list <- df_list[empty_check > 0]
          
          for(i in seq_along(df_list)){
            combined_df <- cbind(combined_df, df_list[[i]])
          }
          combined_df <- combined_df[-c(1)]
          
          return(combined_df)
    }
      
  # Function to Extract Correlated Pairs
    extract_correlated_pairs <- function(data, corr_threshold) {
          # Make Correlation Matrix
            cor_matrix <- cor(data)
          
          # Find indices of highly correlated pairs
            high_cor <- which(abs(cor_matrix) > corr_threshold, arr.ind = TRUE)
          
          # Filter out the diagonal (self-correlations) and duplicates
            pairs <- data.frame(row = rownames(cor_matrix)[high_cor[, 1]],
                                column = rownames(cor_matrix)[high_cor[, 2]],
                                cor_value = cor_matrix[high_cor])
            pairs <- subset(pairs, row != column)
            pairs <- pairs[!duplicated(t(apply(pairs[, 1:2], 1, sort))), ]
            
            return(pairs)
        }
      
  # Function to Perform ERGM Correlation Checks
    ergm_correlation_check <- function(data, corr_threshold) {
        
        # Extract Correlated Pairs
          correlated_pairs <- extract_correlated_pairs(data, corr_threshold)
        
        # Plot Correlations
          p_cor <- recordPlot(pairs.panels(data, pch=21, las=1))
        
          warning("Outputted pairs have high levels of multicolinearity.")
        
        # Create Output List
          output_list <- list(correlated_pairs = correlated_pairs,
                              correlation_plots = p_cor)
          
          return(output_list)
      }
      
  # Function to Check Degree Correlations on Multiple Networks
    ergm_correlation_check_iterative <- function(meta_network_list, corr_threshold) {
      
      # Village Output List
        network_list <- vector("list", length(meta_network_list))
        names(network_list) <- paste0("meta_network_", seq_along(meta_network_list))
      
        for(i in seq_along(network_list)){
          
          # Compute Degree for Each Network
          # To generalize, give list names, or request named network list
            data <- degree_comparer(meta_network_list[[i]])
          
          # Check Correlations - fix correlation checks function to return the plot object
            correlation_checks <- ergm_correlation_check(data, corr_threshold)
            
            network_list[[i]] <- correlation_checks
          }
        
      return(network_list)
    }
  
# ATTRIBUTE(S) CHECK
  
  # Function to Make Pareto Plot of Node Attributes to Determine Inclusion in ERGM
    pareto_plot <- function(data, plot_attribute){
      
      # Filter to Village for Plotting and Compute Cumulative Percentage
        df <- data %>%
          filter(attribute == plot_attribute) %>%
          arrange(desc(proportion)) %>%
          mutate(cumulative_percent = cumsum(proportion))
      
      # Make Pareto Plot
        if (max(df$proportion) > 0.5) {
          plot <- suppressWarnings(ggplot(df, aes(x = reorder(values, -proportion), y = proportion)) +
            geom_bar(stat = "identity", fill = "skyblue") +
            geom_line(aes(y = cumulative_percent / max(cumulative_percent), group = 1),
                      color = "red", size = 1) +
            geom_point(aes(y = cumulative_percent / max(cumulative_percent)), color = "red", size = 2) +
            scale_y_continuous(name = "Proportion", limits = c(0, 1),
                               sec.axis = sec_axis(trans = ~ ., 
                                                   name = "Cumulative Percent", 
                                                   labels = scales::percent)) +
            labs(x = "Attribute Category", y = "Proportion of Individuals in Category",
                 title = paste0(plot_attribute)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)))
        } else {
          plot <- suppressWarnings(ggplot(df, aes(x = reorder(values, -proportion), y = proportion)) +
            geom_bar(stat = "identity", fill = "skyblue") +
            geom_line(aes(y = cumulative_percent / max(cumulative_percent) * 0.5, group = 1),
                      color = "red", size = 1) +
            geom_point(aes(y = cumulative_percent / max(cumulative_percent) * 0.5), color = "red", size = 2) +
            scale_y_continuous(name = "Proportion", limits = c(0, 0.5),
                               sec.axis = sec_axis(trans = ~ . * 2, 
                                                   name = "Cumulative Percent", 
                                                   labels = scales::percent)) +
            labs(x = "Attribute Category", y = "Proportion of Individuals in Category",
                 title = paste0(plot_attribute)) +
            theme_bw())
        }
      
      # Return Plot
        return(plot)
        
        rm(df)
    }
    
  # Function to Perform ERGM Attribute Checks
    ergm_attributes_check <- function(predicted_network, attributes_list, attribute_min_threshold, attribute_max_threshold) {
      
      # Identify All Attributes Contained in the Network
        identified_attributes <- names(predicted_network$val[[1]])
      
      # Checking that the attributes identified in attributes_list are present
        found_attributes <- identified_attributes[(identified_attributes %in% attributes_list)]
        if(length(found_attributes) == length(attributes_list)){
          
          # Isolating Attribute Values of Found Attributes
            attribute_values <- vector('list', length(found_attributes))
            names(attribute_values) <- found_attributes
            for (i in seq_along(found_attributes)){
              attribute_values[[i]] <- network::get.vertex.attribute(predicted_network, found_attributes[[i]])
            }
          
          # Create data.frame of Attributes for Evaluation
              attributes_df <- data.frame(attribute_values)
            }else{
              print("Not All Attributes Found")
              break  
            }
      
      # Creating Proportions
        attributes_outputs <- vector('list', length(found_attributes))
        for (i in seq_along(found_attributes)){
          attribute_table <- table(attributes_df[[i]])
          values_df <- data.frame(values = names(attribute_table), 
                                  counts = as.integer(attribute_table))
          
          values_df$proportion <- values_df$counts/sum(values_df$counts)
          
          values_df <- cbind(rep(found_attributes[[i]], nrow(values_df)), values_df)
          colnames(values_df)[[1]] <- c('attribute')
          
          attributes_outputs[[i]] <- values_df
        }
        
      # Stacking to Make Rule Determinations
        attributes_outputs <- do.call('rbind', attributes_outputs)
      
      # Apply Threshold
        attribute_scores <- attributes_outputs[(attributes_outputs$proportion >= attribute_max_threshold | attributes_outputs$proportion < attribute_min_threshold), ]
        attribute_scores <- unique(attribute_scores[, c("attribute", "values", "proportion")])
      
      # Call Pareto Plot
        pareto_plots <- vector('list', length(found_attributes))
        names(pareto_plots) <- found_attributes
        for (i in seq_along(found_attributes)) {
          pareto_plots[[i]] <- pareto_plot(attributes_outputs, found_attributes[[i]])
        }
      
      # Combine Outputs in One List
        output_list <- list(attribute_scores = attribute_scores,
                            pareto_plots = pareto_plots)
      
      # Return Output List
        return(output_list)
    }
    
  # Function to Perform ERGM Attribute Checks on Multiple Networks
    ergm_attributes_check_iterative <- function(meta_network_list, attributes_list = NULL, min_threshold = 0.2, max_threshold = 0.8) {
      output_list <- vector("list", length(meta_network_list))
      names(output_list) <- names(meta_network_list)
      for(i in seq_along(meta_network_list)){
        if(is.null(attributes_list) == TRUE){
          attributes <- names(meta_network_list[[i]][[1]]$val[[1]])
          attributes <- attributes[attributes != "na"]
          attributes <- attributes[attributes != "vertex.names"]
        }else{
          attributes <- attributes_list
        }
        
        results <- vector("list", length(meta_network_list[[i]]))
        names(results) <- names(meta_network_list[[i]])
        
        for(j in seq_along(meta_network_list[[i]])){
          results[[j]] <- ergm_attributes_check(meta_network_list[[i]][[j]], attributes, min_threshold, max_threshold)
        }
        output_list[[i]] <- results
      }
      return(output_list)
    }
