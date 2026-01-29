# Script to Reduce Test Networks to Largest Component
# Created by Tyler M. Barrett, Ph.D on December 19, 2025

# Load Packages
  library(tidyverse)
  library(igraph)
  library(ggraph)

# Set Input and Output Directories
  input_dir <- "C:/Users/tmbar/OneDrive - Duke University/Active Projects/diffusion_sim/SocFeedback_Simulator_Windows/test_data/synthetic_village_m_graphml"
  output_dir <- "C:/Users/tmbar/OneDrive - Duke University/Active Projects/diffusion_sim/SocFeedback_Simulator_Windows/test_data/synthetic_village_m_graphml_largest_component"
  
# Get a List of All graphml Files
  graphml_files <- list.files(input_dir, pattern = "\\.graphml$", full.names = TRUE)

# Set Plotting Parameters
  n_plots <- 5  # Number of random networks to plot for checking
  set.seed(123)  # For reproducible random sampling
  
# Process Each File
  process_results <- map_dfr(graphml_files, function(file_path) {
    
    # Extract filename for progress tracking
    filename <- basename(file_path)
    
    # Read the graph (directed)
    g <- read_graph(file_path, format = "graphml")
    
    # Get the largest weakly connected component
    # mode = "weak" treats directed edges as undirected for finding components
    components <- components(g, mode = "weak")
    largest_component <- which.max(components$csize)
    
    # Extract vertices in largest component
    vertices_in_largest <- which(components$membership == largest_component)
    g_largest <- induced_subgraph(g, vertices_in_largest)
    
    # Create output path
    output_path <- file.path(output_dir, filename)
    
    # Write the largest component graph
    write_graph(g_largest, output_path, format = "graphml")
    
    # Return summary statistics
    tibble(
      file = filename,
      original_nodes = vcount(g),
      original_edges = ecount(g),
      largest_comp_nodes = vcount(g_largest),
      largest_comp_edges = ecount(g_largest),
      proportion_nodes = vcount(g_largest) / vcount(g)
    )
  }, .progress = TRUE)
  
  # Display summary
  cat("\n=== Processing Summary ===\n")
  cat("Total files processed:", nrow(process_results), "\n\n")
  
  # Show statistics
  cat("Network size statistics:\n")
  print(summary(process_results %>% select(original_nodes:proportion_nodes)))
  
  # Save detailed results
  write_csv(process_results, file.path(output_dir, "processing_log.csv"))
  cat("\nDetailed log saved to:", file.path(output_dir, "processing_log.csv"), "\n")
  
  # Plot random sample of networks for visual checking
  cat("\n=== Creating plots for visual checking ===\n")
  
  # Create plots directory
  plots_output_dir <- file.path(output_dir, "network_plots")
  if (!dir.exists(plots_output_dir)) {
    dir.create(plots_output_dir, recursive = TRUE)
  }
  
  # Get all filenames
  all_files <- process_results$file
  
  # Sample files to plot (or use all if fewer than n_plots)
  files_to_plot <- if (length(all_files) <= n_plots) {
    all_files
  } else {
    sample(all_files, n_plots)
  }
  
  # Create plots
  walk(files_to_plot, function(filename) {
    # Read the reduced network
    g <- read_graph(file.path(output_dir, filename), format = "graphml")
    
    # Create plot filename
    plot_filename <- str_replace(filename, "\\.graphml$", ".png")
    plot_path <- file.path(plots_output_dir, plot_filename)
    
    # Set layout based on network size
    if (vcount(g) < 100) {
      layout_type <- "fr"
    } else if (vcount(g) < 500) {
      layout_type <- "drl"
    } else {
      layout_type <- "lgl"
    }
    
    # Create ggraph plot
    p <- ggraph(g, layout = layout_type) +
      geom_edge_link(
        arrow = arrow(length = unit(2, 'mm'), type = "closed"),
        end_cap = circle(2, 'mm'),
        alpha = 0.3,
        color = "gray50",
        width = 0.3
      ) +
      geom_node_point(
        color = "navy",
        fill = "skyblue",
        shape = 21,
        size = 2,
        stroke = 0.5
      ) +
      labs(
        title = filename,
        subtitle = paste0(vcount(g), " nodes, ", ecount(g), " edges")
      ) +
      theme_graph(
        base_family = "sans",
        background = "white"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10)
      )
    
    # Save plot
    ggsave(
      filename = plot_path,
      plot = p,
      width = 8,
      height = 8,
      dpi = 100,
      bg = "white"
    )
    
    cat("  Plotted:", filename, "\n")
  })
  
  cat("\nPlots saved to:", plots_output_dir, "\n")
  cat("Number of plots created:", length(files_to_plot), "\n")