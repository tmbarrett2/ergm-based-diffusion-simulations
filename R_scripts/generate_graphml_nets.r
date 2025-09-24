# Julia/R Integration Tests
# Tyler Barrett & Jonathan H. Morgan, Ph.D.
# 24 September 2025 

################
#   PACKAGES   #
################

library(igraph)

########################################
#   IMPORT DATA & WRITE graphml FILES  #
########################################

# Import Network List
  setwd("/workspace/data/sim_test_data/")
  synthetic_nets_village_m <- readRDS("synthetic_nets_village_m.rds")

# Write graphml files
  setwd("/workspace/data/sim_test_data/synthetic_village_m_graphml")
  for(i in 1:length(synthetic_nets_village_m)){
    igraph::write_graph(synthetic_nets_village_m[[i]], file = paste0("net_", i, ".graphml"), format = "graphml")
  }

#############################
#   JULIA EXECUTABLE TETS   #
#############################


#############
#   TESTS   #
#############

# Test
  test_net <- igraph::read_graph("WeakCore1_3.2_2.net", format = c("pajek"))
  test_net

  setwd("/workspace/data/sim_test_data/synthetic_village_m_graphml")
  igraph::write_graph(test_net, file = "WeakCore1_3.2_2.graphml", format = "graphml")
