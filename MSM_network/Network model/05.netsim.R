##
## Network simulation for ART-Net Data
## v1: 2018-08
##

## Packages ##
rm(list = ls())
suppressMessages(library("EpiModelHIV"))


## Inputs ##
city <- Sys.getenv("CITY")
if (city == "A") {
  city_name <- "Atlanta"
} else {
  city_name <- "San Francisco"
}


## Load Data ##
fn <- paste("data/artnet.NetEst", gsub(" ", "", city_name), "rda", sep = ".")
est <- readRDS(file = fn)

# Init network sim
nw <- list()
for (i in 1:3) {
  x <- est[[i]]
  nw[[i]] <- simulate(x$fit, basis = x$fit$newnetwork,
                      control = control.simulate.ergm(MCMC.burnin = 2e5))
}
## Dynamic network sim

sim_network <- function(est, nsteps) {
  

  
  # Dynamic time loop
  for (at in 1:nsteps) {
    # Main #
    deg_dist_casl <- as.numeric(summary(nw[[2]] ~ sociality(base = 0), at = at))
    nw[[1]] <- set.vertex.attribute(nw[[1]], attrname = "deg.casl", value = deg_dist_casl)
    nw[[1]] <- suppressWarnings(simulate(nw[[1]],
                                         formation = est[[1]]$formation,
                                         dissolution = est[[1]]$coef.diss$dissolution,
                                         coef.form = est[[1]]$coef.form,
                                         coef.diss = est[[1]]$coef.diss$coef.crude,
                                         time.start = at,
                                         time.slices = 1,
                                         time.offset = 0,
                                         monitor = "all",
                                         output = "networkDynamic"))
    
    deg_dist_main <- as.numeric(summary(nw[[1]] ~ sociality(base = 0), at = at))
    nw[[2]] <- set.vertex.attribute(nw[[2]], attrname = "deg.main", value = deg_dist_main)
    
    # Casual #
    nw[[2]] <- suppressWarnings(simulate(nw[[2]],
                                         formation = est[[2]]$formation,
                                         dissolution = est[[2]]$coef.diss$dissolution,
                                         coef.form = est[[2]]$coef.form,
                                         coef.diss = est[[2]]$coef.diss$coef.crude,
                                         time.start = at,
                                         time.slices = 1,
                                         time.offset = 0,
                                         monitor = "all",
                                         output = "networkDynamic"))
    
    deg_dist_main <- as.numeric(summary(nw[[1]] ~ sociality(base = 0), at = at))
    deg_dist_casl <- as.numeric(summary(nw[[2]] ~ sociality(base = 0), at = at))
    deg_dist_tot <- pmin(deg_dist_main + deg_dist_casl, 3)
    nw[[3]] <- set.vertex.attribute(nw[[3]], attrname = "deg.tot", value = deg_dist_tot)
    
    # One-Off #
    nw[[3]] <- suppressWarnings(simulate(nw[[3]],
                                         formation = est[[3]]$formation,
                                         dissolution = est[[3]]$coef.diss$dissolution,
                                         coef.form = est[[3]]$coef.form,
                                         coef.diss = est[[3]]$coef.diss$coef.crude,
                                         time.start = at,
                                         time.slices = 1,
                                         time.offset = 0,
                                         monitor = "all",
                                         output = "networkDynamic"))
    
    cat("\n Step ", at, "/", nsteps)
  }
  
  return(nw)
}
out <- sim_network(est, nsteps = 12)

fns <- strsplit(fn, "[.]")[[1]]
fn.new <- paste(fns[1], "NetSim", fns[3], "rda", sep = ".")
saveRDS(out, file = fn.new)


# create influence sample


# Initialize an empty data frame with the required columns
edge_df <- data.frame(node1 = integer(), node2 = integer(), starting = numeric(), ending = numeric())

for (n in 1:3) {
  num_edges <- network.edgecount(out[[n]])
  
  for (i in 1:num_edges) {
    # Get the nodes connected by the edge
    node1 <- out[[n]][["mel"]][[i]]$inl
    node2 <- out[[n]][["mel"]][[i]]$outl
    
    # Get the activation periods for the edge
    activity_periods <- get.edge.activity(out[[n]], e = i)
    
    # Add each activation period as a new row to the data frame
    edge_df <- rbind(edge_df, data.frame(
      node1 = node1,
      node2 = node2,
      starting = activity_periods[[1]][1, 1],
      ending = activity_periods[[1]][1, 2]
    ))
  }
}



edge_df$starting <- ifelse(edge_df$starting == -Inf, 0, edge_df$starting)
edge_df$ending <- ifelse(edge_df$ending == Inf, 12, edge_df$ending)


# Number of initial nodes to select
num_initial_nodes <- 100

# Randomly select 100 nodes from the network as initial nodes
initial_nodes <- sample(c(1:1000), num_initial_nodes)

# List to store the final influence samples for all initial nodes
final_influence_samples <- list()

# Iterate over each of the 100 initial nodes
for (i in 1:num_initial_nodes) {
  # Get the current initial node
  initial_node <- initial_nodes[i]
  
  # Initialize the influence sample list for this initial node
  influence_sample <- list()
  influence_sample[[1]] <- initial_node
  
  # Iterate through each time step to find influenced nodes
  for (t in 0:12) {
    # Get the nodes influenced in the previous time step
    current_nodes <- influence_sample[[length(influence_sample)]]
    
    # Find edges active at the current time step that involve any of the current nodes
    new_connections <- unique(c(
      edge_df$node1[edge_df$starting == t & edge_df$node2 %in% current_nodes],
      edge_df$node2[edge_df$starting == t & edge_df$node1 %in% current_nodes]
    ))
    
    # Combine with the current nodes to form the new influenced set
    new_nodes <- unique(c(current_nodes, new_connections))
    
    # Add the new set of influenced nodes to the influence sample list
    influence_sample[[length(influence_sample) + 1]] <- new_nodes
  }
  
  # Get the final influence sample as a unique set of nodes
  final_influence_samples[[i]] <- unique(unlist(influence_sample))
}





# Extend simulation to 120 weeks
out <- sim_network(est, nsteps = 120)

# Save the simulated network
fns <- strsplit(fn, "[.]")[[1]]
fn.new <- paste(fns[1], "NetSim", fns[3], "rda", sep = ".")
saveRDS(out, file = fn.new)

# Initialize an empty data frame with the required columns
edge_df <- data.frame(node1 = integer(), node2 = integer(), starting = numeric(), ending = numeric())

# Extract activation information from the simulation
for (n in 1:3) {
  num_edges <- network.edgecount(out[[n]])
  
  for (i in 1:num_edges) {
    # Get the nodes connected by the edge
    node1 <- out[[n]][["mel"]][[i]]$inl
    node2 <- out[[n]][["mel"]][[i]]$outl
    
    # Get the activation periods for the edge
    activity_periods <- get.edge.activity(out[[n]], e = i)
    
    # Add each activation period as a new row to the data frame
    edge_df <- rbind(edge_df, data.frame(
      node1 = node1,
      node2 = node2,
      starting = activity_periods[[1]][1, 1],
      ending = activity_periods[[1]][1, 2]
    ))
  }
}

# Handle infinite start and end times
edge_df$starting <- ifelse(edge_df$starting == -Inf, 0, edge_df$starting)
edge_df$ending <- ifelse(edge_df$ending == Inf, 120, edge_df$ending)

# Number of total final influence samples and nodes per interval
num_intervals <- 10  # 120 weeks / 12 weeks
nodes_per_interval <- 300
total_samples <- num_intervals * nodes_per_interval  # 3000 total samples
num_collections <- 50  # Number of influence sample collections
a=0
# List to store all influence sample collections
all_collections <- list()

# Iterate to create 50 influence sample collections
for (collection_idx in 1:num_collections) {
  # List to store the final influence samples for the current collection
  all_influence_samples <- list()
  
  # Iterate over each 12-week interval to construct influence samples
  for (interval in 1:num_intervals) {
    # Select 300 random initial nodes for this interval
    initial_nodes <- sample(c(1:1000), nodes_per_interval)
    
    for (i in 1:nodes_per_interval) {
      # Get the current initial node
      initial_node <- initial_nodes[i]
      
      # Initialize the influence sample list for this initial node
      influence_sample <- list() 
      influence_sample[[1]] <- initial_node
      
      # Iterate through each time step within the current 12-week interval
      start_time <- (interval - 1) * 12
      end_time <- interval * 12
      
      for (t in start_time:(end_time - 1)) {
        # Get the nodes influenced in the previous time step
        current_nodes <- influence_sample[[length(influence_sample)]]
        
        # Find edges active at the current time step that involve any of the current nodes
        new_connections <- unique(c(
          edge_df$node1[edge_df$starting == t & edge_df$node2 %in% current_nodes],
          edge_df$node2[edge_df$starting == t & edge_df$node1 %in% current_nodes]
        ))
        
        # Combine with the current nodes to form the new influenced set
        new_nodes <- unique(c(current_nodes, new_connections))
        
        # Add the new set of influenced nodes to the influence sample list
        influence_sample[[length(influence_sample) + 1]] <- new_nodes
      }
      
      # Get the final influence sample as a unique set of nodes
      all_influence_samples[[length(all_influence_samples) + 1]] <- unique(unlist(influence_sample))
    }
  }
  
  # Store the collection of influence samples
  all_collections[[collection_idx]] <- all_influence_samples
  a=a+1
  print(a)
}

# Load the required library for JSON handling
library(jsonlite)

# Save the `all_collections` as a JSON file
write_json(all_collections, path = "influence_samples_msm.json")
