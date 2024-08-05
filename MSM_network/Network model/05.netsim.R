
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


## Dynamic network sim

sim_network <- function(est, nsteps = 52*5) {

  # Init network sim
  nw <- list()
  for (i in 1:3) {
    x <- est[[i]]
    nw[[i]] <- simulate(x$fit, basis = x$fit$newnetwork,
                        control = control.simulate.ergm(MCMC.burnin = 2e5))
  }

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
out <- sim_network(est, nsteps = 260)

fns <- strsplit(fn, "[.]")[[1]]
fn.new <- paste(fns[1], "NetSim", fns[3], "rda", sep = ".")
saveRDS(out, file = fn.new)


library('network')
library(igraph)
library("networkDynamic")
library("ndtv")

dynBeach<-out[[1]]

vertex_size_factor <- 0.05  


plot(dynBeach, 
     vertex.cex = vertex_size_factor,  
     edge.col = "gray")   



static_network <- as.network(dynBeach, matrix.type = "adjacency")
adj_matrix <- as.matrix.network(static_network)
plot(static_network, 
     vertex.cex = vertex_size_factor, 
     edge.col = "gray")
# Save the adjacency matrix to a CSV file
write.csv(adj_matrix, file = "MSM_network.csv", row.names = TRUE)


compute.animation(dynBeach, animation.mode = "kamadakawai",
                  slice.par=list(start=0, end=260, interval=1, 
                                 aggregate.dur=1, rule='any'))

render.d3movie(dynBeach, usearrows = F, 
               displaylabels = F, label=dynBeach %v% "media",
               bg="#ffffff", vertex.border="#333333",
               edge.col = '#55555599',
               vertex.tooltip = paste("<b>Name:</b>", (dynBeach %v% "media") , "<br>",
                                      "<b>Type:</b>", (dynBeach %v% "type.label")),
               edge.tooltip = paste("<b>Edge type:</b>", (dynBeach %e% "type"), "<br>", 
                                    "<b>Edge weight:</b>", (dynBeach %e% "weight" ) ),
               launchBrowser=T, filename="Media-Network-Dynamic.html",
               render.par=list(tween.frames = 30, show.time = F),
               plot.par=list(mar=c(0,0,0,0)))
