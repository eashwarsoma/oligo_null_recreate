library(ggplot2)
library(reshape)
library(scales)
library(dplyr)

#####Figure 1####
#Function for base model
base_model <- function (double_time, t_init, t_fin) {
  r <- log(2)/double_time
  t_vals <- t_init:t_fin
  N_vals <- exp(r*t_vals)
  base.tab <- cbind(t_vals, N_vals)
  return(base.tab)
}

base.tab <- base_model(double_time = 150, 
                      t_init = 0, 
                      t_fin = 1300)


plot(base.tab)

#Function for perturbed model
pert_model <- function (double_time, int.frac, cell_kill, base.tab) {
  
  r <- log(2)/double_time
  
  #Gets the index of the base model to find intervention
  ind <- ceiling(dim(base.tab)[1] * int.frac)
  
  #Gets the cell number at the intervention
  N_start <- base.tab[ind, 2]
  
  #Gets the time against which to plot the cell growth in perturbed model
  t_plot <- base.tab[ind:dim(base.tab)[1], 1]
  
  #Gets the end points of the t_plot
  min <- range(base.tab[ind:dim(base.tab)[1], 1])[1]
  max <- range(base.tab[ind:dim(base.tab)[1], 1])[2]
  
  #Gets the actual time duration to grow this model
  t_grow <- 0:(max-min)
  
  #Finds the N values for the growth duration time
  N_vals <- (N_start-cell_kill) * exp(r*t_grow)
  
  #Combines this with the actual t plot values for base model
  pert.tab <- cbind(t_plot, N_vals)
  return(pert.tab)
  
}

pert.375 <- pert_model(150, .375, 5, base.tab)
pert.5 <- pert_model(150, .5, 5, base.tab)
pert.625 <- pert_model(150, .625, 5, base.tab)


colnames(base.tab) <- c("t_vals", "base_mod")
colnames(pert.375) <- c("t_vals", "pert.375")
colnames(pert.5) <- c("t_vals", "pert.5")
colnames(pert.625) <- c("t_vals", "pert.625")



#combining the data
full.dat <- merge(as.data.frame(base.tab), 
                  as.data.frame(pert.375), by = "t_vals", all = TRUE)
full.dat <- merge(full.dat, 
                  as.data.frame(pert.5), by = "t_vals", all = TRUE)
full.dat <- merge(full.dat, 
                  as.data.frame(pert.625), by = "t_vals", all = TRUE)

#metling the data
full.dat.melt <- melt(full.dat, id = "t_vals")

ggplot(data = full.dat.melt, 
       aes(x = t_vals, y = value, 
           color = variable)) + geom_point() + 
  labs(title="Effect of Treatment (killing 5 cells) at Different Points of Treatment",
        x ="Time Steps", y = "Number of Cells", color = "Treatment History") +
  scale_x_continuous(limits = c(400, 1275)) +
  scale_y_continuous(limits = c(0, 150)) +
  scale_color_manual(labels = c("no treatment", "treat at 3/8 of course", 
                                "treat at 4/8 of course", "treat at 5/8 of course"),
                     values = c("black", "blue", "green", "red")) +
  geom_hline(aes(yintercept = 100))
  

#####FIgure 2 looking at the isocline lines####
#Function to create a table showing delta t base
param_vis <- function(Nd_min, Nd_max, Nc_min, Nc_max, double_time, step) {
  #getting growth rate
  r <- log(2)/double_time
  
  #Doing things in log for scale sake
  #Getting range of Nd values
  Nd_range <- seq(log(Nd_min), log(Nd_max), step)
  
  #Getting range of Nc values
  Nc_range <- seq(log(Nc_min), log(Nc_max), step)
  
  #Creating a table of Nd and Nc value combinations
  N_tab <- expand.grid(Nd_range, Nc_range)
  
  #Getting the indexes at which Nc exceeds Nd
  inds <- which(N_tab[,2] > N_tab[,1])
  
  #Removing those values
  N_tab <- N_tab[-inds, ]
  
  #Now actually calcuating the associated survival diff
  t_diff <- log(exp(N_tab[,1])/(exp(N_tab[,1])-exp(N_tab[,2]))) * (1/r)
  
  #Combing everything into one big table
  fin.tab <- cbind(r, exp(N_tab[,1]), exp(N_tab[,2]), t_diff)
  
  return(as.data.frame(fin.tab))
  
  
}

doub.100 <- param_vis(10^9, 10^13, 10^9, 10^13, 100, .01)
colnames(doub.100) <- c("r", "Nd", "Nc", "td")
doub.200 <- param_vis(10^9, 10^13, 10^9, 10^13, 200, .01)
colnames(doub.200) <- c("r", "Nd", "Nc", "td")
doub.300 <- param_vis(10^9, 10^13, 10^9, 10^13, 300, .01)
colnames(doub.300) <- c("r", "Nd", "Nc", "td")
doub.400 <- param_vis(10^9, 10^13, 10^9, 10^13, 400, .01)
colnames(doub.400) <- c("r", "Nd", "Nc", "td")

all.data.grad <- rbind(doub.100, doub.200, doub.300, doub.400)

#Converting to factor
all.data.grad$r <- as.factor(all.data.grad$r)
levels(all.data.grad$r) <- c("400 days", "300 days", "200 days", "100 days")

plot <- ggplot(data = all.data.grad, aes(x = Nd, y = Nc)) + 
  facet_wrap(~ r, ncol = 2, scales = "free") +
  geom_tile(aes(fill = td)) +
  scale_x_continuous(trans='log10', limits = c(1e9, 1e13)) +
  scale_y_continuous(trans='log10', limits = c(1e9, 1e13)) +
  scale_fill_gradientn(colours = heat.colors(3), trans = "log10", 
                      limits=c(0.02, 1600), breaks = c(.1, 1, 10, 100, 1000), 
                      oob=squish)

#slope for isoclines
iso.cl <- function(double_time, iso_time, Nd_min, Nd_max, step) {
  r <- log(2)/double_time
  t <- iso_time
  m <- exp(r*t)
  slop <- -1/m
  
  Nd_ran <- seq(log(Nd_min), log(Nd_max), step)
  
  Nc_ran <- exp(Nd_ran) + slop*exp(Nd_ran)
  
  iso.lin <- cbind(exp(Nd_ran), Nc_ran)
  
}

line.400 <- cbind(as.data.frame(iso.cl(400, 100, 10^9, 10^13, .01)), "400 days")
line.300 <- cbind(as.data.frame(iso.cl(300, 100, 10^9, 10^13, .01)), "300 days")
line.200 <- cbind(as.data.frame(iso.cl(200, 100, 10^9, 10^13, .01)), "200 days")
line.100 <- cbind(as.data.frame(iso.cl(100, 100, 10^9, 10^13, .01)), "100 days")
colnames(line.400) <- c("x", "y", "r")
colnames(line.300) <- c("x", "y", "r")
colnames(line.200) <- c("x", "y", "r")
colnames(line.100) <- c("x", "y", "r")

lins <- rbind(line.400, line.300, line.200, line.100)

plot + geom_line(data = lins, aes(x = x, y = y)) + 
  labs(title="Possible Parameter Values for Trials",
       x ="Number of Cells at Presentation", y = "Number of Cells Killed", 
       fill = "Survival Advantage") 



