library(jsonlite)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(scales) # to access break formatting functions
library(gridExtra)
library(grid)



# Figure 1
{
  # Active set method (solved to optimality)
  active_set_data = read.csv("figure_1_active_set.csv") %>%
    arrange(E,time_factor,iteration) %>%
    mutate(T = as.factor(T),E = as.factor(E)) %>%
    
    # Get the cumulative time 
    group_by(E,time_factor) %>%
    mutate(cum_time = cumsum(time)) %>% 

    # Extract the approximation gap
    group_by(time_factor,E) %>%
    mutate(min_obj_val = min(obj_val)) %>%
    mutate(approximation = obj_val / min_obj_val) %>% 
    
    # Identify the times that active set method gets within 10%, 1%,
    # and 0.1% of optimal
    group_by(time_factor,T,E,approximation,cum_time) %>%
    expand(approximation_level = c(1.1,1.01,1.001)) %>%
    filter(approximation < approximation_level + 1e-6) %>%  
    group_by(E,time_factor,T,approximation_level) %>%
    summarize(time = min(cum_time)) %>%
    ungroup() %>%
    mutate(gap = round(1000*(approximation_level-1))/ 1000) %>%
    
    # Add the name of the method and remove unnecessary columns
    mutate(method = "Active set method") %>%
    select(-approximation_level)


  # Barrier method
  barrier_data = read.csv("figure_1_barrier.csv") %>%
    select(E, time_factor,T,obj_val,y_nonzero,time,zeta_nonzero,Method, gap) %>%
    mutate(method = ifelse(Method == 0, "Primal simplex", ifelse(Method == 1, "Dual simplex",  ifelse(Method == 2, "Barrier", "Automatic"  )))) %>%
    select(-Method, -zeta_nonzero, -y_nonzero, -obj_val) %>%
    mutate(E = as.factor(E),
           T = as.factor(T))
  
  # Combine the data and format it
  data = rbind(active_set_data, barrier_data ) %>%
    mutate(gap = ifelse(abs(gap - 0.001) < 1e-7, "0.1%", ifelse(abs(gap - 0.01) < 1e-7, "1%", "10%")))
  
  plot_figure_1 = 
    data %>% 
    ggplot(aes(x=gap,y=time,fill=method)) +
    geom_col(position = "dodge") + 
    facet_grid(.~ T) + 
    theme_bw(base_size = 16) + 
    xlab("Approximation quality") + 
    ylab("Time (seconds)") +
    theme(          legend.title = element_blank(),
                    legend.position=c(0.5,.85),
                    legend.spacing.y = unit(0.0, "mm"), 
                    legend.key.size = unit(0.4, "cm"),
                    legend.text=element_text(size=11),
                    legend.background = element_blank(),
                    legend.box.background = element_rect(colour = "black"),
                    legend.spacing.x = unit(.1, 'cm'),
                    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )
  
  png(file="../output/figure_1.png",
      width = 12, height = 5, units = 'in', res = 300)
  grid.arrange(plot_figure_1)
  dev.off()
}





# Size of LPs (recorded from Gurobi output)
{
  # Experiments for (D) with E=5 factories and time_factor \in {2,4,6,8,10}
  # Num constraints does not include nonnegativity constraints
  df = data.frame(
    time_factor = c(2,4,6,8,10),
    num_constraints = c(63162,248418,555786,985266,1536858),
    num_variables = c(29100,113484,253164,448140,698412)
  )
}


# Markovian (not currently in paper)
{
  plot_markov = read.csv("markov.csv") %>%
    mutate(E = as.factor(E)) %>%
    filter(E %in% c(10,20)) %>%
    ggplot(aes(x=k,y=obj_val,color=E)) +
    geom_line() +
    scale_x_continuous(breaks=seq(1,10)) + 
    geom_point() + 
    theme_bw(base_size = 16) +
    theme(panel.grid.minor.x = element_blank()) + 
    xlab("L") + 
    ylab("Objective value") + 
    facet_grid(.~T)
  
  png(file="../output/plot_markov.png",
      width = 8, height = 5, units = 'in', res = 300)
  grid.arrange(plot_markov)
  dev.off()
}

