
library(jsonlite)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(scales) # to access break formatting functions
library(gridExtra)
library(grid)


# Gurobi solvers gap (Figure appendix)
{
  
  
  data_activeset_haihao = read.csv("activeset_dec_3_2024.csv")
  data_dual_haihao = read.csv("dual_dec_3_2024.csv")
  #dual_data = read.csv("dual_E5_012.csv") %>%
  #  filter(Method != 2) %>%
  #  rbind(read.csv("dual_E5.csv") %>% mutate(Method = 2)) %>%
  #  filter(time_factor %in% c(2,4,6,8,10))
  
  #data = dual_data %>%
  data = data_dual_haihao %>%
    select(E, time_factor,T,obj_val,y_nonzero,time,zeta_nonzero,Method) %>%
    mutate(method = ifelse(Method == 0, "Primal simplex", ifelse(Method == 1, "Dual simplex",  ifelse(Method == 2, "Barrier",  ifelse(Method == 6, "PDLP", "Automatic"  ))))) %>%
    select(-Method) %>%
    mutate(E = as.factor(E),
           T = as.factor(T)) %>%
    group_by(T,time_factor,E,method) %>%
    mutate(min_obj_val = min(obj_val)) %>%
    mutate(approximation = obj_val / min_obj_val)
  
  
  
  data_bar =  data %>%      
    filter(E == 5, method != "New Algorithm - No Removal") %>%
    group_by(E,time_factor,T,obj_val,y_nonzero,time,zeta_nonzero,method,min_obj_val,approximation) %>% 
    expand(approximation_level = c(1.1,1.01,1.001)) %>%
    filter(approximation < approximation_level + 1e-6) %>%  
    group_by(E,time_factor,T,method,approximation_level) %>%
    summarize(time = min(time)) %>%
    ungroup() %>% 
    filter(method == "New Algorithm - Removal" | approximation_level == 1.01) %>% 
    mutate(approximation_level = ifelse(method == "New Algorithm - Removal", approximation_level, 1.00)) %>%
    mutate(approximation_level = as.factor(approximation_level)) %>%
    filter(time_factor %in% c(2,3,4,6,8,10))
  
  
  
  plot_gurobi = 
    data %>% 
    mutate(x_val = method) %>%
    mutate(time = min(time,1000)) %>%
    mutate(reached_limit = time == 1000) %>%
    ggplot(aes(x=x_val,y=time,fill=method)) +
    geom_col() + 
    facet_grid(.  ~ T ) + 
    theme_bw(base_size = 16) + 
    #scale_y_continuous(limits=c(0,1000)) +
    xlab("") +
    ylab("Time (seconds)") +
    theme(          legend.title = element_blank(),
                    legend.position=c(0.1,.85),
                    legend.spacing.y = unit(0.0, "mm"), 
                    legend.key.size = unit(0.4, "cm"),
                    legend.text=element_text(size=11),
                    legend.background = element_blank(),
                    legend.box.background = element_rect(colour = "black"),
                    legend.spacing.x = unit(.1, 'cm'),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  png(file="../output/plot_gurobi.png",
      width = 12, height = 5, units = 'in', res = 300)
  grid.arrange(plot_gurobi)
  dev.off()
}
