library(jsonlite)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(scales) # to access break formatting functions
library(gridExtra)
library(grid)

# Changing E (versus time)
{
  new_data = read.csv("figure_2.csv") %>%
    arrange(E,time_factor,removal,iteration) %>%
    group_by(E,time_factor,removal) %>%
    mutate(cum_time = cumsum(time)) %>% 
    mutate(time = cum_time) %>%
    mutate(obj_val = round(obj_val)) %>%
    mutate(T = as.factor(T),E = as.factor(E)) %>%
    filter(E %in% c(10,20,30,40,50))
  
  plot_fixed_T = new_data %>%
    ggplot(aes(x=time,y=obj_val)) +
    geom_step() + 
    theme_bw(base_size = 16)+
    facet_grid(. ~ E) +
    xlab("Time (seconds)") + 
    ylab("Objective value") +
    theme(          legend.title = element_blank(),
                    legend.position=c(0.5,.85),
                    legend.spacing.y = unit(0.0, "mm"), 
                    legend.key.size = unit(0.4, "cm"),
                    legend.text=element_text(size=11),
                    legend.background = element_blank(),
                    legend.box.background = element_rect(colour = "black"),
                    legend.spacing.x = unit(.1, 'cm'),
    )
  png(file="../output/figure_2.png",
      width = 12, height = 5, units = 'in', res = 300)
  grid.arrange(plot_fixed_T)
  dev.off()
  
}
