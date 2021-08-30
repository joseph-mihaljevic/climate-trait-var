##################################################
## Project: Climate trait var 
## Script purpose: Figures from published data files
## Date: July 2021
## Author: JR Mihaljevic
##################################################
library(tidyverse)
library(grid)
library(gridExtra)
library(mgcv)
library(qgam)
library(viridis)

options(dplyr.width = Inf, dplyr.print_max = 1e9)
##################################################
##################################################

# Read in the data from Shocket et al. 2020 eLife
## Identify data files:
csv_list = system("ls ./Six-Viruses-Temp/Trait_Fits/*", intern = TRUE)
csv_list = csv_list[str_detect(csv_list, ".csv")]
n_csv = length(csv_list)

## Read in the files to data frames:
mosq_df_list = vector("list", n_csv)
for(i in 1:n_csv){
    mosq_df_list[[i]] = read_csv(csv_list[i])
}

## Need data from:
# A1F - mosq_df_list[[2]] - CHECKED AND NOT RELEVANT...
# A4A - mosq_df_list[[3]] - OK
# A9A - mosq_df_list[[6]] - OK
# A5A - mosq_df_list[[7]] - 

raw_data_list = vector("list", 4)

## FILTER EACH DF APPROPRIATELY:
### A1F - mosq_df_list[[2]] - NOT RELEVANT
# raw_data_list[[1]] = 
#     # Egg viability, Cx. theileri
#     mosq_df_list[[2]] %>%
#     select(Trait = trait.name, Temp = T, 
#            Host = host.code, Value = trait) %>%
#     filter(Host == "Cthe")

### A4A - mosq_df_list[[3]] 
raw_data_list[[1]] = 
    # MDR, Cx. pipiens
    mosq_df_list[[3]] %>%
    select(Trait = trait.name, Temp = T, 
           Host = host.code, Path = paras.code,
           Value = trait) %>%
    mutate(Value = 1/Value) %>% # to get MDR
    filter(Host == "Cpip")

### A9A - mosq_df_list[[6]] 
raw_data_list[[2]] = 
    # c, WNV, Cx. pipiens
    mosq_df_list[[6]] %>%
    select(Trait = trait.name, Temp = T, 
           Host = host.code, Path = paras.code,
           Value = trait) %>%
    filter(Trait == "c",
           Host == "Cpip",
           Path == "WNV")

### A5A - mosq_df_list[[7]]
raw_data_list[[3]] = 
    # lifespan, Cx. pipiens
    mosq_df_list[[7]] %>%
    select(Trait = trait.name, Temp = T, 
           Host = host.code, Path = paras.code,
           Value = trait) %>%
    #mutate(Value = 1/Value) %>% # to get mu
    filter(Host == "Cpip")

##################################################
##################################################
# Bring in the Tadpole data (Altman et al 2016)
tadpole_df = 
    read_csv("./Altman/Altman_tadpoles_data.csv") %>%
    select(Value = AdjPropEncysted, 
           Temp = PerfTemp)

raw_data_list[[4]] = tadpole_df
##################################################
##################################################
# Calculate the CV per temp:

cv_func = function(x){
    cv = sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

cv_data_list = vector("list", 4)

for(i in 1:4){
    
    temp_df = 
        raw_data_list[[i]] %>%
        group_by(Temp) %>%
        mutate(n_per_temp = n()) %>%
        filter(n_per_temp >= 3) %>%
        summarize(CV = cv_func(Value),
                  Range = max(Value) - min(Value))
    
    cv_data_list[[i]] = temp_df
    
    plot(temp_df$CV ~ temp_df$Temp)
    plot(temp_df$Range ~ temp_df$Temp)
    
}


3##################################################
##################################################
# Generate predictions

qseq = c(0.025,0.5,0.975)

pred_df_list = vector("list", length = length(raw_data_list))

for(i in 1:length(pred_df_list)){
    
    temp_df = raw_data_list[[i]]
    
    # Run the quantile regression gam
    n_k=4
    gam_1 = mqgam(list(Value~s(Temp, bs="cs", k=n_k), 
                       ~s(Temp, bs="cs", k=n_k)),
                  data = temp_df,
                  qu=qseq)
    # New df for predictions:
    new_df = data.frame(
        Temp = seq(min(temp_df$Temp),max(temp_df$Temp),length.out=100)
    )
    
    # Create predictions
    pred_df_list[[i]] = data.frame(
        med = qdo(gam_1, 0.5, predict, newdata = new_df),
        lo = qdo(gam_1, 0.025, predict, newdata = new_df),
        hi = qdo(gam_1, 0.975, predict, newdata = new_df),
        Temp = new_df$Temp
    )
    
}

##################################################
##################################################
# Plot the output

plot_fit_list = vector("list", length = 4)
plot_cv_list = vector("list", length = 4)
y_labs = c(
    "Development rate, MDR",
    "Infection efficiency, c",
    "Lifespan (days)",
    "Proportion encysted parasites"
)

for(i in 1:4){
    
    temp_df = pred_df_list[[i]]
    
    temp_df$lo[which(temp_df$lo < 0)] = 0
    
    if(i==2 | i == 4){
        temp_df$hi[which(temp_df$hi > 1.0)] = 1.0
    }
    
    plot_fit_list[[i]] =
        ggplot(temp_df) + 
        geom_ribbon(aes(x = Temp, ymin = lo, ymax = hi), 
                    fill = "#D3C7F4") +
        geom_path(aes(x = Temp, y = med), 
                  col = viridis(45, option = "G")[10]) +
        geom_point(data = raw_data_list[[i]],
                   aes(x = Temp, y = Value), alpha = 0.4) + 
        #scale_x_continuous(limits = c(10, 40)) +
        coord_cartesian(ylim=c(0,NA)) +
        labs(x="", y=y_labs[i]) +
        #scale_y_continuous(limits = c(0,NA)) +
        theme_classic() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_text(color="black",size=9),
            axis.text.x = element_text(color="black",size=10),
            axis.text.y = element_text(color="black",size=10,
                                       angle = 90, hjust = 0.5)
        )
    
    if(i < 4){
        plot_fit_list[[i]] = 
            plot_fit_list[[i]] + 
            scale_x_continuous(limits=c(15,35))
    }else{
        plot_fit_list[[i]] = 
            plot_fit_list[[i]] + 
            scale_x_continuous(limits=c(12,28))
    }
    
    plot_cv_list[[i]] =
        ggplot(cv_data_list[[i]]) + 
        geom_point(aes(x = Temp, y = Range), size = 2) + 
        geom_path(aes(x = Temp, y = Range)) + 
        #coord_cartesian(ylim=c(0,NA)) +
        labs(x="", y="Trait Value Range") +
        theme_classic() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_text(color="black",size=9),
            axis.text.x = element_text(color="black",size=10),
            axis.text.y = element_text(color="black",size=10,
                                       angle = 90, hjust = 0.5)
        )
    
    if(i < 4){
        plot_cv_list[[i]] = 
            plot_cv_list[[i]] + 
            scale_x_continuous(limits=c(15,35))
    }else{
        plot_cv_list[[i]] = 
            plot_cv_list[[i]] + 
            scale_x_continuous(limits=c(12,28))
    }
    
}

##################################################
##################################################
# Combine the plots:


fit_plots = arrangeGrob(
    plot_fit_list[[3]], plot_fit_list[[2]], plot_fit_list[[4]],
    ncol=1,
    bottom=textGrob(expression("Temperature ("*degree*"C)"),
                    y=unit(1,"lines"))
    )

cv_plots = arrangeGrob(
    plot_cv_list[[3]], plot_cv_list[[2]], plot_cv_list[[4]],
    ncol=1,
    bottom=textGrob(expression("Temperature ("*degree*"C)"),
                    y=unit(1,"lines"))
)

plot_combo = arrangeGrob(
    fit_plots, cv_plots, ncol = 2
)

quartz(height=6, width = 4.5)
grid.arrange(plot_combo)

ggsave(
    plot = plot_combo,
    file = "trait_data.png",
    device = "png",
    height = 6,
    width = 4.5
)
