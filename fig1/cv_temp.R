##################################################
## Project: Climate trait var 
## Script purpose: concept figure generation 
## Date: June 2021
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

# HELPER FUNCS:

# Generate data from the mean and CV
gamma_mean = function(mu, sigma){
    # gamma mu = shape/rate
    # gamma sigma^2 = shape/(rate^2)
    # rate = mu/(sigma^2)
    # shape = rate * mu
    
    rate = mu / (sigma^2)
    shape = rate * mu
    return(list(shape = shape, rate = rate))
}

# Draw data from a gamma with mean and cv
data_func = function(x, n_samp){
    
    sigma_temp = x[1]*x[2]
    gp =  gamma_mean(x[1],sigma_temp)
    
    data_temp = rgamma(n_samp, shape=gp$shape, rate=gp$rate)
    
    return(data_temp)
    
}

##################################################
##################################################
# Generate different functional relationships
# Trait var (CV) vs. temp
# CV = sigma / mu
# Sigma = CV*mu

##################################################
##################################################

## 1. U-shaped CV: 
##    f(x) = a(x - b)^2 + c
##    (b,c) = vertex

## MEAN:
# Gaussian function:
# f(x) = a exp(- (x-mu)^2 / (2sigma^2) ) #a = height

# Generate a mean relationship:
peak = 4
mu = 20
sigma = 6

temp_seq= seq(10,30,by=2)
mean_seq = peak * exp(- (temp_seq-mu)^2 / (2*sigma^2) )

plot(mean_seq~temp_seq, ylim=c(0,peak))

# CV

a_u = 0.003
b_u = mu
c_u = 0.002

cv_seq = a_u*(temp_seq - b_u)^2 + c_u

plot(cv_seq~temp_seq)

##################################################
##################################################
# Generate Data from the mean and CV
n_samp = 40

duh = apply(cbind(mean_seq,cv_seq), 1, 
            data_func, n_samp = n_samp, simplify = FALSE)

data_df = data.frame(
    trait = unlist(duh),
    temp = rep(temp_seq, each = n_samp)
)

data_cv_df = data.frame(
    cv = cv_seq,
    temp = temp_seq
)
# Generate predictions
qseq = c(0.025,0.5,0.975)
gam_1 = mqgam(list(trait~s(temp, bs="cs"), ~s(temp, bs="cs")),
              data = data_df,
              qu=qseq)

# Create predictions:
new_df = data.frame(
    temp = seq(10,30,length.out=100)
)

pred_df = data.frame(
    med = qdo(gam_1, 0.5, predict, newdata = new_df),
    lo = qdo(gam_1, 0.025, predict, newdata = new_df),
    hi = qdo(gam_1, 0.975, predict, newdata = new_df),
    temp = seq(10,30,length.out=100)
)

u_plot_fit =
    ggplot(pred_df) + 
    geom_ribbon(aes(x = temp, ymin = lo, ymax = hi), 
                fill = "#AAC9D2") +
    geom_path(aes(x = temp, y = med), 
              col = viridis(45, option = "G")[20]) +
    geom_point(data = data_df,
               aes(x = temp, y = trait), alpha = 0.4) + 
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

u_plot_cv = 
    ggplot(data_cv_df) +
    geom_point(aes(x=temp,y=cv), size = 1.5) +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

##################################################
##################################################

## 2. Gaussian-shaped: 
##    f(x) = a exp(- (x-mu)^2 / (2sigma^2) ) 
##    a = height

a_g = 0.2
mu_g = mu
sigma_g = 6

cv_seq = a_g*exp(-(temp_seq - mu_g)^2 / (2*sigma_g^2))

plot(cv_seq~temp_seq)

##################################################
##################################################
# Generate Data from the mean and CV
n_samp = 40

duh = apply(cbind(mean_seq,cv_seq), 1, 
            data_func, n_samp = n_samp, simplify = FALSE)

data_df = data.frame(
    trait = unlist(duh),
    temp = rep(temp_seq, each = n_samp)
)

data_cv_df = data.frame(
    cv = cv_seq,
    temp = temp_seq
)

# Generate predictions
qseq = c(0.025,0.5,0.975)
gam_1 = mqgam(list(trait~s(temp, bs="cs"), ~s(temp, bs="cs")),
              data = data_df,
              qu=qseq)

# Create predictions:
new_df = data.frame(
    temp = seq(10,30,length.out=100)
)

pred_df = data.frame(
    med = qdo(gam_1, 0.5, predict, newdata = new_df),
    lo = qdo(gam_1, 0.025, predict, newdata = new_df),
    hi = qdo(gam_1, 0.975, predict, newdata = new_df),
    temp = seq(10,30,length.out=100)
)

g_plot_fit =
    ggplot(pred_df) + 
    geom_ribbon(aes(x = temp, ymin = lo, ymax = hi), 
                fill = "#AAC9D2") +
    geom_path(aes(x = temp, y = med), 
              col = viridis(45, option = "G")[20]) +
    geom_point(data = data_df,
               aes(x = temp, y = trait), alpha = 0.4) + 
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

g_plot_cv = 
    ggplot(data_cv_df) +
    geom_point(aes(x=temp,y=cv), size = 1.5) +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

##################################################
##################################################

## 3. Increasing CV: 
##    f(x) = a*exp(-(b*x))

# Redo the mean:
mu = 14
sigma = 9
mean_seq = peak * exp(- (temp_seq-mu)^2 / (2*sigma^2) )
plot(mean_seq~temp_seq)

# Decreasing CV
a_i = 0.02
b_i = -0.18

cv_seq = a_i * exp(- b_i * (temp_seq-10))

plot(cv_seq~temp_seq)

##################################################
##################################################
# Generate Data from the mean and CV
n_samp = 40

duh = apply(cbind(mean_seq,cv_seq), 1, 
            data_func, n_samp = n_samp, simplify = FALSE)

data_df = data.frame(
    trait = unlist(duh),
    temp = rep(temp_seq, each = n_samp)
)

data_cv_df = data.frame(
    cv = cv_seq,
    temp = temp_seq
)

# Generate predictions
qseq = c(0.025,0.5,0.975)
gam_1 = mqgam(list(trait~s(temp, bs="cs"), ~s(temp, bs="cs")),
              data = data_df,
              qu=qseq)

# Create predictions:
new_df = data.frame(
    temp = seq(10,30,length.out=100)
)

pred_df = data.frame(
    med = qdo(gam_1, 0.5, predict, newdata = new_df),
    lo = qdo(gam_1, 0.025, predict, newdata = new_df),
    hi = qdo(gam_1, 0.975, predict, newdata = new_df),
    temp = seq(10,30,length.out=100)
)

i_plot_fit =
    ggplot(pred_df) + 
    geom_ribbon(aes(x = temp, ymin = lo, ymax = hi), 
                fill = "#AAC9D2") +
    geom_path(aes(x = temp, y = med), 
              col = viridis(45, option = "G")[20]) +
    geom_point(data = data_df,
               aes(x = temp, y = trait), alpha = 0.4) + 
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

i_plot_cv = 
    ggplot(data_cv_df) +
    geom_point(aes(x=temp,y=cv), size = 1.5) +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

##################################################
##################################################

## 3. Decreasing CV: 
##    f(x) = a*exp(-(b*x))

# Redo the mean:
mu = 14
sigma = 9
mean_seq = peak * exp(- (temp_seq-mu)^2 / (2*sigma^2) )
plot(mean_seq~temp_seq)

# Decreasing CV
a_d = 0.3
b_d = 0.08

cv_seq = a_d * exp(- b_d * (temp_seq-10))

plot(cv_seq~temp_seq)

##################################################
##################################################
# Generate Data from the mean and CV
n_samp = 40

duh = apply(cbind(mean_seq,cv_seq), 1, 
            data_func, n_samp = n_samp, simplify = FALSE)

data_df = data.frame(
    trait = unlist(duh),
    temp = rep(temp_seq, each = n_samp)
)

data_cv_df = data.frame(
    cv = cv_seq,
    temp = temp_seq
)

# Generate predictions
qseq = c(0.025,0.5,0.975)
gam_1 = mqgam(list(trait~s(temp, bs="cs"), ~s(temp, bs="cs")),
              data = data_df,
              qu=qseq)

# Create predictions:
new_df = data.frame(
    temp = seq(10,30,length.out=100)
)

pred_df = data.frame(
    med = qdo(gam_1, 0.5, predict, newdata = new_df),
    lo = qdo(gam_1, 0.025, predict, newdata = new_df),
    hi = qdo(gam_1, 0.975, predict, newdata = new_df),
    temp = seq(10,30,length.out=100)
)

d_plot_fit =
    ggplot(pred_df) + 
    geom_ribbon(aes(x = temp, ymin = lo, ymax = hi), 
                fill = "#AAC9D2") +
    geom_path(aes(x = temp, y = med), 
              col = viridis(45, option = "G")[20]) +
    geom_point(data = data_df,
               aes(x = temp, y = trait), alpha = 0.4) + 
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

d_plot_cv = 
    ggplot(data_cv_df) +
    geom_point(aes(x=temp,y=cv), size = 1.5) +
    theme_classic() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

##################################################
##################################################
##################################################
##################################################

# Combine plots:

fit_plots = arrangeGrob(
    u_plot_fit,g_plot_fit,d_plot_fit,i_plot_fit,
    ncol=1, nrow = 4,
    bottom=textGrob(expression("Temperature ("*degree*"C)"),
                    y=unit(1,"lines")),
    top=textGrob("Trait Value",rot=0,hjust=0.5,
                  y=unit(0.5,"lines"))
    )

#quartz(height=8, width = 4.5)
#grid.arrange(fit_plots)

cv_plots = arrangeGrob(
    u_plot_cv,g_plot_cv,d_plot_cv,i_plot_cv,
    ncol=1, nrow = 4,
    bottom=textGrob(expression("Temperature ("*degree*"C)"),
                    y=unit(1,"lines")),
    top=textGrob("Coefficient of Variation",rot=0,hjust=0.5,
                 y=unit(0.5,"lines"))
)

#quartz(height=10, width = 4.5)
#grid.arrange(cv_plots)

plot_combo = arrangeGrob(
    fit_plots,cv_plots,
    ncol=2,nrow=1
)

quartz(height=7, width = 5)
grid.arrange(plot_combo)

ggsave(
    plot = plot_combo,
    file = "trait_cv_plot.png",
    device = "png",
    height = 7.5,
    width = 5
)
