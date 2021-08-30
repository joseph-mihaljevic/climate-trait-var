
#-------------------------------------------------------
#-------------------------------------------------------
library(deSolve)
library(tidyverse)
library(grid)
library(gridExtra)
library(rpart)

set.seed(5)
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------

# RELATIONSHIP BETWEEN CV and Temperature:
T_seq = seq(0, 45, by = 0.1)

# Equation: 
# qT(T-T_min)(T_max - T)^(1/2)

q_beta = 10^-3.7;
T_min_beta = 12;
T_max_beta = 34;

# CV gaussian? 
a_CV = 1.2
b_CV = 15
c_CV = 8

# Gaussian
CV_t = a_CV * exp( -0.5 * ((T_seq - b_CV)/c_CV)^2 )

plot(CV_t ~ T_seq, type = "l",
     xlab = "Temp (C)", ylab = expression("C"["(T)"]))

#-------------------------------------------------------
#-------------------------------------------------------
# Create a realistic temperature sequence:

time_seq = c(1:90)
T_temp = 40 * exp(-0.5 * (time_seq + 10 - 55)^2 / 800) + rnorm(length(time_seq), 0, 2.5)
plot(T_temp ~ time_seq, type = "l")

# Fit a piece-wise constant (regression tree)
piece_mod = rpart(T_temp ~ time_seq,
                  control=rpart.control(maxdepth=5, cp=.0001))

# Quick check on results:
T_piece = predict(piece_mod, data.frame(time_seq = time_seq))
plot(T_temp ~ time_seq, type = "l")
lines(T_piece ~ time_seq, type = "l")


CV_piece = a_CV * exp( -0.5 * ((T_piece - b_CV)/c_CV)^2 )
plot(CV_piece ~ time_seq, type = "l") 
mean_CV = round(mean(CV_piece),2)
 
#-------------------------------------------------------
#-------------------------------------------------------
# Quick ggplot of the CV:

CV_df = 
  data.frame(CV_f = rep(mean_CV, length(time_seq)),
             CV_t = CV_piece, 
             time = time_seq) %>%
  pivot_longer(CV_f:CV_t, names_to = "CV_type",
               values_to = "CV") %>%
  mutate(CV_type = case_when(
    CV_type == "CV_f" ~ "Fixed C", 
    CV_type == "CV_t" ~ "Temperature-Dependent C",
    TRUE ~ "NA"
  ))

#expression("Temp."*degree*"C")
CV_plot = 
  ggplot(CV_df) +
  annotate("segment", x = min(time_seq), xend = max(time_seq),
           y = mean_CV, yend = mean_CV, linetype = 2) +
  geom_path(aes(x = time, y = CV), color = "blue", size = 1.5) +
  scale_y_continuous(limits = c(0, 1.5)) +
  scale_x_continuous(limits = c(0,max(time_seq)), breaks = c(0,30,60,90)) +
  theme_classic() +
  labs(x = "Day of Annual Epizootic", y = "Coef. of Variation, C") +
  facet_wrap(~CV_type, nrow = 2) +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 13),
    strip.text = element_text(color = "black", size = 12)

  )

#-------------------------------------------------------
#-------------------------------------------------------
# CV vs TEMP

cv_temp_df = data.frame(CV_t, T_seq)

CV_temp_plot =
  ggplot(cv_temp_df %>% filter(T_seq >=10)) + 
  geom_path(aes(x = T_seq, y = CV_t), color = "blue", size = 1.5) +
  scale_x_continuous(limits = c(5,45)) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme_classic() +
  labs(x = expression("Temp. ("*degree*"C)"), y = "Coef. of Var., C") +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 13),
    strip.text = element_text(color = "black", size = 12)
  )

time_temp_df = data.frame(time_seq, T_temp, T_piece)

Temp_time_plot =
  ggplot(time_temp_df) + 
  geom_path(aes(x = time_seq, y = T_temp), color = "skyblue2", size = 1) +
  geom_path(aes(x = time_seq, y = T_piece), color = "blue", size = 1.5) +
  scale_y_continuous(limits = c(5, 45), breaks = c(10,20,30,40)) +
  scale_x_continuous(limits = c(0,max(time_seq)), breaks = c(0,30,60,90)) +
  theme_classic() +
  labs(x = "Time (days)", y = expression("Temp. ("*degree*"C)")) +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 13),
    strip.text = element_text(color = "black", size = 12)
  )

plot_CV_temp_combo = 
  arrangeGrob(CV_temp_plot, Temp_time_plot, nrow = 2)

ggsave(
  "CV_temp_combo.png",
  plot_CV_temp_combo,
  device = "png",
  dpi = 300,
  height = 3.7, width = 3
)


#-------------------------------------------------------
#-------------------------------------------------------
# Function to solve the ODE system

forced_SIR_ode = function(t, y, params){
  with(as.list(c(y, params)), {
    
    # storage
    ## We need to have a vector of size equal to the number of ODEs in the system
    dydt = rep(0, 3)
    
    # equations:
    dydt[1] = - beta * y[1] * y[4] * ((y[1]/S_zero)^(CV^2))
    dydt[2] =  beta * y[1] * y[4] * ((y[1]/S_zero)^(CV^2)) - (gamma) * y[2]
    dydt[3] = (gamma) * y[2] - mu * y[3]
    dydt[4] = psi*y[3] - phi*y[4]

    return(list(dydt))
    
  })
}

#-------------------------------------------------------
#-------------------------------------------------------

# Define a function that will solve the ODE over one day 
# (i.e. one "step") at a time

solve_ODE_1step = function(inits, params){
  
  # Time to solve ODE, one day at a time:
  time_ode = seq(0, 1, by = 0.1)
 
  # Run your ODE with the inputs to the function:
  out_temp = rk(y = inits, 
                 times = time_ode, 
                 func = forced_SIR_ode, 
                 method="rk45f", # Runge-Kutta 4-5 method
                 parms = params,
                rtol = 1e-10,
                atol = 1e-10)
  
  # Object 'out' is a matrix
  # Specify the column names
  # Then convert to data.frame object
  colnames(out_temp) = c("time", "S", "E", "I", "P")
  out_temp = data.frame(out_temp)
  
  # Return the data frame
  return(out_temp)
  
}

#-------------------------------------------------------
#-------------------------------------------------------

# Allow CV to be fixed or vary with temperature following the 
# thermal performance curve generated above. 

CV_func = function(i, CV_0 = 0.6, CV_type = "temp"){
  
  if(CV_type == "fixed"){
    CV_temp = CV_0
  }
  
  if(CV_type == "temp"){
    
    # From the piece_wise function:
    CV_temp = as.numeric(CV_piece[i])
    
    # offset_days = 100
    # T_temp = 45 * exp(-0.5 * (i + offset_days - 180)^2 / 5000) - 10 + rnorm(1, 0, 2.5)
    # # GAUSSIAN:
    # CV_temp =a_CV * exp( -0.5 * ((T_temp - b_CV)/c_CV)^2 )
    
    if(CV_temp < 0 | is.na(CV_temp)){
      CV_temp = 0.0
    }
    
  }
  
  return(CV_temp)
  
}

#-------------------------------------------------------
#-------------------------------------------------------

# Allow CV to be fixed or vary with temperature following the 
# thermal performance curve generated above. 

beta_func = function(i, beta_0 = 0.05, beta_type = "temp"){
  
  if(beta_type == "fixed"){
    beta_temp = beta_0
  }
  
  if(beta_type == "temp"){
    offset_days = 100
    T_temp = 45 * exp(-0.5 * (i + offset_days - 180)^2 / 5000) - 10 + rnorm(1, 0, 2.5)
    beta_temp = q_beta * T_temp * (T_temp - T_min_beta) * (T_max_beta - T_temp)^(1/2);
    
    if(beta_temp < 0 | is.na(beta_temp)){
      beta_temp = 0.0
    }
    
  }
  
  return(beta_temp)
  
}

#-------------------------------------------------------
#-------------------------------------------------------

# Define a function that will solve the model for all 
# days within the annual epizootic and output the final state


intra_annual_func = function(t_max, params, beta_type, CV_type){
  
  
  for(i in 1:t_max){
    
    if(i == 1){ # If this is the first time looping...
      # Set your initial values manually:
      S0 = params$S_zero
      E0 = 0
      I0 = 0
      P0 = params$P_zero
    }else{ # otherwise...
      # What am I doing here??
      S0 = df_temp[nrow(df_temp), 2]
      E0 = df_temp[nrow(df_temp), 3]
      I0 = df_temp[nrow(df_temp), 4]
      P0 = df_temp[nrow(df_temp), 5]
    }
    
    inits_temp = c(S0, E0, I0, P0)
    
    # Calculate your \beta and CV value at a specific value of t
    beta_temp = beta_func(i, beta_0 = 0.005, beta_type = beta_type)
    #CV_temp = 0.6
    CV_temp = CV_func(i, CV_0 = mean_CV, CV_type = CV_type) #0.5
    #print(CV_temp)
    
    # Store your params:
    ODE_params_temp = c(CV = CV_temp,
                        beta = beta_temp,
                        S_zero = params$S_zero,
                        gamma = params$gamma,
                        mu = params$mu,
                        psi = params$psi,
                        phi = params$phi)
    
    # Solve over one day
    df_temp = solve_ODE_1step(inits = inits_temp,
                              params = ODE_params_temp)
    
    # What is the rest of this doing??
    time_ode_fix = seq(i, i+1, by = .1)
    df_temp$time = time_ode_fix
    
    if(i == 1){
      out = df_temp
    }else{
      out = rbind(out, df_temp)
    }
    
  }
  
  out$R_inf = 1 - out$S/out$S[1]
  
  return(list(end_epi = out[nrow(out),],
              R_inf = data.frame(time = out$time-1,
                                 R_inf = out$R_inf))
  )
  
}

#-------------------------------------------------------
#-------------------------------------------------------

# Define a function that will solve the host demography 
# between years (winter season)

inter_annual_func = function(inits, params){
  
  next_gen = list()
  
  next_gen$S = params$lambda * inits$S# * exp((rnorm(1, 0, 0.0001)))
  next_gen$P = params$omega * inits$P
    
  return(next_gen)
  
}

#-------------------------------------------------------
#-------------------------------------------------------
# Run with the two CV types (fixed or temperature-dependent)

these_CV_types = c("fixed", "temp")

# How many years do you want to solve?
n_year = 50

# STORAGE:
annual_store =
  array(0, 
        dim = c(n_year+1, 2, 2),
        dimnames = list(
          Gen = 1:(n_year+1),
          Class = c("S", "P"),
          CV = c("Fixed C", "Temperature-Dependent C")
        ))

# RUN THE LOOP
for(k in 1:2){
  
  # INITIAL PARAMETER VALUES:
  S_zero_init = 18.1 # 
  P_zero_init = 2.1e-2 # 
  
  params = list(
    
    t_max = max(time_seq),
    
    # EPIZOOTIC
    gamma = 1/5, # incubation
    mu = 1/5, # infectious mortality
    psi = 1.5, # shedding (pathogen/day/host)
    phi = 0.3, # pathogen decay
    CV = NULL,
    
    S_zero = NULL,
    P_zero = NULL,
    
    # DEMOGRAPHY
    lambda = 5.5,
    omega = 0.15 #0.15
  )
  
  for(i in 1:n_year){
    
    if(i == 1){
      params$S_zero = S_zero_init
      params$P_zero = P_zero_init
      
      annual_store[i, 1, k] = S_zero_init
      annual_store[i, 2, k] = P_zero_init
      
    }else{
      params$S_zero = next_inits$S
      params$R_zero = next_inits$P
    }
    
    end_epizootic = 
      intra_annual_func(params$t_max, params, 
                        beta_type = "fixed", 
                        CV_type = these_CV_types[k])
    
    # inits_temp = list(
    #   fracI = 1 - (end_epizootic$S / params$S_zero),
    #   oldR = params$R_zero
    # )
    
    next_inits = 
      inter_annual_func(end_epizootic$end_epi, params)
    
    annual_store[i+1, 1, k] = next_inits$S
    annual_store[i+1, 2, k] = next_inits$P
    
  }
  
  
}

  
#-------------------------------------------------------
#-------------------------------------------------------

# Plot the dynamics. 
{
  multi_gen_df = 
    data.frame(ftable(annual_store)) %>%
    mutate(Density = Freq,
           Gen = as.numeric(Gen))%>%
    select(-Freq)
  
  mean_dens_df = 
    multi_gen_df %>%
    filter(Class == "S") %>%
    group_by(CV) %>%
    summarize(mean_Density = mean(Density, na.rm=TRUE))
  #mean_dens_df
  
  multi_gen_plot = 
    ggplot(multi_gen_df %>% filter(Class == "S")) +
    annotate("segment", x = 0, xend = n_year,
             y = 0, yend = 0, linetype = 2) +
    geom_point(aes(x = Gen, y = (Density)), 
               shape = 19, color = "red", size = 1.5) +
    geom_path(aes(x = Gen, y = (Density)), 
              color = "red", size = 1.1) +
    theme_classic() +
    facet_wrap(~CV, nrow = 2) +
    labs(x = "Generation", y = "Host Density") +
    theme(
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 13),
      strip.text = element_text(color = "black", size = 12)
    )
  
  multi_gen_plot
  
  multi_gen_plot_virus = 
    ggplot(multi_gen_df %>% filter(Class == "P")) +
    annotate("segment", x = 0, xend = n_year,
             y = 0, yend = 0, linetype = 2) +
    geom_point(aes(x = Gen, y = Density), 
               shape = 19, color = "red", size = 1.5) +
    geom_path(aes(x = Gen, y = Density), 
              color = "red", size = 1.1) +
    theme_classic() +
    facet_wrap(~CV, nrow = 2) +
    labs(x = "Generation", y = "Pathogen Density") +
    theme(
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 13),
      strip.text = element_text(color = "black", size = 12)
    )
  
  multi_gen_plot_virus
}



plot_combo = 
  arrangeGrob(CV_plot, multi_gen_plot,
              layout_matrix = matrix(c(1,2,2,
                                       1,2,2), 
                                     nrow = 2, byrow = T))

ggsave(
  "multigen_combo_plot.png",
  plot_combo,
  device = "png",
  dpi = 300,
  height = 4.5, width = 10
)


#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------

# Look at R_inf over time as function of initial Host density

S_zero_seq = ((seq((0.1), (50), length.out = 30)))
S_zero_df = data.frame(
  S_zero = paste0(1:length(S_zero_seq)),
  S_zero_lab = factor(S_zero_seq, levels = S_zero_seq)
)

# Storage:
# rinf_fixed_list = vector("list", length(S_zero_seq))
# rinf_temp_list = vector("list", length(S_zero_seq))
rinf_fixed_vec = vector("numeric", length(S_zero_seq))
rinf_temp_vec = vector("numeric", length(S_zero_seq))

params = list(
  
  t_max = max(time_seq),
  
  # EPIZOOTIC
  gamma = 1/5, # incubation
  mu = 1/5, # infectious mortality
  psi = 1.5, # shedding (pathogen/day/host)
  phi = 0.3, # pathogen decay
  CV = NULL,
  
  S_zero = NULL,
  P_zero = 1.0,#2.1e-2,
  
  # DEMOGRAPHY
  lambda = 5.5,
  omega = 0.15 #0.15
)

# RUN THE LOOP
for(k in 1:length(S_zero_seq)){
  
  # INITIAL PARAMETER VALUES:
  S_zero_init = S_zero_seq[k]

  params$S_zero = S_zero_init
  
  for(j in 1:2){
    epi_traj = 
      intra_annual_func(params$t_max, params, 
                        beta_type = "fixed", 
                        CV_type = these_CV_types[j])
    
    if(j==1){
      #rinf_fixed_list[[k]] = epi_traj$R_inf
      rinf_fixed_vec[k] = epi_traj$end_epi$R_inf
    }else{
      #rinf_temp_list[[k]] = epi_traj$R_inf
      rinf_temp_vec[k] = epi_traj$end_epi$R_inf
    }
    
    rm(epi_traj)
  }
  
}

# rinf_fixed_df =
#   map_dfr(
#     rinf_fixed_list,
#     ~.,
#     .id = "S_zero"
#   ) %>%
#   left_join(y=S_zero_df, by = "S_zero")

# ggplot(rinf_fixed_df) +
#   geom_path(aes(x = time, y = R_inf, color=S_zero_lab)) +
#   scale_color_viridis_d("S(0)")

# rinf_temp_df =
#   map_dfr(
#     rinf_temp_list,
#     ~.,
#     .id = "S_zero"
#   ) %>%
#   left_join(y=S_zero_df, by = "S_zero")

# ggplot(rinf_temp_df) +
#   geom_path(aes(x = time, y = R_inf, color=S_zero_lab)) +
#   scale_color_viridis_d("S(0)")

# rinf_full_df = 
#   rbind(rinf_fixed_df, rinf_temp_df)
# rinf_full_df$Type = rep(c("Fixed C", "Temperature-Dependent C"),
#                         each = nrow(rinf_temp_df))

# ggplot(rinf_full_df) +
#   geom_path(aes(x = time, y = R_inf, color=S_zero_lab)) +
#   scale_color_viridis_d("S(0)") +
#   facet_wrap(~Type, ncol=1, nrow=2)

r_full_df = data.frame(
  R_inf = c(rinf_fixed_vec,rinf_temp_vec),
  S_zero = rep(S_zero_seq, 2),
  Type = rep(c("Fixed C", "Temperature-Dependent C"),
             each = length((rinf_fixed_vec)))
)

ggplot(r_full_df) +
  geom_path(aes(x = S_zero, y = R_inf, linetype=Type)) +
  theme_classic() +
  labs(x = "Host Density, S(0)", y = "Cumulative Frac. Infected") +
  scale_linetype("") +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 13),
    strip.text = element_text(color = "black", size = 12),
    legend.position = c(0.7,0.2)
  )


