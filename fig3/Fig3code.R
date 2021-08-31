library(ggplot2)
library(tidyverse)
library(deSolve)
library(scales)
library(cowplot)

SIRmod <- function(times, State, Pars) {

     with(as.list(c(State, Pars)), {

      dS <-   -nu * (S/S0)^V * S * I
      dE <-    nu * (S/S0)^V * S * I - delta * E 
      dI <-    (delta * E) 
	return(list(c(dS,  dE, dI)))
	})
 }


   ## Transmission parameters

   nu = .0000001
   V = 5
   V2 = 0.25
   delta = 0.0001
   S0 = 1000000

 parsSIR    <- c(nu = nu, V = V, delta = delta, S0 = S0)
 parsSIR_lv    <- c(nu = nu, V = V2, delta = delta, S0 = S0)


   times <- seq(0, 356, 1)

outSIR_highV <- as.data.frame(ode(func = SIRmod, y = c(S = S0-1, I = 1,  E = 0),  parms = parsSIR, times = times))

outSIR_highV$fi <- nu * (outSIR_highV$S/S0)^V

outSIR_lowV <- as.data.frame(ode(func = SIRmod, y = c(S = S0-1, I = 1,  E = 0), parms = parsSIR_lv, times = times))

outSIR_lowV$fi <- nu * (outSIR_lowV$S/S0)^V2

outSIR_highV$Variation <- "High"
outSIR_lowV$Variation <- "Low"

out <- rbind(outSIR_highV, outSIR_lowV)

names(out)[2:5] <- c("Susceptible", "Infectious", "Exposed", "fi")

outL <-   out%>% gather(key = Compartment, value = Number, -c(time, Variation)) %>% data.frame()
outL$Compartment <- factor(outL$Compartment, levels = c("Susceptible", "Infectious", "Exposed", "fi"))


scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

FI_plot <- outL[outL$Compartment == "fi" & outL$time %in% c(0:365),] %>%
     ggplot(aes(x = time, y = Number, group = Variation, color = Variation)) +  ylab(expression(paste("Force of infection ", bar(beta), "[", "S"[t],"/S"[0],"]"^C^2, sep = ""))) + # ylab("nu[S(t)/S(0)]^C2") + 
     xlab("Time, days") + geom_line(size =1) + scale_x_continuous(breaks = seq(0,365, 50)) +
    scale_y_continuous(label=scientific_10)+ scale_color_manual(values=c("#FF0000","#0E4C92"))+
     theme_bw() + theme(text = element_text(size = 20), strip.background = element_blank())
FI_plot


SIR_grant_plot <- outL[!(outL$Compartment %in% c("Exposed","fi")),]  %>%     ggplot(aes(x = time, y = Number, group = Variation, color = Variation)) + ylab("Number of hosts") + 
     xlab("Time, days") + geom_line(size =1) + facet_grid(~Compartment) + scale_color_manual(values=c("#FF0000","#0E4C92"))+
 scale_y_continuous(label=scientific_10)+
     theme_bw() + theme(text = element_text(size = 20), strip.background = element_blank())
SIR_grant_plot


Figure_1_cartoon <- plot_grid(FI_plot, SIR_grant_plot, ncol = 1)

ggsave("Figure_1_cartoon.png", Figure_1_cartoon, height =8, width = 8)




ggsave("SIR_grant_plot.pdf", SIR_grant_plot, width = 8,  height = 3)








   par(mfrow = c(3, 2))
   plot(outSIR[,2], x = outSIR[,1], type = "l", main = "S")
   plot(outSIR[,3], x = outSIR[,1], type = "l", main = "I")
   plot(outSIR[,4], x = outSIR[,1], type = "l", main = "R")
   plot(outSIR[,5], x = outSIR[,1], type = "l", main = "P")
   plot(outSIR[,2] + outSIR[,3] + outSIR[,4], x = outSIR[,1], type = "l", main = "N")
