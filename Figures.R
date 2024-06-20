# This script used R 4.3.1. 
# Author: Yu Zhou (zhouy.work@gmail.com)

# set working folder 
setwd("C:\\Users\\yz2872\\OneDrive - Cornell University\\ClarkU\\Paper3\\2_Constrained_results")

# Package names
packages <- c("R.matlab", "Hmisc", "matrixStats", "gghalves", "GGally", "ggplot2", "viridis", "MASS", "ismev", "dplyr", "broom", "multcompView")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# load packages
library(R.matlab)
library(Hmisc)
library(matrixStats)
library(gghalves)
library(GGally)
library(ggplot2)
library(viridis)
library(MASS)
library(ismev)
library(dplyr)
library(broom)
library(multcompView)

# define all forest groups 
fttppd_all <- c('HighProd_DouglasFir', 'LowProd_DouglasFir', 'LowProd_PonderosaPine','LowProd_FirSpruceMHemlock');
# define parameters 
para.names <- c("Emax",  "τwood", "αAGwood",  "αwood","αCWD", "τCWD",  "τsnag", "Q10", "kslow", "SFFlitter")

### plot parameters (currently using) ######

for (i_fttppd in 1:length(fttppd_all)) {
  fttppd <- fttppd_all[i_fttppd];
  # read from HDA0 (constrained by all pools at the same time; as benchmark)
  data <- readMat(paste0(fttppd, "_paras_posterior_HD_0.mat", sep = ""))
  parameters.keep0 <- data$para.posterior
  
  # read from HDA1 (constrained by live biomass)
  data <- readMat(paste0(fttppd, "_paras_posterior_HD_1.mat", sep = ""))
  parameters.keep1 <- data$para.posterior
  
  # read from HDA2 (constrained by live and then dead biomass)
  data <- readMat(paste0(fttppd, "_paras_posterior_HD_2.mat", sep = ""))
  parameters.keep2 <- data$para.posterior
  
  # read from HDA3 (constrained by live, then dead biomass, and SOC at the end)
  data <- readMat(paste0(fttppd, "_paras_posterior_HD_3.mat", sep = ""))
  parameters.keep3 <- data$para.posterior
  
  # define complied para data frame
  para <- data.frame(matrix(data = NA, ncol = 12, nrow = (ncol(parameters.keep1) + ncol(parameters.keep2) + ncol(parameters.keep3) + ncol(parameters.keep0))))
  colnames(para) <- c("HDA", "Emax",  "τwood", "αAGwood",  "αwood",
                      "αCWD", "τCWD",  "τsnag", "Q10", "kslow", "SFFlitter", "Foresttype")
  # add posterior parameters
  para[, 2:11] <- t(cbind(parameters.keep1, parameters.keep2, parameters.keep3, parameters.keep0))
  # add forest group index
  para$Foresttype <- i_fttppd
  
  # add HDA index
  para$HDA[1:ncol(parameters.keep1)] <- 1
  para$HDA[(ncol(parameters.keep1) + 1): (ncol(parameters.keep2)+ ncol(parameters.keep1))] <- 2
  para$HDA[(ncol(parameters.keep1) + ncol(parameters.keep2)+ 1): (ncol(parameters.keep3)+ ncol(parameters.keep2)+ ncol(parameters.keep1))] <- 3
  para$HDA[(ncol(parameters.keep3) + ncol(parameters.keep1) + ncol(parameters.keep2)+ 1): (ncol(parameters.keep0) + ncol(parameters.keep3)+ ncol(parameters.keep2)+ ncol(parameters.keep1))] <- 0
  
  para$HDA <- as.factor(para$HDA)
  para$Foresttype <- as.factor(para$Foresttype)
  
  if (i_fttppd == 1) {
    all.para <- para
  } else {
    all.para <- rbind(all.para, para)
  }
}

# convert % to 100%
all.para$αAGwood <- all.para$αAGwood *100
all.para$αwood <- all.para$αwood *100
all.para$αCWD <- all.para$αCWD *100

para.names.outrank <- c("Emax",  "τwood", "αAGwood",  "αwood",
                "αCWD", "τCWD",  "τsnag", "SFFlitter", "Q10", "kslow")
para.min <- c(0.2, 50, 20, 20, 10, 1, 1, 0.1, 1, 0.01)
para.max <- c(3, 200, 90, 80, 80, 30, 40, 3, 2, 0.3)
para.def <- c(0.65,   75,   75,   33, 50, 10, 10, 1, 1.5, 0.2)
xlab2 <- bquote(para.label);
para.label <- c("g[.(i)]")

# For figure 4 in main text
for (i_para in 1:length(para.def)) { 
  eval(parse(text = paste("all.para$current <- all.para$", para.names.outrank[i_para], sep = "")))
  if (i_para == 2 | i_para == 6 | i_para ==7) {
    all.para$current <- log10(all.para$current)
  }
  
  p <- ggplot(all.para, aes(x=Foresttype, y=current)) +
    geom_half_violin(data = all.para %>% filter(HDA == "1"),
                     alpha = 0.5, fill = "#0c84c6", trim = TRUE, colour = "grey20", linewidth = 0.3,
                     side = "l") +
    geom_half_violin( data = all.para %>% filter(HDA == "2"),
                      alpha = 0.5, fill = "#ffbd66", trim = TRUE, colour = "grey20", linewidth = 0.3,
                      side = "l") +
    geom_half_violin(data = all.para %>% filter(HDA == "3"),
                     alpha = 0.5, fill = "#f74d4d", trim = TRUE, colour = "grey20", linewidth = 0.3,
                     side = "r") +
    # geom_half_violin(data = all.para %>% filter(HDA == "0"),
    #                  alpha = 0.5, fill = "grey20", trim = TRUE, colour = "grey20", linewidth = 0.3,
    #                  side = "r") +
    geom_hline(yintercept = para.def[i_para], linetype = "dashed", colour = "grey20", linewidth = 0.3) +
    xlab("Forest group type") +
    scale_x_discrete(breaks=c("1","2","3", "4"),
                     labels=c("HP Df", "LP Df", "LP Pp", "LP F/S/MH")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 14),axis.text = element_text(size = 14, colour = "black"))  
  
  if (i_para ==1 ) { p <- p + ylab(bquote(E[max]~(g~C~MJ^-1))) + ggtitle('(a)')+ theme(plot.title = element_text( hjust = 0.01, vjust = - 10, size = 14))}
  if (i_para ==2 ) { p <- p + ylab(bquote(τ[wood]~('log10(year)'))) + ggtitle('(b)')+ theme(plot.title = element_text( hjust = 0.01, vjust = - 10, size = 14))}
  if (i_para ==3 ) { p <- p + ylab(bquote(α[AGwood]~('% of wood'))) + ggtitle('(c)')+ theme(plot.title = element_text( hjust = 0.01, vjust = - 10, size = 14))}
  if (i_para ==4 ) { p <- p + ylab(bquote(α[wood]~('% of NPP'))) + ggtitle('(d)')+ theme(plot.title = element_text( hjust = 0.01, vjust = - 10, size = 14))}
  if (i_para ==5 ) { p <- p + ylab(bquote(α[CWD]~('%'))) + ggtitle('(e)')+ theme(plot.title = element_text( hjust = 0.01, vjust = - 10, size = 14))}
  if (i_para ==6 ) { p <- p + ylab(bquote(τ[CWD]~('log10(year)'))) + ggtitle('(f)')+ theme(plot.title = element_text( hjust = 0.01, vjust = - 10, size = 14))}
  if (i_para ==7 ) { p <- p + ylab(bquote(τ[snag]~('log10(year)'))) + ggtitle('(g)')+ theme(plot.title = element_text( hjust = 0.01, vjust = - 10, size = 14))}
  if (i_para ==8 ) { p <- p + ylab(bquote(S[FF~litter])) + ggtitle('(h)')+ theme(plot.title = element_text( hjust = 0.01, vjust = - 10, size = 14))}
  if (i_para ==9 ) { p <- p + ylab(bquote(Q[10])) + ggtitle('(i)')+ theme(plot.title = element_text( hjust = 0.01, vjust = - 10, size = 14))}
  if (i_para ==10 ) { p <- p + ylab(bquote(k[slow]~(year^-1))) + ggtitle('(j)')+ theme(plot.title = element_text( hjust = 0.01, vjust = - 10, size = 14))}

  ggsave(paste0("Fig4_", para.names.outrank[i_para], ".tiff"), units="in", width=6, height=3, dpi=400)
}



lower_density_plot <- function(data, mapping, N=100, ...){
  
  get_density <- function(x, y, n ) {
    dens <- MASS::kde2d(x = x, y = y, n = n)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  X <- eval_data_col(data, mapping$x)
  Y <- eval_data_col(data, mapping$y)
  # Remove rows with missing or infinite values
  finite_idx <- is.finite(X) & is.finite(Y)
  X <- X[finite_idx]
  Y <- Y[finite_idx]
  
  data <- data[finite_idx,]
  
  data$density <- get_density(x=X, y=Y, n=N)
  
  p <- ggplot(data, mapping) +
    geom_point(aes(colour=density), ...) +
    scale_color_viridis()      
  p
}
# For figure S2 in the supplement 
# parameter correlations in HDA1
para.currHDA <- all.para[all.para$HDA == 1, c("Foresttype","Emax","τwood","αAGwood","αwood")]
para.currHDA$τwood <- log10(para.currHDA$τwood)
para.currHDA$Foresttype <- factor(para.currHDA$Foresttype,
                                  levels = c(1, 2, 3, 4),
                                  labels = c('HP Df',
                                             'LP Df',
                                             'LP Pp',
                                             'LP F/S/MH'))
p <- ggpairs(para.currHDA, columns =  c("Emax","τwood","αAGwood","αwood"),
        upper = list(continuous = "cor"),
        diag = list(continuous = "densityDiag"),
        lower = list(continuous = lower_density_plot),
        mapping = ggplot2::aes(color = Foresttype)) +
  theme_bw()
ggsave(paste0("FigS2_HDA1", "paracorr", ".tiff"), p, units="in", width=6, height=6, dpi=300)

# For figure S3 in the supplement 
# parameter correlations in HDA2
para.currHDA <- all.para[all.para$HDA == 2, c("Foresttype","Emax","τwood","αAGwood","αwood","αCWD", "τCWD",  "τsnag", "SFFlitter", "Q10", "kslow")]
para.currHDA$τwood <- log10(para.currHDA$τwood)
para.currHDA$τsnag <- log10(para.currHDA$τsnag)
para.currHDA$τCWD <- log10(para.currHDA$τCWD)
para.currHDA$Foresttype <- factor(para.currHDA$Foresttype,
                                  levels = c(1, 2, 3, 4),
                                  labels = c('HP Df',
                                             'LP Df',
                                             'LP Pp',
                                             'LP F/S/MH'))
p <- ggpairs(para.currHDA, columns =  c("Emax","τwood","αAGwood","αwood","αCWD", "τCWD",  "τsnag", "SFFlitter", "Q10", "kslow"),
        upper = list(continuous = "cor"),
        diag = list(continuous = "densityDiag"),
        lower = list(continuous = lower_density_plot),
        mapping = ggplot2::aes(color = Foresttype)) +
  theme_bw()
ggsave(paste0("FigS3_HDA2", "paracorr", ".tiff"), p, units="in", width=16, height=16, dpi=300)

# For figure S4 in the supplement 
# parameter correlations in HDA3
para.currHDA <- all.para[all.para$HDA == 3, c("Foresttype","Emax","τwood","αAGwood","αwood","αCWD", "τCWD",  "τsnag", "SFFlitter", "Q10", "kslow")]
para.currHDA$τwood <- log10(para.currHDA$τwood)
para.currHDA$τsnag <- log10(para.currHDA$τsnag)
para.currHDA$τCWD <- log10(para.currHDA$τCWD)
para.currHDA$Foresttype <- factor(para.currHDA$Foresttype,
                                  levels = c(1, 2, 3, 4),
                                  labels = c('HP Df',
                                             'LP Df',
                                             'LP Pp',
                                             'LP F/S/MH'))
p <- ggpairs(para.currHDA, columns =  c("Emax","τwood","αAGwood","αwood","αCWD", "τCWD",  "τsnag", "SFFlitter", "Q10", "kslow"),
             upper = list(continuous = "cor"),
             diag = list(continuous = "densityDiag"),
             lower = list(continuous = lower_density_plot),
             mapping = ggplot2::aes(color = Foresttype)) +
  theme_bw()
ggsave(paste0("FigS4_HDA3", "paracorr", ".tiff"), p, units="in", width=16, height=16, dpi=300)


#### get parameter correlations 
# correlations and plots
cor.mtest <- function(x, y, method) {
  p.mat <- matrix(NA, ncol(x), ncol(y))
  for (i in 1:ncol(x)) {
    for (j in 1:ncol(y)) {
      tmp <- cor.test(x[, i], y[, j], method = method,exact=FALSE)
      p.mat[i, j] <- tmp$p.value
    }
  }
  colnames(p.mat) <- colnames(y)
  rownames(p.mat) <- colnames(x)
  p.mat
}
plot_corpara <- function(i_HDA, corr.p.para, corr.para){
  pos <- which(corr.p.para >= 0.01)
  corr.para[pos] <- NA
  corr.para[corr.para == 1] <- NA
  # 
  # corr.para[!is.na(corr.para)] <- 1
  # corr.para.sum <- rowSums(corr.para, na.rm = T, dims = 2)
  corr.para.neg <- corr.para.pos <- array(c(NA), dim = c(10,10,4))
  corr.para.neg[corr.para < 0] <- -1
  corr.para.pos[corr.para > 0] <- 1
  corr.para.negsum <- rowSums(corr.para.neg, na.rm = T, dims = 2)
  corr.para.possum <- rowSums(corr.para.pos, na.rm = T, dims = 2)
  pos <- which(corr.para.negsum <= -3)
  corr.para.sum[pos] <- corr.para.negsum[pos]
  pos <- which(corr.para.possum >= 3)
  corr.para.sum[pos] <- corr.para.possum[pos]
  
  corr.df <- data.frame(matrix(data = NA, nrow = 10*10, ncol = 3))
  colnames(corr.df) <- c("x","y","r")
  para.names  <- c("Emax",  "τwood", "αAGwood",  "αwood",
                   "αCWD", "τCWD",  "τsnag", "Q10", "kslow", "SFFlitter")
  i.row <- 1
  for (i in 1:10) {
    for (j in 1:10) {
      corr.df$x[i.row] <- para.names[i]
      corr.df$y[i.row] <- para.names[j]
      if (i == j) {
        corr.df$r[i.row] <- corr.para.sum[i,j] <- 0
      } else {
        corr.df$r[i.row] <- corr.para.sum[i,j]
      }
      i.row <- i.row + 1
    }
  }
  corr.df$r[corr.df$r == 0] <- NA
  corr.df$ord.x <- factor(corr.df$x, ordered=TRUE, levels = c("Emax",  "τwood", "αAGwood",  "αwood",
                                                              "αCWD", "τCWD",  "τsnag", "SFFlitter", "Q10", "kslow"))
  corr.df$ord.y <- factor(corr.df$y, ordered=TRUE, levels = rev(c("Emax",  "τwood", "αAGwood",  "αwood",
                                                                  "αCWD", "τCWD",  "τsnag", "SFFlitter", "Q10", "kslow")))
  # Get upper-left triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  corr.df$x[corr.df$ord.x == "Emax"] <- corr.df$y[corr.df$ord.y == "Emax"] <- 1
  corr.df$x[corr.df$ord.x == "τwood"] <- corr.df$y[corr.df$ord.y == "τwood"] <- 2
  corr.df$x[corr.df$ord.x == "αAGwood"] <- corr.df$y[corr.df$ord.y == "αAGwood"] <- 3
  corr.df$x[corr.df$ord.x == "αwood"] <- corr.df$y[corr.df$ord.y == "αwood"] <- 4
  corr.df$x[corr.df$ord.x == "αCWD"] <- corr.df$y[corr.df$ord.y == "αCWD"] <- 5
  corr.df$x[corr.df$ord.x == "τCWD"] <- corr.df$y[corr.df$ord.y == "τCWD"] <- 6
  corr.df$x[corr.df$ord.x == "τsnag"] <- corr.df$y[corr.df$ord.y == "τsnag"] <- 7
  corr.df$x[corr.df$ord.x == "SFFlitter"] <- corr.df$y[corr.df$ord.y == "SFFlitter"] <- 8
  corr.df$x[corr.df$ord.x == "Q10"] <- corr.df$y[corr.df$ord.y == "Q10"] <- 9
  corr.df$x[corr.df$ord.x == "kslow"] <- corr.df$y[corr.df$ord.y == "kslow"] <- 10
  
  corr.df$r[is.na(corr.df$r)] <- 0
  for (i in 1:10) {
    for (j in i:10) {
      if (i == 1&j==1) {
        lower_triangle_df <- corr.df[corr.df$x == i & corr.df$y == j,]
      } else {
        lower_triangle_df <- rbind(lower_triangle_df,  corr.df[corr.df$x == i & corr.df$y == j,])
      }
    }
  }
  lower_triangle_df$r[lower_triangle_df$r <3 & lower_triangle_df$r >-3] <- NA
  ggplot(lower_triangle_df,aes(x = ord.x, y = ord.y, fill = r)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "#075AFF",
                         mid = "white",
                         high = "red3",
                         na.value = "white") +
    geom_text(aes(label = r), color = "black", size = 4) +
    coord_fixed() +
    guides(fill = guide_colourbar(title = "")) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 12), # Adjust size of x-axis labels
          axis.text.y = element_text(size = 12),  # Adjust size of y-axis labels 
          axis.title = element_blank(),
          panel.grid.major = element_blank(),   # Remove major grid lines
          panel.grid.minor = element_blank(),   # Remove minor grid lines
          axis.line = element_blank(), # Remove axis lines)   # Remove minor grid lines) 
          panel.border = element_blank())   + guides(fill = FALSE)  # Remove legend         
  
  ggsave(paste0("Fig_S5_HDA",i_HDA, "_para_correlation.tiff"), units="in", width=6, height=6, dpi=400)
}

para.gev <- matrix(data = NA, nrow = 12, ncol = 10) #forest type * HDA; paras
# for HDA1
rm(corr.p.para, corr.para)
corr.para <- array(c(NA), dim = c(10,10,4))
corr.p.para <- array(c(NA), dim = c(10,10,4))
corr.para.sum <- array(c(NA), dim = c(10,10))
for (i_fttppd in 1:length(fttppd_all)) {
  curr_para <- all.para[all.para$HDA == 1 & all.para$Foresttype == i_fttppd, ]
  pos <- which(is.na(curr_para$Emax))
  curr_para <- curr_para[-pos,] 
  curr_para <- (curr_para[,para.names])
  corr.para[1:4,1:4,i_fttppd] <- cor((curr_para), method = "pearson")[1:4,1:4]
  corr.p.para[1:4,1:4,i_fttppd] <- cor.mtest(curr_para, curr_para, method = "pearson")[1:4,1:4]
  for (i_para in 1:length(para.names)) {
    eval(parse(text = paste("tmp <- curr_para$", para.names[i_para], sep = "")))
    mu <- gev.fit(tmp,type = c("mle"), show = F)$mle[1]
    para.gev[(i_fttppd-1)*3+1,i_para] <- mu
  }
}
plot_corpara(i_HDA = 1, corr.p.para, corr.para)
# for HDA2
rm(corr.p.para, corr.para)
corr.para <- array(c(NA), dim = c(10,10,4))
corr.p.para <- array(c(NA), dim = c(10,10,4))
corr.para.sum <- array(c(NA), dim = c(10,10))
for (i_fttppd in 1:length(fttppd_all)) {
  curr_para <- all.para[all.para$HDA == 2 & all.para$Foresttype == i_fttppd, ]
  pos <- which(is.na(curr_para$Emax))
  curr_para <- curr_para[-pos,] 
  curr_para <- (curr_para[,para.names])
  corr.para[1:9,1:9,i_fttppd] <- cor((curr_para), method = "pearson")[1:9,1:9]
  corr.p.para[1:9,1:9,i_fttppd] <- cor.mtest(curr_para, curr_para, method = "pearson")[1:9,1:9]
  for (i_para in 1:length(para.names)) {
    eval(parse(text = paste("tmp <- curr_para$", para.names[i_para], sep = "")))
    mu <- gev.fit(tmp,type = c("mle"), show = F)$mle[1]
    # mu <- gevFit(tmp,type = c("mle"))@fit$par.ests[2]
    para.gev[(i_fttppd-1)*3+2,i_para] <- mu
  }
}
plot_corpara(i_HDA = 2, corr.p.para, corr.para)

# For figure S5 in the supplement 
# for HDA3
rm(corr.p.para, corr.para)
corr.para <- array(c(NA), dim = c(10,10,4))
corr.p.para <- array(c(NA), dim = c(10,10,4))
corr.para.sum <- array(c(NA), dim = c(10,10))
for (i_fttppd in 1:length(fttppd_all)) {
  curr_para <- all.para[all.para$HDA == 3 & all.para$Foresttype == i_fttppd, ]
  pos <- which(is.na(curr_para$Emax))
  curr_para <- curr_para[-pos,] 
  curr_para <- (curr_para[,para.names])
  corr.para[1:10,1:10,i_fttppd] <- cor((curr_para), method = "pearson")[1:10,1:10]
  corr.p.para[1:10,1:10,i_fttppd] <- cor.mtest(curr_para, curr_para, method = "pearson")[1:10,1:10]
  for (i_para in 1:length(para.names)) {
    eval(parse(text = paste("tmp <- curr_para$", para.names[i_para], sep = "")))
    mu <- gev.fit(tmp,type = c("mle"), show = F)$mle[1]
    # mu <- gevFit(tmp,type = c("mle"))@fit$par.ests[2]
    para.gev[(i_fttppd-1)*3+3,i_para] <- mu
  }
}
plot_corpara(i_HDA = 3, corr.p.para, corr.para)

# For Table 2 in the main text
colnames(para.gev) <- para.names
rownames(para.gev) <- c("HP Df HDA1", "HP Df HDA2", "LP Df HDA3", "LP Df HDA1", "LP Df HDA2", "LP Df HDA3",
                        "LP Pp HDA1", "LP Pp HDA2", "LP Pp HDA3", "LP F/S/MH HDA1", "LP F/S/MHf HDA2", "LP F/S/MH HDA3")
write.csv(para.gev, file = "Table2_gev_para.csv")

### plot RMSE of carbon pools (currently using) ######
# carbon pools for each age class
for (i_fttppd in 1:length(fttppd_all)) {
  fttppd <- fttppd_all[i_fttppd];
  # observation from FIA
  data.obs <- readMat(paste0(fttppd, "_carbonpools.mat", sep=""))
  obs.mean <- data.obs$carbonpools.obs[,1:5] / 1000 # convert to kg
  obs.sd <- data.obs$carbonpools.obs[,6:10]/ 1000
  for (i_HD in 1:3) {
    data <- readMat(paste0(fttppd,"_modeled_C_stock_40_HD_", i_HD,".mat", sep = ""))
    curr_rmse <- data.frame(matrix(data = NA, ncol = 7, nrow = 500))
    colnames(curr_rmse) <- c("HDA","Foresttype", "AGB.rmse",  "CWD.rmse", "FFlitter.rmse",  "Snag.rmse","SOC.rmse")
    for (i in 1:500) {
      curr_rmse$AGB.rmse[i] <- sqrt(mean((obs.mean[,1] - data$mod.AGB.40[,i]/ 1000)^2, na.rm = T))
      curr_rmse$CWD.rmse[i] <- sqrt(mean((obs.mean[,2] - data$mod.CWD.40[,i]/ 1000)^2, na.rm = T))
      curr_rmse$FFlitter.rmse[i] <- sqrt(mean((obs.mean[,3] - data$mod.FFlitter.40[,i]/ 1000)^2, na.rm = T))
      curr_rmse$Snag.rmse[i] <- sqrt(mean((obs.mean[,4] - data$mod.Snag.40[,i]/ 1000)^2, na.rm = T))
      curr_rmse$SOC.rmse[i] <- sqrt(mean((obs.mean[,5] - data$mod.SOC.40[,i]/ 1000)^2, na.rm = T))
      curr_rmse$HDA[i] <-  i_HD
      curr_rmse$Foresttype[i] <- i_fttppd
    }
    if (i_fttppd == 1 & i_HD == 1) {
      Cstock_rmse <- curr_rmse
    } else {
      Cstock_rmse <- rbind(Cstock_rmse, curr_rmse)
    }
  }
}

# just for display
for (i_tmp in 1:4) {
  tmp_foresttype <- i_tmp
  tmp_hda <- 3
  # mean(Cstock_rmse$AGB.rmse[Cstock_rmse$HDA == tmp_hda & Cstock_rmse$Foresttype == tmp_foresttype])
  # mean(Cstock_rmse$Snag.rmse[Cstock_rmse$HDA == tmp_hda & Cstock_rmse$Foresttype == tmp_foresttype])
  # mean(Cstock_rmse$CWD.rmse[Cstock_rmse$HDA == tmp_hda & Cstock_rmse$Foresttype == tmp_foresttype])
  # print(mean(Cstock_rmse$FFlitter.rmse[Cstock_rmse$HDA == tmp_hda & Cstock_rmse$Foresttype == tmp_foresttype]))
  print(mean(Cstock_rmse$SOC.rmse[Cstock_rmse$HDA == tmp_hda & Cstock_rmse$Foresttype == tmp_foresttype]))
}

Cstock_rmse$HDA <- as.factor(Cstock_rmse$HDA)
Cstock_rmse$Foresttype <- as.factor(Cstock_rmse$Foresttype)
C.stocks.name <- c("AGB", "CWD", "FFlitter","Snag", "SOC")

# Function to perform ANOVA and Tukey HSD, then generate significance letters
generate_letters <- function(df) {
  ## test of mean
  # anova <- aov(current ~ HDA, data = df)
  # tukey <- TukeyHSD(anova)
  # cld <- multcompLetters4(anova, tukey)
  # Tk <- group_by(df, HDA) %>%
  #   summarise(mean=mean(current), quant = quantile(current, probs = 0.75)) %>%
  #   arrange(desc(mean))
  # # extracting the compact letter display and adding to the Tk table
  # cld <- as.data.frame.list(cld$HDA)
  # Tk$cld <- cld$Letters
  # # Extract the letters and ensure they are properly named
  # letters_df <- as.data.frame(cld$Letters)
  # letters_df$HDA <- rownames(letters_df)
  # letters_df$Foresttype <- unique(df$Foresttype)
  # rownames(letters_df) <- NULL
  # return(letters_df)
  
  ## test of distribution 
  wilcox_test_12 <- wilcox.test(current ~ HDA, data = df %>% filter(HDA %in% c(1, 2)))
  wilcox_test_13 <- wilcox.test(current ~ HDA, data = df %>% filter(HDA %in% c(1, 3)))
  wilcox_test_23 <- wilcox.test(current ~ HDA, data = df %>% filter(HDA %in% c(2, 3)))
  
  # Extract p-values and adjust for multiple comparisons using Bonferroni correction
  p_values <- c(wilcox_test_12$p.value, wilcox_test_13$p.value, wilcox_test_23$p.value)
  names(p_values) <- c("1-2", "1-3", "2-3")
  p_adjusted <- p.adjust(p_values, method = "bonferroni")
  
  dif3 <- p_adjusted < 0.05
  dif3L <- multcompLetters(dif3)
  # Create a matrix to store pairwise comparison results
  # comparison_matrix <- matrix(FALSE, nrow = 3, ncol = 1)
  # colnames(comparison_matrix) <- rownames(comparison_matrix) <- c("1", "2", "3")
  # 
  # # Fill the matrix with results of pairwise comparisons
  # comparison_matrix["1", "2"] <- p_adjusted["1-2"] < 0.05
  # comparison_matrix["1", "3"] <- p_adjusted["1-3"] < 0.05
  # comparison_matrix["2", "3"] <- p_adjusted["2-3"] < 0.05
  # 
  # Generate compact letter display based on the pairwise comparisons
  # cld <- multcompLetters(comparison_matrix)
  Tk <- group_by(df, HDA) %>%
    summarise(mean=mean(current), quant = quantile(current, probs = 0.75)) %>%
    arrange(desc(mean))
  
  # extracting the compact letter display and adding to the Tk table
  # cld <- as.data.frame.list(cld$HDA)
  Tk$cld <- dif3L$Letters
  # Extract the letters and ensure they are properly named
  letters_df <- as.data.frame(dif3L$Letters)
  letters_df$HDA <- rownames(letters_df)
  letters_df$Foresttype <- unique(df$Foresttype)
  rownames(letters_df) <- NULL
  return(letters_df)
}

# For figure 3 in the main text
for (i_stock in 1:5) {
  eval(parse(text = paste("Cstock_rmse$current <- Cstock_rmse$", C.stocks.name[i_stock], sep = "")))
  df <- Cstock_rmse

  # Apply the function to each forest type group and get letters
  pairwise_letters <- df %>%
    group_by(Foresttype) %>%
    do(generate_letters(.)) %>%
    ungroup()
  
  # Ensure the column names are correct
  colnames(pairwise_letters) <- c("Letter", "HDA", "Foresttype")
  
  # Check if 'Letter' column exists
  if (!"Letter" %in% colnames(pairwise_letters)) {
    stop("Letter column not found in pairwise_letters")
  }
  
  # Check the pairwise_letters data frame
  print(pairwise_letters)
  # Merge the letters with the original data for plotting
  df_letters <- df %>%
    left_join(pairwise_letters, by = c("Foresttype", "HDA"))
  
  # Custom colors for HDA
  custom_colors <- c("1" = "#0c84c6", "2" = "#ffbd66", "3" = "#f74d4d")
  
  p <- ggplot(df, aes(x = Foresttype, y = current, fill = HDA)) +
    geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.5,color = 'grey20', linewidth = 0.2, outlier.shape = NA, width = 0.9) +
    scale_fill_manual(values = custom_colors, labels = c("HDA1", "HDA2", "HDA3")) +  # Set custom colors and alpha
    xlab("Forest group type") + ylab(bquote(RMSE~(kg~C~m^-2)))  + 
    scale_x_discrete( labels=c("HP Df", "LP Df", "LP Pp", "LP F/S/MH")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 12),axis.text = element_text(size = 12, colour = "black"),legend.position = "none")
  p <- p + geom_text(data = df_letters, aes(x = Foresttype, y = 0, label = Letter, size = 5, color = HDA),
                     position = position_dodge(width = 0.75), vjust = 1) +scale_color_manual(values = custom_colors)

  if (i_stock ==1 ) { p <- p + ggtitle('(a) AGB') + ylim(0, 4) +
    theme(plot.title = element_text( hjust = 0.05, vjust = -10, size = 14))}
  if (i_stock ==2 ) { p <- p + ggtitle('(c) CWD') + ylim(0, 8) +
    theme(plot.title = element_text( hjust = 0.05, vjust = -10, size = 14))}
  if (i_stock ==3 ) { p <- p + ggtitle('(d) FF litter') + ylim(0, 6) +
    theme(plot.title = element_text( hjust = 0.05, vjust = -10, size = 14))}
  if (i_stock ==4 ) { p <- p + ggtitle('(b) Snag') + ylim(0, 6) +
    theme(plot.title = element_text( hjust = 0.05, vjust = -10, size = 14))}
  if (i_stock ==5 ) { p <- p + ggtitle('(e) SOC') + ylim(0, 20) +
    theme(plot.title = element_text( hjust =0.05, vjust = -10, size = 14))}
    # Add significant letters to the plot

  
  ggsave(paste0("Fig3_", C.stocks.name[i_stock], "_rmse.tiff"),plot = p, units="in", width=4, height=4, dpi=400)
  
  if (i_stock == 5) {
    p <- ggplot(df, aes(x = Foresttype, y = current, fill = HDA)) +
      geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.5,color = 'grey20', linewidth = 0.2, outlier.shape = NA, width = 0.9) +
      scale_fill_manual(values = custom_colors, labels = c("HDA1", "HDA2", "HDA3")) +  # Set custom colors and alpha
      xlab("Forest group type") + ylab(bquote(RMSE~(kg~C~m^-2)))  + 
      scale_x_discrete( labels=c("HP Df", "LP Df", "LP Pp", "LP F/S/MH")) +
      theme_bw() + theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank()) + 
      theme(text = element_text(size = 12),axis.text = element_text(size = 12, colour = "black"),legend.position = "none")
    # Add geom_text with a color aesthetic to create a legend
    p <- p + geom_text(data = df_letters, aes(x = Inf, y = 0, label = Letter, size = 5, color = HDA)) +
      theme(legend.position = "right", legend.text = element_text(size = 16),   # Increase legend text size
            legend.key.size = unit(3, "lines"))   # Increase size of legend key)
    ggsave(paste0("Fig3_", "_rmse_for_legend.tiff"),plot = p, units="in", width=4, height=4, dpi=400)
    
  }
  
}

### plot age-carbon pools (for supplement) ######
# For figure S1 in the supplement 
# plot carbon pools
for (i_fttppd in 1:length(fttppd_all)) {
  fttppd <- fttppd_all[i_fttppd];
  # FIA observation
  data <- readMat(paste0(fttppd, "_carbonpools.mat", sep=""))
  obs.mean <- data$carbonpools.obs[,1:5] / 1000
  obs.sd <- data$carbonpools.obs[,6:10]/ 1000
  
  carbonpool <- data.frame(matrix(data = NA, ncol = (5*5+2), nrow = (40*3)))
  #provide column names
  colnames(carbonpool) <- c("HDA", "Standage", "AGB", "CWD",   "FFlitter","Snag", "SOC",
                            "AGBmax", "CWDmax",  "FFlittermax", "Snagmax", "SOCmax",
                            "AGBmin",  "CWDmin",  "FFlittermin","Snagmin", "SOCmin",
                            "AGBobs",  "CWDobs",  "FFlitterobs", "Snagobs","SOCobs",
                            "AGBobssd",   "CWDobssd",  "FFlitterobssd","Snagobssd", "SOCobssd")
  carbonpool[1:40, 18:22] <- obs.mean
  carbonpool[1:40, 23:27] <- obs.sd
  carbonpool$HDA[1:40] <- 1
  carbonpool$HDA[41:80] <- 2
  carbonpool$HDA[81:120] <- 3
  carbonpool$Standage[1:40] <- seq(from = 2.5, to = 200, by = 5)
  carbonpool$Standage[41:80] <- seq(from = 2.5, to = 200, by = 5)
  carbonpool$Standage[81:120] <- seq(from = 2.5, to = 200, by = 5)
  # HDA 1
  data <- readMat(paste0(fttppd,"_modeled_C_stock_40_HD_", 1,".mat", sep = ""))
  carbonpool$AGB[carbonpool$HDA ==1] <- rowMeans(data$mod.AGB.40)/1000
  carbonpool$AGBmax[carbonpool$HDA ==1] <- rowMaxs(data$mod.AGB.40)/1000
  carbonpool$AGBmin[carbonpool$HDA ==1] <- rowMins(data$mod.AGB.40)/1000
  carbonpool$CWD[carbonpool$HDA ==1] <- rowMeans(data$mod.CWD.40)/1000
  carbonpool$CWDmax[carbonpool$HDA ==1] <- rowMaxs(data$mod.CWD.40)/1000
  carbonpool$CWDmin[carbonpool$HDA ==1] <- rowMins(data$mod.CWD.40)/1000
  carbonpool$FFlitter[carbonpool$HDA ==1] <- rowMeans(data$mod.FFlitter.40)/1000
  carbonpool$FFlittermax[carbonpool$HDA ==1] <- rowMaxs(data$mod.FFlitter.40)/1000
  carbonpool$FFlittermin[carbonpool$HDA ==1] <- rowMins(data$mod.FFlitter.40)/1000
  carbonpool$Snag[carbonpool$HDA ==1] <- rowMeans(data$mod.Snag.40)/1000
  carbonpool$Snagmax[carbonpool$HDA ==1] <- rowMaxs(data$mod.Snag.40)/1000
  carbonpool$Snagmin[carbonpool$HDA ==1] <- rowMins(data$mod.Snag.40)/1000
  carbonpool$SOC[carbonpool$HDA ==1] <- rowMeans(data$mod.SOC.40)/1000
  carbonpool$SOCmax[carbonpool$HDA ==1] <- rowMaxs(data$mod.SOC.40)/1000
  carbonpool$SOCmin[carbonpool$HDA ==1] <- rowMins(data$mod.SOC.40)/1000
  # HDA 2
  data <- readMat(paste0(fttppd,"_modeled_C_stock_40_HD_", 2,".mat", sep = ""))
  carbonpool$AGB[carbonpool$HDA ==2] <- rowMeans(data$mod.AGB.40)/1000
  carbonpool$AGBmax[carbonpool$HDA ==2] <- rowMaxs(data$mod.AGB.40)/1000
  carbonpool$AGBmin[carbonpool$HDA ==2] <- rowMins(data$mod.AGB.40)/1000
  carbonpool$CWD[carbonpool$HDA ==2] <- rowMeans(data$mod.CWD.40)/1000
  carbonpool$CWDmax[carbonpool$HDA ==2] <- rowMaxs(data$mod.CWD.40)/1000
  carbonpool$CWDmin[carbonpool$HDA ==2] <- rowMins(data$mod.CWD.40)/1000
  carbonpool$FFlitter[carbonpool$HDA ==2] <- rowMeans(data$mod.FFlitter.40)/1000
  carbonpool$FFlittermax[carbonpool$HDA ==2] <- rowMaxs(data$mod.FFlitter.40)/1000
  carbonpool$FFlittermin[carbonpool$HDA ==2] <- rowMins(data$mod.FFlitter.40)/1000
  carbonpool$Snag[carbonpool$HDA ==2] <- rowMeans(data$mod.Snag.40)/1000
  carbonpool$Snagmax[carbonpool$HDA ==2] <- rowMaxs(data$mod.Snag.40)/1000
  carbonpool$Snagmin[carbonpool$HDA ==2] <- rowMins(data$mod.Snag.40)/1000
  carbonpool$SOC[carbonpool$HDA ==2] <- rowMeans(data$mod.SOC.40)/1000
  carbonpool$SOCmax[carbonpool$HDA ==2] <- rowMaxs(data$mod.SOC.40)/1000
  carbonpool$SOCmin[carbonpool$HDA ==2] <- rowMins(data$mod.SOC.40)/1000
  # HDA 3
  data <- readMat(paste0(fttppd,"_modeled_C_stock_40_HD_", 3,".mat", sep = ""))
  carbonpool$AGB[carbonpool$HDA ==3] <- rowMeans(data$mod.AGB.40)/1000
  carbonpool$AGBmax[carbonpool$HDA ==3] <- rowMaxs(data$mod.AGB.40)/1000
  carbonpool$AGBmin[carbonpool$HDA ==3] <- rowMins(data$mod.AGB.40)/1000
  carbonpool$CWD[carbonpool$HDA ==3] <- rowMeans(data$mod.CWD.40)/1000
  carbonpool$CWDmax[carbonpool$HDA ==3] <- rowMaxs(data$mod.CWD.40)/1000
  carbonpool$CWDmin[carbonpool$HDA ==3] <- rowMins(data$mod.CWD.40)/1000
  carbonpool$FFlitter[carbonpool$HDA ==3] <- rowMeans(data$mod.FFlitter.40)/1000
  carbonpool$FFlittermax[carbonpool$HDA ==3] <- rowMaxs(data$mod.FFlitter.40)/1000
  carbonpool$FFlittermin[carbonpool$HDA ==3] <- rowMins(data$mod.FFlitter.40)/1000
  carbonpool$Snag[carbonpool$HDA ==3] <- rowMeans(data$mod.Snag.40)/1000
  carbonpool$Snagmax[carbonpool$HDA ==3] <- rowMaxs(data$mod.Snag.40)/1000
  carbonpool$Snagmin[carbonpool$HDA ==3] <- rowMins(data$mod.Snag.40)/1000
  carbonpool$SOC[carbonpool$HDA ==3] <- rowMeans(data$mod.SOC.40)/1000
  carbonpool$SOCmax[carbonpool$HDA ==3] <- rowMaxs(data$mod.SOC.40)/1000
  carbonpool$SOCmin[carbonpool$HDA ==3] <- rowMins(data$mod.SOC.40)/1000
  
  carbonpool[carbonpool < 0] <- 0
  
  carbonpool$HDA <- as.factor(carbonpool$HDA)
  custom_colors <- c("1" = "#0c84c6", "2" = "#ffbd66", "3" = "#f74d4d")
  
  # write.csv(carbonpool, file = paste0(fttppd,"_carbonpools.csv"))
  ggplot(data = carbonpool, aes(x = Standage, group = HDA)) + 
    geom_line(aes(y = AGB, color = HDA), size = 1) + geom_point(aes(y = AGB, color = HDA)) +
    geom_point(aes(y = AGBobs), color = 'grey20') +
    geom_errorbar(aes(ymin = AGBobs-AGBobssd, ymax = AGBobs+AGBobssd), width=.2, color = 'grey20')+
    scale_color_manual(values=custom_colors) +
    geom_ribbon(aes(y = AGB, ymin = AGBmin, ymax = AGBmax, fill = HDA), alpha = .2) +
    scale_fill_manual(values=custom_colors) +
    xlab("Years") + ylab(bquote(AGB~(kg~C~m^-2)))  + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    ggtitle(paste0("(a", i_fttppd, ")", sep = "")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 14),axis.text = element_text(size = 14, colour = "black")) + theme(plot.title = element_text( hjust = 0.05, vjust = - 8, size = 14)) 
  ggsave(paste0("FigS1_", fttppd, "_1AGB.tiff"), units="in", width=5, height=4, dpi=500)
  
  ggplot(data = carbonpool, aes(x = Standage, group = HDA)) + 
    geom_line(aes(y = Snag, color = HDA), size = 1) + geom_point(aes(y = Snag, color = HDA)) +
    geom_point(aes(y = Snagobs), color = 'grey20') +
    geom_errorbar(aes(ymin = Snagobs-Snagobssd, ymax = Snagobs+Snagobssd), width=.2, color = 'grey20')+
    scale_color_manual(values=custom_colors) +
    geom_ribbon(aes(y = Snag, ymin = Snagmin, ymax = Snagmax, fill = HDA), alpha = .2) +
    scale_fill_manual(values=custom_colors) +
    xlab("Years") + ylab(bquote(Snag~(kg~C~m^-2)))  + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    ggtitle(paste0("(b", i_fttppd, ")", sep = "")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 14),axis.text = element_text(size = 14, colour = "black")) + theme(plot.title = element_text( hjust = 0.05, vjust = - 8, size = 14)) 
  ggsave(paste0("FigS1_", fttppd, "_2Snag.tiff"), units="in", width=5, height=4, dpi=500)
  # Coarse woody debris
  ggplot(data = carbonpool, aes(x = Standage, group = HDA)) + 
    geom_line(aes(y = CWD, color = HDA), size = 1) + geom_point(aes(y = CWD, color = HDA)) +
    geom_point(aes(y = CWDobs), color = 'grey20') +
    geom_errorbar(aes(ymin = CWDobs-CWDobssd, ymax = CWDobs+CWDobssd), width=.2, color = 'grey20')+
    scale_color_manual(values=custom_colors) +
    geom_ribbon(aes(y = CWD, ymin = CWDmin, ymax = CWDmax, fill = HDA), alpha = .2) +
    scale_fill_manual(values=custom_colors) +
    xlab("Years") + ylab(bquote(CWD~(kg~C~m^-2)))  + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    ggtitle(paste0("(c", i_fttppd, ")", sep = "")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 14),axis.text = element_text(size = 14, colour = "black")) + theme(plot.title = element_text( hjust = 0.05, vjust = - 8, size = 14)) 
  ggsave(paste0("FigS1_", fttppd, "_3CWD.tiff"), units="in", width=5, height=4, dpi=500)
  
  ggplot(data = carbonpool, aes(x = Standage, group = HDA)) + 
    geom_line(aes(y = FFlitter, color = HDA), size = 1) + geom_point(aes(y = FFlitter, color = HDA)) +
    geom_point(aes(y = FFlitterobs), color = 'grey20') +
    geom_errorbar(aes(ymin = FFlitterobs-FFlitterobssd, ymax = FFlitterobs+FFlitterobssd), width=.2, color = 'grey20')+
    scale_color_manual(values=custom_colors) +
    geom_ribbon(aes(y = FFlitter, ymin = FFlittermin, ymax = FFlittermax, fill = HDA), alpha = .2) +
    scale_fill_manual(values=custom_colors) +
    xlab("Years") + ylab(bquote(FF~litter~(kg~C~m^-2)))  + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    ggtitle(paste0("(d", i_fttppd, ")", sep = "")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 14),axis.text = element_text(size = 14, colour = "black")) + theme(plot.title = element_text( hjust = 0.05, vjust = - 8, size = 14)) 
  ggsave(paste0("FigS1_", fttppd, "_4FFlitter.tiff"), units="in", width=5, height=4, dpi=500)
  
  ggplot(data = carbonpool, aes(x = Standage, group = HDA)) + 
    geom_line(aes(y = SOC, color = HDA), size = 1) + geom_point(aes(y = SOC, color = HDA)) +
    geom_point(aes(y = SOCobs), color = 'grey20') +
    geom_errorbar(aes(ymin = SOCobs-SOCobssd, ymax = SOCobs+SOCobssd), width=.2, color = 'grey20')+
    scale_color_manual(values=custom_colors) +
    geom_ribbon(aes(y = SOC, ymin = SOCmin, ymax = SOCmax, fill = HDA), alpha = .2) +
    scale_fill_manual(values=custom_colors) +
    xlab("Years") + ylab(bquote(SOC~(kg~C~m^-2)))  + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    ggtitle(paste0("(e", i_fttppd, ")", sep = "")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 14),axis.text = element_text(size = 14, colour = "black")) + theme(plot.title = element_text( hjust = 0.05, vjust = - 8, size = 14)) 
  ggsave(paste0("FigS1_", fttppd, "_5SOC.tiff"), units="in", width=5, height=4, dpi=500)
}

### plot age-carbon fluxes under different HDA (currently using) ######
rm(Cfluxes)
labels <- c("HP Df", "LP Df", "LP Pp", "LP F/S/MH")
for (i_fttppd in 1:length(fttppd_all)) {
  fttppd <- fttppd_all[i_fttppd];

  for (i_HD in 1:3) {
    data <- readMat(paste0(fttppd,"_modeled_C_HD_", i_HD,".mat", sep = ""))
    curr_flux <- data.frame(matrix(data = NA, ncol = 12, nrow = 200)) #500*200
    colnames(curr_flux) <- c("HDA","Foresttype","Standage", "NPP",  "Rh", "NEP", "NPPmin",  "Rhmin", "NEPmin", "NPPmax",  "Rhmax", "NEPmax")

    for (i_age in 1:200) {
      curr_flux$NPP[i_age] <- median(data$mod.NPP[i_age,],na.rm = T)
      curr_flux$Rh[i_age] <-  median(data$mod.Rh[i_age,],na.rm = T)
      curr_flux$NEP[i_age] <- median(data$mod.NEP[i_age,],na.rm = T)
      curr_flux$NPPmin[i_age] <- quantile(data$mod.NPP[i_age,],na.rm = T,probs = c(0.25))
      curr_flux$Rhmin[i_age] <-  quantile(data$mod.Rh[i_age,],na.rm = T,probs = c(0.25))
      curr_flux$NEPmin[i_age] <- quantile(data$mod.NEP[i_age,],na.rm = T,probs = c(0.25))
      curr_flux$NPPmax[i_age] <- quantile(data$mod.NPP[i_age,],na.rm = T,probs = c(0.75))
      curr_flux$Rhmax[i_age] <-  quantile(data$mod.Rh[i_age,],na.rm = T,probs = c(0.75))
      curr_flux$NEPmax[i_age] <- quantile(data$mod.NEP[i_age,],na.rm = T,probs = c(0.75))
      curr_flux$HDA[i_age] <-  i_HD
      curr_flux$Foresttype[i_age] <- i_fttppd
      curr_flux$Standage[i_age] <- i_age
    }
    if (i_fttppd == 1 & i_HD == 1) {
      Cfluxes <- curr_flux
    } else {
      Cfluxes <- rbind(Cfluxes, curr_flux)
    }
  }
}
Cfluxes$HDA <- factor(Cfluxes$HDA, levels = c(1, 2, 3))
Cfluxes$Foresttype <- as.factor(Cfluxes$Foresttype)

for (i_tmp in 1:4) {
  tmp_foresttype <- i_tmp
  tmp_hda <- 1
  print(mean(Cfluxes$NEP[Cfluxes$HDA == tmp_hda & Cfluxes$Foresttype == tmp_foresttype & Cfluxes$Standage >= 170]))
}

# Custom colors for HDA
custom_colors <- c("1" = "#0c84c6", "2" = "#ffbd66", "3" = "#f74d4d")

# For figure 5 in the main text
for (i in 1:4) {
  fttppd <- fttppd_all[i]
  df <- Cfluxes[Cfluxes$Foresttype == i,]
  
  ggplot(data = df, aes(x = Standage, group = HDA)) +  
    geom_line(aes(y = NEP, color = HDA), size = 1) + scale_color_manual(values=custom_colors) +
    geom_ribbon(aes(y = NEP, ymin =NEPmin, ymax = NEPmax, fill = HDA), alpha = .2) +
    scale_fill_manual(values=custom_colors) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black") +
    xlab("Years") + ylab(bquote(NEP~(g~C~m^-2)))  + ylim(-2000, 2000) +
    # scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    ggtitle(paste0("(c", i, ") ", labels[i], sep = "")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 14),axis.text = element_text(size = 14, colour = "black")) + theme(plot.title = element_text( hjust = 0.05, vjust = - 8, size = 14)) 
  
  ggsave(paste0("Fig5_", fttppd, "_NEP.tiff"), units="in", width=5, height=4, dpi=400)
  
  ggplot(data = df, aes(x = Standage, group = HDA)) +  
    geom_line(aes(y = NPP, color = HDA), size = 1) + scale_color_manual(values=custom_colors) +
    geom_ribbon(aes(y = NPP, ymin = NPPmin, ymax = NPPmax, fill = HDA), alpha = .2) +
    scale_fill_manual(values=custom_colors) +
    xlab("Years") +  ylab(bquote(NPP~(g~C~m^-2)))   +  ylim(0, 6500) +
    # scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    ggtitle(paste0("(a", i, ") ", labels[i], sep = "")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 14),axis.text = element_text(size = 14, colour = "black")) + theme(plot.title = element_text( hjust = 0.05, vjust = - 8, size = 14)) 
  
  ggsave(paste0("Fig5_", fttppd, "_NPP.tiff"), units="in", width=5, height=4, dpi=400)
  
  ggplot(data = df, aes(x = Standage, group = HDA)) +  
    geom_line(aes(y = Rh, color = HDA), size = 1) + scale_color_manual(values=custom_colors) +
    geom_ribbon(aes(y = Rh, ymin = Rhmin, ymax = Rhmax, fill = HDA), alpha = .2) +
    scale_fill_manual(values=custom_colors) +
    xlab("Years") +  ylab(bquote(Rh~(g~C~m^-2)))   + ylim(0, 6500) +
    # scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    ggtitle(paste0("(b", i, ") ", labels[i], sep = "")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 14),axis.text = element_text(size = 14, colour = "black")) + theme(plot.title = element_text( hjust = 0.05, vjust = - 8, size = 14)) 
  
  ggsave(paste0("Fig5_", fttppd, "_Rh.tiff"), units="in", width=5, height=4, dpi=400)
  
}
