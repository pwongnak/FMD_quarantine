library(stringr)
library(ggplot2)
library(ggpubr)
###########################################################################
fileList <- list.files("mainResults_2023_2")
marketList <- fileList[grep("marketList", fileList)]
###########################################################################
settings_market <- do.call(rbind.data.frame, str_split(marketList, "_"))[,c(3,5,7,8)]
colnames(settings_market) <- c("Imm", "r", "vac", "nWeek")
settings_market$nWeek <- str_extract(settings_market$nWeek, "[0-9]")
settings_market <- as.data.frame(apply(settings_market, 2,as.numeric))
###########################################################################
SUMMARY <- list()
for(i in 1:nrow(settings_market)){
  dat_all <- read.csv(paste0("mainResults_2023_2/", marketList[i]))
  Summary <- NULL
  for(j in 1:max(dat_all$rep)){
    dat <- dat_all[dat_all$rep == j,]
    dat <- dat[-which(dat$week %in% 1:5 ),]
    
    flagS <- which(str_detect(colnames(dat), "S"))
    flagE <- which(str_detect(colnames(dat), "E"))
    flagIA <- which(str_detect(colnames(dat), "IA"))
    flagIS <- which(str_detect(colnames(dat), "Is"))
    flagR <- which(str_detect(colnames(dat), "R"))

    dat$S <- dat[,flagS]
    dat$E <- dat[,flagE] 
    dat$IA <- rowSums(dat[,flagIA])
    dat$IS <- dat[,flagIS]
    dat$R <- dat[,flagR]

    risk <- (dat$IA)/((dat$S + dat$E + dat$IA + dat$R))
    prot <- (dat$R)/((dat$S + dat$E + dat$IA + dat$R))
    erad  <- (dat$IS)/((dat$S + dat$E + dat$IA + dat$R + dat$IS))

    weekNotSafe <- sum(!(dat$E == 0 & dat$IA == 0))/nrow(dat)

    summary <- data.frame("risk" = risk,
                          "prot" = prot,
                          "erad" = erad,
                          "weekNotSafe" = weekNotSafe)
    Summary <- rbind(Summary, summary)
}

SUMMARY[[i]] <- as.data.frame(t(apply(Summary,2,quantile, c(0.025, 0.5, 0.975))))
}
SUMMARY[[1]]

Risk <- do.call(rbind.data.frame, (lapply(SUMMARY, "[", 1,)))
Prot <- do.call(rbind.data.frame, (lapply(SUMMARY, "[", 2,)))
Erad <- do.call(rbind.data.frame, (lapply(SUMMARY, "[", 3,)))
WeekNotSafe <- do.call(rbind.data.frame, (lapply(SUMMARY, "[", 4,)))

summary_all <- data.frame(settings_market, Risk, Prot, Erad, WeekNotSafe)
colnames(summary_all)[5:16] <- c("risk_low", "risk_med", "risk_up",
                                 "prot_low", "prot_med", "prot_up",
                                 "erad_low", "erad_med", "erad_up",
                                 "weeknotsafe_low", "weeknotsafe_med", "weeknotsafe_up")

write.csv(summary_all, "summarize_results_new.csv", row.names = F)


library(ggplot2)
colnames(summary_all)

data_plot <- summary_all[summary_all$r %in% c(0,0.5,1) & summary_all$vac %in% c(0,0.5,1),]
data_plot$r <- as.factor(data_plot$r)

data_plot$nWeek <- as.factor(data_plot$nWeek)
levels(data_plot$nWeek) <- c("t[max]==1~week", "t[max]==2~weeks", 
                             "t[max]==3~weeks", "t[max]==4~weeks")

data_plot$vac <- as.factor(data_plot$vac)
levels(data_plot$vac) <-   c("e[v]==0.0", "e[v]==0.5", 
                             "e[v]==1.0")

name <- "Reduction factor for\nbetween-batch \ntransmission"

G1 <- ggplot(data_plot, aes(x = Imm, y = risk_med, gr = r)) +
  geom_ribbon(aes(ymin = risk_low, ymax = risk_up, fill = r), alpha = 0.15) +
  geom_line(aes(col = r), alpha = 0.75) +
  geom_point(aes(col = r, shape = r), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#537FE7", "#183A1D", "#CD0404"), name =name) +
  scale_fill_manual(values = c("#537FE7", "#183A1D", "#CD0404"), name =name) +
  scale_shape_manual(values = c(16, 17, 15), name =name) +
  facet_grid(vac~nWeek, labeller = label_parsed) +
  ylim(0,1) +
  theme_bw() +
  geom_hline(yintercept = c(0,0.5,1) , linetype = "dashed", linewidth = 0.2) +
  xlab(bquote("Proportion of non-infected cattle with immunity at arrival, P"[R(0)])) +
  ylab(bquote("Risk of exported cattle being infected, P"[cattle]))  +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  ggtitle(bquote("Risk of exported cattle being infected, P"[cattle]))

G2 <- ggplot(data_plot, aes(x = Imm, y = weeknotsafe_med, gr = r)) +
  geom_ribbon(aes(ymin = weeknotsafe_low, ymax = weeknotsafe_up, fill = r), alpha = 0.15) +
  geom_line(aes(col = r), alpha = 0.75) +
  geom_point(aes(col = r, shape = r), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#537FE7", "#183A1D", "#CD0404"), name =name) +
  scale_fill_manual(values = c("#537FE7", "#183A1D", "#CD0404"), name =name) +
  scale_shape_manual(values = c(16, 17, 15), name =name) +
  facet_grid(vac~nWeek, labeller = label_parsed) +
  ylim(0,1) +
  theme_bw() +
  geom_hline(yintercept = c(0,0.5,1) , linetype = "dashed", linewidth = 0.2) +
  xlab(bquote("Proportion of non-infected cattle with immunity at arrival, P"[R(0)])) +
  ylab(bquote("Risk of weeks exporting infected cattle, P"[week]))  +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  ggtitle(bquote("Risk of weeks exporting infected cattle, P"[week]))

G3 <- ggplot(data_plot, aes(x = Imm, y = erad_med, gr = r)) +
  geom_ribbon(aes(ymin = erad_low, ymax = erad_up, fill = r), alpha = 0.15) +
  geom_line(aes(col = r), alpha = 0.75) +
  geom_point(aes(col = r, shape = r), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#537FE7", "#183A1D", "#CD0404"), name =name) +
  scale_fill_manual(values = c("#537FE7", "#183A1D", "#CD0404"), name =name) +
  scale_shape_manual(values = c(16, 17, 15), name =name) +
  facet_grid(vac~nWeek, labeller = label_parsed) +
  ylim(0,1) +
  theme_bw() +
  geom_hline(yintercept = c(0,0.5,1) , linetype = "dashed", linewidth = 0.2) +
  xlab(bquote("Proportion of non-infected cattle with immunity at arrival, P"[R(0)])) +
  ylab(bquote("Risk of removing symptomatic cattle before release, P"[remove]))  +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  ggtitle(bquote("Risk of removing symptomatic cattle before release, P"[remove]))


G4 <- ggplot(data_plot, aes(x = Imm, y = prot_med, gr = r)) +
  geom_ribbon(aes(ymin = prot_low, ymax = prot_up, fill = r), alpha = 0.15) +
  geom_line(aes(col = r), alpha = 0.75) +
  geom_point(aes(col = r, shape = r), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#537FE7", "#183A1D", "#CD0404"), name =name) +
  scale_fill_manual(values = c("#537FE7", "#183A1D", "#CD0404"), name =name) +
  scale_shape_manual(values = c(16, 17, 15), name =name) +
  facet_grid(vac~nWeek, labeller = label_parsed) +
  ylim(0,1) +
  theme_bw() +
  geom_hline(yintercept = c(0,0.5,1) , linetype = "dashed", linewidth = 0.2) +
  xlab(bquote("Proportion of non-infected cattle with immunity at arrival, P"[R(0)])) +
  ylab(bquote("Probability of releasing immunized cattle, P"[immunized]))  +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  ggtitle(bquote("Probability of releasing immunized cattle, P"[immunized]))

png("Fig4.png", width = 11, height = 11, units = "in", res = 400)
ggarrange(G1, G2, ncol = 1, nrow = 2, labels = "AUTO")
dev.off()

png("Fig5.png", width = 11, height = 11, units = "in", res = 400)
ggarrange(G3, G4, ncol = 1, nrow = 2, labels = "AUTO")
dev.off()

#####################################################################################################################
data_plot2 <- summary_all[summary_all$vac %in% c(0,0.5,1),]

data_plot2$r <- as.numeric(data_plot2$r)

data_plot2$nWeek <- as.factor(data_plot2$nWeek)
levels(data_plot2$nWeek) <- c("t[max]==1~week", "t[max]==2~weeks", 
                             "t[max]==3~weeks", "t[max]==4~weeks")

data_plot2$vac <- as.factor(data_plot2$vac)
levels(data_plot2$vac) <-   c("e[v]==0.0", "e[v]==0.5", 
                             "e[v]==1.0")

S1 <- ggplot(data_plot2, aes(x = Imm, y = r)) +
  geom_tile(aes(fill = risk_med)) +
  facet_grid(vac~nWeek, labeller = label_parsed) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high="red",
                         midpoint = 0.25,
                       limits = c(0,0.5),
                       name = "Values") +
  xlab(bquote("Proportion of non-infected cattle with immunity at arrival, P"[R(0)])) +
  ylab(bquote("Reduction factor for between-batch  transmission, r"))  +
  ggtitle(bquote("Risk of exported cattle being infected, P"[cattle])) +
  theme(strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) 

S2 <- ggplot(data_plot2, aes(x = Imm, y = r)) +
  geom_tile(aes(fill = weeknotsafe_med)) +
  facet_grid(vac~nWeek, labeller = label_parsed) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high="red",
                       midpoint = 0.5,
                       limits = c(0,1),
                       name = "Values") +
  xlab(bquote("Proportion of non-infected cattle with immunity at arrival, P"[R(0)])) +
  ylab(bquote("Reduction factor for between-batch  transmission, r"))  +
  ggtitle(bquote("Risk of weeks exporting infected cattle, P"[week])) +
theme(strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) 

S3 <- ggplot(data_plot2, aes(x = Imm, y = r)) +
  geom_tile(aes(fill = erad_med)) +
  facet_grid(vac~nWeek, labeller = label_parsed) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high="red",
                       midpoint = 0.15,
                       limits = c(0,0.3),
                       name = "Values") +
  xlab(bquote("Proportion of non-infected cattle with immunity at arrival, P"[R(0)])) +
  ylab(bquote("Reduction factor for between-batch  transmission, r"))  +
  ggtitle(bquote("Risk of removing symptomatic cattle before release, P"[remove])) +
theme(strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) 


S4 <- ggplot(data_plot2, aes(x = Imm, y = r)) +
  geom_tile(aes(fill = prot_med)) +
  facet_grid(vac~nWeek, labeller = label_parsed) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high="blue",
                       midpoint = 0.5,
                       limits = c(0,1),
                       name = "Values") +
  xlab(bquote("Proportion of non-infected cattle with immunity at arrival, P"[R(0)])) +
  ylab(bquote("Reduction factor for between-batch  transmission, r"))  +
  ggtitle(bquote("Probability of releasing immunized cattle, P"[immunized])) +
theme(strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) 

png("FigS3.png", width = 11, height = 11, units = "in", res = 400)
ggarrange(S1, S2, ncol = 1, nrow = 2, labels = "AUTO")
dev.off()

png("FigS4.png", width = 11, height = 11, units = "in", res = 400)
ggarrange(S3, S4, ncol = 1, nrow = 2, labels = "AUTO")
dev.off()
#####################################################################################################################
summary_all
