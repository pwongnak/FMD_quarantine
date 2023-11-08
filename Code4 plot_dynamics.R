library(stringr)
###########################################################################
fileList <- list.files("mainResults_2023_2")
outputALL <- fileList[grep("outputALL", fileList)]
###########################################################################
settings_market <- do.call(rbind.data.frame, str_split(outputALL, "_"))[,c(3,5,7,8)]
colnames(settings_market) <- c("Imm", "r", "vac", "nWeek")
settings_market$nWeek <- str_extract(settings_market$nWeek, "[0-9]")
settings_market <- as.data.frame(apply(settings_market, 2,as.numeric))
settings_market <- settings_market[settings_market$nWeek == 3 &
                  settings_market$r %in% c(0, 0.5, 1) &
                  settings_market$vac %in% c(0, 0.5, 1) &
                  settings_market$Imm %in% c(0, 0.5, 1),]
rowname_list <- as.numeric(row.names(settings_market))
###########################################################################
SUMMARY <- NULL

for(i in 1:nrow(settings_market)){
  dat_all <- read.csv(paste0("mainResults_2023_2/", outputALL[rowname_list[i]]))
  I_ind <- grepl("I", colnames(dat_all))
  dat_all$I_all <- rowSums(dat_all[,I_ind])
  dat_all$I_prop <- dat_all$I_all/dat_all$N

  med_data <- aggregate(list("I_prop" = dat_all$I_prop), by = list("day" = dat_all$day), FUN = median)
  med_data$r <- settings_market$r[i]
  med_data$vac <- settings_market$vac[i]
  med_data$Imm <- settings_market$Imm[i]
  
  SUMMARY <- rbind(SUMMARY, med_data)
}


SUMMARY$vac <- as.factor(SUMMARY$vac)
levels(SUMMARY$vac) <-   c("e[v]==0.0", "e[v]==0.5", 
                             "e[v]==1.0")

SUMMARY$Imm <- as.factor(SUMMARY$Imm)
levels(SUMMARY$Imm) <-   c("P[R(0)]==0.0", "P[R(0)]==0.5", 
                           "P[R(0)]==1.0")

G <- ggplot(SUMMARY, aes(x = day, y = I_prop, col = as.factor(r))) +
  geom_line(linewidth = 0.7) +
  facet_grid(vac ~ Imm, labeller = label_parsed) +
  theme_bw() +
  scale_color_manual(values = rev(c("#EB455F", "#F49D1A", "#227C70")), 
                     name ="Reduction factor for\nbetween-batch \ntransmission") +
  xlab("Time step (days)") +
  ylab("Proportion of infectious cattle within the quarantine center") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  geom_vline(xintercept = 5*7, col = "red", linetype = "dashed", linewidth = 0.5)


png("Fig S1.png", width = 12, height = 6, units = "in", res = 400)
G
dev.off()
