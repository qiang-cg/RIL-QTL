setwd("E:\\2018-2024-博士\\RIL相关分析\\RIL426_final_results\\表型数据分析\\231016")
library(ggplot2)
library(cowplot)
library(patchwork)
library(svglite)
library(RColorBrewer)
library(plotrix)

RIL <- read.csv("RIL_426_lines_pheno.csv", head=T)
par <- read.csv("231016-RIL_NR_parental_phenotype.csv", head=T)
par$CD <- par$CD*0.1

t.test(data=par, AWL~ID)

n_phe <- subset(par, ID=="NIV")
r_phe <- subset(par, ID=="RUF")

par_stat_lt <- list()
for(tr in names(par)[-1]){
  n <- subset(n_phe, select=tr)
  r <- subset(r_phe, select=tr)
  test <- t.test(r, n)
  par_stat_lt[[tr]] <- data.frame(n_indn=length(which(n[,1]!="NA")), 
                                  r_indn=length(which(r[,1]!="NA")),
                                  mean_n=test$estimate[1], mean_r=test$estimate[2],
                                  sd_n=sd(n[,1], na.rm=T), sd_r=sd(r[,1], na.rm=T),
                                  t.statistic=test$statistic, p.value=test$p.value)
}
par_stat_df <- do.call(rbind,par_stat_lt)
par_stat_df[,3:7] <- round(par_stat_df[,3:7],2)
apply(par_stat_df,1,function(x){
  sig <- ifelse(x[8] > 0.05, "", 
           ifelse(x[8] <= 0.05 & x[8] > 0.01, "*",
             ifelse(x[8] <= 0.01 & x[8] > 0.001, "**",
               ifelse(x[8] <= 0.001, "***"))))
  return(sig)
}) -> sig_lev
par_stat_df$sig_level <- sig_lev

write.csv(par_stat_df, file="231217_par_statistic.csv", quote=F)

#########
pheno_name <- c("ANL","FH","PE","AWL","CD","CH","CL","FLA","FLL","FLW","GL",
                "GWE","GWI","PL","PS","SN")
RIL_pheno <- RIL[,pheno_name]

######
shp.test <- apply(RIL_pheno, 2, function(x){
  st <- shapiro.test(x)
  W_sta <- st$statistic
  P_value <- st$p.value
  st_res <- c(W=W_sta, P=P_value)
  return(st_res)
})
round(shp.test, 4)

#########
###1. histogram plot for phenotypes

count_max <- max(hist(RIL_pheno$ANL, breaks=20, plot=F)$counts)

sd_list <- list()
for(p in pheno_name){  
  n.se <- sd(n_phe[,p], na.rm=T)
  r.se <- sd(r_phe[,p], na.rm=T)
  n.mean <- mean(n_phe[,p], na.rm=T)
  r.mean <- mean(r_phe[,p], na.rm=T)
  sd_list[[p]] <- round(data.frame(N_SD=n.se, N_AVG=n.mean, R_SD=r.se, R_AVG=r.mean),4)
}  

sz=12
lysz=0.8
N_col="#990000" 
R_col="#003399"
##########################
df <- sd_list$ANL
anl <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=ANL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=ct_max-ct_max/10, yend=ct_max+ct_max/10, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=ct_max, yend=ct_max, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=ct_max-ct_max/10, yend=ct_max+ct_max/10, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=ct_max, yend=ct_max, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$ANL), max(RIL_pheno$ANL), length.out = 4),1)) + 
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

ct_max <- max(hist(RIL_pheno$FH, breaks=20, plot=F)$counts)
df <- sd_list$FH
fh <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=FH), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=75, yend=85, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=80, yend=80, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=75, yend=85, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=80, yend=80, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$FH), max(RIL_pheno$FH), length.out = 4),0)) + 
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$PE
pe <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=PE), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=82, yend=96, linewidth=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=89, yend=89, linewidth=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=82, yend=96, linewidth=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=89, yend=89, linewidth=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$PE), max(RIL_pheno$PE), length.out = 4),0)) + 
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$AWL
awl <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=AWL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=62, yend=74, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=68, yend=68, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=54, yend=66, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=60, yend=60, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$AWL), max(RIL_pheno$AWL), length.out = 4), 1)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$CD
cd <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=CD), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=45, yend=55, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=50, yend=50, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=50, yend=60, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=55, yend=55, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$CD), max(RIL_pheno$CD), length.out = 4), 2)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$CH
ch <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=CH), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=95, yend=115, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=105, yend=105, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=110, yend=130, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=120, yend=120, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$CH, na.rm=T), max(RIL_pheno$CH, na.rm=T), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$CL
cl <-  ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=CL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=50, yend=60, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=55, yend=55, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=55, yend=65, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=60, yend=60, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$CL), max(RIL_pheno$CL), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$FLA
fla <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=FLA), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=70, yend=84, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=77, yend=77, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=70, yend=84, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=77, yend=77, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$FLA), max(RIL_pheno$FLA), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$FLL
fll <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=FLL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=55, yend=65, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=60, yend=60, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=45, yend=55, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=50, yend=50, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$FLL), max(RIL_pheno$FLL), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$FLW
flw <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=FLW), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=100, yend=120, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=110, yend=110, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=90, yend=110, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=100, yend=100, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$FLW), max(RIL_pheno$FLW), length.out = 4), 1)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$GL
gl <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=GL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=55, yend=67, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=61, yend=61, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=50, yend=62, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=56, yend=56, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$GL, na.rm=T), max(RIL_pheno$GL, na.rm=T), length.out = 4), 1)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$GWE
gwe <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=GWE), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=60, yend=70, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=65, yend=65, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=50, yend=60, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=55, yend=55, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$GWE,na.rm=T), max(RIL_pheno$GWE,na.rm=T), length.out = 4), 1)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$GWI
gwi <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=GWI), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=60, yend=70, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=65, yend=65, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=55, yend=65, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=60, yend=60, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$GWI, na.rm=T), max(RIL_pheno$GWI, na.rm=T), length.out = 4), 2)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$PL
pl <-  ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=PL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=50, yend=60, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=55, yend=55, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=60, yend=70, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=65, yend=65, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$PL), max(RIL_pheno$PL), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$PS
ps <-  ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=PS), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=75, yend=89, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=82, yend=82, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=76, yend=89, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=82, yend=82, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$PS), max(RIL_pheno$PS), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$SN
sn <- ggplot() + 
  geom_histogram(data=RIL_pheno, aes(x=SN), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=70, yend=84, size=lysz, colour=N_col) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=77, yend=77, size=lysz, colour=N_col) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=60, yend=74, size=lysz, colour=R_col) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=67, yend=67, size=lysz, colour=R_col) +
  scale_x_continuous(breaks = round(seq(min(RIL_pheno$SN), max(RIL_pheno$SN), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

(anl/fh/pe/awl/cd/ch/cl/fla/fll/flw/gl/gwe/gwi/pl/ps/sn) +
  plot_layout(ncol = 4, nrow=4,  heights = rep(1,16))
ggsave("231017_RIL_pheno_distri_2.svg", width=10, height=7)

svglite("231017_RIL_pheno_distri_2.svg", width=10, height=7)
(anl/fh/pe/awl/cd/ch/cl/fla/fll/flw/gl/gwe/gwi/pl/ps/sn) +
  plot_layout(ncol = 4, nrow=4,  heights = rep(1,16))
dev.off()

svglite("250309_RIL_pheno_distri_2.svg", width=10, height=7)
(anl/fh/pe/awl/cd/ch/cl/fla/fll/flw/gl/gwe/gwi/pl/ps/sn) +
  plot_layout(ncol = 4, nrow=4,  heights = rep(1,16))
dev.off()



##################
###2. correleation plot
library(ggcorrplot)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(reshape2)

colnames <- c("ANL","AWL","CD","CH","CL","FH","FLA","FLL","FLW","GL",
              "GWE","GWI","PAF","PAS","PAT","PE","PL","PS","SN")
pheno_cor <- RIL[,colnames]
res <- cor(pheno_cor, use="complete.obs", method="pearson")
cor.sig <- cor_pmat(pheno_cor, method="pearson", alternative = "two.sided")

corr_c_sig <- ifelse(cor.sig > 0.05, "",
                     ifelse(cor.sig <= 0.05 & cor.sig > 0.01, "*",
                            ifelse(cor.sig <= 0.01 & cor.sig > 0.001, "**", "***")))

svglite(file="250310_pheno_corr_heatmap.svg", height=8, width=8.5)
ComplexHeatmap::pheatmap(res, col=color_pal2, show_row_dend=FALSE, 
        row_names_side = "left", cluster_rows = FALSE,  cluster_cols = FALSE,
        border_color=NA,  
        treeheight_col=100, fontsize=14, annotation_legend=FALSE, 
        display_numbers=corr_c_sig, fontsize_number=10)
dev.off()

cor.sig.v <- as.vector(cor.sig) 
length(which(cor.sig.v < 0.05 & cor.sig.v != 0))
res.v <- as.vector(res)

#############################################
#######3. day length and flowering time
dl <- read.csv("Lingshui_day_length_2019.csv", head=T)
names(dl)[1] <- "Day"
dl$type <- "A"
dl$type[which(dl$Length>13)] <- "B"
dl$type <- as.factor(dl$type)

#pheno <- read.csv("426_pheno.csv", head=T)
ril_fh <- RIL$FH
ril_fhd <- ril_fh + 122
ril_fhd <- sort(ril_fhd) 
ril_fhd_df <- data.frame(Day=rle(ril_fhd)$values, ind=rle(ril_fhd)$lengths)

f2_pheno <- read.csv("f2pheno.csv", head=T)
f2_fh <- f2_pheno$FDH
f2_fh <- sort(f2_fh)
f2_fh_df <- data.frame(Day=rle(f2_fh)$values, ind=rle(f2_fh)$lengths)

###
f2_par <- read.csv("F2_parent_phe.csv", head=T)
f2_par$fhd <- f2_par$FH - 37

f2_niv <- f2_par %>% 
          subset(Parent=="N" & fhd != 'NA', select=c("fhd")) %>%
          arrange(fhd)     
f2_niv_fhd <- data.frame(Day=rle(f2_niv$fhd)$value, ind=rle(f2_niv$fhd)$length)

f2_niv <- f2_par %>% 
          subset(Parent=="N" & fhd != 'NA', select=c("fhd")) %>%
          arrange(fhd)     
f2_niv_fhd <- data.frame(Day=rle(f2_niv$fhd)$value, ind=rle(f2_niv$fhd)$length)

f2_ruf <- f2_par %>% 
          subset(Parent=="R" & fhd != 'NA', select=c("fhd")) %>%
          arrange(fhd)     
f2_ruf_fhd <- data.frame(Day=rle(f2_ruf$fhd)$value, ind=rle(f2_ruf$fhd)$length)

###
ril_par <- read.csv("231016-RIL_NR亲本表型.csv", head=T)
ril_par$fhd <- ril_par$FH + 122

ril_niv <- ril_par %>% 
           subset(ID=="NIV" & fhd != 'NA', select=c("fhd")) %>%
           arrange(fhd)
ril_niv_fhd <- data.frame(Day=rle(ril_niv$fhd)$value, ind=rle(ril_niv$fhd)$length)

ril_ruf <- ril_par %>% 
           subset(ID=="RUF" & fhd != 'NA', select=c("fhd")) %>%
           arrange(fhd)
ril_ruf_fhd <- data.frame(Day=rle(ril_ruf$fhd)$value, ind=rle(ril_ruf$fhd)$length)

N_col="#990000" 
R_col="#003399"

###
ldp <- ggplot() + 
  scale_y_continuous(expand=c(0,0), limits=c(0,36), name=" ",
                    # breaks=c(10,20,30), labels=c(10,20,30),
                     sec.axis=sec_axis(~./1.5, name=" ")) +
  scale_x_continuous(expand=c(0,0), limits=c(0,365), breaks=grep("15", dl$X2019), 
                     labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
                              "Sep","Oct","Nov","Dec"), name=" ") + 
  theme_bw() +
  geom_area(data=dl, aes(x=Day, y=Length*1.5), colour="black", fill="white") +
  #scale_fill_manual(values=c("gray90","white")) + 
  geom_col(data=ril_fhd_df, aes(x=Day, y=ind), fill="gray60", alpha=0.8) +
  geom_col(data=ril_ruf_fhd, aes(x=Day, y=ind), fill=R_col, alpha=0.6) +
  geom_col(data=ril_niv_fhd, aes(x=Day, y=ind), fill=N_col, alpha=0.6) +
  geom_col(data=f2_fh_df, aes(x=Day, y=ind), fill="gray60", alpha=1) +
  geom_col(data=f2_ruf_fhd, aes(x=Day, y=ind), fill=R_col, alpha=0.6) +
  geom_col(data=f2_niv_fhd, aes(x=Day, y=ind), fill=N_col, alpha=0.6) +
  theme(panel.background=element_rect(fill="gray30"), panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), legend.position="none", 
        axis.title=element_text(size=16, family="Times"), 
        axis.text=element_text(size=12, family="Times")) +
  annotate("segment", x=122, xend=122, y=5, yend=0, linetype="dashed", 
           colour="#008B8B", size=0.5) + #RIL播种
  annotate("segment", x=328, xend=328, y=5, yend=0,  linetype="dashed",
           colour="#8B5A2B", size=0.5)  #F2播种
 
library(svglite)  
svglite("250111_flwoering time and day length-3.svg", width=8, height=5)
ldp
dev.off()



library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,"GnBu")

##########
f2_fh <- f2_pheno$FD
ril_fh <- RIL$FH

ks.test(f2_fh, ril_fh)


###################
###4. first heading of parental lines under Spring and Fall

mean(subset(f2_par, Parent=="N")$FH, na.rm=T)
sd(subset(f2_par, Parent=="N")$FH, na.rm=T)

mean(subset(ril_par, ID=="NIV")$FH, na.rm=T)
sd(subset(ril_par, ID=="NIV")$FH, na.rm=T)

###
mean(subset(f2_par, Parent=="R")$FH, na.rm=T)
sd(subset(f2_par, Parent=="R")$FH, na.rm=T)

mean(subset(ril_par, ID=="RUF")$FH, na.rm=T)
sd(subset(ril_par, ID=="RUF")$FH, na.rm=T)

wilcox.test(subset(f2_par, Parent=="R")$FH, subset(ril_par, ID=="RUF")$FH)
wilcox.test(subset(f2_par, Parent=="N")$FH, subset(ril_par, ID=="NIV")$FH)

f2_fh_r <- subset(f2_fh_r, Parent=="R")$FH
f2_fh_n <- subset(f2_par, Parent=="N")$FH

ril_fh_r <- subset(ril_par, ID=="RUF")$FH
ril_fh_n <- subset(ril_par, ID=="NIV")$FH

mean_al <- c(mean(f2_fh_r, na.rm=T), mean(ril_fh_r, na.rm=T),
             mean(f2_fh_n, na.rm=T), mean(ril_fh_n, na.rm=T))
sd_al <- c(sd(f2_fh_r, na.rm=T), sd(ril_fh_r, na.rm=T), 
           sd(f2_fh_n, na.rm=T), sd(ril_fh_n, na.rm=T))
fh_df <- data.frame(type=c("RUF_SD", "RUF_LD", "NIV_SD", "NIV_LD"),
                    cond=c("SD", "LD", "SD", "LD"),
                    SP=c("RUF","RUF","NIV","NIV"),
                    FH=mean_al, FH_sd=sd_al)
fh_df$type <- factor(fh_df$type, levels=c("RUF_SD", "RUF_LD", "NIV_SD", "NIV_LD"))
fh_df$cond <- factor(fh_df$cond, levels=c("SD","LD"))


df_al <- rbind(data.frame(TP="NIV_SD", cond="SD", SP="NIV", FH=f2_fh_n),
               data.frame(TP="RUF_SD", cond="SD", SP="RUF", FH=f2_fh_r),
               data.frame(TP="NIV_LD", cond="LD", SP="NIV", FH=ril_fh_n),
               data.frame(TP="RUF_LD", cond="LD", SP="RUF", FH=ril_fh_r)
              ) 
df_al$cond <- factor(df_al$cond, levels=c("SD","LD"))


stat_test <- compare_means(FH~cond, data=df_al, method = "t.test", group.by="SP")
stat_test <- stat_test %>% 
  mutate(y.position =c(90, 140))
stat_test[2,8] <- "***"

svglite(file="240501_par_LD_SD.svg", width=5, height=6)
ggplot(df_al, aes(x = SP, y = FH, fill = cond)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(0.7), width = 0.7) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", position = position_dodge(0.8), width = 0.2) +
  geom_jitter(aes(color = cond), 
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8),
              size = 3, shape = 21, stroke = 0.3, fill = "black") +
  scale_fill_manual(values = c("SD" = "gray40", "LD" = "gray70")) +
  scale_color_manual(values = c("SD" = "black", "LD" = "black")) +
  scale_y_continuous(limits = c(0, 150), breaks = c(50,100,150), expand = c(0, 0)) +
  theme_classic() + 
  theme(legend.position = "none", text=element_text(size=14))
dev.off()


####################
###5. barplot of perennial capacity
na <- names(RIL)
LH1 <- data.frame(LH=RIL[,20], Trait=rep("PA for the first stage", nrow(RIL)))
LH2 <- data.frame(LH=RIL[,21], Trait=rep("PA for the second stage", nrow(RIL)))
LH3 <- data.frame(LH=RIL[,22], Trait=rep("PA for the third stage", nrow(RIL)))

LH_n <- rbind(LH1, LH2, LH3)
LH_n$Trait <- factor(LH_n$Trait, levels=c("PA for the first stage","PA for the second stage",
                                          "PA for the third stage"))
LH_tn <- as.data.frame(table(LH_n))

display.brewer.pal(9, "Greys")
col=brewer.pal(9, "Greys")[c(4,6,8)]
library(svglite)
svglite("250508_426_life_history_dist.svg", width=8, height=6)
ggplot(LH_tn, aes(x=LH, y=Freq)) +
  geom_bar(stat="identity", aes(fill=Trait), position=position_dodge(0.8), width=0.8) +
  theme_bw() +
  theme(text=element_text(family="Times New Roman"), 
        panel.grid=element_blank(), legend.position=c(0.65,0.85),
        legend.title=element_blank(), legend.text=element_blank(),
        axis.title=element_text(size=20), axis.text=element_text(size=16)) +
        ylab("No. individuals") + xlab("Perennating ability") +
        scale_fill_manual(values=c(col))
dev.off()

