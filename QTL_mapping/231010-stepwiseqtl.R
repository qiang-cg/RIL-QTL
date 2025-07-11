library(qtl)
library(ggplot2)
library(stringr)

###1.ANL
load("RIL426.ANL.sc2_perms100.RData")
ANL_perms <- c(ANL.sc2.perm1,ANL.sc2.perm2,ANL.sc2.perm3,ANL.sc2.perm4,ANL.sc2.perm5,
               ANL.sc2.perm6,ANL.sc2.perm7,ANL.sc2.perm8,ANL.sc2.perm9,ANL.sc2.perm10)
summary(ANL_perms)
ANL_pen <- calc.penalties(ANL_perms)

ANL.sw1 <- stepwiseqtl(cross, pheno.col=12, max.qtl=15, method="hk", penalties=ANL_pen)
ANL.fit <- fitqtl(cross, pheno.col=12, qtl=ANL.sw1, method="hk", get.ests=TRUE,
  formula=y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q10 + Q11 + Q12 + Q13 + Q14)
summary(ANL.fit)

###2.FH
FH_perms_files <- list.files()[grep("RIL426.FH.sc2_perms100", list.files())]
for(f in FH_perms_files){load(f)}
FH_perms <- c(FH.FH.sc2.perm1,FH.sc2.perm2,FH.sc2.perm3,FH.sc2.perm4,FH.sc2.perm5,
                 FH.sc2.perm6,FH.sc2.perm7,FH.sc2.perm8,FH.sc2.perm9,FH.sc2.perm10)
summary(FH_perms)
FH.pen <- calc.penalties(FH_perms)

FH.sw1 <- stepwiseqtl(cross, pheno.col=1, max.qtl=15, method="hk", penalties=FH.pen)
FH.fit <- fitqtl(cross, pheno.col=1, qtl=FH.sw1, method="hk", get.ests=TRUE, 
  formula=y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q10 + Q11
  + Q5:Q7 + Q6:Q7 + Q7:Q8)
summary(FH.fit)

###3.PE
load("RIL426.PE.sc2_perms100.RData")
PE_perms <- c(PE.sc2.perm1,PE.sc2.perm2,PE.sc2.perm3,PE.sc2.perm4,PE.sc2.perm5,
              PE.sc2.perm6,PE.sc2.perm7,PE.sc2.perm8,PE.sc2.perm9,PE.sc2.perm10)
summary(PE_perms)
PE.pen <- calc.penalties(PE_perms)

PE.sw1 <- stepwiseqtl(cross, pheno.col=9, max.qtl=15, method="hk", penalties=PE.pen)
PE.fit <- fitqtl(cross, pheno.col=9, qtl=PE.sw1, method="hk", get.ests=TRUE, 
  formula=y ~ Q1 + Q2 + Q3 + Q4 + Q5)
summary(PE.fit)

###4.AWL
load("RIL426.AWL.sc2_perms100.RData")
AWL_perms <- c(AWL.sc2.perm1,AWL.sc2.perm2,AWL.sc2.perm3,AWL.sc2.perm4,AWL.sc2.perm5,
               AWL.sc2.perm6,AWL.sc2.perm7,AWL.sc2.perm8,AWL.sc2.perm9,AWL.sc2.perm10)
summary(AWL_perms)
AWL.pen <- calc.penalties(AWL_perms)

AWL.sw1 <- stepwiseqtl(cross, pheno.col=11, max.qtl=15, method="hk", penalties=AWL.pen)
AWL.fit <- fitqtl(cross, pheno.col=11, qtl=AWL.sw1, method="hk", get.ests=TRUE,
                  formula=y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q10 + Q11 + 
                  Q12 + Q3:Q7 + Q8:Q9 + Q10:Q11 + Q4:Q12)
summary(AWL.fit)

###5.CD
load("RIL426.CD.sc2_perms100.RData")
CD_perms <- c(CD.sc2.perm1,CD.sc2.perm2,CD.sc2.perm3,CD.sc2.perm4,CD.sc2.perm5,
              CD.sc2.perm6,CD.sc2.perm7,CD.sc2.perm8,CD.sc2.perm9,CD.sc2.perm10)
summary(CD_perms)
CD.pen <- calc.penalties(CD_perms)

CD.sw1 <- stepwiseqtl(cross, pheno.col=4, max.qtl=6, method="hk", penalties=CD.pen)

###6.CH
load("RIL426.CH.sc2_perms100.RData")
CH_perms <- c(CH.sc2.perm1,CH.sc2.perm2,CH.sc2.perm3,CH.sc2.perm4,CH.sc2.perm5,
              CH.sc2.perm6,CH.sc2.perm7,CH.sc2.perm8,CH.sc2.perm9,CH.sc2.perm10)
summary(CH_perms)

CH.pen <- calc.penalties(CH_perms)
CH.sw <- stepwiseqtl(cross, pheno.col=3, max.qtl=10, method="hk", penalties=CH.pen)
CH.fit <- fitqtl(cross, pheno.col=3, qtl=CH.sw, formula=y ~ Q1 + Q2 + Q3 + Q2:Q3, 
  method="hk", get.ests=TRUE)
summary(CH.fit)

###7.CL
load("RIL426.CL.sc2_perms200.RData")
CL_perms <- c(CL.sc2.perm1,CL.sc2.perm2,CL.sc2.perm3,CL.sc2.perm4,CL.sc2.perm5,
              CL.sc2.perm6,CL.sc2.perm7,CL.sc2.perm8,CL.sc2.perm9,CL.sc2.perm10)
summary(CL_perms)
CL.pen <- calc.penalties(CL_perms)

names(cross$pheno)
CL.sw <- stepwiseqtl(cross, pheno.col=2, max.qtl=15, method="hk", penalties=CL.pen)
CL.fit <- fitqtl(cross, pheno.col=2, qtl=CL.sw, 
  formula=y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q10 + Q11 + Q12 + Q3:Q8, 
  method="hk", get.ests=TRUE)
summary(CL.fit)

###8.FLA
load("RIL426.FLA.sc2_perms100.RData")
FLA_perms <- c(FLA.sc2.perm1,FLA.sc2.perm2,FLA.sc2.perm3,FLA.sc2.perm4,FLA.sc2.perm5,
               FLA.sc2.perm6,FLA.sc2.perm7,FLA.sc2.perm8,FLA.sc2.perm9,FLA.sc2.perm10)
summary(FLA_perms)
FLA.pen <- calc.penalties(FLA_perms)

names(cross$phen)
FLA.sw <- stepwiseqtl(cross, pheno.col=5, max.qtl=15, method="hk", penalties=FLA.pen)
FLA.fit <- fitqtl(cross, pheno.col=5, qtl=FLA.sw, formula=attributes(FLA.sw)$formula, 
  method="hk", get.ests=TRUE)
summary(FLA.fit)

###9.FLL
load("RIL426.FLL.sc2_perms100.RData")
FLL_perms <- c(FLL.sc2.perm1,FLL.sc2.perm2,FLL.sc2.perm3,FLL.sc2.perm4,FLL.sc2.perm5,
               FLL.sc2.perm6,FLL.sc2.perm7,FLL.sc2.perm8,FLL.sc2.perm9,FLL.sc2.perm10)
summary(FLL_perms)
FLL.pen <- calc.penalties(FLL_perms)
names(cross$phen)
FLL.sw <- stepwiseqtl(cross, pheno.col=6, max.qtl=15, method="hk", penalties=FLL.pen)
FLL.fit <- fitqtl(cross, pheno.col=6, qtl=FLL.sw, formula=attributes(FLL.sw)$formula, 
  method="hk", get.ests=TRUE)
summary(FLL.fit)

###10.FLW
load("RIL426.FLW.sc2_perms100.RData")
FLW_perms <- c(FLW.sc2.perm1,FLW.sc2.perm2,FLW.sc2.perm3,FLW.sc2.perm4,FLW.sc2.perm5,
               FLW.sc2.perm6,FLW.sc2.perm7,FLW.sc2.perm8,FLW.sc2.perm9,FLW.sc2.perm10)
summary(FLW_perms)
FLW.pen <- calc.penalties(FLW_perms)
names(cross$phen)
FLW.sw <- stepwiseqtl(cross, pheno.col=7, max.qtl=15, method="hk", penalties=FLW.pen)
FLW.fit <- fitqtl(cross, pheno.col=7, qtl=FLW.sw, formula=attributes(FLW.sw)$formula, 
  method="hk", get.ests=TRUE)
summary(FLW.fit)

###11.GL
load("RIL426.GL.sc2_perms100.RData")
GL_perms <- c(GL.sc2.perm1,GL.sc2.perm2,GL.sc2.perm3,GL.sc2.perm4,GL.sc2.perm5,
              GL.sc2.perm6,GL.sc2.perm7,GL.sc2.perm8,GL.sc2.perm9,GL.sc2.perm10)
GL.pen <- calc.penalties(GL_perms)
names(cross$phen)
GL.sw <- stepwiseqtl(cross, pheno.col=23, max.qtl=15, method="hk", penalties=GL.pen)
GL.fit <- fitqtl(cross, pheno.col=23, qtl=GL.sw, formula=attributes(GL.sw)$formula, 
  method="hk", get.ests=TRUE)
summary(GL.fit)

###12.GWE
load("RIL426.GW30.sc2_perms100.RData")
GWE_perms <- c(GW30.sc2.perm1,GW30.sc2.perm2,GW30.sc2.perm3,GW30.sc2.perm4,GW30.sc2.perm5,
              GW30.sc2.perm6,GW30.sc2.perm7,GW30.sc2.perm8,GW30.sc2.perm9,GW30.sc2.perm10)
GWE_pen <- calc.penalties(GWI_perms)

GWE_sw1 <- stepwiseqtl(cross, pheno.col=22, max.qtl=10, method="hk", penalties=GWE_pen)
GWE_fit <- fitqtl(cross, pheno.col=22, qtl=GWI_sw1, method="hk", get.ests=TRUE,
  formula=y ~ Q1 + Q2)
summary(GWI_fit)

###13.GWI
load("RIL426.GWI.sc2_perms100.RData")
GWI_perms <- c(GWI.sc2.perm1,GWI.sc2.perm2,GWI.sc2.perm3,GWI.sc2.perm4,GWI.sc2.perm5,
               GWI.sc2.perm6,GWI.sc2.perm7,GWI.sc2.perm8,GWI.sc2.perm9,GWI.sc2.perm10)
summary(GWI_perms)
names(cross$pheno)
GWI.pen <- calc.penalties(GWI_perms)
GWI_sw1 <- stepwiseqtl(cross, pheno.col=24, max.qtl=12, method="hk", penalties=GWI_pen)
GWI.fit <- fitqtl(cross, pheno.col=24, qtl=GWI_sw1, formula=attributes(GWI_sw1)$formula, 
  method="hk", get.ests=TRUE)
summary(GWI.fit)

###14.PL
load("RIL426.PL.sc2_perms100.RData")
PL_perms <- c(PL.sc2.perm1,PL.sc2.perm2,PL.sc2.perm3,PL.sc2.perm4,PL.sc2.perm5,
              PL.sc2.perm6,PL.sc2.perm7,PL.sc2.perm8,PL.sc2.perm9,PL.sc2.perm10)
summary(PL_perms)
PL_pen <- calc.penalties(PL_perms)
names(cross$pheno)
PL_sw1 <- stepwiseqtl(cross, pheno.col=10, max.qtl=15, method="hk", penalties=PL_pen)
PL.fit <- fitqtl(cross, pheno.col=10, qtl=PL_sw1, formula=attributes(PL_sw1)$formula, 
  method="hk", get.ests=TRUE)
summary(PL.fit)

###15.PS
load("RIL426.PS.sc2_perms100.RData")
PS_perms <- c(PS.sc2.perm1,PS.sc2.perm2,PS.sc2.perm3,PS.sc2.perm4,PS.sc2.perm5,
              PS.sc2.perm6,PS.sc2.perm7,PS.sc2.perm8,PS.sc2.perm9,PS.sc2.perm10)
summary(PS_perms)
PS_pen <- calc.penalties(PS_perms)
names(cross$pheno)
PS_sw1 <- stepwiseqtl(cross, pheno.col=8, max.qtl=5, method="hk", penalties=PS_pen)
PS.fit <- fitqtl(cross, pheno.col=8, qtl=PS_sw1, formula=attributes(PS_sw1)$formula, 
  method="hk", get.ests=TRUE)
summary(PS.fit)

###16.SN
load("RIL426.SN.sc2_perms100.RData")
SN_perms <- c(SN.sc2.perm1,SN.sc2.perm2,SN.sc2.perm3,SN.sc2.perm4,SN.sc2.perm5,
              SN.sc2.perm6,SN.sc2.perm7,SN.sc2.perm8,SN.sc2.perm9,SN.sc2.perm10)
summary(SN_perms)
names(cross$pheno)
SN_pen <- calc.penalties(SN_perms)
SN_sw1 <- stepwiseqtl(cross, pheno.col=13, max.qtl=8, method="hk", penalties=SN_pen)
SN.fit <- fitqtl(cross, pheno.col=13, qtl=SN_sw1, formula=attributes(SN_sw1)$formula, 
  method="hk", get.ests=TRUE)
summary(SN.fit)

###17.PC_1
load("RIL426.LH_1.sc2_perms100.RData")
PC_1_perms <- c(PC_1.sc2.perm1,PC_1.sc2.perm2,PC_1.sc2.perm3,PC_1.sc2.perm4,PC_1.sc2.perm5,
                PC_1.sc2.perm6,PC_1.sc2.perm7,PC_1.sc2.perm8,PC_1.sc2.perm9,PC_1.sc2.perm10)
summary(PC_1_perms)
PC_1.pen <- calc.penalties(PC_1_perms)
names(cross$pheno)

PC_1.sw1 <- stepwiseqtl(cross, pheno.col=19, max.qtl=10, method="hk", penalties=PC_1.pen)
PC_1.fit <- fitqtl(cross, pheno.col=19, qtl=PC_1.sw1, method="hk", get.ests=TRUE,
                   formula=attributes(PC_1.sw1)$formula)
summary(PC_1.fit)

###18.PC_2
load("RIL426.PC_2.sc2_perms100.RData")
PC_2_perms <- c(PC_2.sc2.perm1,PC_2.sc2.perm2,PC_2.sc2.perm3,PC_2.sc2.perm4,PC_2.sc2.perm5,
                PC_2.sc2.perm6,PC_2.sc2.perm7,PC_2.sc2.perm8,PC_2.sc2.perm9,PC_2.sc2.perm10)
summary(PC_2_perms)
PC_2.pen <- calc.penalties(PC_2_perms)

PC_2.sw1 <- stepwiseqtl(cross, pheno.col=20, max.qtl=10, method="hk", penalties=PC_2.pen)
PC_2.fit <- fitqtl(cross, pheno.col=20, qtl=PC_2.sw1, method="hk", get.ests=TRUE,
                   formula=attributes(PC_2.sw1)$formula)
summary(PC_2.fit)

###22.PC_2
load("RIL426.PC_3.sc2_perms100.RData")
PC_3_perms <- c(PC_3.sc2.perm1,PC_3.sc2.perm2,PC_3.sc2.perm3,PC_3.sc2.perm4,PC_3.sc2.perm5,
                PC_3.sc2.perm6,PC_3.sc2.perm7,PC_3.sc2.perm8,PC_3.sc2.perm9,PC_3.sc2.perm10)
summary(PC_3_perms)
PC_3.pen <- calc.penalties(PC_3_perms)

PC_3.sw1 <- stepwiseqtl(cross, pheno.col=21, max.qtl=5, method="hk", penalties=PC_3.pen)
PC_3.fit <- fitqtl(cross, pheno.col=21, qtl=PC_3.sw1, method="hk", get.ests=TRUE,
                   formula=attributes(PC_3.sw1)$formula)
summary(PC_3.fit)

####################
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,"Set1")
col <- brewer.pal(9,"Set1")[1:4]

PC_1.lod <- attr(PC_1.sw1, "lodprofile")
PC_1 <- PC_1.lod[["4@13.9"]

PC_2.lod <- attr(PC_2.sw1, "lodprofile")
PC_2 <- PC_2.lod[["4@14.3"]]

PC_3.lod <- attr(PC_3.sw1, "lodprofile")
PC_3 <- PC_3.lod[["4@19.2"]]

p <- ggplot() + 
  geom_line(data=PC_1, aes(x=pos, y=lod), color=col[1], size=1) +
  geom_line(data=PC_2, aes(x=pos, y=lod), color=col[2], size=1) +
  geom_line(data=PC_3, aes(x=pos, y=lod), color=col[3], size=1) +
  annotate("rect", xmin=13.4, xmax=14.88, ymin=6, ymax=7.5, alpha=1, fill=col[1]) + 
  annotate("rect", xmin=12.4, xmax=14.61, ymin=4, ymax=5.5, alpha=1, fill=col[2]) + 
  annotate("rect", xmin=9.73, xmax=20.41, ymin=2, ymax=3.5, alpha=1, fill=col[3]) +
  annotate("rect", xmin=12.4, xmax=14.88, ymin=0, ymax=1.5, fill=col[4]) +
  geom_vline(xintercept=c(12.4, 14.88), linetype="dashed", colour="gray60") + 
  xlab("Chromosome 4 (cM)") + ylab("LOD value") + 
  theme_classic() +
  theme(axis.title=element_text(family="serif", size=18), 
        axis.text=element_text(family="serif", size=16),
        plot.margin=unit(c(4,10,4,4), "pt"))

library(svglite)
svglite(file="231217_PA4_plot.svg", height=5, width=9) 
p
dev.off()

###################################
qtl_info_extract <- function(rqtl, fit){
    left_ran <- c()
	right_ran <- c()
	pos <- c()
	genet_left <- c()
	Q <- c()
	start_cM <- c(); end_cM <- c(); peak_cM <- c()
	for(i in 1:rqtl$n.qtl){
	    interval<-lodint(rqtl, qtl.index=i, drop=1.5, expandtomarkers=TRUE)
        mar <- row.names(interval)
        left_ran <- c(left_ran, unlist(strsplit(mar[1], split="_"))[2])
        right_ran <- c(right_ran, unlist(strsplit(mar[3], split="_"))[3])
		pos <- c(pos, paste(unlist(strsplit(mar[2], split="_"))[2:3], collapse="_"))
		Q <- c(Q, paste("Q",i,sep=""))
		start_cM <- c(start_cM, interval[1,2]); end_cM <- c(end_cM, interval[3,2])
        peak_cM <- c(peak_cM, interval[2,2])
	}
    left_ran <- as.numeric(left_ran)
    right_ran <- as.numeric(right_ran)
    int_size <- right_ran - left_ran	
    add <- fit$ests[seq(2,2*rqtl$n.qtl+1,2)]
	dom <- fit$ests[seq(3,2*rqtl$n.qtl+1,2)]
	add <- round(add, 4)
	dom <- round(dom, 4)
	start_cM <- round(start_cM, 2)
    end_cM <- round(end_cM, 2)
    peak_cM <- round(peak_cM, 2)
	if(qtl$n.qtl == 1){
	    LOD <- fit$result.full[1,4]
        PVE <- fit$result.full[1,5]
    }else{
        PVE <- c(t(fit$result.drop[1:rqtl$n.qtl,4]))
	    LOD <- c(t(fit$result.drop[1:rqtl$n.qtl,3]))
    } 
    LOD <- round(LOD, 2)
    PVE <- round(PVE, 4)	
	qtl_int <- data.frame(Q=Q, chr=rqtl$chr, start_cM=start_cM, end_cM=end_cM, gent_peak=peak_cM, start=left_ran, end=right_ran, peak_pos=pos,  
	                      int_size=int_size, LOD=LOD, PVE=PVE, ADD=add, DOM=dom)
	###information of epstasis
	dim <- attributes(fit$result.drop)$dim
	if(dim[1]>rqtl$n.qtl){
		qtl_stat <- as.data.frame(fit$result.drop[1:dim[1],])
	    epi <- qtl_stat[(qtl$n.qtl+1):dim[1],]
		#epi %>% round(2)
		epi$LOD <- round(epi$LOD,2); epi$"%var" <- round(epi$"%var",4); epi$"F value" <- round(epi$"F value",2)
	}else{epi <- "none epistasis detected"}
	qtl_info <- list()
	qtl_info[["qtl interval info"]] <- qtl_int
	qtl_info[["qtl epistasis"]] <- epi
	return(qtl_info)
}	

#######################
FH.lod <- attr(FH.sw1, "lodprofile")

library(dplyr)
library(svglite)

svglite(file="231215.RIL_FH.svg", width=10, height=5)
plotLodProfile(FH.sw1, qtl.labels=F, showallchr=T, bandcol="gray70",
               incl.markers=F, ylab="LOD value")
abline(h=4.11, col="red")
dev.off()

FH.sw1.lod <- attr(FH.sw1, "lodprofile")
#FH_sw1_lod_df <- do.call(rbind, FH.sw1.lod)
FH_sw1_lod_df <- bind_rows(
  lapply(names(FH.sw1.lod), function(name) {
    df <- FH.sw1.lod[[name]]
    df$region <- name
    return(df)
  }),
  .id = "id"
)
head(FH_sw1_lod_df)

FH_sw1_lod_df$chr <- factor(FH_sw1_lod_df$chr, levels=seq(1:12), 
                            labels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                                     "chr7", "chr8", "chr9", "chr10", "chr11", "chr12"))

map <- read.table("E:\\2018-2024-博士\\RIL相关分析\\RIL426_final_results\\QTL分布图\\Genetic_Map.txt", head=F)
names(map) <- c("bin", "genetic_pos")
bin_reg <- as.data.frame(matrix(unlist(strsplit(map$bin, split="_")), byrow=T, ncol=3))
names(bin_reg) <- c("chr", "start", "end")
map <- cbind(map, bin_reg)
map$chr <- as.factor(map$chr)
map$phy_pos <- (as.numeric(map$end) + as.numeric(map$start))/2

map$chr <- factor(map$chr, levels=c("chr1","chr2","chr3","chr4","chr5","chr6",
                                    "chr7","chr8","chr9","chr10","chr11","chr12"))
chr_len <- map %>% group_by(chr) %>% summarise(len=max(genetic_pos)+10) %>%
           mutate(chr_pos=cumsum(len))
chr_len$chr_pos <- chr_len$chr_pos - chr_len$len

map_cum <- chr_len %>% left_join(map, ., by="chr") %>% mutate(pos_cum=genetic_pos+chr_pos)

x_axis <- map_cum %>% group_by(chr) %>% summarize(center=(max(pos_cum)+min(pos_cum))/2)

FH_sw1_lod_df_cum <- chr_len %>% left_join(FH_sw1_lod_df, ., by="chr") %>% mutate(pos_cum=pos+chr_pos)
FH_sw1_lod_df_cum

col <- c("1@17.6"="#3914AF", "2@119.3"="#A68900", "3@44.4"="#3914AF", "6@22.1"="#A68900", "6@39.2"="#A68900", "7@42.8"="#3914AF", 
         "7@99.1"="#3914AF", "8@22.3"="#A68900", "8@61.9"="#A68900", "9@48.8"="#3914AF", "10@28.9"="#A68900")

qtl <- read.table("E:\\2018-2024-博士\\RIL相关分析\\RIL426_final_results\\QTL分布图\\QTL_pos.txt", head=T, sep="\t")
FH_RIL_qtl <- subset(qtl, Traits=="First heading")
FH_RIL_qtl <- FH_RIL_qtl[c(2,3,4,5,7:11),]

FH_RIL_qtl_cum <- chr_len %>% left_join(FH_RIL_qtl, ., by="chr") %>% 
  mutate(start_cum=start_cM+chr_pos, end_cum=end_cM+chr_pos, peak_cum=peak_cM+chr_pos)

FH_RIL_qtl_cum <- edit(FH_RIL_qtl_cum)

RIL_FH_plot <- ggplot(data=FH_sw1_lod_df_cum, aes(x=pos_cum, y=lod, color=region)) +
  #annotate("rect", xmin=FH_RIL_qtl_cum$start_cum, xmax=FH_RIL_qtl_cum$end_cum, 
  #         ymin=0, ymax=Inf, fill="#FF1300", alpha=0.5) +
  #annotate("segment", x=FH_RIL_qtl_cum$peak_cum, xend=FH_RIL_qtl_cum$peak_cum, 
  #          y=0, yend=Inf, color = "#cb3b3b", linewidth=1) + 
  geom_line(linewidth=1.2) +
  geom_hline(yintercept=4.11, col="red") +
  scale_colour_manual(values=col) + 
  scale_x_continuous(limits=c(0, 1401.6), breaks=x_axis$center, labels=x_axis$chr) + 
  theme_classic() +
  theme(legend.position="none", axis.line=element_line(linewidth=1.1))
  
svglite(file="250629_RIL_FH_stw.svg", width=10, height=3.5)
RIL_FH_plot
dev.off()

########
FH.cim <- cim(cross, pheno.col=1, n.marcovar=11, method="hk", map.function="kosambi")

svglite(file="231215.RIL_FH_cim.svg", width=10, height=5)
plot(FH.cim,qtl.labels=F, showallchr=T, bandcol="gray70",
               incl.markers=F, ylab="LOD value")
abline(h=4.11, col="red")
dev.off()

FH.cim$chr <- factor(FH.cim$chr, levels=1:12, 
                     labels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                              "chr7", "chr8", "chr9", "chr10", "chr11", "chr12"))

FH.cim_cum <- chr_len %>% left_join(FH.cim, ., by="chr") %>% mutate(pos_cum=pos+chr_pos)

RIL_FH_cim_plot <- ggplot(data=FH.cim_cum, aes(x=pos_cum, y=lod, color=chr)) +
  geom_line(linewidth=1.2) +
  scale_colour_manual(values=rep(c("#3914AF", "#A68900"),6)) + 
  scale_x_continuous(limits=c(0, 1401.6), breaks=x_axis$center, labels=x_axis$chr) + 
  geom_hline(yintercept=4.11, col="red") +
  theme_classic() +
  theme(legend.position="none", axis.line=element_line(linewidth=1.1))

svglite(file="250507_RIL_FH_cim.svg", width=10, height=3.5)
RIL_FH_cim_plot
dev.off()



