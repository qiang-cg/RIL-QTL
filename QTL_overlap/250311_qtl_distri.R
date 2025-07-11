library(ggplot2)
library(svglite)
library(dplyr)
library(cowplot)

map <- read.table("Genetic_Map.txt", head=F)
names(map) <- c("bin", "genetic_pos")
bin_reg <- as.data.frame(matrix(unlist(strsplit(map$bin, split="_")), byrow=T, ncol=3))
names(bin_reg) <- c("chr", "start", "end")
map <- cbind(map, bin_reg)
map$chr <- as.factor(map$chr)

qtl <- read.table("QTL_pos.txt", head=T, sep="\t")
qtl <- subset(qtl, Cart != "CR")
qtl$Cart <- factor(qtl$Cart, levels=c("RI","FR","LH"), labels=c("RR","HR","PA"))

qtl <- arrange(qtl, Traits)
n=1
for(t in unique(qtl$Traits)){
  qtl[qtl$Traits==t, 2] <- n
  n=n+1
}

#########################
###1.QTL distribution horizontal 

map$chr <- factor(map$chr, levels=c("chr1","chr2","chr3","chr4","chr5","chr6",
                                    "chr7","chr8","chr9","chr10","chr11","chr12"))
chr_len <- map %>% group_by(chr) %>% summarise(len=max(genetic_pos)+10) %>%
           mutate(chr_pos=cumsum(len))
chr_len$chr_pos <- chr_len$chr_pos - chr_len$len

map_cum <- chr_len %>% left_join(map, ., by="chr") %>% mutate(pos_cum=genetic_pos+chr_pos)

qtl_cum <- chr_len %>% left_join(qtl, ., by="chr") %>% 
  mutate(start_cum=start_cM+chr_pos, end_cum=end_cM+chr_pos, peak_cum=peak_cM+chr_pos)

qtl_cum <- arrange(qtl_cum, Traits)
#qtl_cum$Traits <- factor(qtl$Traits, levels=sort(unique(qtl_cum$Traits), decreasing=TRUE))

chro_height <- 0.5; gap <- 0.8+0.5+0.5; height=0.8
qtl_cum$he_s <- chro_height + gap
qtl_cum$he_e <- chro_height + gap + height

x_axis <- map_cum %>% group_by(chr) %>% summarize(center=(max(pos_cum)+min(pos_cum))/2)

n=0
for(t in rev(unique(qtl_cum$Traits))){
  row_id <- which(qtl_cum$Traits==t)
  qtl_cum$he_s[row_id] <- qtl_cum$he_s[row_id] + gap*n
  qtl_cum$he_e[row_id] <- qtl_cum$he_e[row_id] + gap*n + height
  n <- n + 1
}

qtl_cum$region_size <- qtl_cum$end_cum - qtl_cum$start_cum
for(q in 1:length(qtl_cum$region_size)){
  if (qtl_cum$region_size[q] < 2){
    mid_pos <- (qtl_cum$start_cum[q] + qtl_cum$end_cum[q])/2 
    qtl_cum$start_cum[q] <- mid_pos - 1
    qtl_cum$end_cum[q]   <- mid_pos + 1
  }
}

y_axis_center <- unique((qtl_cum$he_s+qtl_cum$he_e)/2)

hs <- data.frame(chr=c("chr1", "chr1", "chr3", "chr6", "chr6", "chr7", "chr8", "chr9"), 
                 start=c(11.00, 111.51, 107.92, 16.37, 35.4, 96.25, 0.56, 32.09), 
                 end=c(30.54, 132.64, 125.39, 24.98, 39.94, 99.10, 23.55, 66.77))	
hs <- chr_len %>% left_join(hs, ., by="chr") %>% 
  mutate(start_cum=start+chr_pos, end_cum=end+chr_pos)

trait_name <- c("Anther length", "Awn length", "Culm habit", "Culm length", "First heading",
                "Flag leaf attitude", "Flag leaf length", "Flag leaf width", "Grain length", "Grain weight",
                "Grain width", "Perennating ability of the first stage", "Perennating ability of the second stage",
                "Perennating ability of the third stage", "Panicle exsertion", "Panicle length", "Panicle shape", 
                "Spikelet number")

p <- ggplot() +
  annotate("segment", x=map_cum$pos_cum, xend=map_cum$pos_cum, y=0, yend=0.8, color="black") +
  annotate("rect", xmin=hs$start_cum, xmax=hs$end_cum, ymin=0.8, ymax=Inf, 
           fill="gray", alpha=.6)+
  annotate("rect", xmin=qtl_cum$start_cum, xmax=qtl_cum$end_cum, 
           ymin=qtl_cum$he_s, ymax=qtl_cum$he_e, fill="black") + 
  scale_y_continuous(breaks=y_axis_center, 
                     #labels=unique(qtl_cum$Traits),
                     labels=trait_name,
                     expand=c(0,0), limits=c(0,36)) +
  scale_x_continuous(breaks=x_axis$center, labels=x_axis$chr) + 
  theme_bw() + 
  theme(axis.title=element_blank(), 
        panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(),
        #panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        text=element_text(size=18, family="Times new roman"))
  
svglite(file="250504_QTL_distri.svg", width=15, height=6)
p
dev.off()


#################
###2. QTL distribution vertical
###The fixed distance between chromosome should be drawn with the 'gap'.
p <- ggplot() +
     annotate("segment", x=-1, xend=-1, y=0, yend=-max(map$genetic_pos)) + 
     theme_bw() +
     theme(axis.text=element_blank(), panel.grid=element_blank(), 
           panel.border=element_blank(), axis.title=element_blank()) 
              
n=0
#gap=20; 
bar_width=1; bar_dis=0.2; chro_width=0.8  ###
chro_end_pos = 0; chro_gap = 3
for(c in unique(map$chr)[1:6]){
  map_data <- subset(map, chr==c)
  qtl_data <- subset(qtl, chr==c)
  qtl_data <- qtl_data[order(qtl_data$Cart, qtl_data$start_cM,qtl_data$Trait_ID),]
  #chrom_pos <- gap*n
  chrom_pos <- chro_end_pos + chro_gap
  qtl_data$x_right <- rep(chrom_pos+chro_width+0.5, nrow(qtl_data))
  qtl_data$x_left  <- qtl_data$x_right+bar_width
  qtl_data$layer <- 0
  for(i in 2:nrow(qtl_data)){
    start_i <- qtl_data[i,6];end_i <- qtl_data[i,7]
    start_p <- qtl_data[1:(i-1),6];end_p <- qtl_data[1:(i-1),7]
    ly <- 0
    ly_c <- c()
    for(j in 1:length(start_p)){
      if(end_i >= start_p[j] & end_p[j] >= start_i){ly_c <- c(ly_c, qtl_data$layer[j]);ly <- max(ly_c)+1}
    }
    qtl_data$layer[i] <- ly
    #if("TRUE" %in% c(start_i <= end_p)){
    #   qtl_data[i,12:13] <- qtl_data[i-1,12:13]+(bar_width+bar_dis)}
  }
  qtl_data[,13:14] <- qtl_data[,13:14]+qtl_data$layer*(bar_width+bar_dis)
  chro_end_pos <- max(qtl_data$x_left)
  p <- p + 
       annotate("segment", x=chrom_pos, xend=chrom_pos+chro_width, y=-map_data$genetic_pos, yend=-map_data$genetic_pos) +
       annotate("segment", x=chrom_pos+(chro_width/2), xend=chrom_pos+(chro_width/2), y=0, yend=-max(map_data$genetic_pos)) +
       annotate("rect", xmin=qtl_data$x_right, xmax=qtl_data$x_left, 
                ymin=-qtl_data$start_cM, ymax=-qtl_data$end_cM, fill="gray45", alpha=0.7) + 
       #annotate("text", x=qtl_data$x_right+0.25, y=-qtl_data$peak_cM, label=as.character(qtl_data$Trait_ID), 
       #         colour = "black", size=8, family="Times New Roman", fontface="bold")       
       #annotate("segment", x=qtl_data$x_right, xend=qtl_data$x_left, 
       #          y=-qtl_data$peak_cM, yend=-qtl_data$peak_cM, colour="white") + 
       annotate("text", x=qtl_data$x_right+bar_width/2, y=-(qtl_data$start_cM+qtl_data$end_cM)/2,
                label=as.character(qtl_data$Trait_ID), colour = "black", 
                size=6, family="Times New Roman", fontface="bold")      
  n=n+1
}

x_left_max <- max(qtl_data$x_left)

svglite("250311_qtl_dist_chr1-6.svg", width=21, height=6)
p
dev.off()


###chr7-chr12###
#library(dplyr)
#map$chr <- factor(map$chr, levels=c("chr1","chr2","chr3","chr4","chr5","chr6",
#                                    "chr7","chr8","chr9","chr10","chr11","chr12"))

p <- ggplot() +
     annotate("segment", x=-1, xend=-1, y=0, yend=-max(map$genetic_pos)) + 
     theme_bw() +
     theme(axis.text=element_blank(), panel.grid=element_blank(), 
           panel.border=element_blank(), axis.title=element_blank()) 
              
n=0
#gap=18; bar_width=1; bar_dis=0.2  ###
chro_end_pos = 0; 
#chro_gap = 4
for(c in unique(map$chr)[7:12]){
  map_data <- subset(map, chr==c)
  qtl_data <- subset(qtl, chr==c)
  qtl_data <- qtl_data[order(qtl_data$Cart, qtl_data$start_cM,qtl_data$Trait_ID),]
  #chrom_pos <- gap*n
  chrom_pos <- chro_end_pos + chro_gap
  qtl_data$x_right <- rep(chrom_pos+chro_width+0.5, nrow(qtl_data))
  qtl_data$x_left  <- qtl_data$x_right+bar_width
  qtl_data$layer <- 0
  for(i in 2:nrow(qtl_data)){
    start_i <- qtl_data[i,6];end_i <- qtl_data[i,7]
    start_p <- qtl_data[1:(i-1),6];end_p <- qtl_data[1:(i-1),7]
    ly <- 0
    for(j in 1:length(start_p)){
      if(end_i >= start_p[j] & end_p[j] >= start_i){ly <- qtl_data$layer[j]+1}
    }
    qtl_data$layer[i] <- ly
    #if("TRUE" %in% c(start_i <= end_p)){
    #   qtl_data[i,12:13] <- qtl_data[i-1,12:13]+(bar_width+bar_dis)}
  }
  qtl_data[,13:14] <- qtl_data[,13:14]+qtl_data$layer*(bar_width+bar_dis)
  chro_end_pos <- max(qtl_data$x_left) 
  p <- p + 
       annotate("segment", x=chrom_pos, xend=chrom_pos+chro_width, y=-map_data$genetic_pos, yend=-map_data$genetic_pos) +
       annotate("segment", x=chrom_pos+(chro_width/2), xend=chrom_pos+(chro_width/2), y=0, yend=-max(map_data$genetic_pos)) +
       annotate("rect", xmin=qtl_data$x_right, xmax=qtl_data$x_left, 
                ymin=-qtl_data$start_cM, ymax=-qtl_data$end_cM, fill="gray45", alpha=0.7) + 
       #annotate("text", x=qtl_data$x_right+0.25, y=-qtl_data$peak_cM, label=as.character(qtl_data$Trait_ID), 
       #         colour = "black", size=8, family="Times New Roman", fontface="bold")       
       annotate("text", x=qtl_data$x_right+bar_width/2, y=-(qtl_data$start_cM+qtl_data$end_cM)/2,
                label=as.character(qtl_data$Trait_ID), colour = "black", 
                size=6, family="Times New Roman", fontface="bold")
  n=n+1
}

p <- p+annotate("rect", xmin=x_left_max-1, xmax=x_left_max, 
                ymin=0, ymax=0, fill="white", alpha=0.7)

svglite("250311_qtl_dist_chr7-12.svg", width=21, height=6)
p
dev.off()

#####################
###3. overlap index
qtl <- qtl[,-ncol(qtl)]

tr <- unique(qtl$Traits)
tr_ab <- combn(tr,2)

QTL_OP <- list()
Trait_op=c(); overlap_num=c();op_index=c()
for(i in 1:ncol(tr_ab)){
  ta <- tr_ab[1,i]; tb <- tr_ab[2,i]
  Ta <- subset(qtl, Traits==ta)
  Tb <- subset(qtl, Traits==tb)
  Tab_num <- 0
  overlap_qtl <- list()
  op <- paste(ta,"-",tb)
  for(q in 1:nrow(Ta)){
    TaC <- Ta$chr[q]; TaS <- Ta$start_cM[q]; TaE <- Ta$end_cM[q]
    ol <- subset(Tb, chr==TaC & start_cM<=TaE & end_cM>=TaS)      
    Tab_num <- Tab_num + nrow(ol)
    if(nrow(ol)>0){overlap_qtl[[q]] <- cbind(Ta[q,],ol)}
  }
  QTL_OP[[op]] <- do.call(rbind, overlap_qtl)
  Trait_op=c(Trait_op,op); overlap_num <- c(overlap_num,Tab_num)
  index=Tab_num/(nrow(Ta)+nrow(Tb)-Tab_num)
  op_index=c(op_index,index)
}  

OP_num <- data.frame(Trait_op=Trait_op, overlap_num=overlap_num, overlap_index=op_index)
write.csv(OP_num, file="240207.qtl_overlap.csv", quote=F, row.names=F)

trait_p <- as.data.frame(matrix(unlist(strsplit(OP_num$Trait_op, split=" - ")),ncol=2, byrow=T))
names(trait_p) <- c("trait1", "trait2")
OP_num_n <- cbind(trait_p, OP_num) 

### overlap index heatmap
tmp <- unique(c(OP_num_n$trait1, OP_num_n$trait2))

idx_df <- list()
for(t in tmp){
  df_tm <- subset(OP_num_n, trait1==t | trait2==t)  
  diff <- which(df_tm$trait1 != t)
  df_tm[diff, "trait2"] <- df_tm[diff, "trait1"]
  df_tm[diff, "trait1"] <- t
  idx_df[[t]] <- df_tm
}
idx_df <- do.call(rbind, idx_df)
df <- data.frame(trait1=tmp, trait2=tmp, Trait_op=0, overlap_num=0, overlap_index=1)
idx_df <- rbind(idx_df, df)

#trait_order <- c("CH","BLSC","LC","GWE","GWI","AWL","SC","FLW","SN","FLL","FH","CL",
                 "PL","ANL","PE","GL","FLA","PS","AWC","PAT","PAF","PAS")

trait_order <- tmp

idx_df$trait1 <- factor(idx_df$trait1, levels=rev(trait_order))
idx_df$trait2 <- factor(idx_df$trait2, levels=trait_order)

#############

idx_p <- ggplot(idx_df, aes(x=trait2, y=trait1, fill=overlap_index)) +
  geom_tile() +
  scale_fill_gradient2(low = 'gray',high ='blue') +
  annotate("text", x=idx_df$trait1, y=idx_df$trait2, label=round(idx_df$overlap_index,2),
           size=4) +
  theme_bw() + 
  theme(panel.border = element_blank(),
        axis.ticks=element_blank(), 
        axis.title=element_blank(),
        axis.text=element_text(size=14, family="Times"),
        axis.text.x=element_text(angle=90, vjust=0.2))
  
svglite(idx_p, file="250312_overlap_iddex_all.svg", height=7, width=9)
idx_p
dev.off()

ppi=300
png(paste("240515_overlap_iddex_all.png", sep=""), height = 7*ppi, width = 8.5*ppi, res = 300)
idx_p
dev.off()

###########
###calculate mean overlap index for each trait
tr <- unique(c(QTL_corr$Var2, QTL_corr$Var1))
corr_sta <- list()
for(t in tr){
  tr_corr <- subset(QTL_corr, Var2==t | Var1 == t)
  corr_mean <- mean(abs(tr_corr$corr_value))
  ovlp_mean <- mean(tr_corr$overlap_index)
  corr_sta[[t]] <- data.frame(trait=t, corr_mean=corr_mean, overlap_mean=ovlp_mean)
}
corr_sta <- do.call(rbind, corr_sta)

corr_sta[,2:3] <- round(corr_sta[,2:3], 3)
#write.csv(corr_sta, file="240207_overlap_idx_avg.csv", row.names=F, quote=F)
#####

###random sample for overlap index
qtl <- read.table("240207_QTL_pos_abrr_trait_name.txt", head=T)

qtl_sim <- qtl[,c(1,4:7)]
qtl_sim$size <- qtl_sim$end_cM - qtl_sim$start_cM

overlap_permutation <- list()
for(p in 1:1000){
  samp <- sample(1:nrow(map), 133, replace=TRUE)
  pos_s <- map[samp,]
  qtl_s <- list()
  for(i in 1:length(samp)){
    qtl_info <- qtl_sim[i, 1:2]
    chr_s <- pos_s$chr[i]
    chr_max <- max(subset(map, chr==chr_s)$genetic_pos)
    genet_pos_s <- pos_s$genetic_pos[i]
    qtl_size <- qtl_sim$size[i]
    if((genet_pos_s + qtl_size) <= chr_max){
      qtl_s[[i]] <- cbind(qtl_info, data.frame(chr=chr_s, start=genet_pos_s, end=genet_pos_s+qtl_size))
    }else{
      qtl_s[[i]] <- cbind(qtl_info, data.frame(chr=chr_s, start=genet_pos_s-qtl_size, end=genet_pos_s))
    }
  }
  qtl_s_df <- do.call(rbind, qtl_s)
  op_index <- overlap_index(qtl_s_df)

  tr <- unique(qtl_s_df$Traits)
  corr_sta_s <- list()
  for(t in tr){
    tr_corr <- subset(op_index, Trait1==t | Trait2==t)
    ovlp_mean <- mean(tr_corr$overlap_index)
    corr_sta_s[[t]] <- data.frame(trait=t, overlap_mean=ovlp_mean)
  }
  corr_sta_df_s <- do.call(rbind, corr_sta_s)
  overlap_permutation[[p]] <- corr_sta_df_s
}
permutation <- do.call(cbind, overlap_permutation)

permutation <- permutation[,seq(2,2000,2)]

th_al <- round(apply(permutation, 1, function(x){quantile(x, 0.95)}),3)
th_al <- as.data.frame(th_al)

corr_sta[,2:3] <- round(corr_sta[,2:3], 3)

corr_sta <- merge(corr_sta, th_al, by="row.names")

#write.csv(corr_sta, file="240221_overlap_idx_avg.csv", row.names=F, quote=F)

###overlap p-value estimate
ANL_rand <- permutation["ANL",]
mean(ANL_rand >= 0.076)

FH_rand <- permutation["FH",] 
mean(FH_rand >= 0.103)

p_value <- c()
for (tr in rownames(permutation)){
  sim_res <- permutation[tr,]
  om <- corr_sta[corr_sta$trait == tr, 4]
  p <- mean(sim_res >= om)
  p_value <- c(p_value, p) 
}  

p_value_df <- data.frame(trait = rownames(permutation), p_value = p_value)
corr_sta <- merge(corr_sta, p_value_df, by="trait")

write.csv(corr_sta, file="250126_overlap_idx_avg.csv", row.names=F, quote=F)

###random test plot
CL_rand <- t(permutation["CL", ])[,1]
CL_rand_plot <- ggplot() +
  geom_density(fill="gray60", aes(x=CL_rand), alpha=0.2, color="black") +
  xlab("Culm length") + ylab("Density") + 
  scale_y_continuous(expand = c (0, 0)) +
  scale_x_continuous(limits=c(0,0.12), breaks=c(0.03,0.06, 0.09, 0.12), labels=c(0.03,0.06, 0.09, 0.12)) +
  geom_vline(xintercept=quantile(CL_rand, 0.95), color="red") +
  geom_vline(xintercept=0.103, color="blue") +
  theme_classic() +
  theme(axis.title=element_text(size=14, family="Times New Roman"),
        axis.text=element_text(size=12, family="Times New Roman"))

FH_rand <- t(permutation["FH", ])[,1]
FH_rand_plot <- ggplot() +
  geom_density(fill="gray60", aes(x=FH_rand), alpha=0.2, color="black") +
  xlab("First heading") + ylab("") + 
  scale_y_continuous (expand = c (0, 0)) +
  scale_x_continuous(limits=c(0,0.12), breaks=c(0.03,0.06, 0.09, 0.12), labels=c(0.03,0.06, 0.09, 0.12)) +
  geom_vline(xintercept=quantile(FH_rand, 0.95), color="red") +
  geom_vline(xintercept=0.103, color="blue") +
  theme_classic() +
  theme(axis.title=element_text(size=16, family="Times New Roman"),
        axis.text=element_text(size=14, family="Times New Roman"))

FLA_rand <- t(permutation["FLA", ])[,1]
FLA_rand_plot <- ggplot() +
  geom_density(fill="gray60", aes(x=FLA_rand), alpha=0.2, color="black") +
  xlab("Flag leaf attitude") + ylab("") + 
  scale_y_continuous (expand = c (0, 0)) +
  scale_x_continuous(limits=c(0,0.12), breaks=c(0.03,0.06, 0.09, 0.12), labels=c(0.03,0.06, 0.09, 0.12)) +
  geom_vline(xintercept=quantile(FLA_rand, 0.95), color="red") +
  geom_vline(xintercept=0.103, color="blue") +
  theme_classic() +
  theme(axis.title=element_text(size=14, family="Times New Roman"),
        axis.text=element_text(size=12, family="Times New Roman"))

GL_rand <- t(permutation["GL", ])[,1]
GL_rand_plot <- ggplot() +
  geom_density(fill="gray60", aes(x=GL_rand), alpha=0.2, color="black") +
  xlab("Grain length") + ylab("Density") + 
  scale_y_continuous (expand = c (0, 0)) +
  scale_x_continuous(limits=c(0,0.12), breaks=c(0.03,0.06, 0.09, 0.12), labels=c(0.03,0.06, 0.09, 0.12)) +
  geom_vline(xintercept=quantile(GL_rand, 0.95), color="red") +
  geom_vline(xintercept=0.103, color="blue") +
  theme_classic() +
  theme(axis.title=element_text(size=14, family="Times New Roman"),
        axis.text=element_text(size=12, family="Times New Roman"))

PAF_rand <- t(permutation["PAF", ])[,1]
PAF_rand_plot <- ggplot() +
  geom_density(fill="gray60", aes(x=PAF_rand), alpha=0.2, color="black") +
  xlab("PA for the first stage") + ylab("") + 
  scale_y_continuous (expand = c (0, 0)) +
  scale_x_continuous(limits=c(0,0.12), breaks=c(0.03,0.06, 0.09, 0.12), labels=c(0.03,0.06, 0.09, 0.12)) +
  geom_vline(xintercept=quantile(PAF_rand, 0.95), color="red") +
  geom_vline(xintercept=0.103, color="blue") +
  theme_classic() +
  theme(axis.title=element_text(size=14, family="Times New Roman"),
        axis.text=element_text(size=12, family="Times New Roman"))

PAS_rand <- t(permutation["PAS", ])[,1]
PAS_rand_plot <- ggplot() +
  geom_density(fill="gray60", aes(x=PAS_rand), alpha=0.2, color="black") +
  xlab("PA for the second stage") + ylab("") + 
  scale_y_continuous (expand = c (0, 0)) +
  scale_x_continuous(limits=c(0,0.12), breaks=c(0.03,0.06, 0.09, 0.12), labels=c(0.03,0.06, 0.09, 0.12)) +
  geom_vline(xintercept=quantile(PAS_rand, 0.95), color="red") +
  geom_vline(xintercept=0.103, color="blue") +
  theme_classic() +
  theme(axis.title=element_text(size=14, family="Times New Roman"),
        axis.text=element_text(size=12, family="Times New Roman"))

svglite(file="250313_overlap_rand_distr.svg", height=5, width=9)
plot_grid(CL_rand_plot, FH_rand_plot, FLA_rand_plot, GL_rand_plot, PAF_rand_plot, PAS_rand_plot, nrow=2, ncol=3)
dev.off()


###overlap_index function
overlap_index <- function(qtl_df){
  tr <- unique(qtl_df$Traits)
  tr_ab <- combn(tr,2)
  QTL_OP <- list()
  Trait1 <- c(); Trait2 <- c()
  overlap_num=c();op_index=c()
  for(i in 1:ncol(tr_ab)){
    ta <- tr_ab[1,i]; tb <- tr_ab[2,i]
    Ta <- subset(qtl_df, Traits==ta)
    Tb <- subset(qtl_df, Traits==tb)
    Tab_num <- 0
    for(q in 1:nrow(Ta)){
      TaC <- Ta$chr[q]; TaS <- Ta$start[q]; TaE <- Ta$end[q]
      ol <- subset(Tb, chr==TaC & start<=TaE & end>=TaS)      
      Tab_num <- Tab_num + nrow(ol)
    }
    Trait1=c(Trait1,ta); Trait2 <- c(Trait2,tb)
    overlap_num <- c(overlap_num,Tab_num)
    index=Tab_num/(nrow(Ta)+nrow(Tb)-Tab_num)
    op_index=c(op_index,index)
  }  
  OP_num <- data.frame(Trait1=Trait1, Trait2=Trait2, overlap_num=overlap_num, overlap_index=op_index)
  return(OP_num)
}


###################
#####Hotspots######
library(dplyr)

qtl$chr <- factor(qtl$chr, levels=c("chr1","chr2","chr3","chr4","chr5","chr6",
                                    "chr7","chr8","chr9","chr10","chr11","chr12"))
qtl_od <- qtl  %>% arrange(chr, start_cM)
qtl_od$size <- qtl_od$end_cM-qtl_od$start_cM
qtl_od <- subset(qtl_od, size <= 20)

hst_l <- list()
h <- 1
for(c in unique(qtl_od$chr)){
  qtl_chr <- qtl_od %>% subset(chr==c)
  s <- 1; e <- 1; 
  for(i in 2:nrow(qtl_chr)){
    start_p <- qtl_chr$start_cM[s:i-1]; end_p <- qtl_chr$end_cM[s:i-1]
    start_i <- qtl_chr$start_cM[i]; end_i <- qtl_chr$end_cM[i]
    if(TRUE %in% c(start_p <= end_i) & TRUE %in% c(end_p >= start_i)){
      e <- i
    }else{
      if(e-s >= 2){
        block <- qtl_chr[s:e,]
        hst <- paste("hotspot", h, sep="")
        hst_l[[hst]] <- block
        h <- h+1; s <- i; e <-i
      }else{s <- i; e <-i}
    } 
    if(i==nrow(qtl_chr) & e-s >= 2){
      block <- qtl_chr[s:e,]
      hst <- paste("hotspot", h, sep="")
      hst_l[[hst]] <- block 
      h <- h+1 
    }
  }
}

for(i in 1:length(hst_l)){
  cat(names(hst_l)[i], file="240215_hotspot.tsv", sep="\n", append=TRUE)
  write.table(hst_l[[i]], file="240215_hotspot.tsv", sep="\t", append=TRUE, row.names=FALSE, quote=F)
}

for(i in 1:length(hst_l)){
  vec <- paste(paste("q",hst_l[[i]]$QTL, sep=""), collapse=", ")
  print(vec)
}

#################
qtl_number <- c()
for(c in unique(qtl_od$chr)){
  qtl_n <- qtl_od %>% subset(chr==c) %>% nrow()
  qtl_number <- c(qtl_number, qtl_n)
}  
chr <- unique(qtl_od$chr)

chr_length <- c(43270923,35937250,36413819,35502694,29958434,31248787,29697621,28443022,
                23012720,23207287,29021106,27531856)
all_length <- 373245519

chr_len_r <- chr_length/all_length

chisq_chr_len <- chisq.test(qtl_number, p=chr_len_r, B=10000)
chisq_chr_len$stdres

genet_len <- map %>% group_by(chr) %>% summarise(chr_genet_len=max(genetic_pos))
genet_len_r <- genet_len$chr_genet_len/sum(genet_len$chr_genet_len)

chrsq_genet_len_r <- chisq.test(qtl_number, p=genet_len_r, B=10000)
chrsq_genet_len_r$stdres

gene_1 <- c(3421,2737,3050,2135,1968,1896,1830,1602,1321,1230,1386,1367)
gene_2 <- c(1247,994,1031,805,719,859,696,666,525,609,653,532)
gene_3 <- c(93,71,76,55,35,64,40,44,39,38,50,48)
gene_4 <- c(211,162,203,145,133,122,101,104,79,95,100,82)
chr_gene_numb <- gene_1+gene_2+gene_3+gene_4

chr_gene_r <- chr_gene_numb/sum(chr_gene_numb)

chrsq_gene_r <- chisq.test(qtl_number, p=chr_gene_r, B=10000)

###
D <- c(); p <- c()
for(c in unique(qtl_od$chr)){
  qtl_genet_pos <- qtl_od %>% subset(chr==c) %>% subset(select="peak_cM")
  qtl_genet_pos <- c(qtl_genet_pos$peak_cM)
  genet_map <- map %>% subset(chr==c)
  max_pos <- max(genet_map$genetic_pos)
  ks <- ks.test(qtl_genet_pos, punif, 0, max_pos, simulate.p.value = TRUE, B = 10000)
  D <- c(D, ks$statistic[[1]])
  p <- c(p, ks$p.value)
}


chr <- c("chr1","chr2","chr3","chr4","chr5","chr6",
         "chr7","chr8","chr9","chr10","chr11","chr12")

cluster_test <- data.frame(chr=chr, size=chr_length, size_percent=round(chr_len_r,3), 
                           size_stdres=round(chisq_chr_len$stdres,3), 
                           genes=chr_gene_numb, gene_percent=round(chr_gene_r,3), 
                           genes_stdres=round(chrsq_gene_r$stdres,3),
                           QTLs=qtl_number, QTLs_percent=round(qtl_number/sum(qtl_number),3),
                           D=round(D,2), ks.p=round(p,3))

write.table(cluster_test, file="240215_cluster_test.tsv", quote=F, row.names=F, sep="\t")
