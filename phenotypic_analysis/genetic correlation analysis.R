library(sommer)
library(rrBLUP)
library(ggplot2) 
library(reshape2) 
library(svglite)
library(dplyr)

pheno <- read.csv("426_pheno.csv")
geno  <- read.csv("231023_RIL_geno.csv", row.names=1)
markp <- read.csv("231023_RIL_marker_position.csv")

names(pheno)[1] <- "id"

pheno$idd <- pheno$id; pheno$ide <- pheno$id
row.names(pheno) <- pheno$id

geno[geno=="a"] <- -1
geno[geno=="h"] <- 0
geno[geno=="b"] <- 1

rn <- row.names(geno)
geno <- as.data.frame(lapply(geno, as.numeric))
row.names(geno) <- rn
geno <- as.matrix(geno)


# 1) 计算 K（基因组关系矩阵）
K <- A.mat(as.matrix(geno))  # rrBLUP::A.mat

# 2) 准备表型
pheno2 <- pheno
rownames(pheno2) <- pheno2$line
phe <- sort(names(pheno)[-c(1,21,22)])

tra_comb <- as.data.frame(combn(phe, 2))

results <- data.frame(trait1 = character(),
                      trait2 = character(),
                      rg = numeric(),
                      se = numeric(),
                      p = numeric(),
                      stringsAsFactors = FALSE)
for(i in 1:ncol(tra_comb)){
  tra_comb_i <- tra_comb[,i]
  tra1 <- tra_comb_i[1]; tra2 <- tra_comb_i[2]
  form <- as.formula(paste0("cbind(", tra1, ", ", tra2, ") ~ 1"))
  fit_i <- mmer(form,
            random = ~ vsr(id, Gu = K),
            rcov = ~ vsr(units),
            data = pheno2, dateWarning=FALSE)

  out <- get_rg_p(fit_i)
  results <- rbind(results, data.frame(trait1 = out$traits[1],
                                       trait2 = out$traits[2],
                                       rg = out$rg,
                                       se = out$se,
                                       p = out$p,
                                       stringsAsFactors = FALSE))
}


get_rg_p <- function(fit){
  sv_names <- names(fit$sigmaVector)
  idx_var1 <- which(grepl("^u:id\\.", sv_names) & grepl("-FH$|FH-FH", sv_names) & grepl("FH-FH|FH-.*FH", sv_names) ) 
  G <- fit$sigma$`u:id`
  trait_names <- colnames(G)
  name_v1 <- paste0("u:id.", trait_names[1], "-", trait_names[1])
  name_cov<- paste0("u:id.", trait_names[1], "-", trait_names[2])
  name_v2 <- paste0("u:id.", trait_names[2], "-", trait_names[2])
  sv_idx <- match(c(name_v1, name_cov, name_v2), sv_names)
  if(any(is.na(sv_idx))){
    stop("无法在 fit$sigmaVector 中找到对应的遗传分量名字。请检查 fit 对象的 sigmaVector 名称。")
  }
  
  v1  <- fit$sigmaVector[ sv_idx[1] ]
  cov <- fit$sigmaVector[ sv_idx[2] ]
  v2  <- fit$sigmaVector[ sv_idx[3] ]
  rg  <- cov / sqrt(v1 * v2)
  
  V_theta_full <- fit$sigmaSE
  V_theta <- V_theta_full[ sv_idx, sv_idx, drop = FALSE ]
  
  denom <- sqrt(v1 * v2)
  g1 <- -0.5 * cov / (v1 * denom)  
  g2 <- 1 / denom                   
  g3 <- -0.5 * cov / (v2 * denom)   
  g <- matrix(c(g1, g2, g3), ncol = 1)
  
  var_rg <- as.numeric(t(g) %*% V_theta %*% g)
  
  if(var_rg < 0 && var_rg > -1e-12) var_rg <- 0
  if(var_rg < 0){
    warning("估计的 rg 方差为负，可能是 asymptotic var-cov 数值问题；设置 SE = NA。")
    se_rg <- NA
    pval <- NA
  } else {
    se_rg <- sqrt(var_rg)
    z <- rg / se_rg
    pval <- 2 * (1 - pnorm(abs(z)))
  }
  
  return(list(rg = rg, se = se_rg, p = pval, v1 = v1, cov = cov, v2 = v2, traits = trait_names))
}


###
mat <- acast(results, trait1 ~ trait2, value.var = "rg")
mat_full <- mat
mat_full[lower.tri(mat_full)] <- t(mat_full)[lower.tri(mat_full)]

###
df_plot <- melt(mat_full, varnames = c("Trait1", "Trait2"), value.name = "rg")

ggplot(df_plot, aes(x = Trait1, y = Trait2, fill = rg)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limit = c(-1,1), 
                       name = "Genetic\ncorrelation") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggplot(df_plot, aes(x = Trait1, y = Trait2, fill = rg)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rg)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limit = c(-1,1), 
                       name = "Genetic\ncorrelation") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())



trait_al <- sort(unique(c(results$trait1, results$trait2)))

results2 <- results
results2$trait1 <- factor(results2$trait1, levels=rev(trait_al))
results2$trait2 <- factor(results2$trait2, levels=trait_al)

results3 <- results2[,c(2,1,3,4,5)]
names(results3)[c(1,2)] <- c("trait1","trait2")
res_add <- data.frame(trait1=phe, trait2=phe, rg=1, se=0, p=0)

sig <-  apply(results3, 1, function(x){
	  #value <- round(as.numeric(x[3]),2)
        pvalue <- as.numeric(x[5])
	  if(pvalue > 0.05){tx <- " "}
	  if(pvalue <= 0.05 & pvalue > 0.01){tx <- "*"}
	  if(pvalue <= 0.01 & pvalue > 0.001){tx <- "**"}
	  if(pvalue <= 0.001){tx <- "***"}
	  return(tx)
	})
results3$sig <- sig

write.table(results3, file="251007_genetic_correlation.tsv", sep="\t", quote=F, row.names=F)

length(which(results3$p <= 0.05))

results_n <- rbind(results2, results3, res_add)

results_n <- arrange(results_n, "trait1", "trait1")

p <- ggplot(results_n, aes(x=trait2, y=trait1, fill=rg)) +
  geom_tile() +
  scale_fill_gradient2(low="#B2182B", mid="white", high="#2166AC", limits=c(-1, 1)) +
  annotate("text", x=results2$trait1, y=results2$trait2, label=round(results2$rg,2),
           size=4) +
  annotate("text", x=results3$trait1, y=results3$trait2, label=results3$sig,
           size=4) +
  theme_bw() + 
  theme(panel.border = element_blank(),
        axis.ticks=element_blank(), 
        axis.title=element_blank(),
        axis.text=element_text(size=14, family="Times"),
        axis.text.x=element_text(angle=90, vjust=0.2))

svglite(file="250928_genetic_correlation.svg", height=8, width=9)
p
dev.off()

###
pheno_cor <- read.csv("250928_phenotypic_correlation.csv", head=T)
overlap_idx <- read.csv("240207.qtl_overlap_corr_all.csv", head=T)
overlap_idx <- arrange(overlap_idx, Var1, Var2)
overlap_idx <- overlap_idx[,-1]
names(overlap_idx)[1:2] <- c("trait1","trait2")
 
cor_all <- full_join(overlap_idx, results, by=c("trait1","trait2"))

cor.test(cor_all$corr_value, cor_all$rg)
cor.test(cor_all$overlap_index, abs(cor_all$rg))

plot(cor_all$overlap_index, cor_all$rg)

write.csv(results, file="250928_genetic_correlation.csv", quote=F, row.names=F)

######
library(ggcorrplot)

colnames <- c("ANL","AWL","CD","CH","CL","FH","FLA","FLL","FLW","GL",
              "GWE","GWI","PC_I","PC_II","PC_III","PE","PL","PS","SN")
pheno_cor <- pheno[,colnames]
cor_res <- cor(pheno_cor, use="complete.obs", method="pearson")
cor_long <- melt(cor_res, varnames = c("trait1", "trait2"), value.name = "cor")

cor.sig <- cor_pmat(pheno_cor, method="pearson", alternative = "two.sided")
corr_c_sig <- ifelse(cor.sig > 0.05, "",
                     ifelse(cor.sig <= 0.05 & cor.sig > 0.01, "*",
                            ifelse(cor.sig <= 0.01 & cor.sig > 0.001, "**", "***")))
cor_sig_long <- melt(corr_c_sig, varnames = c("trait1", "trait2"), value.name = "cor_sig")

rg_df <- results
rg_df$rg_sig <- ifelse(rg_df$p > 0.05, "",
                       ifelse(rg_df$p <= 0.05 & rg_df$p > 0.01, "*",
                              ifelse(rg_df$p <= 0.01 & rg_df$p > 0.001, "**", "***")))



plot_df <- merge(rg_df, cor_long, by=c("trait1", "trait2"), all.x=TRUE)
plot_df <- merge(plot_df, cor_sig_long, by=c("trait1", "trait2"), all.x=TRUE)

df <- data.frame(trait1=colnames, trait2=colnames, rg=rep(1,19), se=rep(0,19),p=rep(0,19),
                 rg_sig=rep("",19), cor=rep(1,19),cor_sig=rep("",19))

plot_df <- rbind(plot_df, df)

plot_df <- arrange(plot_df, trait1, trait2)
plot_df$trait1 <- factor(plot_df$trait1, level=rev(sort(unique(plot_df$trait1))))
plot_df$trait2 <- factor(plot_df$trait2, level=unique(plot_df$trait2))

p <- ggplot() +
  geom_tile(data=plot_df, aes(x=trait2, y=trait1, fill=rg)) +
  geom_tile(data=plot_df, aes(x=trait1, y=trait2, fill=cor)) +
  scale_fill_gradient2(low="#B2182B", mid="white", high="#2166AC", limits=c(-1, 1)) +
  annotate("text", x=plot_df$trait2, y=plot_df$trait1, label=plot_df$rg_sig,
           size=4) +
  annotate("text", x=plot_df$trait1, y=plot_df$trait2, label=plot_df$cor_sig,
           size=4) +
  theme_bw() + 
  theme(panel.border = element_blank(),
        axis.ticks=element_blank(), 
        axis.title=element_blank(),
        axis.text=element_text(size=14, family="Times"),
        axis.text.x=element_text(angle=90, vjust=0.2))

svglite(file="251016_genetic_correlation.svg", height=8, width=9)
p
dev.off()







