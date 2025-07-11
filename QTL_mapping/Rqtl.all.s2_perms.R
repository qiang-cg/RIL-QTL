library(qtl)
library(snow)

setwd("/mnt/ge-jbod/qiangchenggen/datapool/RIL/snpbinner_2/F6_426_2nd_sdf_p1e-7/bin_5k/rqtl_sc2_perms")
load("cross.cg.RData")

cross<-read.cross(format ="csvr", file="bin_5k.csv", estimate.map=F, genotypes=c("a", "h", "b"),
                  alleles=c("a", "b"),na.strings=c("-", "NA"))
cross<-convert2bcsft(cross, BC.gen=0, F.gen=6)
newmap=est.map(cross, map.function="kosambi", tol=0.01,
               maxit=1000, n.cluster=4)
cross= replace.map(cross, newmap)
cross<-calc.genoprob(cross, map.function="kosambi")

###FH
set.seed(10001)
FH.sc2.perm1 <- scantwo(cross, pheno.col=1, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(10002)
FH.sc2.perm2 <- scantwo(cross, pheno.col=1, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(10003)
FH.sc2.perm3<- scantwo(cross, pheno.col=1, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(10004)
FH.sc2.perm4 <- scantwo(cross, pheno.col=1, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(10005)
FH.sc2.perm5 <- scantwo(cross, pheno.col=1, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(10006)
FH.sc2.perm6 <- scantwo(cross, pheno.col=1, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(10007)
FH.sc2.perm7 <- scantwo(cross, pheno.col=1, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(10008)
FH.sc2.perm8 <- scantwo(cross, pheno.col=1, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(10009)
FH.sc2.perm9 <- scantwo(cross, pheno.col=1, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(10010)
FH.sc2.perm10 <- scantwo(cross, pheno.col=1, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.FH.sc2_perms100.RData")

###CL
set.seed(11001)
CL.sc2.perm1 <- scantwo(cross, pheno.col=2, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(11002)
CL.sc2.perm2 <- scantwo(cross, pheno.col=2, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(11003)
CL.sc2.perm3<- scantwo(cross, pheno.col=2, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(11004)
CL.sc2.perm4 <- scantwo(cross, pheno.col=2, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(11005)
CL.sc2.perm5 <- scantwo(cross, pheno.col=2, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(11006)
CL.sc2.perm6 <- scantwo(cross, pheno.col=2, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(11007)
CL.sc2.perm7 <- scantwo(cross, pheno.col=2, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(11008)
CL.sc2.perm8 <- scantwo(cross, pheno.col=2, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(11009)
CL.sc2.perm9 <- scantwo(cross, pheno.col=2, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(11010)
CL.sc2.perm10 <- scantwo(cross, pheno.col=2, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.CL.sc2_perms100.RData")

###CH
set.seed(12001)
CH.sc2.perm1 <- scantwo(cross, pheno.col=3, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(12002)
CH.sc2.perm2 <- scantwo(cross, pheno.col=3, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(12003)
CH.sc2.perm3<- scantwo(cross, pheno.col=3, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(12004)
CH.sc2.perm4 <- scantwo(cross, pheno.col=3, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(12005)
CH.sc2.perm5 <- scantwo(cross, pheno.col=3, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(12006)
CH.sc2.perm6 <- scantwo(cross, pheno.col=3, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(12007)
CH.sc2.perm7 <- scantwo(cross, pheno.col=3, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(1208)
CH.sc2.perm8 <- scantwo(cross, pheno.col=3, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(12009)
CH.sc2.perm9 <- scantwo(cross, pheno.col=3, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(12010)
CH.sc2.perm10 <- scantwo(cross, pheno.col=3, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.CH.sc2_perms100.RData")

###CD
set.seed(13001)
CD.sc2.perm1 <- scantwo(cross, pheno.col=4, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(13002)
CD.sc2.perm2 <- scantwo(cross, pheno.col=4, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(13003)
CD.sc2.perm3<- scantwo(cross, pheno.col=4, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(13004)
CD.sc2.perm4 <- scantwo(cross, pheno.col=4, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(13005)
CD.sc2.perm5 <- scantwo(cross, pheno.col=4, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(13006)
CD.sc2.perm6 <- scantwo(cross, pheno.col=4, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(13007)
CD.sc2.perm7 <- scantwo(cross, pheno.col=4, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(13008)
CD.sc2.perm8 <- scantwo(cross, pheno.col=4, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(13009)
CD.sc2.perm9 <- scantwo(cross, pheno.col=4, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(13010)
CD.sc2.perm10 <- scantwo(cross, pheno.col=4, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.CD.sc2_perms100.RData")

###FLA
set.seed(14001)
FLA.sc2.perm1 <- scantwo(cross, pheno.col=5, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(14002)
FLA.sc2.perm2 <- scantwo(cross, pheno.col=5, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(14003)
FLA.sc2.perm3<- scantwo(cross, pheno.col=5, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(14004)
FLA.sc2.perm4 <- scantwo(cross, pheno.col=5, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(14005)
FLA.sc2.perm5 <- scantwo(cross, pheno.col=5, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(14006)
FLA.sc2.perm6 <- scantwo(cross, pheno.col=5, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(14007)
FLA.sc2.perm7 <- scantwo(cross, pheno.col=5, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(14008)
FLA.sc2.perm8 <- scantwo(cross, pheno.col=5, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(14009)
FLA.sc2.perm9 <- scantwo(cross, pheno.col=5, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(14010)
FLA.sc2.perm10 <- scantwo(cross, pheno.col=5, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.FLA.sc2_perms100.RData")

###FLL
set.seed(15001)
FLL.sc2.perm1 <- scantwo(cross, pheno.col=6, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(15002)
FLL.sc2.perm2 <- scantwo(cross, pheno.col=6, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(15003)
FLL.sc2.perm3<- scantwo(cross, pheno.col=6, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(15004)
FLL.sc2.perm4 <- scantwo(cross, pheno.col=6, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(15005)
FLL.sc2.perm5 <- scantwo(cross, pheno.col=6, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(15006)
FLL.sc2.perm6 <- scantwo(cross, pheno.col=6, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(15007)
FLL.sc2.perm7 <- scantwo(cross, pheno.col=6, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(15008)
FLL.sc2.perm8 <- scantwo(cross, pheno.col=6, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(15009)
FLL.sc2.perm9 <- scantwo(cross, pheno.col=6, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(15010)
FLL.sc2.perm10 <- scantwo(cross, pheno.col=6, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.FLL.sc2_perms100.RData")

###FLW
set.seed(16001)
FLW.sc2.perm1 <- scantwo(cross, pheno.col=7, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(16002)
FLW.sc2.perm2 <- scantwo(cross, pheno.col=7, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(16003)
FLW.sc2.perm3<- scantwo(cross, pheno.col=7, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(16004)
FLW.sc2.perm4 <- scantwo(cross, pheno.col=7, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(16005)
FLW.sc2.perm5 <- scantwo(cross, pheno.col=7, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(16006)
FLW.sc2.perm6 <- scantwo(cross, pheno.col=7, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(16007)
FLW.sc2.perm7 <- scantwo(cross, pheno.col=7, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(16008)
FLW.sc2.perm8 <- scantwo(cross, pheno.col=7, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(16009)
FLW.sc2.perm9 <- scantwo(cross, pheno.col=7, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(16010)
FLW.sc2.perm10 <- scantwo(cross, pheno.col=7, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.FLW.sc2_perms100.RData")

###PS
set.seed(17001)
PS.sc2.perm1 <- scantwo(cross, pheno.col=8, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(17002)
PS.sc2.perm2 <- scantwo(cross, pheno.col=8, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(17003)
PS.sc2.perm3<- scantwo(cross, pheno.col=8, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(17004)
PS.sc2.perm4 <- scantwo(cross, pheno.col=8, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(17005)
PS.sc2.perm5 <- scantwo(cross, pheno.col=8, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(17006)
PS.sc2.perm6 <- scantwo(cross, pheno.col=8, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(17007)
PS.sc2.perm7 <- scantwo(cross, pheno.col=8, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(17008)
PS.sc2.perm8 <- scantwo(cross, pheno.col=8, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(17009)
PS.sc2.perm9 <- scantwo(cross, pheno.col=8, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(17010)
PS.sc2.perm10 <- scantwo(cross, pheno.col=8, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.PS.sc2_perms100.RData")

###PE
set.seed(18001)
PE.sc2.perm1 <- scantwo(cross, pheno.col=9, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(18002)
PE.sc2.perm2 <- scantwo(cross, pheno.col=9, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(18003)
PE.sc2.perm3<- scantwo(cross, pheno.col=9, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(18004)
PE.sc2.perm4 <- scantwo(cross, pheno.col=9, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(18005)
PE.sc2.perm5 <- scantwo(cross, pheno.col=9, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(18006)
PE.sc2.perm6 <- scantwo(cross, pheno.col=9, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(18007)
PE.sc2.perm7 <- scantwo(cross, pheno.col=9, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(18008)
PE.sc2.perm8 <- scantwo(cross, pheno.col=9, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(18009)
PE.sc2.perm9 <- scantwo(cross, pheno.col=9, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(18010)
PE.sc2.perm10 <- scantwo(cross, pheno.col=9, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.PE.sc2_perms100.RData")

###PL
set.seed(19001)
PL.sc2.perm1 <- scantwo(cross, pheno.col=10, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(19002)
PL.sc2.perm2 <- scantwo(cross, pheno.col=10, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(19003)
PL.sc2.perm3<- scantwo(cross, pheno.col=10, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(19004)
PL.sc2.perm4 <- scantwo(cross, pheno.col=10, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(19005)
PL.sc2.perm5 <- scantwo(cross, pheno.col=10, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(19006)
PL.sc2.perm6 <- scantwo(cross, pheno.col=10, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(19007)
PL.sc2.perm7 <- scantwo(cross, pheno.col=10, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(19008)
PL.sc2.perm8 <- scantwo(cross, pheno.col=10, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(19009)
PL.sc2.perm9 <- scantwo(cross, pheno.col=10, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(19010)
PL.sc2.perm10 <- scantwo(cross, pheno.col=10, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.PL.sc2_perm100.RData")

###AWL
set.seed(20001)
AWL.sc2.perm1 <- scantwo(cross, pheno.col=11, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(20002)
AWL.sc2.perm2 <- scantwo(cross, pheno.col=11, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(20003)
AWL.sc2.perm3<- scantwo(cross, pheno.col=11, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(20004)
AWL.sc2.perm4 <- scantwo(cross, pheno.col=11, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(20005)
AWL.sc2.perm5 <- scantwo(cross, pheno.col=11, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(20006)
AWL.sc2.perm6 <- scantwo(cross, pheno.col=11, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(20007)
AWL.sc2.perm7 <- scantwo(cross, pheno.col=11, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(20008)
AWL.sc2.perm8 <- scantwo(cross, pheno.col=11, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(20009)
AWL.sc2.perm9 <- scantwo(cross, pheno.col=11, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(20010)
AWL.sc2.perm10 <- scantwo(cross, pheno.col=11, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.AWL.sc2_perm100.RData")

###ANL
set.seed(21001)
ANL.sc2.perm1 <- scantwo(cross, pheno.col=12, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(21002)
ANL.sc2.perm2 <- scantwo(cross, pheno.col=12, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(21003)
ANL.sc2.perm3<- scantwo(cross, pheno.col=12, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(21004)
ANL.sc2.perm4 <- scantwo(cross, pheno.col=12, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(21005)
ANL.sc2.perm5 <- scantwo(cross, pheno.col=12, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(21006)
ANL.sc2.perm6 <- scantwo(cross, pheno.col=12, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(21007)
ANL.sc2.perm7 <- scantwo(cross, pheno.col=12, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(21008)
ANL.sc2.perm8 <- scantwo(cross, pheno.col=12, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(21009)
ANL.sc2.perm9 <- scantwo(cross, pheno.col=12, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(21010)
ANL.sc2.perm10 <- scantwo(cross, pheno.col=12, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.ANL.sc2_perm100.RData")

###SN
set.seed(22001)
SN.sc2.perm1 <- scantwo(cross, pheno.col=13, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(22002)
SN.sc2.perm2 <- scantwo(cross, pheno.col=13, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(22003)
SN.sc2.perm3 <- scantwo(cross, pheno.col=13, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(22004)
SN.sc2.perm4 <- scantwo(cross, pheno.col=13, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(22005)
SN.sc2.perm5 <- scantwo(cross, pheno.col=13, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(22006)
SN.sc2.perm6 <- scantwo(cross, pheno.col=13, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(22007)
SN.sc2.perm7 <- scantwo(cross, pheno.col=13, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(22008)
SN.sc2.perm8 <- scantwo(cross, pheno.col=13, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(22009)
SN.sc2.perm9 <- scantwo(cross, pheno.col=13, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(22010)
SN.sc2.perm10 <- scantwo(cross, pheno.col=13, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.SN.sc2_perm100.RData")

###PC_1
set.seed(23001)
PC_1.sc2.perm1 <- scantwo(cross, pheno.col=14, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(23002)
PC_1.sc2.perm2 <- scantwo(cross, pheno.col=14, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(23003)
PC_1.sc2.perm3<- scantwo(cross, pheno.col=14, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(23004)
PC_1.sc2.perm4 <- scantwo(cross, pheno.col=14, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(23005)
PC_1.sc2.perm5 <- scantwo(cross, pheno.col=14, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(23006)
PC_1.sc2.perm6 <- scantwo(cross, pheno.col=14, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(23007)
PC_1.sc2.perm7 <- scantwo(cross, pheno.col=14, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(23008)
PC_1.sc2.perm8 <- scantwo(cross, pheno.col=14, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(23009)
PC_1.sc2.perm9 <- scantwo(cross, pheno.col=14, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(23010)
PC_1.sc2.perm10 <- scantwo(cross, pheno.col=14, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.PC_1.sc2_perm100.RData")

###PC_2
set.seed(24001)
PC_2.sc2.perm1 <- scantwo(cross, pheno.col=15, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(24002)
PC_2.sc2.perm2 <- scantwo(cross, pheno.col=15, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(24003)
PC_2.sc2.perm3<- scantwo(cross, pheno.col=15, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(24004)
PC_2.sc2.perm4 <- scantwo(cross, pheno.col=15, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(24005)
PC_2.sc2.perm5 <- scantwo(cross, pheno.col=15, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(24006)
PC_2.sc2.perm6 <- scantwo(cross, pheno.col=15, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(24007)
PC_2.sc2.perm7 <- scantwo(cross, pheno.col=15, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(24008)
PC_2.sc2.perm8 <- scantwo(cross, pheno.col=15, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(24009)
PC_2.sc2.perm9 <- scantwo(cross, pheno.col=15, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(24010)
PC_2.sc2.perm10 <- scantwo(cross, pheno.col=15, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.PC_2.sc2_perm100.RData")

###PC_3
set.seed(25001)
PC_3.sc2.perm1 <- scantwo(cross, pheno.col=16, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(25002)
PC_3.sc2.perm2 <- scantwo(cross, pheno.col=16, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(25003)
PC_3.sc2.perm3<- scantwo(cross, pheno.col=16, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(25004)
PC_3.sc2.perm4 <- scantwo(cross, pheno.col=16, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(25005)
PC_3.sc2.perm5 <- scantwo(cross, pheno.col=16, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(25006)
PC_3.sc2.perm6 <- scantwo(cross, pheno.col=16, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(25007)
PC_3.sc2.perm7 <- scantwo(cross, pheno.col=16, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(25008)
PC_3.sc2.perm8 <- scantwo(cross, pheno.col=16, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(25009)
PC_3.sc2.perm9 <- scantwo(cross, pheno.col=16, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(25010)
PC_3.sc2.perm10 <- scantwo(cross, pheno.col=16, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.PC_3.sc2_perm100.RData")

###GWE
set.seed(26001)
GWE.sc2.perm1 <- scantwo(cross, pheno.col=17, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(26002)
GWE.sc2.perm2 <- scantwo(cross, pheno.col=17, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(26003)
GWE.sc2.perm3<- scantwo(cross, pheno.col=17, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(26004)
GWE.sc2.perm4 <- scantwo(cross, pheno.col=17, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(26005)
GWE.sc2.perm5 <- scantwo(cross, pheno.col=17, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(26006)
GWE.sc2.perm6 <- scantwo(cross, pheno.col=17, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(26007)
GWE.sc2.perm7 <- scantwo(cross, pheno.col=17, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(26008)
GWE.sc2.perm8 <- scantwo(cross, pheno.col=17, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(26009)
GWE.sc2.perm9 <- scantwo(cross, pheno.col=17, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(26010)
GWE.sc2.perm10 <- scantwo(cross, pheno.col=17, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL426.GWE.sc2_perm100.RData")

###GL
set.seed(27001)
GL.sc2.perm1 <- scantwo(cross, pheno.col=18, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(27002)
GL.sc2.perm2 <- scantwo(cross, pheno.col=18, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(27003)
GL.sc2.perm3<- scantwo(cross, pheno.col=18, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(27004)
GL.sc2.perm4 <- scantwo(cross, pheno.col=18, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(27005)
GL.sc2.perm5 <- scantwo(cross, pheno.col=18, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(27006)
GL.sc2.perm6 <- scantwo(cross, pheno.col=18, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(27007)
GL.sc2.perm7 <- scantwo(cross, pheno.col=18, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(27008)
GL.sc2.perm8 <- scantwo(cross, pheno.col=18, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(27009)
GL.sc2.perm9 <- scantwo(cross, pheno.col=18, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(27010)
GL.sc2.perm10 <- scantwo(cross, pheno.col=18, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL427.GL.sc2_perm100.RData")

###GWI
set.seed(28001)
GWI.sc2.perm1 <- scantwo(cross, pheno.col=19, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(28002)
GWI.sc2.perm2 <- scantwo(cross, pheno.col=19, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(28003)
GWI.sc2.perm3 <- scantwo(cross, pheno.col=19, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(28004)
GWI.sc2.perm4 <- scantwo(cross, pheno.col=19, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(28005)
GWI.sc2.perm5 <- scantwo(cross, pheno.col=19, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(28006)
GWI.sc2.perm6 <- scantwo(cross, pheno.col=19, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(28007)
GWI.sc2.perm7 <- scantwo(cross, pheno.col=19, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(28008)
GWI.sc2.perm8 <- scantwo(cross, pheno.col=19, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(28009)
GWI.sc2.perm9 <- scantwo(cross, pheno.col=19, model="normal", method="hk", n.perm=100, n.cluster=10)

set.seed(28010)
GWI.sc2.perm10 <- scantwo(cross, pheno.col=19, model="normal", method="hk", n.perm=100, n.cluster=10)

save.image("RIL428.GWI.sc2_perm100.RData")
