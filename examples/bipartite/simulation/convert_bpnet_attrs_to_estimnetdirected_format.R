#!/usr/bin/env Rscript
#
# convert_bpnet_attrs_to_estimnetdirected_format.R
#
# Read binary, categorical, and contiuous attributes in BPNet format
# (where attributes for mode A, P, and both (AP) nodes are in separate files)
# and write in cwd in EStimNEtDirected format, where they are in the same
# file (stil separate for binary, continuous, categorical types though)
# but with NA for nodes in the mode to which they are not applicable.
#
# Writes output to files in cwd (WARNING: overwrites):
#
# binattr_all.txt  catattr_all.txt  conattr_all.txt
#
# ADS 30 May 2022
#

binattrA <- read.table('binattrA.txt', header=T)
binattrP <- read.table('binattrP.txt', header=T)
binattrAP <- read.table('binattrAP.txt', header=T)
conattrA <- read.table('conattrA.txt', header=T)
conattrP <- read.table('conattrP.txt', header=T)
conattrAP <- read.table('conattrAP.txt', header=T)
catattrA <- read.table('catattrA.txt', header=T)
catattrP <- read.table('catattrP.txt', header=T)
catattrAP <- read.table('catattrAP.txt', header=T)

num_A <- nrow(binattrA)
num_P <- nrow(binattrP)
N <- num_A + num_P
stopifnot(nrow(binattrAP) == N)

binattr <- data.frame(binattrA = c(binattrA$binattrA,  rep(NA, num_P)), binattrP = c(rep(NA, num_A), binattrP$binattrP),  binattrAP = binattrAP$binattrAP)
catattr <- data.frame(catattrA = c(catattrA$catattrA,  rep(NA, num_P)), catattrP = c(rep(NA, num_A), catattrP$catattrP),  catattrAP = catattrAP$catattrAP)
conattr <- data.frame(conattrA = c(conattrA$conattrA,  rep(NA, num_P)), conattrP = c(rep(NA, num_A), conattrP$conattrP),  conattrAP = conattrAP$conattrAP)

summary(binattr)
summary(catattr)
summary(conattr)

write.table(binattr, 'binattr_all.txt', col.names=T,row.names=F, quote=F)
write.table(catattr, 'catattr_all.txt', col.names=T,row.names=F, quote=F)
write.table(conattr, 'conattr_all.txt', col.names=T,row.names=F, quote=F)
