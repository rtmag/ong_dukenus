h2az = read.table(pipe("grep -v -P '#|genes' h2azOnly_20bp_2kb_QNORM.rmat"),sep="\t")

ach2az = read.table(pipe("grep -v -P '#|genes' ach2azOnly_20bp_2kb_CPM.rmat"),sep="\t")
#########################################################################################################################
# 95% confidence interval
h2_DRB=colMeans(h2az)[1:200]
h2_ACTD=colMeans(h2az)[201:400]
h2_DMSO=colMeans(h2az)[401:600]

ac_DRB=colMeans(ach2az)[1:200]
ac_ACTD=colMeans(ach2az)[201:400]
ac_DMSO=colMeans(ach2az)[401:600]

h2_DRB_conf=apply(h2az[,1:200],2,function(x) t.test(x)$conf.int)
h2_ACTD_conf=apply(h2az[,201:400],2,function(x) t.test(x)$conf.int)
h2_DMSO_conf=apply(h2az[,401:600],2,function(x) t.test(x)$conf.int)

ac_DRB_conf=apply(ach2az[,1:200],2,function(x) t.test(x)$conf.int)
ac_ACTD_conf=apply(ach2az[,201:400],2,function(x) t.test(x)$conf.int)
ac_DMSO_conf=apply(ach2az[,401:600],2,function(x) t.test(x)$conf.int)

#
