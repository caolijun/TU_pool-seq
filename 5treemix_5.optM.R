library(OptM)
library(RColorBrewer)
library(R.utils)
source("/work/share/kdy_caolijun/software/treemix-1.13/src/plotting_funcs.R")
linear = optM("./10mig_100repeat") ##文件夹名
#linear = optM("./10mig_100repeat",method = "linear") ##文件夹名
plot_optM(linear, method = "Evanno", plot = F, pdf = "OptM.pdf")
#plot_optM(linear, method = "linear", plot = F, pdf = "OptM_linear.pdf")

# 绘制Treemix graph
#for (i in 1:10){
#  pdf(paste('migration_1m_', i, '.pdf', sep = ""), width = 8, height = 7)
#  plot_tree(paste('./10mig_100repeat/mite_tree_1mig_bootstrap_', i, sep=""))
#  plot_resid(paste('./10mig_100repeat/mite_tree_1mig_bootstrap_', i, sep=""), 'popnames.txt')
#  dev.off()
#}
#with lowest lk
#pdf(paste('migration_1m_', 3, '.pdf', sep = ""), width = 8, height = 7)
#plot_tree(paste('./opt1mig_1000repeat/mite_tree_1mig_bootstrap_', 3, sep=""))
#plot_resid(paste('./opt1mig_1000repeat/mite_tree_1mig_bootstrap_', 3, sep=""), 'popnames.txt')
#dev.off()
