source("/work/share/kdy_caolijun/software/treemix-1.13/src/plotting_funcs.R")
# Get the variation explained by each model
data <- matrix(nrow=1000, ncol=1) 
for (i in 1:1000) {
data[i,] <- get_f(paste0("0mig_SCCD_1000repeat/mite_tree_1mig_bootstrap_",i))
}
mean(data[,1])
median(data[,1])
