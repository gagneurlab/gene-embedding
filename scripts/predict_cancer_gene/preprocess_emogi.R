

library(data.table)

logit <- function(x){
    return(log(x) - log(1 - x))
}


dt <- fread(snakemake@input[[1]])
dt <- dt[, .(gene_id = ID, label, pred = Mean_Pred)]
dt[, pred := logit(pred * 0.99 + 0.005)]

fwrite(dt,
       snakemake@output[[1]], 
       sep = '\t')
