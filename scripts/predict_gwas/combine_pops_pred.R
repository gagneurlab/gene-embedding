

library(data.table)

tables <- lapply(snakemake@input, fread)
table <- rbindlist(tables)

fwrite(table, snakemake@output[['pred']])
