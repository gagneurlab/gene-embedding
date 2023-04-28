# load packages
library(ggplot2)
library(data.table)
library(dplyr)
library(cowplot)
library(khroma)

# get path for metric table
metric_table_path <- snakemake@input[["metric_table_path"]] 
# get path for plots
box_plot_path <- snakemake@output[["box_plot_path"]]
box_rdata_path <- snakemake@output[["box_rdata_path"]]

#bar_plot_path <- snakemake@output[["bar_plot_path"]]
#bar_rdata_path <- snakemake@output[["bar_rdata_path"]]
# get min performance table path
min_performance_path <- snakemake@input[["min_performance_path"]]

label_dict <- snakemake@params[["emb_labels"]]
box_order <- snakemake@params[["box_order"]]

# load metric table
table_list <- lapply(metric_table_path, function(x){
  fread(x)
})

metric_table <- do.call(rbind, table_list)

# load min performance table
min_perf_df <- fread(min_performance_path)

# load select trait list
select_traits <- snakemake@params[["select_traits"]]

# filter for select traits if any
if (length(select_traits) > 0){
  metric_table <- filter(metric_table, disease %in% as.vector(select_traits))
  min_perf_df <- filter(min_perf_df, disease %in% as.vector(select_traits))
}


min_perf_df$label <- toupper(min_perf_df$label)

# New facet label names 
disease_labels <- min_perf_df$label
names(disease_labels) <- min_perf_df$disease

d_labeller <- as_labeller(disease_labels)


# new emb names in df
for(combi in names(label_dict)){
  metric_table$combination <- replace(metric_table$combination,
                                      metric_table$combination == combi,
                                      label_dict[[combi]])
}

# new emb names in order
for(i in seq(1, length(box_order))){
  box_order[[i]] <-label_dict[[box_order[[i]]]]
}


# merge au_prc values and faceting information for facet order to work
plot_df <- merge(metric_table, min_perf_df, by='disease')

# order according to box_order
plot_df$combination <- factor(plot_df$combination,
                            levels = as.vector(box_order),
                            ordered = TRUE)

# order according to number of genes per disease 
plot_df$disease <- factor(plot_df$disease,
                            levels = min_perf_df$disease,
                            ordered = TRUE)

# min performance values as factor
# plot_df$min_auprc <- factor(plot_df$min_auprc,
#                            levels = min_perf_df$disease,
#                            ordered = TRUE)

# create box plot
p <- ggplot(plot_df,  aes(y=au_prc, fill=combination)) +
  geom_boxplot(width=5) +
  geom_hline(aes(yintercept = min_auprc), color = 'black', linetype = 2) + 
  facet_wrap(~ disease, nrow=1, labeller=d_labeller) + 
  # scale_fill_mediumcontrast() +
  theme_cowplot(font_size = 18) + 
  background_grid(major = 'y') +
  theme(legend.pos = 'bottom',axis.ticks.x = element_blank(), axis.text.x = element_blank() ) + 
  guides(fill = guide_legend(nrow = 3)) +
  labs(fill='Embedding') +
  ylab("Average precision") +
  xlab(NULL)

# save plot
ggsave(filename=box_plot_path,
       plot=p, 
       width = 12,
       height = 8,
       dpi=400)

# save plot as Rdata object
saveRDS(p, file = box_rdata_path)

################## create subsetted plot for poster #######################

plot_df <- plot_df[grep("string", combination, ignore.case = T)]

p <- ggplot(plot_df,  aes(y=au_prc, fill=combination)) +
  geom_boxplot(width=5) +
  geom_hline(aes(yintercept = min_auprc), color = 'black', linetype = 2) + 
  facet_wrap(~ disease, nrow=1, labeller=d_labeller) + 
  scale_fill_mediumcontrast() +
  theme_cowplot(font_size = 18) + 
  background_grid(major = 'y') +
  theme(legend.pos = 'bottom',axis.ticks.x = element_blank(), axis.text.x = element_blank() ) + 
  guides(fill = guide_legend(nrow = 2)) +
  labs(fill='Embedding') +
  ylab("Average precision") +
  xlab(NULL)

# save plot
ggsave(filename=paste0(box_plot_path, "string_subset.png"),
       plot=p, 
       width = 7,
       height = 8,
       dpi=400)

# save plot as Rdata object
saveRDS(p, file = box_rdata_path)








##########################################################################


# create bar plot 
# metric_table <- filter(metric_table, fold == max(metric_table$fold))

#p <- ggplot(plot_df,  aes(x=combination, y=au_prc, fill=combination)) +
 # geom_col(position = 'dodge') + 
#  geom_hline(aes(yintercept = min_auprc), color = 'black', linetype = 2) + 
#  facet_wrap(~ disease, nrow=1, labeller=d_labeller) + 
#  scale_fill_mediumcontrast() +
#  theme_cowplot() + 
#  background_grid() +
#  theme(legend.pos = 'bottom',axis.ticks.x = element_blank(), axis.text.x = element_blank() ) + 
#  guides(fill = guide_legend(nrow = 2)) +
#  labs(fill='Embedding') +
#  ylab("auPRC")

# save plot
# ggsave(filename=bar_plot_path,
#       plot=p, 
#       width = 12,
#       height = 8,
#       dpi=400)

# save plot as Rdata object
# saveRDS(p, file = bar_rdata_path)
