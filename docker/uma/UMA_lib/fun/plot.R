#!/usr/bin/r

# file: plot.R
# description: data plotting functions
# last update: 23-09-24

plotMetrics <- function(df, cut_offs){
  return(df %>%
           melt() %>%
           ggplot(aes(x = ID, y = value, fill = variable)) +
            geom_bar(stat = "identity") +
            geom_hline(data = cut_offs, aes(yintercept = value), linetype = "dashed") +
            labs(fill = "Metric", x = "Sample", y = "Value") + 
            facet_grid(variable ~., scales = "free_y") +
            scale_fill_paletteer_d("nationalparkcolors::Acadia") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)))
}
