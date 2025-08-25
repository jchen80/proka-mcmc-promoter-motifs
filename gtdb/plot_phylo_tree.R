library(ggtree)
library(treeio)
library(tidyverse)

#tree <- read.tree("nap/paxdb/bacteria/Bacteria_paxdb_mass_molar_nap_TF_abundance_metadata_ogt_all_final_tree.nwk")
tree <- read.tree("gtdb/bacteria_tree.nwk")
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- str_trim(tree$tip.label)

nap_data <- read_csv("nap/bacteria/bacteria_chromatin_abundance.csv")
#nap_data$Species <- gsub("=", " ", nap_data$Species)
#nap_data$Species <- gsub("_", " ", nap_data$Species)
#nap_data$Species <- gsub("\\(", " ", nap_data$Species)
#nap_data$Species <- gsub("\\)", " ", nap_data$Species)
#nap_data$Species <- str_trim(nap_data$Species)

# Basic tree plot to get data frame with tip positions
p <- ggtree(tree) + theme_tree2()
tree_data <- p$data

# Join NAP data to tree data (only tips)
tip_data <- tree_data %>%
  filter(isTip) %>%
  left_join(nap_data, by = c("label" = "species_reformatted")) %>%
  filter(!is.na(nap_noFerritin))

max_x <- max(tree_data$x)
scale_factor <- 0.2
offset <- 0.02 * max_x
vertical_label_offset <- 0.35  # vertical offset for species labels

# Create a data frame for tip labels with vertical offset
tip_labels_offset <- tip_data %>%
  mutate(y_label = y + vertical_label_offset)  # shift species labels up

# Now build the plot
p2 <- ggtree(tree) + 
  theme_tree2() +
  
  # Bars at original y positions
  geom_segment(data = tip_data,
               aes(x = x + offset, 
                   xend = x + offset + nap_noFerritin * max_x * scale_factor,
                   y = y, yend = y),
               color = "steelblue", linewidth = 3) +
  
  # NAP numeric labels aligned with bars
  geom_text(data = tip_data,
            aes(x = x + offset + nap_noFerritin * max_x * scale_factor + 0.01 * max_x,
                y = y, label = sprintf("%.2f", nap_noFerritin)),
            size = 3, hjust = 0) +
  
  # Species labels shifted up vertically
  geom_text(data = tip_labels_offset,
            aes(x = x, y = y_label, label = label),
            size = 3, hjust = 1)  # hjust=1 aligns labels to right of tip

print(p2)


