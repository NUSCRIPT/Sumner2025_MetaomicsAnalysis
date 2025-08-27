library(ggtree)
library(tidyverse)
library(ggsci)

tree <- read.tree("/Users/jacksumner/Library/CloudStorage/CloudMounter-QuestBAS_manuscript/bins/gtdbtk.bac120.decorated.tree")

ggplot(tree, aes(x, y)) + geom_tree(layout="circular") + theme_tree()

ggtree(tree,layout="circular") + 
  geom_tiplab(aes(angle=angle),size=.1,color='blue')

my_bins <- tree$tip.label[grep(pattern = "PT_", tree$tip.label,invert = T)]
reduced_tree <- tidytree::drop.tip(tree, my_bins)

bin_classification <- read_tsv("/Users/jacksumner/Library/CloudStorage/CloudMounter-QuestBAS_manuscript/bins/gtdbtk.bac120.summary.tsv")
bin_classification<- bin_classification %>% 
  separate(classification,into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep=";") %>% 
  mutate(label = user_genome)

reduced_tree <- left_join(reduced_tree,bin_classification, by = "label")

ggtree(reduced_tree,layout="circular",ladderize = T) + 
  geom_tiplab(aes(label=species, color=phylum),show.legend=F,
              hjust = -.1,
              offset = .25,
              size=2.5,
              align =T) + 
  geom_tippoint(shape=21, fill="grey", size=.5) +
  geom_hilight(node=118, fill="brown", alpha=.6,to.bottom = T) +
  # geom_text(aes(label=node), hjust=-.3, size=1, color="red") 
  scale_color_aaas() +
  theme_tree() +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(55,55,55,55), "mm"))

gheatmap(last_plot(), bin_classification %>% 
           select(label, phylum) %>%
           column_to_rownames("label"), 
         offset=0, 
         width=0.1, 
         colnames=FALSE, legend_title="phylum")  +
  scale_fill_aaas() +
  theme(legend.key.size = unit(.01, "in"),
        legend.box.spacing = unit(1, "in"),
        legend.title = element_text(size = 3),
        legend.text = element_text(size = 3),
        legend.position = "right",
        # legend.background = element_blank(),
        legend.box.margin = margin(20, 0, 0,0),
        legend.spacing = unit(10, "pt")
  )


ggsave("Desktop/bin_tree.pdf", units = "in", width= 9.60, height=9.60)
