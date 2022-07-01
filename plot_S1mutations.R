library(tidyverse)
library(treeio)
library(ggtree)
library(ggrepel)
library(scales)
library(RColorBrewer)

options(ignore.negative.edge=TRUE)

nwk_file <- "nextstrain_ncov_open_global_all-time_timetree.nwk"
tsv_file <- "nextstrain_ncov_open_global_all-time_metadata.tsv"

origin_strain <- "Wuhan-Hu-1/2019"

clade2name <- tibble(
  clade =  c("20I (Alpha, V1)", "20H (Beta, V2)",
               "21A (Delta)", "21I (Delta)", "21J (Delta)",
               "21K (Omicron)", "21L (Omicron)",
               "22C (Omicron)",
             "22A (Omicron)", "22B (Omicron)"
               ),
  name = c("Alpha", "Beta", "Delta", "Delta", "Delta",
          rep("Omicron", 3), "Omicron (BA.4)",
          "Omicron (BA.5)")
)

nstr <- read.newick(nwk_file)
df <- read_tsv(tsv_file) %>% 
  select(strain, date, clade_membership, region,
         S1_mutations, emerging_lineage)

df <- left_join(
  df, clade2name, 
  by = c("clade_membership" = "clade")
)

df <- df %>% 
  mutate(name = ifelse(strain == origin_strain, "Origine", name)) %>% 
  replace_na(list(S1_mutations = 0, name = "Other variants")) %>% 
  mutate(name = factor(name, levels = c("Origine", "Other variants",
                                      "Alpha", "Beta", "Delta",
                                      "Omicron", "Omicron (BA.4)",
                                      "Omicron (BA.5)")))

dd <- full_join(
  nstr,
  df, by = c("label" = "strain")
)


# annotation data.frame for ggrepel
dfa <- df %>% 
  filter(name %in%
           c("Origine", "Alpha", "Beta", "Delta", "Omicron")) %>% 
  arrange(date) %>% 
  group_by(name) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  mutate(
    nudgey = ifelse(name == "Omicron", -8, 5)
  ) 


## PLOT

ggs <- ggtree(
  dd, 
  aes(color=name), 
  yscale = "S1_mutations", 
  size = 0.2,
  mrsd="2022-06-04",
  show.legend = FALSE
  ) + 
  geom_tippoint(aes(colour=name)) +
  theme(legend.position = c(0.1, 0.9))+
  scale_size_identity() +
  scale_y_continuous(name = "Spike S1 mutations count",
                        expand = c(0.02, 0.6)) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
                      name = "",
                      breaks = levels(df$name), na.value = "lightgrey") +
  geom_text_repel(data = dfa,
                  aes(ggtree::Date2decimal(date), y = S1_mutations,
                      label = name, group = 1, colour = name),
                      nudge_y = dfa$nudgey,fontface = "bold",
             size = 5, show.legend = F, hjust = "right",
             direction = "y", segment.curvature = 0.15,
             segment.ncp = 4,
             segment.angle = 30
                  ) +
  scale_x_continuous(
    expand = c(0.01,0.05),
    n.breaks = 3,
    labels = number_format(accuracy = 1, big.mark =""))  +
  # guides(color = guide_legend(override.aes = list(size = 3),
  #                             direction = "vertical", ncol = 2)) +
  theme_minimal()

ggs
