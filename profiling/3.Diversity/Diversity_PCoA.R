# alpha and beta diversity


library("tidyverse")
library("phyloseq")
library("vegan")



load("cache/ps_core5t.rda")

## plot shannon diversity

ps_core5t


sample_data <- sample_data(ps_core5t)


rch_4085 <- estimate_richness(ps_core5t, measures = c("Shannon"))
                              
                              
rch_4085_rhz <- rch_4085 %>%
  rownames_to_column(var="Sample_ID") %>%
  left_join(sample_data) %>%
  filter(sample_type == "RHZ")
                                


colors = c("#743282", "#dc9200")


shannon_plot <- ggplot(rch_4085_rhz, aes(x = state, y = Shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = state), size=1, width=0.2) +
  #stat_compare_means(aes(group = state), label = "p.format") +
  scale_color_manual(values=colors) +
  ylab("Shannon Diversity") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

shannon_plot





#### Constrained ordination, 4085 ASVs ####



## https://rdrr.io/rforge/vegan/man/adonis.html
# convert to relative abundance, aka Total Sum Normalization [TSS] and log transform ASV counts
ps <- transform_sample_counts(ps_core5t, function(x) log(x / sum(x) +0.001))
sdat <- data.frame(sample_data(ps))
colnames(sdat)[2] <- "location"
sdat <- filter(sdat, sample_type != "SOL" & location != "Salema Myler")
# exclude Bulk Soil
ps <- subset_samples(ps, !sample_type %in% c("SOL") & name != "Salema Myler")

## get distance matrix
dm <- distance(ps, method = "wunifrac")


### replace missing values with median
#dm[is.na(dm)] <- median(dm, na.rm = TRUE)
#save(dm, file='cache/distance_matrix_wunifrac.rda')


## calculate bray-custis distance matrix
bray.cap.whole <- capscale(as.dist(dm) ~ state + location + genotype, data = sdat, add = T, na.action = na.exclude)
#save(bray.cap.whole, file='cache/bray.cap.whole.rda')


## Permutation ANOVA
permanova <- adonis(dm~state + location + genotype,data = sdat,add = T)

#save(permanova, file = "cache/permanova.rda")

caps <- rownames_to_column(data.frame(scores(bray.cap.whole)$sites), var = "Sample_ID")

bray.cap.whole.axes <- left_join(sdat, caps)



percent_explained <- bray.cap.whole$CCA$eig / sum(bray.cap.whole$CCA$eig) * 100
percent_explained

colors = c("#743282", "#dc9200")

bray.cap.whole.axes <- bray.cap.whole.axes %>%
  mutate(species = ifelse(genotype %in% c("B73", "Mo17"), "Maize", "Soybean"))

bray.cap.whole.axes$genotype <- factor(bray.cap.whole.axes$genotype, levels = c("B73", "BTx623", "Mo17", "K037"))

ord_plot <- ggplot(bray.cap.whole.axes, aes(x = CAP1*10, y = CAP2*10, color = state)) +
  geom_point(size=2) +
  facet_wrap(~genotype, nrow = 2) +
  scale_color_manual(values=colors) +
  labs(x = "Constrained PCo1 (40.3%)", y = "Constrained PCo2 (17.1%)") +
  theme_bw() +
  theme(panel.spacing = unit(1, "lines"))

ord_plot


#unique(bray.cap.whole.axes$location)


### bulk soil



ps <- transform_sample_counts(ps_core5t, function(x) log(x / sum(x) +0.001))
ps <- subset_samples(ps, !sample_type %in% c("RHZ") & name != "Salema Myler")
sdat <- data.frame(sample_data(ps))
colnames(sdat)[2] <- "location"

## get distance matrix
dm <- distance(ps, method = "wunifrac")


### replace missing values with median
#dm[is.na(dm)] <- median(dm, na.rm = TRUE)
#save(dm, file='cache/distance_matrix_wunifrac.rda')


## calculate bray-custis distance matrix
bray.cap.whole <- capscale(as.dist(dm) ~ state + location, data = sdat, add = T, na.action = na.exclude)
#save(bray.cap.whole, file='cache/bray.cap.whole.rda')


## Permutation ANOVA
permanova <- adonis(dm~state + location,data = sdat,add = T)

#save(permanova, file = "cache/permanova.rda")

caps <- rownames_to_column(data.frame(scores(bray.cap.whole)$sites), var = "Sample_ID")

bray.cap.whole.axes <- left_join(sdat, caps)



percent_explained <- bray.cap.whole$CCA$eig / sum(bray.cap.whole$CCA$eig) * 100
percent_explained

colors = c("#743282", "#dc9200")

bray.cap.whole.axes <- bray.cap.whole.axes %>%
  mutate(species = ifelse(genotype %in% c("B73", "Mo17"), "Maize", "Soybean"))

bray.cap.whole.axes$genotype <- factor(bray.cap.whole.axes$genotype, levels = c("B73", "BTx623", "Mo17", "K037"))

ord_plot <- ggplot(bray.cap.whole.axes, aes(x = CAP1*10, y = CAP2*10, color = state)) +
  geom_point(size=2) +
  #facet_wrap(~genotype, nrow = 2) +
  scale_color_manual(values=colors) +
  labs(x = "Constrained PCo1 (40.3%)", y = "Constrained PCo2 (17.1%)") +
  theme_bw() +
  theme(panel.spacing = unit(1, "lines"))

ord_plot



