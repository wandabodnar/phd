# ==============================================================================
# STUDY: Fish assemblage diversity and habitat condition in the Thames Estuary
# AUTHOR: Wanda Bodnar
# ==============================================================================


# --- Load required libraries ---
library(iNEXT)        
library(glmmTMB)      
library(DHARMa)      
library(tidyverse)    
library(emmeans)      
library(vegan)        
library(indicspecies)  
library(pheatmap)     
library(eulerr)
library(ggrepel)
library(car)


# --- Data cleaning & preparation ---
setwd("C:/Users/bodna/Dropbox/Wanda Bodnar/Data/Practice/Stat")

species_raw <- read.csv("data3.csv")
salinity_raw <- read.csv("salinity.csv")

# Clean salinity: average the 3 measurements per location/season
salinity_clean <- salinity_raw %>%
  group_by(Location, Season) %>%
  summarise(Mean_Salinity = mean(ppt), 
            SD_Salinity = sd(ppt), .groups = "drop")

# Pool species counts by location/season/condition
species_cols <- colnames(species_raw)[4:ncol(species_raw)]
pooled_counts <- species_raw %>%
  group_by(Location, Season, Condition) %>%
  summarise(across(all_of(species_cols), sum), .groups = "drop") %>%
  mutate(Assemblage_ID = paste(Location, Season, sep = "_"))



# --- Alpha diversity (iNEXT Framework) ---
# Prepare matrix (Species = rows, Site-Seasons = columns)
inext_input <- t(pooled_counts %>% select(all_of(species_cols)))
colnames(inext_input) <- pooled_counts$Assemblage_ID

# Calculate Asymptotic estimates for q = 0, 1, 2
iNEXT_out <- iNEXT(inext_input, q = c(0, 1, 2), datatype = "abundance", nboot = 50)

# Factor levels for longitudinal plotting (Upstream -> Downstream)
location_levels <- c("Teddington", "Richmond", "Kew Bridge", "Chiswick", 
                     "Putney", "Fulham", "Vauxhall", "Battersea", 
                     "Billingsgate", "Greenwich", "Woolwich", "Barking")

# Hill Numbers 
ggplot(asy_est, aes(x = Location, y = Estimator, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", alpha = 0.8) +
  facet_grid(Diversity ~ Season, scales = "free_y") +
  scale_fill_manual(values = c("Degraded" = "#d7191c", "Created" = "#2c7bb6", 
                               "Restored" = "#fdae61", "Natural" = "#1a9641")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Fish diversity metrics (Hill Numbers)",
       x = "Sampling location (Upstream to Downstream)",
       y = "Asymptotic diversity estimate")

# Merge diversity estimates with environmental data
asy_est <- iNEXT_out$AsyEst %>%
  separate(Assemblage, into = c("Location", "Season"), sep = "_") %>%
  left_join(unique(pooled_counts[, c("Location", "Condition")]), by = "Location") %>%
  left_join(salinity_clean, by = c("Location", "Season")) %>%
  mutate(
    Location = factor(Location, levels = location_levels),
    Condition = factor(Condition, levels = c("Degraded", "Created", "Restored", "Natural")),
    Diversity = factor(Diversity, levels = c("Species richness", "Shannon diversity", "Simpson diversity")))

# Alpha diversity vs. salinity
ggplot(asy_est, aes(x = Mean_Salinity, y = Estimator, color = Condition, fill = Condition)) +
  geom_point(size = 3, alpha = 0.8, shape = 21, color = "black") +
  geom_smooth(data = subset(asy_est, Condition != "Natural"), 
              method = "lm", se = TRUE, alpha = 0.1) +
  facet_wrap(~Diversity, scales = "free_y") +
  scale_color_manual(values = c("Degraded" = "#d7191c", "Created" = "#2c7bb6", 
                                "Restored" = "#fdae61", "Natural" = "#1a9641")) +
  scale_fill_manual(values = c("Degraded" = "#d7191c", "Created" = "#2c7bb6", 
                               "Restored" = "#fdae61", "Natural" = "#1a9641")) +
  theme_bw() +
  labs(title = "Diversity boost across salinity gradient", x = "Mean salinity (ppt)", y = "Hill number estimate")


# --- Inferential Statistics: GLMMs ---
# Filter: exclude natural habitat from modeling due to lack of replication
modeling_df <- asy_est %>% filter(Condition != "Natural")

# q0: Species richness (Negative Binomial)
m0 <- glmmTMB(Estimator ~ Condition + Season + Mean_Salinity + (1|Location), 
              data = modeling_df %>% filter(Diversity == "Species richness"),
              family = nbinom2)
m0

# q1: Shannon diversity (Gamma with Log Link)
m1 <- glmmTMB(Estimator ~ Condition + Season + Mean_Salinity + (1|Location), 
              data = modeling_df %>% filter(Diversity == "Shannon diversity"),
              family = Gamma(link = "log"))
m1

# q2: Simpson diversity (Gamma with Log Link)
m2 <- glmmTMB(Estimator ~ Condition + Season + Mean_Salinity + (1|Location), 
              data = modeling_df %>% filter(Diversity == "Simpson diversity"),
              family = Gamma(link = "log"))
m2

# Function to extract GLMM results for your text
get_glmm_stats <- function(model, model_name) {
  res <- car::Anova(model, type = "III") # Type III for models with interactions/covariates
  print(paste("--- Stats for", model_name, "---"))
  print(res)
}

# Apply
get_glmm_stats(m0, "Species Richness (q=0)")
get_glmm_stats(m1, "Shannon Diversity (q=1)")
get_glmm_stats(m2, "Simpson Diversity (q=2)")


# DHARMa 
# For q0 
res_m0 <- simulateResiduals(m0)
plot(res_m0) 

# For q1 
res_m1 <- simulateResiduals(m1)
plot(res_m1)

# For q2
res_m2 <- simulateResiduals(m2)
plot(res_m2)

#  Post-hoc 
# Run these only if the models show significant effects
emmeans(m0, pairwise ~ Condition, type = "response")
emmeans(m1, pairwise ~ Condition, type = "response")
emmeans(m2, pairwise ~ Condition, type = "response")



# --- Multivariate community composition (vegan) ---
spec_matrix <- pooled_counts %>% 
  select(all_of(species_cols))

meta_final <- pooled_counts %>% 
  left_join(salinity_clean, by = c("Location", "Season"))

# Transformations & dissimilarity
spec_hel <- decostand(spec_matrix, method = "hellinger")
dist_bray <- vegdist(spec_hel, method = "bray")

# PERMANOVA & PERMDISP
perm_res <- adonis2(dist_bray ~ Condition + Season + Mean_Salinity, data = meta_final)

# Extracting from your PERMANOVA
print(perm_res) 

disp_res <- betadisper(dist_bray, meta_final$Condition)
anova(disp_res) # Check for homogeneity of dispersion


# Running NMDS (k=2 dimensions)
nmds_res <- metaMDS(dist_bray, k = 2, trymax = 100)

# Extract site scores for ggplot
nmds_scores <- as.data.frame(scores(nmds_res, display = "sites"))
nmds_scores <- cbind(nmds_scores, meta_final) # Merge with your metadata

# Plot per condition
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Condition, shape = Season)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(fill = Condition), geom = "polygon", alpha = 0.1, level = 0.95) +
  scale_color_manual(values = c("Degraded" = "#d7191c", "Created" = "#2c7bb6", 
                                "Restored" = "#fdae61", "Natural" = "#1a9641")) +
  scale_fill_manual(values = c("Degraded" = "#d7191c", "Created" = "#2c7bb6", 
                               "Restored" = "#fdae61", "Natural" = "#1a9641")) +
  theme_bw() +
  labs(title = "NMDS of fish assemblage composition",
       subtitle = paste("Stress =", round(nmds_res$stress, 3)),
       caption = "Bray-Curtis distance on Hellinger-transformed abundance") +
  theme(panel.grid = element_blank())

# Plot per locaion too
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Condition, shape = Season)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(aes(fill = Condition), geom = "polygon", alpha = 0.1, level = 0.95) +
  geom_text_repel(aes(label = Location), size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("Degraded" = "#d7191c", "Created" = "#2c7bb6", 
                                "Restored" = "#fdae61", "Natural" = "#1a9641")) +
  scale_fill_manual(values = c("Degraded" = "#d7191c", "Created" = "#2c7bb6", 
                               "Restored" = "#fdae61", "Natural" = "#1a9641")) +
  theme_bw() +
  labs(title = "NMDS: Fish assemblage composition",
       subtitle = paste("Stress =", round(nmds_res$stress, 3)),
       caption = "Bray-Curtis distance on Hellinger-transformed abundance")
# Sites closer together are more similar

# (Jaccard / Presence-Absence) 
spec_pa <- decostand(spec_matrix, method = "pa")

dist_jaccard <- vegdist(spec_pa, method = "jaccard", binary = TRUE)

perm_jaccard <- adonis2(dist_jaccard ~ Condition + Season + Mean_Salinity, data = meta_final)
print(perm_jaccard)

nmds_jaccard <- metaMDS(dist_jaccard, k = 2, trymax = 100)

jaccard_scores <- as.data.frame(scores(nmds_jaccard, display = "sites"))
jaccard_scores <- cbind(jaccard_scores, meta_final)

ggplot(jaccard_scores, aes(x = NMDS1, y = NMDS2, color = Condition, shape = Season)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(fill = Condition), geom = "polygon", alpha = 0.1) +
  geom_text_repel(aes(label = Location), size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("Degraded" = "#d7191c", "Created" = "#2c7bb6", 
                                "Restored" = "#fdae61", "Natural" = "#1a9641")) +
  scale_fill_manual(values = c("Degraded" = "#d7191c", "Created" = "#2c7bb6", 
                               "Restored" = "#fdae61", "Natural" = "#1a9641")) +
  theme_bw() +
  labs(title = "NMDS: Species turnover",
       subtitle = paste("Stress =", round(nmds_jaccard$stress, 3)),
       caption = "Binary Jaccard distance on presence-absence data")
# Sites closer have consisent set of species 


# --- ISA (IndVal) ---
# Running ISA based on Location to identify site-specific indicators
inv <- multipatt(spec_hel, meta_final$Location, func = "IndVal.g", control = how(nperm = 999))
summary(inv)

# Extract species with p < 0.05 to use for filtering the heatmap
sig_species <- inv$sign %>%
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  #filter(p.value < 0.05) %>%
  pull(Species)

meta_only <- meta_final %>% 
  select(Location, Season, Condition, Mean_Salinity)

# Combine the clean metadata with the Hellinger matrix
heatmap_data <- cbind(meta_only, spec_hel)

heatmap_avg <- heatmap_data %>%
  group_by(Location, Season) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = -c(Location, Season, Mean_Salinity), 
               names_to = "Species", 
               values_to = "Abundance") %>%
  #filter(Species %in% sig_species) %>%
  mutate(Location = factor(Location, levels = location_levels),
         Species = reorder(Species, Abundance, sum))

# Plot by RA
ggplot(heatmap_avg, aes(x = Location, y = Species, fill = Abundance)) +
  geom_tile(color = "white") +
  facet_wrap(~Season) +  
  scale_fill_gradient(low = "white", high = "darkred", na.value = "white") +
  theme_minimal() +
  labs(title = "Indicator species ",
       subtitle = "Mean Hellinger-transformed abundance",
       fill = "Relative\nabundance",
       x = "Location") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, face = "italic"), 
        strip.background = element_rect(fill = "gray95"),     
        panel.grid = element_blank())


# Salinity gradient 
heatmap_avg2 <- heatmap_avg %>%
  mutate(Location = factor(Location, levels = location_levels),
         Loc_Numeric = as.numeric(Location))

species_geo_order <- heatmap_avg2 %>%
  group_by(Species) %>%
  summarise(
    Centroid = weighted.mean(Loc_Numeric, w = Abundance + 0.0001)) %>%
  arrange(Centroid) %>% 
  pull(Species)

heatmap_avg2$Species <- factor(heatmap_avg2$Species, levels = species_geo_order)

# Plot by gradient
ggplot(heatmap_avg2, aes(x = Location, y = Species, fill = Abundance)) +
  geom_tile(color = "white", linewidth = 0.1) +
  facet_wrap(~Season) +  
  scale_fill_gradient(low = "#f7f7f7", high = "#b2182b", na.value = "white") + 
  theme_minimal() +
  labs(title = "Indicator species ",
       subtitle = "Mean Hellinger-transformed abundance",
       fill = "Relative\nabundance",
       x = "Location") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 8, face = "italic"), 
        strip.background = element_rect(fill = "gray95"),     
        panel.grid = element_blank())

# Condition
inv_cond <- multipatt(spec_hel, meta_final$Condition, func = "IndVal.g")



# --- Global and Partial dbRDA ---

# Global dbRDA (The combined influence of all variables)
final_dbRDA <- dbrda(spec_hel ~ Condition + Season + Mean_Salinity, data = meta_final)

plot(final_dbRDA, type = "n", main = "dbRDA: Global model")
text(final_dbRDA, display = "sites", labels = meta_final$Location, 
     col = "black", cex = 0.7)
text(final_dbRDA, display = "bp", col = "blue", cex = 0.9, font = 2)

# Partial dbRDA (habitat effect, controlling for Salinity)
# This addresses RQ4
dbRDA_partial <- dbrda(spec_hel ~ Condition + Season + Condition(Mean_Salinity), 
                       data = meta_final)
plot(dbRDA_partial, type = "n", main = "dbRDA: Partial model")
text(dbRDA_partial, display = "sites", labels = meta_final$Location, 
     col = "black", cex = 0.7)
text(dbRDA_partial, display = "bp", col = "blue", cex = 0.9, font = 2)

# Significance test for Habitat/Season after removing Salinity noise
anova(dbRDA_partial, by = "terms") 


# --- Variance partitioning ---
var_part <- varpart(spec_hel, ~ Condition, ~ Mean_Salinity, data = meta_final)

# Variance Partitioning
summary(var_part)

# Accessing the exact percentages for your Venn Diagram:
cat("Unique Habitat:", var_part$part$indfract$Adj.R.square[1] * 100, "%\n")
cat("Unique Salinity:", var_part$part$indfract$Adj.R.square[2] * 100, "%\n")
cat("Shared:", var_part$part$indfract$Adj.R.square[3] * 100, "%\n")

# Extract individual fractions correctly:
# [1] = Unique Habitat, [2] = Unique Salinity, [3] = Shared, [4] = Residuals
vals <- var_part$part$indfract

fit <- euler(c(
  "Habitat" = pmax(0, vals$Adj.R.square[1]),
  "Salinity" = pmax(0, vals$Adj.R.square[2]),
  "Habitat&Salinity" = pmax(0, vals$Adj.R.square[3])))

plot(fit, 
     quantities = list(
       type = "percent", 
       font = 2, 
       cex = 0.8,
       labels = paste0(round(vals$Adj.R.square[1:3] * 100, 1), "%")),
     fills = list(fill = c("#fdae61", "#2c7bb6"), alpha = 0.6), # Orange and Blue
     edges = list(lty = 1),
     labels = list(font = 2, cex = 1),
     main = list(label = "Variance partitioning: Habitat vs Salinity", cex = 1.2))



# --- Summary table ---

# Alpha diversity (GLMMs)
alpha_summary <- rbind(
  as.data.frame(car::Anova(m0, type = "III")) %>% mutate(Model = "Richness (q0)"),
  as.data.frame(car::Anova(m1, type = "III")) %>% mutate(Model = "Shannon (q1)"),
  as.data.frame(car::Anova(m2, type = "III")) %>% mutate(Model = "Simpson (q2)")) %>% 
  rownames_to_column("Term")

# Composition (PERMANOVA)
perm_res <- adonis2(dist_bray ~ Condition + Season + Mean_Salinity, 
                    data = meta_final, by = "terms")

perm_jaccard <- adonis2(dist_jaccard ~ Condition + Season + Mean_Salinity, 
                        data = meta_final, by = "terms")

comp_summary <- data.frame(
  Analysis = c("Structure (Bray)", "Structure (Bray)", "Turnover (Jaccard)", "Turnover (Jaccard)"),
  Term = c("Condition", "Season", "Condition", "Season"),
  F_Stat = c(perm_res["Condition", "F"], perm_res["Season", "F"], 
             perm_jaccard["Condition", "F"], perm_jaccard["Season", "F"]),
  P_Value = c(perm_res["Condition", "Pr(>F)"], perm_res["Season", "Pr(>F)"], 
              perm_jaccard["Condition", "Pr(>F)"], perm_jaccard["Season", "Pr(>F)"]))

print(comp_summary)

# Environmental
env_summary <- data.frame(
  Fraction = c("Unique Habitat", "Unique Salinity", "Shared"),
  Adj_R2_Percent = round(vals$Adj.R.square[1:3] * 100, 2))

# All
list(Alpha_Stats = alpha_summary, Composition_Stats = comp_summary, Environment_Stats = env_summary)
