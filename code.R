# ------ Setup and Preparation ------

# install required packages
#install.packages(openxlsx)
#install.packages(readxl)
#install.packages(writexl)
#install.packages(coda)
#install.packages(gemtc)
#install.packages(meta)
#install.packages(netmeta)
#install.packages(ggplot2)
#install.packages(rjags)

# Load required packages
library(openxlsx)
library(readxl)
library(writexl)
library(coda)
library(gemtc)
library(meta)
library(netmeta)
library(ggplot2)
library(rjags)

# Set working directory
setwd("D:/Desktop/Breast_Cancer_Anthracyclines/4. Data_Extraction/Input_Data/Breast_Cancer")

# Load dataset
data <- readxl::read_excel("breast_cancer.xlsx", sheet = "sheet1" )

# ------ Treatment Coding Standardization ------
# Define treatment codes and descriptions
treatments <- data.frame(
  id = 1:9,
  description = c("DC", "EC", "ET", "E-T", "EC+Bev", "DBev", 
                  "EC-Tra", "EC+Lap", "DV")
)

# ------ Network Construction and Visualization ------
# Create network with standardized descriptions
network <- mtc.network(data.ab = data, 
                       description = "BC_network",
                       treatments = treatments)

# Network summary
summary(network)

# Generate network plot (High resolution for publication)
tiff("BC_network.tiff", height = 5000, width = 5000, res = 600)
plot(network, 
     use.description = TRUE, 
     vertex.label.cex = 1.5, 
     vertex.size = 12, 
     vertex.shape = "circle", 
     vertex.label.color = "darkblue", 
     vertex.label.dist = 2, 
     vertex.label.degree = -pi/3, 
     vertex.color = "blue", 
     dynamic.edge.width = TRUE,
     edge.color = "black", 
     vertex.label.font = 2)
dev.off()

# ------ Consistency Model (Random Effects) ------
# Option A: With treatment standardization
model.ran <- mtc.model(network, 
                       type = "consistency",
                       likelihood = "binom", 
                       link = "logit",
                       linearModel = "random", 
                       n.chain = 3,  # 3-4 chains recommended for convergence assessment
                       dic = TRUE)


# Run MCMC simulation
result.ran <- mtc.run(model.ran, 
                      n.adapt = 10000, 
                      n.iter = 100000, 
                      thin = 10)

# ------ Model Diagnostics ------
# Convergence assessment
gelman.diag(result.ran)

# 1. Convergence plot
tiff("convergence_plot.tiff", height = 4000, width = 4000, res = 600)
gelman.plot(result.ran, use.description = TRUE, digits = 2)
dev.off()

# 2. Trace and density plots
tiff("trace_density_plot.tiff", height = 4000, width = 4000, res = 600)
plot(result.ran, use.description = TRUE, digits = 2)
dev.off()

# ------ Forest Plot ------
tiff("forest_plot.tiff", height = 5000, width = 4000, res = 600)
forest_plot <- gemtc::forest(relative.effect(result.ran, t1 = "TAC_H"), 
                             use.description = TRUE)
plot(forest_plot)
dev.off()

# ------ Inconsistency and Heterogeneity Assessment ------
# 1. Node-splitting analysis
mtc.nodesplit.comparisons(network)  # Identify closed loops
result.ns <- mtc.nodesplit(network)
summary.ns <- summary(result.ns)
plot(summary.ns, digits = 3)

# 2. Unrelated mean effects (UME) model for inconsistency
model_ume <- mtc.model(network, 
                       type = "ume", 
                       n.chain = 4, 
                       likelihood = "binom",
                       link = "logit", 
                       linearModel = "random")
result_ume <- mtc.run(model_ume,
                      n.adapt = 10000, 
                      n.iter = 200000, 
                      thin = 10)
summary(result_ume)  # Compare DIC and I² statistics

# 3. ANOVA-type heterogeneity test
result.anohe <- mtc.anohe(network, 
                          n.adapt = 10000, 
                          n.iter = 100000, 
                          thin = 1)
plot(summary(result.anohe))

# ------ Treatment Ranking ------
# League table (odds ratios)
league_table <- as.data.frame(round(exp(relative.effect.table(result.ran)), 2))
write.csv(league_table, "league_table.csv")

# Ranking probabilities
ranks <- rank.probability(result.ran, preferredDirection = -1)  # -1 for harmful outcomes
write.csv(ranks, "ranking_probabilities.csv")

# SUCRA values
sucra_values <- sucra(ranks)
write.csv(sucra_values, "SUCRA.csv")

# Rank quantiles
rank_quantiles <- rank.quantiles(ranks)
write.csv(rank_quantiles, "rank_quantiles.csv")

# ------ Ranking Visualization ------
# Cumulative ranking plot
tiff("cumulative_ranking.tiff", height = 5000, width = 9000, res = 600)
plot(ranks,
     xlab = 'Treatment', 
     ylab = 'Cumulative Probability',
     col = c("#db6968", "#4d97cd", "#99cbeb", "#459943", 
             "#fdc58f", "#e8c559", "#a3d393", "#f8984e"))
dev.off()

# Non-cumulative ranking plot (bar chart)
tiff("noncumulative_ranking.tiff", height = 5000, width = 9000, res = 600)
plot(ranks, 
     xlab = 'Treatment', 
     ylim = c(0, 0.3),
     beside = TRUE,
     col = c("#db6968", "#4d97cd", "#99cbeb", "#459943", 
             "#fdc58f", "#e8c559", "#a3d393", "#f8984e"))
dev.off()

# ------ Custom Ranking Line Plot ------
# Prepare ranking data for ggplot
ranking_data <- read.csv("ranking_data.csv")  # Pre-formatted ranking data
# Expected format: columns = rank, prob, treatment

# Create line plot
ranking_plot <- ggplot(ranking_data, aes(x = rank, y = prob, colour = treatment)) +
  geom_line(linewidth = 2) +
  geom_point(size = 5) +
  labs(x = "Rank", y = "Probability") +
  scale_color_manual(values = c("#db6968", "#4d97cd", "#99cbeb", "#459943",
                                "#fdc58f", "#e8c559", "#a3d393", "#f8984e")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 28),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 35),
    legend.text = element_text(size = 20)
  )

# Export high-resolution plot
tiff("ranking_line_plot.tiff", height = 5000, width = 7000, res = 600)
print(ranking_plot)
dev.off()

