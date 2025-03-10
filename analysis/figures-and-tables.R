### Script Header ===============================================================
# Topic: Ichthyoplankton Phenology Manuscript Figures
# Author: Kathryn Chen
# Contact: ksc005@ucsd.edu
# Date: November 2024

# Objective: To consolidate the code for producing figures for the manuscript. First run the sections "Set Up" and "Set Common Figure Theme," then each of the following sections can be run independently.
### ============================================================================

# Set Up =======================================================================
## location declaration
here::i_am("analysis/figures-and-tables.R")

## load libraries needed to make figures
library(plyr) # for data manipulation
library(tidyverse) # for data manipulation and figures
library(reshape2) # for data manipulation
library(ggpubr) # for multi-panel figures
library(ggrepel) # for text repel in figures
library(AICcmodavg) # for environmental variable model selection
library(stargazer) # for data tables
library(gt) # for data tables
library(gtsummary) # for data tables
library(readxl) # for reading in Excel sheets
library(mgcv) # for GAMs
library(tidymv) # for GAM plotting
library(lmtest) # for testing autocorrelation in linear models
library(here)

# Set Common Figure Theme ======================================================
## Make adjustments to figure style here, including text sizes, linewidths, and fonts.
my.ggplot.theme <- list(
  theme_bw(),
  theme(axis.text = element_text(size = 13, family = "Helvetica"),
        axis.title = element_text(size = 15, family = "Helvetica"),
        axis.ticks = element_line(linewidth = 1.2),

        legend.position = "none",

        panel.border = element_rect(linewidth = 1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        strip.background = element_rect(linewidth = 1.5),
        strip.text = element_text(size = 14, family = "Helvetica"),

        plot.title = element_text(size = 16, family = "Helvetica"))
)

# CTa x Time ===================================================================
## Figure description: One multi-panel figure, with eight panels showing ichthyoplankton central tendency anomaly against time for the entire assemblage, species shifting earlier, species with no long-term linear change, and species shifting later, for CalCOFI (left column) and NH Line (right column).

## Read in data needed to make this figure
cta.time.calcofi <- read.csv(here("data", "fig1.calcofi.csv"))
cta.time.nhl <- read.csv(here("data", "fig1.nhl.csv"))

## Make the four CalCOFI panels
cta.time.1 <- cta.time.calcofi %>%
  filter(group == "assemblage") %>%
  ggplot(aes(x = as.numeric(as.character(Decade)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-30, 30),
                     breaks = c(-15, 0, 15)) +
  scale_x_continuous(limits = c(1950, 2020),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  labs(title = "CalCOFI Assemblage",
       x = "Decade",
       y = "CTa (days)")
cta.time.2 <- cta.time.calcofi %>%
  filter(group == "earlier") %>%
  ggplot(aes(x = as.numeric(as.character(Decade)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-90, 90),
                     breaks = c(-60, 0, 60)) +
  scale_x_continuous(limits = c(1948, 2022),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  labs(title = "",
       x = "Decade",
       y = "CTa (days)")
cta.time.3 <- cta.time.calcofi %>%
  filter(group == "no_change") %>%
  ggplot(aes(x = as.numeric(as.character(Decade)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-30, 30),
                     breaks = c(-15, 0, 15)) +
  scale_x_continuous(limits = c(1950, 2020),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  labs(title = "",
       x = "Decade",
       y = "CTa (days)")
cta.time.4 <- cta.time.calcofi %>%
  filter(group == "later") %>%
  ggplot(aes(x = as.numeric(as.character(Decade)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-60, 60),
                     breaks = c(-30, 0, 30)) +
  scale_x_continuous(limits = c(1950, 2020),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  labs(title = "",
       x = "Decade",
       y = "CTa (days)")

## Make the four NH Line panels
cta.time.5 <- cta.time.nhl %>%
  filter(group == "assemblage") %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-30, 30),
                     breaks = c(-15, 0, 15)) +
  scale_x_continuous(limits = c(.5, 9.5),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "NH Line Assemblage",
       x = "3-y Period",
       y = "CTa (days)")
cta.time.6 <- cta.time.nhl %>%
  filter(group == "earlier") %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-90, 90),
                     breaks = c(-60, 0, 60)) +
  scale_x_continuous(limits = c(.5, 9.5),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "",
       x = "3-y Period",
       y = "CTa (days)")
cta.time.7 <- cta.time.nhl %>%
  filter(group == "no_change") %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-30, 30),
                     breaks = c(-15, 0, 15)) +
  scale_x_continuous(limits = c(.5, 9.5),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "",
       x = "3-y Period",
       y = "CTa (days)")
cta.time.8 <- cta.time.nhl %>%
  filter(group == "later") %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-60, 60),
                     breaks = c(-30, 0, 30)) +
  scale_x_continuous(limits = c(.5, 9.5),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "",
       x = "3-y Period",
       y = "CTa (days)")

## Put the panels together into one figure with panel labels
cta.time <- ggarrange(cta.time.1, cta.time.5, cta.time.2, cta.time.6,
                      cta.time.3, cta.time.7, cta.time.4, cta.time.8,
                      labels = c("(a)", "(b)", "(c)", "(d)",
                                 "(e)", "(f)", "(g)", "(h)"),
                      font.label = list(size = 20, face = "bold"),
                      # hjust = -0.8,
                      vjust = 1.25,
                      ncol = 2, nrow = 4, align = "hv")

## Save figure image with set size
ggsave(filename = "cta.time.png",
       plot = cta.time,
       path = here("figures"), height = 250, width = 180, units = "mm", dpi = 600)

# Boxplots =====================================================================
## Figure description: One multi-panel figure showing side-by-side group boxplots for orders, adult habitats, cross-shore distributions, biogeographic affinities for CalCOFI trend in phenology (top) and NH Line trend in phenology (bottom).

# Read in data needed for this section
design.spp.calcofi <- read.csv("design.spp.calcofi.csv")
design.spp.nhl <- read.csv("design.spp.nhl.csv")

# set colors. these colors are taken from the grafify colorblind-friendly palette "kelly", which can be viewed with library(grafify) and plot_grafify_palette(palette = "kelly")
colorPalette <- data.frame(trait = c("Argentiniformes", "Stomiiformes", "Myctophiformes",
                                     "Pleuronectiformes", "Scombriformes", "P.Scorpaenoidei",
                                     "P.Cottoidei", "Demersal", "Mesopelagic", "Epipelagic",
                                     "Coastal", "Oceanic", "Coastal-Oceanic",
                                     "Wide distribution", "Cool-water", "Warm-water"),
                           color = c("#F3C300", "#F38400", "#848482",
                                     "#E68FAC", "#F99379", "#B3446C",
                                     "#E25822", "#604E97", "#875692", "#F6A600",
                                     "#C2B280", "#0067A5", "#A1CAF1",
                                     "#DCD300", "#8DB600", "#008856"))

# make calcofi boxplots
boxplots.calcofi <- design.spp.calcofi %>%
  select(!c(X, r)) %>%
  melt(id.vars = "trend_in_phenology") %>%
  na.omit()

boxplots.calcofi$value[boxplots.calcofi$value == "Perciformes/Scorpaenoidei"] <-
  "P.Scorpaenoidei" # Change `Perciformes/Suborder` / `Perciformes.Sub` labels to `P.Suborder`

boxplots.calcofi <- boxplots.calcofi %>%
  filter(!(value %in% c("Unfished", "Fished")))

boxplots.calcofi$value = factor(boxplots.calcofi$value, levels = unique(boxplots.calcofi$value))

boxplots.calcofi$color = colorPalette$color[match(boxplots.calcofi$value,
                                                  colorPalette$trait)]

boxplots.calcofi.means <- boxplots.calcofi %>%
  group_by(value) %>%
  reframe(meanTrend = mean(trend_in_phenology))

calcofi.assemblage.mean = mean(design.spp.calcofi$trend_in_phenology)

p1 <-
  ggplot(boxplots.calcofi,
         aes(y = trend_in_phenology,
             x = value)) +
  geom_boxplot(aes(fill = color)) +
  scale_fill_identity() +
  geom_point() +
  geom_point(data = boxplots.calcofi.means,
             aes(x = value, y = meanTrend),
             size = 4, fill = "red", shape = 23) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs(y = expression(paste("Trend in Phenology (d ", y^-1, ")")))

# make nh line boxplots
boxplots.nhl <- design.spp.nhl %>%
  select(!c(X, cor_NHL)) %>%
  melt(id.vars = "trend_NHL_d3y") %>%
  na.omit()

boxplots.nhl$value[boxplots.nhl$value == "Perciformes/Cottoidei"] <-
  "P.Cottoidei" # Change `Perciformes/Suborder` / `Perciformes.Sub` labels to `P.Suborder`

# there's only one epipelagic, one CO, two widely distributed spp, so we'll take those out
boxplots.nhl <- boxplots.nhl %>%
  filter(!(value %in% c("Epipelagic", "Coastal-Oceanic", "Wide distribution",
                        "Unfished", "Fished")))

boxplots.nhl$value = factor(boxplots.nhl$value, levels = c("P.Cottoidei",
                                                           "Myctophiformes",
                                                           "Pleuronectiformes",
                                                           "Demersal", "Mesopelagic",
                                                           "Coastal", "Oceanic", "Cool-water"))

boxplots.nhl$color = colorPalette$color[match(boxplots.nhl$value,
                                              colorPalette$trait)]

boxplots.nhl.means <- boxplots.nhl %>%
  group_by(value) %>%
  reframe(meanTrend = mean(trend_NHL_d3y))

nhl.assemblage.mean = mean(design.spp.nhl$trend_NHL_d3y)

p2 <-
  ggplot(boxplots.nhl,
         aes(y = trend_NHL_d3y,
             x = value)) +
  geom_boxplot(aes(fill = color)) +
  scale_fill_identity() +
  geom_point() +
  geom_point(data = boxplots.nhl.means,
             aes(x = value, y = meanTrend),
             size = 4, fill = "red", shape = 23) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  scale_y_continuous(limits = c(-6, 6),
                     breaks = seq(from = -6, to = 6, by = 2)) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) +
  labs(y = expression(paste("Trend in Phenology (d ", y^-1, ")")))

# Put the boxplots together in one figure with panels labeled
boxplots <- ggarrange(p1, p2,
                      labels = c("(a)", "(b)"),
                      font.label = list(size = 20),
                      # hjust = -1,
                      vjust = 1.25,
                      ncol = 1, nrow = 2, align = "hv")

## Save the trend figure with set size
ggsave(filename = "trendBoxplots.png",
       plot = boxplots,
       path = here("figures"), height = 12, width = 11, units = "in")


# PCA Plots ====================================================================
## Figure description: One multi-panel figure showing the results of PCA on the CalCOFI assemblage (left column) and NH Line assemblage (right column). Panels will show species ordination and eigenvectors of ecological traits.

## Read in data needed to make this figure
pca.markers.calcofi <- read.csv(here("data", "calcofi-pca-markers.csv"))
pca.cors.calcofi <- read.csv(here("data", "calcofi-pca-plots.csv"), row.names = 1)

pca.markers.nhl <- read.csv(here("data", "nhl-pca-markers.csv"))
pca.cors.nhl <- read.csv(here("data", "nhl-pca-plots.csv"), row.names = 1)

## Change `Perciformes/Suborder` labels to `P.Suborder`
row.names(pca.cors.calcofi)[13] <- "P.Scorpaenoidei"
row.names(pca.cors.nhl)[11] <- "P.Cottoidei"

# plot species ordination with markers proportional to magnitude of trend in phenology, and colored by direction of trend (black = negative, white = positive)
pca.plot1.calcofi <-
  ggplot(data = pca.markers.calcofi,
         aes(x = Comp.1, y = Comp.2)) +
  geom_hline(yintercept = 0, colour = "gray65", linetype = 2) +
  geom_vline(xintercept = 0, colour = "gray65", linetype = 2) +
  geom_point(shape = 21, alpha = 0.8,
             aes(size = trend,
                 fill = direction)) +
  scale_fill_manual(values = c("black", "white")) +
  my.ggplot.theme +
  lims(x = c(-6.5, 6.5),
       y = c(-6.5, 6.5)) +
  labs(x = "PC1", y = "PC2", title = "CalCOFI Assemblage")

pca.plot1.nhl <-
  ggplot(data = pca.markers.nhl,
         aes(x = Comp.1, y = Comp.2)) +
  geom_hline(yintercept = 0, colour = "gray65", linetype = 2) +
  geom_vline(xintercept = 0, colour = "gray65", linetype = 2) +
  geom_point(shape = 21,alpha = 0.8,
             aes(size = trend,
                 fill = direction)) +
  scale_fill_manual(values = c("black", "white")) +
  my.ggplot.theme +
  lims(x = c(-8.0, 8.0),
       y = c(-8.0, 8.0)) +
  labs(x = "PC1", y = "PC2", title = "NH Line Assemblage")

# function to draw circles
circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

## taxonomic orders, data frame with arrows coordinates
orders.calcofi <- pca.cors.calcofi[c(11:16),]

arrows.orders.calcofi = data.frame(x1 = rep(0,6),
                                   y1 = rep(0,6),
                                   x2 = orders.calcofi$Comp.1,
                                   y2 = orders.calcofi$Comp.2)

orders.nhl <- pca.cors.nhl[c(10:12),]

arrows.orders.nhl = data.frame(x1 = rep(0,3),
                               y1 = rep(0,3),
                               x2 = orders.nhl$Comp.1,
                               y2 = orders.nhl$Comp.2)

pca.plot2.calcofi <-
  ggplot() +
  geom_path(data = corcir,
            aes(x = x, y = y),
            colour = "gray65") +
  geom_segment(data = arrows.orders.calcofi,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               colour = "gray65") +
  geom_hline(yintercept = 0,
             colour = "gray65",
             linetype = 2) +
  geom_vline(xintercept = 0,
             colour = "gray65",
             linetype = 2) +
  geom_text_repel(data = orders.calcofi,
                  aes(x = Comp.1, y = Comp.2, label = rownames(orders.calcofi)),
                  direction = "both",
                  box.padding = 0.3,
                  label.padding = 0.3,
                  force = 30) +
  my.ggplot.theme +
  lims(x = c(-1.1, 1.1),
       y = c(-1.1, 1.1)) +
  labs(x = "PC1", y = "PC2")

pca.plot2.nhl <-
  ggplot() +
  geom_path(data = corcir,
            aes(x = x, y = y),
            colour = "gray65") +
  geom_segment(data = arrows.orders.nhl,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               colour = "gray65") +
  geom_hline(yintercept = 0,
             colour = "gray65",
             linetype = 2) +
  geom_vline(xintercept = 0,
             colour = "gray65",
             linetype = 2) +
  geom_text_repel(data = orders.nhl,
                  aes(x = Comp.1, y = Comp.2, label = rownames(orders.nhl)),
                  direction = "both",
                  box.padding = 0.05,
                  force = 0.3) +
  my.ggplot.theme +
  lims(x = c(-1.1, 1.1),
       y = c(-1.1, 1.1)) +
  labs(x = "PC1", y = "PC2")

## traits: adult habitat, cross-shore distribution, biogeographic affinity, monthMax
traits1.calcofi <- pca.cors.calcofi[c(1:10),]

arrows.traits1.calcofi = data.frame(x1 = rep(0,10),
                                    y1 = rep(0,10),
                                    x2 = traits1.calcofi$Comp.1,
                                    y2 = traits1.calcofi$Comp.2)

traits1.nhl <- pca.cors.nhl[c(1:9),]

arrows.traits1.nhl = data.frame(x1 = rep(0,9),
                                y1 = rep(0,9),
                                x2 = traits1.nhl$Comp.1,
                                y2 = traits1.nhl$Comp.2)

rownames(traits1.calcofi)[1] <- "MonthMax"
rownames(traits1.nhl)[1] <- "MonthMax"

pca.plot3.calcofi <-
  ggplot() +
  geom_path(data = corcir,
            aes(x = x, y = y),
            colour = "gray65") +
  geom_segment(data = arrows.traits1.calcofi,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               colour = "gray65") +
  geom_hline(yintercept = 0,
             colour = "gray65",
             linetype = 2) +
  geom_vline(xintercept = 0,
             colour = "gray65",
             linetype = 2) +
  geom_text_repel(data = traits1.calcofi,
                  aes(x = Comp.1, y = Comp.2, label = rownames(traits1.calcofi)),
                  direction = "both",
                  box.padding = 0.05,
                  force = 0.3) +
  my.ggplot.theme +
  lims(x = c(-1.1, 1.1),
       y = c(-1.1, 1.1)) +
  labs(x = "PC1", y = "PC2")

pca.plot3.nhl <-
  ggplot() +
  geom_path(data = corcir,
            aes(x = x, y = y),
            colour = "gray65") +
  geom_segment(data = arrows.traits1.nhl,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               colour = "gray65") +
  geom_hline(yintercept = 0,
             colour = "gray65",
             linetype = 2) +
  geom_vline(xintercept = 0,
             colour = "gray65",
             linetype = 2) +
  geom_text_repel(data = traits1.nhl,
                  aes(x = Comp.1, y = Comp.2, label = rownames(traits1.nhl)),
                  direction = "both",
                  box.padding = 0.05,
                  force = 0.3) +
  my.ggplot.theme +
  lims(x = c(-1.1, 1.1),
       y = c(-1.1, 1.1)) +
  labs(x = "PC1", y = "PC2")


pca.plots <- ggarrange(pca.plot1.calcofi, pca.plot1.nhl,
                       pca.plot2.calcofi, pca.plot2.nhl,
                       pca.plot3.calcofi, pca.plot3.nhl,
                       labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
                       font.label = list(size = 20),
                       # hjust = -0.9,
                       vjust = 1.25,
                       ncol = 2, nrow = 3, align = "hv")

ggsave(filename = "pca.plots.png",
       plot = pca.plots,
       path = here("figures"), height = 250, width = 180, units = "mm", dpi = 600)


# PCA Clusters CTa x Time =====================================================
## Figure description: One multi-panel figure showing mean CTa against time for the clusters of species from PCA on the CalCOFI assemblage (left column) and NH Line assemblage (right column).

# read in data needed to make this figure
calcofi.clusters <- read.csv(here("data", "calcofi-cluster-plots.csv"))
nhl.clusters <- read.csv(here("data", "nhl-cluster-plots.csv"))

# attach labels to the decades/periods
decades <- data.frame(decade.number = 1:7,
                      decade.year = seq(from = 1955, to = 2015, by = 10))

calcofi.clusters$Decade.Year <- decades$decade.year[match(calcofi.clusters$Decade,
                                                              decades$decade.number)]

periods <- data.frame(period.number = 1:9,
                      period.year = seq(from = 1998.5, to = 2022.5, by = 3))

nhl.clusters$Period_3y.Year <- periods$period.year[match(nhl.clusters$Period_3y,
                                                         periods$period.number)]


# make the calcofi subplots
calcofi.e <- calcofi.clusters %>%
  filter(cluster == "Epipelagic") %>%
  ggplot(aes(x = as.numeric(as.character(Decade.Year)),
             y = mean_cta_days)) +
  geom_point(size = 2) +
  geom_path(linewidth = 1) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_y_continuous(limits = c(-20, 20)) +
  scale_x_continuous(limits = c(1948, 2022),
                     breaks = c(1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020)) +
  theme_bw() +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  labs(x = "Decade",
       y = "CT Anomaly (days)",
       title = "CalCOFI Assemblage")

calcofi.d <- calcofi.clusters %>%
  filter(cluster == "Demersal") %>%
  ggplot(aes(x = as.numeric(as.character(Decade.Year)),
             y = mean_cta_days)) +
  geom_point(size = 2) +
  geom_path(linewidth = 1) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_y_continuous(limits = c(-20, 20)) +
  scale_x_continuous(limits = c(1948, 2022),
                     breaks = c(1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020)) +
  theme_bw() +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  labs(x = "Decade",
       y = "CT Anomaly (days)",
       title = "")

calcofi.m <- calcofi.clusters %>%
  filter(cluster == "Mesopelagic") %>%
  ggplot(aes(x = as.numeric(as.character(Decade.Year)),
             y = mean_cta_days)) +
  geom_point(size = 2) +
  geom_path(linewidth = 1) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_y_continuous(limits = c(-20, 20)) +
  scale_x_continuous(limits = c(1948, 2022),
                     breaks = c(1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020)) +
  theme_bw() +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  labs(x = "Decade",
       y = "CT Anomaly (days)",
       title = "")


# make the nh line subplots
nhl.e <- nhl.clusters %>%
  filter(cluster == "Epipelagic") %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y.Year)),
             y = mean_cta_days)) +
  geom_point(size = 2) +
  geom_path(linewidth = 1) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_y_continuous(limits = c(-75, 75),
                     breaks = c(-80, -40, 0, 40, 80)) +
  scale_x_continuous(limits = c(1996, 2024),
                     breaks = seq(from = 1997, to = 2026, by = 3)) +
  theme_bw() +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  labs(x = "Decade",
       y = "CT Anomaly (days)",
       title = "NH Line Assemblage")

nhl.d <- nhl.clusters %>%
  filter(cluster == "Demersal") %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y.Year)),
             y = mean_cta_days)) +
  geom_point(size = 2) +
  geom_path(linewidth = 1) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_y_continuous(limits = c(-45, 45),
                     breaks = c(-40, 0, 40)) +
  scale_x_continuous(limits = c(1996, 2024),
                     breaks = seq(from = 1997, to = 2026, by = 3)) +
  theme_bw() +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  labs(x = "Decade",
       y = "CT Anomaly (days)",
       title = "")

nhl.m <- nhl.clusters %>%
  filter(cluster == "Mesopelagic") %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y.Year)),
             y = mean_cta_days)) +
  geom_point(size = 2) +
  geom_path(linewidth = 1) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_y_continuous(limits = c(-75, 75),
                     breaks = c(-80, -40, 0, 40, 80)) +
  scale_x_continuous(limits = c(1996, 2024),
                     breaks = seq(from = 1997, to = 2026, by = 3)) +
  theme_bw() +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  labs(x = "Decade",
       y = "CT Anomaly (days)",
       title = "")

# put together
cluster.plots <- ggarrange(calcofi.e, nhl.e,
                           calcofi.m, nhl.m,
                           calcofi.d, nhl.d,
                           labels = c("(a)", "(b)", "(c)",
                                      "(e)", "(f)", "(g)"),
                           font.label = list(size = 20),
                           # hjust = -0.9,
                           vjust = 1.25,
                           ncol = 2, nrow = 3, align = "hv")

ggsave(filename = "cluster.plots.png",
       plot = cluster.plots,
       path = here("figures"), height = 12, width = 9.6, units = "in")

# Environmental Variable Phenology =============================================
## Figure description:

## seven vars: 3 CalCOFI (T, ZDV, BCUI), 4 NHL (T, BCUI, N.COPE, S.COPE)
t.calcofi <- read.csv(here("data", "calcofi-temp-ct.csv"))
zdv.calcofi <- read.csv(here("data", "calcofi-zdv-ct.csv"))
bcui.calcofi <- read.csv(here("data", "calcofi-bcui-ct.csv"))

t.nhl <- read.csv(here("data", "nhl-temp-ct.csv"))
cope.nhl <- read.csv(here("data", "nhl-cope-ct.csv"))
bcui.nhl <- read.csv(here("data", "nhl-bcui-ct.csv"))

env.1 <- ggplot(t.calcofi,
                aes(x = as.numeric(as.character(decade)),
                    y = ct_anomaly_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = ct_anomaly_days - se_cta_days,
                     ymax = ct_anomaly_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_y_continuous(limits = c(-30, 30)) +
  scale_x_continuous(limits = c(1948, 2022),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  labs(title = "CalCOFI Temperature",
       x = "Decade",
       y = "CTa (days)")
env.2 <- ggplot(bcui.calcofi,
                aes(x = as.numeric(as.character(decade)),
                    y = ct_anomaly_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = ct_anomaly_days - cta_se,
                     ymax = ct_anomaly_days + cta_se),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_x_continuous(limits = c(1948, 2022),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  labs(title = "CalCOFI Upwelling",
       x = "Decade",
       y = "CTa (days)")
env.3 <- ggplot(zdv.calcofi,
                aes(x = as.numeric(as.character(decade)),
                    y = ct_anomaly_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = ct_anomaly_days - cta_se,
                     ymax = ct_anomaly_days + cta_se),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_x_continuous(limits = c(1948, 2022),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  labs(title = "CalCOFI Mesozooplankton",
       x = "Decade",
       y = "CTa (days)")
env.4 <- ggplot(t.nhl,
                aes(x = as.numeric(as.character(Period_3y)),
                    y = ct_anomaly_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = ct_anomaly_days - se_cta_days,
                     ymax = ct_anomaly_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_x_continuous(limits = c(0.1, 9.9),
                     breaks = seq(from = 0.5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  labs(title = "NHL Temperature",
       x = "3-y Period",
       y = "CTa (days)")
env.5 <- ggplot(bcui.nhl,
                aes(x = as.numeric(as.character(Period_3y)),
                    y = ct_anomaly_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = ct_anomaly_days - se_cta_days,
                     ymax = ct_anomaly_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_x_continuous(limits = c(0.1, 9.9),
                     breaks = seq(from = 0.5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  labs(title = "NHL Upwelling",
       x = "3-y Period",
       y = "CTa (days)")
env.6 <-
  ggplot(cope.nhl,
         aes(x = as.numeric(as.character(Period_3y)),
             y = ct_anomaly_days,
             group = Assemblage)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2, aes(lty = Assemblage)) +
  geom_linerange(aes(ymin = ct_anomaly_days - se_cta_days,
                     ymax = ct_anomaly_days + se_cta_days,
                     lty = Assemblage),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_x_continuous(limits = c(0.1, 9.9),
                     breaks = seq(from = 0.5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  theme(legend.position = c(.2, .85),
        legend.background = element_blank(),
        legend.text = element_text(size = 13, family = "Helvetica"),
        legend.title = element_blank()) +
  labs(title = "NHL Copepods",
       x = "3-y Period",
       y = "CTa (days)")

env.vars <- ggarrange(env.1, env.4, env.2, env.5, env.3,
                      env.6,
                      labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
                      font.label = list(size = 20),
                      # hjust = -1,
                      vjust = 1.25,
                      ncol = 2, nrow = 3, align = "hv")

ggsave(filename = "env.vars.png",
       plot = env.vars,
       path = here("figures"), height = 250, width = 180, units = "mm", dpi = 600)


# Environmental x Ichthyoplankton Candidate Models =============================
## Figure description:

## Put env vars CTa and ichthyo CTa together
zdv_ct <- read.csv(here("data", "calcofi-zdv-ct.csv"))
zdv_july_ct <- read.csv(here("data", "calcofi-zdv-july-ct.csv"))

temp_ct <- read.csv(here("data", "calcofi-temp-ct.csv"))
temp_july_ct <- read.csv(here("data", "calcofi-temp-july-ct.csv"))

bcui_ct <- read.csv(here("data", "calcofi-bcui-ct.csv"))
bcui_july_ct <- read.csv(here("data", "calcofi-bcui-july-ct.csv"))

# load in our ichthyoplankton species data
spp_cta <- read.csv(here("data", "calcofi-ct-summary.csv"), check.names = F)

# which spp are year-shifted
spp_start_july <- c("Argentina sialis", "Aristostomias scintillans",
                    "Bathylagus pacificus", "Diogenichthys atlanticus",
                    "Electrona risso",
                    "Hippoglossina stomata", "Leuroglossus stilbius",
                    "Lipolagus ochotensis",
                    "Merluccius productus",
                    "Oxylebius pictus", "Pleuronichthys decurrens",
                    "Protomyctophum crockeri", "Pseudobathylagus milleri",
                    "Scorpaenichthys marmoratus", "Sebastes goodei",
                    "Sebastes jordani", "Sebastes levis",
                    "Sebastes paucispinis", "Stomias atriventer",
                    "Tetragonurus cuvieri")

spp_cta <- spp_cta %>%
  mutate(yearstart = ifelse(scientific_name %in% spp_start_july, "July", "January"))

# let us put together a dataframe with the desired variables all together
cand.mods.calcofi <- spp_cta[,c(1,2,7,9,11)]

cand.mods.calcofi <- cand.mods.calcofi %>%
  dplyr::rename(ichthyo_cta = ct_anomaly_days)

cand.mods.calcofi$temp_cta <- ifelse(cand.mods.calcofi$yearstart == "July",
                          temp_july_ct$ct_anomaly_days[match(cand.mods.calcofi$Decade,
                                                            temp_july_ct$New_Decade)],
                          temp_ct$ct_anomaly_days[match(cand.mods.calcofi$Decade,
                                                       temp_ct$decade
                          )])

cand.mods.calcofi$zdv_cta <- ifelse(cand.mods.calcofi$yearstart == "July",
                          zdv_july_ct$ct_anomaly_days[match(cand.mods.calcofi$Decade,
                                                            zdv_july_ct$New_Decade)],
                          zdv_ct$ct_anomaly_days[match(cand.mods.calcofi$Decade,
                                                       zdv_ct$decade
                          )])

cand.mods.calcofi$bcui_cta <- ifelse(cand.mods.calcofi$yearstart == "July",
                           bcui_july_ct$ct_anomaly_days[match(cand.mods.calcofi$Decade,
                                                              bcui_july_ct$New_Decade)],
                           bcui_ct$ct_anomaly_days[match(cand.mods.calcofi$Decade,
                                                         bcui_ct$decade
                           )])

# load in all our environmental variable CT data
temp <- read.csv(here("data", "nhl-temp-ct.csv"))
temp_oct <- read.csv(here("data", "nhl-temp-oct-ct.csv"))
temp_july <- read.csv(here("data", "nhl-temp-july-ct.csv"))

onlypos_bcui <- read.csv(here("data", "nhl-bcui-ct.csv"))
onlypos_bcui_oct <- read.csv(here("data", "nhl-bcui-oct-ct.csv"))
onlypos_bcui_july <- read.csv(here("data", "nhl-bcui-july-ct.csv"))

ns_cope <- read.csv(here("data", "nhl-cope-ct.csv"))
ns_cope_oct <- read.csv(here("data", "nhl-cope-oct-ct.csv"))
ns_cope_july <- read.csv(here("data", "nhl-cope-july-ct.csv"))

# load in our ichthyoplankton species data
spp_cta <- read.csv(here("data", "nhl-ct-summary.csv"))

cand.mods.nhl <- spp_cta %>%
  dplyr::select(Taxa, Period_3y, ct_anomaly_days, group) %>%
  dplyr::rename(ichthyo_cta = ct_anomaly_days)

# indicate which spp are year-shifted
spp_start_oct <- c("Glyptocephalus zachirus", "Lipolagus ochotensis")
spp_start_july <- c("Hemilepidotus hemilepidotus", "Hexagrammos decagrammus",
                    "Leptocottus armatus", "Liparis fucensis", "Ophiodon elongatus",
                    "Parophrys vetulus", "Scorpaenichthys marmoratus")

cand.mods.nhl <- cand.mods.nhl %>%
  mutate(yearstart = ifelse(Taxa %in% spp_start_oct, "October",
                            ifelse(Taxa %in% spp_start_july, "July", "January")))

# attach shifted env var data to spp. combine into one df. i will add each var one by one for ease of checking work.
cand.mods.nhl$temp_cta <- ifelse(cand.mods.nhl$yearstart == "October",
                           temp_oct$ct_anomaly_days,
                           ifelse(cand.mods.nhl$yearstart == "July",
                                  temp_july$ct_anomaly_days,
                                  temp$ct_anomaly_days))

cand.mods.nhl$onlypos_bcui_cta <- ifelse(cand.mods.nhl$yearstart == "October",
                                    onlypos_bcui_oct$ct_anomaly_days,
                                    ifelse(cand.mods.nhl$yearstart == "July",
                                           onlypos_bcui_july$ct_anomaly_days,
                                           onlypos_bcui$ct_anomaly_days))

cand.mods.nhl$n_cope_cta <- ifelse(cand.mods.nhl$yearstart == "October",
                              ns_cope_oct[ns_cope_oct$Assemblage == "Northern",]$ct_anomaly_days,
                              ifelse(cand.mods.nhl$yearstart == "July",
                                     ns_cope_july[ns_cope_july$Assemblage == "Northern",]$ct_anomaly_days,
                                     ns_cope[ns_cope$Assemblage == "Northern",]$ct_anomaly_days))

cand.mods.nhl$s_cope_cta <- ifelse(cand.mods.nhl$yearstart == "October",
                              ns_cope_oct[ns_cope_oct$Assemblage == "Southern",]$ct_anomaly_days,
                              ifelse(cand.mods.nhl$yearstart == "July",
                                     ns_cope_july[ns_cope_july$Assemblage == "Southern",]$ct_anomaly_days,
                                     ns_cope[ns_cope$Assemblage == "Southern",]$ct_anomaly_days))


# CalCOFI
Modnames.calcofi <- c("T + ZDV + BCUI", "T + ZDV", "T + BCUI", "T",
                      "ZDV + BCUI", "ZDV", "BCUI", "Intercept only")

earlier_spp.calcofi <- cand.mods.calcofi %>%
  filter(group == "earlier")
no_change_spp.calcofi <- cand.mods.calcofi %>%
  filter(group == "no_change")
later_spp.calcofi <- cand.mods.calcofi %>%
  filter(group == "later")

Earlier.Cand.mod.calcofi <- list()
Earlier.Cand.mod.calcofi[[1]] <- lm(ichthyo_cta ~ temp_cta + zdv_cta + bcui_cta,
                                    data = earlier_spp.calcofi)
Earlier.Cand.mod.calcofi[[2]] <- lm(ichthyo_cta ~ temp_cta + zdv_cta,
                                    data = earlier_spp.calcofi)
Earlier.Cand.mod.calcofi[[3]] <- lm(ichthyo_cta ~ temp_cta + bcui_cta,
                                    data = earlier_spp.calcofi)
Earlier.Cand.mod.calcofi[[4]] <- lm(ichthyo_cta ~ temp_cta,
                                    data = earlier_spp.calcofi)
Earlier.Cand.mod.calcofi[[5]] <- lm(ichthyo_cta ~ zdv_cta + bcui_cta,
                                    data = earlier_spp.calcofi)
Earlier.Cand.mod.calcofi[[6]] <- lm(ichthyo_cta ~ zdv_cta,
                                    data = earlier_spp.calcofi)
Earlier.Cand.mod.calcofi[[7]] <- lm(ichthyo_cta ~ bcui_cta,
                                    data = earlier_spp.calcofi)
Earlier.Cand.mod.calcofi[[8]] <- lm(ichthyo_cta ~ 1,
                                    data = earlier_spp.calcofi)

Earlier.table.calcofi <- aictab(cand.set = Earlier.Cand.mod.calcofi,
                                modnames = Modnames.calcofi,
                                second.ord = T)

Earlier.table.calcofi2 <- Earlier.table.calcofi[which(round(
  Earlier.table.calcofi$Delta_AICc, 2) < 2.00),]


NoChange.Cand.mod.calcofi <- list()
NoChange.Cand.mod.calcofi[[1]] <- lm(ichthyo_cta ~ temp_cta + zdv_cta + bcui_cta,
                                     data = no_change_spp.calcofi)
NoChange.Cand.mod.calcofi[[2]] <- lm(ichthyo_cta ~ temp_cta + zdv_cta,
                                     data = no_change_spp.calcofi)
NoChange.Cand.mod.calcofi[[3]] <- lm(ichthyo_cta ~ temp_cta + bcui_cta,
                                     data = no_change_spp.calcofi)
NoChange.Cand.mod.calcofi[[4]] <- lm(ichthyo_cta ~ temp_cta,
                                     data = no_change_spp.calcofi)
NoChange.Cand.mod.calcofi[[5]] <- lm(ichthyo_cta ~ zdv_cta + bcui_cta,
                                     data = no_change_spp.calcofi)
NoChange.Cand.mod.calcofi[[6]] <- lm(ichthyo_cta ~ zdv_cta,
                                     data = no_change_spp.calcofi)
NoChange.Cand.mod.calcofi[[7]] <- lm(ichthyo_cta ~ bcui_cta,
                                     data = no_change_spp.calcofi)
NoChange.Cand.mod.calcofi[[8]] <- lm(ichthyo_cta ~ 1,
                                     data = no_change_spp.calcofi)

NoChange.table.calcofi <- aictab(cand.set = NoChange.Cand.mod.calcofi,
                                 modnames = Modnames.calcofi,
                                 second.ord = T)

NoChange.table.calcofi2 <- NoChange.table.calcofi[which(round(
  NoChange.table.calcofi$Delta_AICc, 2) < 2.00),]


Later.Cand.mod.calcofi <- list()
Later.Cand.mod.calcofi[[1]] <- lm(ichthyo_cta ~ temp_cta + zdv_cta + bcui_cta,
                                  data = later_spp.calcofi)
Later.Cand.mod.calcofi[[2]] <- lm(ichthyo_cta ~ temp_cta + zdv_cta,
                                  data = later_spp.calcofi)
Later.Cand.mod.calcofi[[3]] <- lm(ichthyo_cta ~ temp_cta + bcui_cta,
                                  data = later_spp.calcofi)
Later.Cand.mod.calcofi[[4]] <- lm(ichthyo_cta ~ temp_cta,
                                  data = later_spp.calcofi)
Later.Cand.mod.calcofi[[5]] <- lm(ichthyo_cta ~ zdv_cta + bcui_cta,
                                  data = later_spp.calcofi)
Later.Cand.mod.calcofi[[6]] <- lm(ichthyo_cta ~ zdv_cta,
                                  data = later_spp.calcofi)
Later.Cand.mod.calcofi[[7]] <- lm(ichthyo_cta ~ bcui_cta,
                                  data = later_spp.calcofi)
Later.Cand.mod.calcofi[[8]] <- lm(ichthyo_cta ~ 1,
                                  data = later_spp.calcofi)

Later.table.calcofi <- aictab(cand.set = Later.Cand.mod.calcofi,
                              modnames = Modnames.calcofi,
                              second.ord = T)

Later.table.calcofi2 <- Later.table.calcofi[which(round(
  Later.table.calcofi$Delta_AICc, 2) < 2.00),]


colnames(Earlier.table.calcofi2) <- c("Model", "K", "AICc", "DeltaAICc",
                                      "ModelLik", "AICcWt", "LL", "Cum.Wt")
Earlier.table.calcofi2$`Dependent variable` <- "Earlier phenology group"
Earlier.table.calcofi3 <- Earlier.table.calcofi2[,c(9, 1:4, 6, 8, 7)]

colnames(NoChange.table.calcofi2) <- c("Model", "K", "AICc", "DeltaAICc",
                                       "ModelLik", "AICcWt", "LL", "Cum.Wt")
NoChange.table.calcofi2$`Dependent variable` <- c("No linear change group",
                                                  rep("", times = 3))
NoChange.table.calcofi3 <- NoChange.table.calcofi2[,c(9, 1:4, 6, 8, 7)]

colnames(Later.table.calcofi2) <- c("Model", "K", "AICc", "DeltaAICc",
                                    "ModelLik", "AICcWt", "LL", "Cum.Wt")
Later.table.calcofi2$`Dependent variable` <- "Later phenology group"
Later.table.calcofi3 <- Later.table.calcofi2[,c(9, 1:4, 6, 8, 7)]

Global.table.calcofi <- rbind(Earlier.table.calcofi3,
                              NoChange.table.calcofi3,
                              Later.table.calcofi3)


# Newport Hydrographic Line
Modnames.nhl <- c("T + BCUI + N.cope + S.cope", "T + BCUI + N.cope",
                  "T + BCUI + S.cope", "T + N.cope + S.cope",
                  "T + N.cope", "T + S.cope", "T + BCUI", "T",
                  "BCUI + N.cope + S.cope", "BCUI + N.cope", "BCUI + S.cope",
                  "BCUI", "N.cope + S.cope", "N.cope", "S.cope", "Intercept only")

earlier_spp.nhl <- cand.mods.nhl %>%
  filter(group == "earlier")
no_change_spp.nhl <- cand.mods.nhl %>%
  filter(group == "no_change")
later_spp.nhl <- cand.mods.nhl %>%
  filter(group == "later")

Earlier.Cand.mod.nhl <- list()
Earlier.Cand.mod.nhl[[1]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta +
                                  n_cope_cta + s_cope_cta,
                                data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[2]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta +
                                  n_cope_cta,
                                data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[3]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta +
                                  s_cope_cta,
                                data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[4]] <- lm(ichthyo_cta ~ temp_cta + n_cope_cta + s_cope_cta,
                                data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[5]] <- lm(ichthyo_cta ~ temp_cta + n_cope_cta,
                                data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[6]] <- lm(ichthyo_cta ~ temp_cta + s_cope_cta,
                                data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[7]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta,
                                data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[8]] <- lm(ichthyo_cta ~ temp_cta,
                                data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[9]] <- lm(ichthyo_cta ~ onlypos_bcui_cta + n_cope_cta +
                                  s_cope_cta,
                                data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[10]] <- lm(ichthyo_cta ~ onlypos_bcui_cta + n_cope_cta,
                                 data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[11]] <- lm(ichthyo_cta ~ onlypos_bcui_cta + s_cope_cta,
                                 data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[12]] <- lm(ichthyo_cta ~ onlypos_bcui_cta,
                                 data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[13]] <- lm(ichthyo_cta ~ n_cope_cta + s_cope_cta,
                                 data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[14]] <- lm(ichthyo_cta ~ n_cope_cta,
                                 data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[15]] <- lm(ichthyo_cta ~ s_cope_cta,
                                 data = earlier_spp.nhl)
Earlier.Cand.mod.nhl[[16]] <- lm(ichthyo_cta ~ 1,
                                 data = earlier_spp.nhl)

Earlier.table.nhl <- aictab(cand.set = Earlier.Cand.mod.nhl,
                            modnames = Modnames.nhl,
                            second.ord = T)

Earlier.table.nhl2 <- Earlier.table.nhl[which(round(
  Earlier.table.nhl$Delta_AICc, 2) < 2.00),]

NoChange.Cand.mod.nhl <- list()
NoChange.Cand.mod.nhl[[1]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta +
                                   n_cope_cta + s_cope_cta,
                                 data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[2]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta +
                                   n_cope_cta,
                                 data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[3]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta +
                                   s_cope_cta,
                                 data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[4]] <- lm(ichthyo_cta ~ temp_cta + n_cope_cta + s_cope_cta,
                                 data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[5]] <- lm(ichthyo_cta ~ temp_cta + n_cope_cta,
                                 data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[6]] <- lm(ichthyo_cta ~ temp_cta + s_cope_cta,
                                 data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[7]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta,
                                 data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[8]] <- lm(ichthyo_cta ~ temp_cta,
                                 data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[9]] <- lm(ichthyo_cta ~ onlypos_bcui_cta + n_cope_cta +
                                   s_cope_cta,
                                 data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[10]] <- lm(ichthyo_cta ~ onlypos_bcui_cta + n_cope_cta,
                                  data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[11]] <- lm(ichthyo_cta ~ onlypos_bcui_cta + s_cope_cta,
                                  data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[12]] <- lm(ichthyo_cta ~ onlypos_bcui_cta,
                                  data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[13]] <- lm(ichthyo_cta ~ n_cope_cta + s_cope_cta,
                                  data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[14]] <- lm(ichthyo_cta ~ n_cope_cta,
                                  data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[15]] <- lm(ichthyo_cta ~ s_cope_cta,
                                  data = no_change_spp.nhl)
NoChange.Cand.mod.nhl[[16]] <- lm(ichthyo_cta ~ 1,
                                  data = no_change_spp.nhl)

NoChange.table.nhl <- aictab(cand.set = NoChange.Cand.mod.nhl,
                             modnames = Modnames.nhl,
                             second.ord = T)

NoChange.table.nhl2 <- NoChange.table.nhl[which(round(
  NoChange.table.nhl$Delta_AICc, 2) < 2.00),]

Later.Cand.mod.nhl <- list()
Later.Cand.mod.nhl[[1]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta +
                                n_cope_cta + s_cope_cta,
                              data = later_spp.nhl)
Later.Cand.mod.nhl[[2]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta +
                                n_cope_cta,
                              data = later_spp.nhl)
Later.Cand.mod.nhl[[3]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta +
                                s_cope_cta + lag(ichthyo_cta),
                              data = later_spp.nhl)
Later.Cand.mod.nhl[[4]] <- lm(ichthyo_cta ~ temp_cta + n_cope_cta + s_cope_cta + lag(ichthyo_cta),
                              data = later_spp.nhl)
Later.Cand.mod.nhl[[5]] <- lm(ichthyo_cta ~ temp_cta + n_cope_cta,
                              data = later_spp.nhl)
Later.Cand.mod.nhl[[6]] <- lm(ichthyo_cta ~ temp_cta + s_cope_cta + lag(ichthyo_cta),
                              data = later_spp.nhl)
Later.Cand.mod.nhl[[7]] <- lm(ichthyo_cta ~ temp_cta + onlypos_bcui_cta + lag(ichthyo_cta),
                              data = later_spp.nhl)
Later.Cand.mod.nhl[[8]] <- lm(ichthyo_cta ~ temp_cta + lag(ichthyo_cta),
                              data = later_spp.nhl)
Later.Cand.mod.nhl[[9]] <- lm(ichthyo_cta ~ onlypos_bcui_cta + n_cope_cta +
                                s_cope_cta,
                              data = later_spp.nhl)
Later.Cand.mod.nhl[[10]] <- lm(ichthyo_cta ~ onlypos_bcui_cta + n_cope_cta,
                               data = later_spp.nhl)
Later.Cand.mod.nhl[[11]] <- lm(ichthyo_cta ~ onlypos_bcui_cta + s_cope_cta,
                               data = later_spp.nhl)
Later.Cand.mod.nhl[[12]] <- lm(ichthyo_cta ~ onlypos_bcui_cta,
                               data = later_spp.nhl)
Later.Cand.mod.nhl[[13]] <- lm(ichthyo_cta ~ n_cope_cta + s_cope_cta,
                               data = later_spp.nhl)
Later.Cand.mod.nhl[[14]] <- lm(ichthyo_cta ~ n_cope_cta,
                               data = later_spp.nhl)
Later.Cand.mod.nhl[[15]] <- lm(ichthyo_cta ~ s_cope_cta + lag(ichthyo_cta),
                               data = later_spp.nhl)
Later.Cand.mod.nhl[[16]] <- lm(ichthyo_cta ~ 1,
                               data = later_spp.nhl)

Later.table.nhl <- aictab(cand.set = Later.Cand.mod.nhl,
                          modnames = Modnames.nhl,
                          second.ord = T)

Later.table.nhl2 <- Later.table.nhl[which(round(
  Later.table.nhl$Delta_AICc, 2) < 2.00),]


colnames(Earlier.table.nhl2) <- c("Model", "K", "AICc", "DeltaAICc",
                                  "ModelLik", "AICcWt", "LL", "Cum.Wt")
Earlier.table.nhl2$`Dependent variable` <- c("Earlier phenology group",
                                             "")
Earlier.table.nhl3 <- Earlier.table.nhl2[,c(9, 1:4, 6, 8, 7)]

colnames(NoChange.table.nhl2) <- c("Model", "K", "AICc", "DeltaAICc",
                                   "ModelLik", "AICcWt", "LL", "Cum.Wt")
NoChange.table.nhl2$`Dependent variable` <- c("No linear change group",
                                              rep("", times = 5))
NoChange.table.nhl3 <- NoChange.table.nhl2[,c(9, 1:4, 6, 8, 7)]

colnames(Later.table.nhl2) <- c("Model", "K", "AICc", "DeltaAICc",
                                "ModelLik", "AICcWt", "LL", "Cum.Wt")
Later.table.nhl2$`Dependent variable` <- c("Later phenology group", "")
Later.table.nhl3 <- Later.table.nhl2[,c(9, 1:4, 6, 8, 7)]

Global.table.nhl <- rbind(Earlier.table.nhl3, NoChange.table.nhl3, Later.table.nhl3)

stargazer(Global.table.calcofi, Global.table.nhl,
          type = "html",
          title = c("CalCOFI", "NH Line"),
          summary = F,
          digits = 2, digit.separator = "",
          rownames = F)
# at this point I have copied the html to text editor and adjusted the following:
# style header:
# <style>
#   td {
#     font-size: 19.5px;
#   }
# caption {
#   font-size: 21.5px;
# }
# </style>
#
# width = "850" inside each table style tag
#
# <tr><td colspan="8" style="border-bottom: 0.5px solid black"></td></tr> between each dependent variable group

# 1.24.25 Sample model validation ==========================================
## This section requires running the previous code section.
Earlier.bestmod.calcofi <- lm(ichthyo_cta ~ temp_cta + zdv_cta + bcui_cta,
                              data = earlier_spp.calcofi)
Nochange.bestmod.calcofi <- lm(ichthyo_cta ~ zdv_cta,
                               data = no_change_spp.calcofi)
Later.bestmod.calcofi <- lm(ichthyo_cta ~ temp_cta,
                            data = later_spp.calcofi)

Earlier.bestmod.nhl <- lm(ichthyo_cta ~ n_cope_cta,
                          data = earlier_spp.nhl)
Nochange.bestmod.nhl <- lm(ichthyo_cta ~ n_cope_cta + s_cope_cta,
                           data = no_change_spp.nhl)


mod <- Nochange.bestmod.nhl
E <- rstandard(mod)
F <- fitted(mod)
plot(F,E)

hist(E)

qqnorm(E)

library(olsrr)
ols_coll_diag(Earlier.bestmod.calcofi)

# Environmental x Ichthyoplankton Most Plausible Models ========================
## This section requires running the previous code section.
## Figure description:
Earlier.bestmod.calcofi <- lm(ichthyo_cta ~ temp_cta + zdv_cta + bcui_cta,
                              data = earlier_spp.calcofi)
Nochange.bestmod.calcofi <- lm(ichthyo_cta ~ zdv_cta,
                               data = no_change_spp.calcofi)
Later.bestmod.calcofi <- lm(ichthyo_cta ~ temp_cta,
                            data = later_spp.calcofi)

Earlier.bestmod.nhl <- lm(ichthyo_cta ~ n_cope_cta,
                          data = earlier_spp.nhl)
Nochange.bestmod.nhl <- lm(ichthyo_cta ~ n_cope_cta + s_cope_cta,
                           data = no_change_spp.nhl)
Later.bestmod.nhl <- lm(ichthyo_cta ~ n_cope_cta,
                        data = later_spp.nhl)

# tables with stargazer
stargazer(Earlier.bestmod.calcofi, Nochange.bestmod.calcofi,
          Later.bestmod.calcofi,
          type = "html",
          column.labels = c("Earlier", "No Linear Change", "Later"),
          covariate.labels = c("Temperature CTa", "Zooplankton CTa", "Upwelling CTa", "Intercept"),
          dep.var.caption = "CalCOFI Ichthyoplankton CTa",
          dep.var.labels = "<em> Phenology Change Group: </em>",
          model.numbers = F, # omit model numbers
          omit.stat = c("ser", "f"), # omit F statistic and RSE
          report = "vc*", # omit coefficient standard errors
          notes.label = "", # omit the "Note:"
          title = "")

stargazer(Earlier.bestmod.nhl, Nochange.bestmod.nhl, Later.bestmod.nhl,
          type = "html",
          column.labels = c("Earlier", "No Linear Change", "Later"),
          covariate.labels = c("Northern Copepod CTa", "Southern Copepod CTa",
                               "Intercept"),
          dep.var.caption = "NH Line Ichthyoplankton CTa",
          dep.var.labels = "<em> Phenology Change Group: </em>",
          model.numbers = F, # omit model numbers
          omit.stat = c("F", "ser"), # omit F statistic and Residual Std. Error
          report = "vc*", # omit coefficient standard errors
          notes.label = "", # omit the "Note:"
          title = "")


# at this point I have copied the html to text editor and adjusted the following:
# style header:
# <style>
#   td {
#     font-size: 19.5px;
#   }
# caption {
#   font-size: 21.5px;
# }
# </style>
#
# width = "470" inside each table style tag
#
# bold the titles
# italicize R2 and p and p-values
# put spaces between the stars and p-values


# Climate Mode Plots ===========================================================
## Figure description:

## Read in data needed to make this figure
climate.calcofi <- read.csv(here("data", "climate.calcofi.csv"))
climate.nhl <- read.csv(here("data", "climate.nhl.csv"))


climate.1 <- climate.calcofi %>%
  filter(mode == "enso") %>%
  mutate(phase = factor(phase,
                        levels = c("La Nina", "Neutral", "El Nino"),
                        labels = c("La Niña", "Neutral", "El Niño")),
         group = factor(group,
                        levels = c("earlier",
                                     "no_change",
                                     "later"),
                        labels = c("Earlier Phenology",
                                   "No Linear Change",
                                   "Later Phenology"))) %>%
  ggplot(aes(x = phase,
             y = ct_anomaly_days,
             group = phase)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  facet_grid(.~ group) +
  geom_boxplot(linewidth = 1.2, outlier.size = 2) +
  my.ggplot.theme +
  labs(x = "",
       y = "CTa (days)")


annotation.enso.nhl <- data.frame(
  enso = c(rep(c("La Nina", "Neutral", "El Nino"), times = 3)),
  group = c(rep(c("earlier", "no_change", "later"), each = 3)),
  y = -60,
  label = c("","","",
            "AB", "A", "B",
            "", "", "")
)

annotation.enso.nhl$enso = factor(annotation.enso.nhl$enso,
                                  levels = c("La Nina", "Neutral", "El Nino"),
                                  labels = c("La Niña", "Neutral", "El Niño"))

annotation.enso.nhl$group = factor(annotation.enso.nhl$group,
                                   levels = c("earlier",
                                              "no_change",
                                              "later"),
                                   labels = c("Earlier Phenology",
                                              "No Linear Change",
                                              "Later Phenology"))

climate.2 <- climate.nhl %>%
  filter(mode == "enso") %>%
  mutate(phase = factor(phase,
                        levels = c("L", "N", "E"),
                        labels = c("La Niña", "Neutral", "El Niño")),
         group = factor(group,
                        levels = c("earlier",
                                   "no_change",
                                   "later"),
                        labels = c("Earlier Phenology",
                                   "No Linear Change",
                                   "Later Phenology"))) %>%
  ggplot(aes(x = phase,
             y = ct_anomaly_days,
             group = phase)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  geom_text(data = annotation.enso.nhl,
            aes(x = enso, y = y, group = group, label = label),
            color = "black",
            size = 7,
            fontface = "bold") +
  facet_grid(.~ group) +
  geom_boxplot(linewidth = 1.2, outlier.size = 2) +
  my.ggplot.theme +
  labs(x = "Phase",
       y = "CTa (days)")

climate <- ggarrange(climate.1, # CalCOFI x ENSO,
                     climate.2, # NHL x ENSO
                     labels = c("(a)", "(b)"),
                     font.label = list(size = 20),
                     # hjust = -1,
                     vjust = 1.25,
                     ncol = 1, nrow = 2, align = "hv")

ggsave(filename = "climate.png",
       plot = climate,
       path = here("figures"), height = 200, width = 180, units = "mm", dpi = 600)


# Species Table ===============================================================
## Table description: List of species with their ecological traits and calculated trends.

## Read in data needed to make this table
calcofi.spp <- read.csv(here("data", "calcofi-eco-traits.csv"))
calcofi.common <- read.csv(here("data", "calcofi-common-names.csv"))
calcofi.summary <- read.csv(here("data", "calcofi-ct-summary.csv"))
calcofi.seasons <- read.csv(here("data", "calcofi-seasons.csv"))

nhl.spp <- read_excel(here("data", "nhl-eco-traits.xlsx"))
nhl.summary <- read.csv(here("data", "nhl-ct-summary.csv"))
nhl.seasons <- read.csv(here("data", "nhl-seasons.csv"))
nhl.seasons$Taxa <- gsub("[.]", " ", nhl.seasons$Taxa)
nhl.summary$Taxa <- gsub("[.]", " ", nhl.summary$Taxa)


## Organizing the species table
calcofi.spp2 <- calcofi.spp %>%
  dplyr::select(!c(adult_trophic_level, fishing_status, mean_ct, sd_ct)) %>%
  dplyr::rename(`Scientific Name` = scientific_name,
                Order = order,
                `Adult Habitat` = adult_habitat,
                `Cross-shore Distribution` = cross_shore_distribution,
                `Biogeographic Affinity` = biogeographic_affinity,
                `CalCOFI MonthMax` = month_maximum,
                `CalCOFI r` = r,
                `CalCOFI Trend` = trend_in_phenology) %>%
  mutate(`CalCOFI Trend` = format(round(`CalCOFI Trend` / 10, digits = 2), nsmall = 2), # get trend in d/y
         `CalCOFI r` = format(round(`CalCOFI r`, digits = 2), nsmall = 2),
         `CalCOFI Group` = calcofi.summary$group[match(calcofi.spp$scientific_name,calcofi.summary$scientific_name)],
         `Common Name` = calcofi.common$common_name[match(calcofi.spp$scientific_name,
                                                            calcofi.common$scientific_name)],
         `CalCOFI Season` = calcofi.seasons$season[match(calcofi.spp$scientific_name,
                                                         calcofi.seasons$scientific_name)]) %>%
  mutate(`CalCOFI Group` = ifelse(`CalCOFI Group` == "earlier",
                                  "Earlier",
                                  ifelse(`CalCOFI Group` == "later",
                                         "Later", "No Change")),
         `Adult Habitat` = ifelse(`Adult Habitat` == "Demersal",
                                  "D",
                                  ifelse(`Adult Habitat` == "Mesopelagic",
                                         "M", "E")),
         `Cross-shore Distribution` = ifelse(`Cross-shore Distribution` == "Coastal",
                                             "C",
                                             ifelse(`Cross-shore Distribution` == "Oceanic",
                                                    "O", "CO")),
         `Biogeographic Affinity` = ifelse(`Biogeographic Affinity` == "Cool-water",
                                           "C",
                                           ifelse(`Biogeographic Affinity` == "Warm-water",
                                                  "W", "WD")),
         `CalCOFI Season` = ifelse(`CalCOFI Season` == "fall_winter",
                                   "FW", "SS"))
calcofi.spp2 <- calcofi.spp2[,c(1:8, 16,17,18)]


nhl.spp2 <- nhl.spp %>%
  dplyr::select(!c(season_4_NHL, calculated_yearstart_NHL,
                   calculated_mean_ct_NHL, unshift_ct_NHL_months,
                   sd_ct_NHL_months, sd_ct_NHL_days,
                   adult_trophic_level, fishing_status,
                   group.3333_NHL)) %>%
  dplyr::rename(`Scientific Name` = scientific_name,
                `Common Name` = FB_name,
                Order = order,
                `Adult Habitat` = adult_habitat,
                `Cross-shore Distribution` = cross_shore_distribution,
                `Biogeographic Affinity` = biogeographic_affinity,
                `NH Line MonthMax` = month_maximum_NHL,
                `NH Line r` = cor_NHL,
                `NH Line Trend` = trend_NHL_d3y) %>%
  mutate(`NH Line Trend` = format(round(`NH Line Trend` / 3, digits = 2), nsmall = 2), # get trend in d/y
         `NH Line r` = format(round(`NH Line r`, digits = 2), nsmall = 2),
         `NH Line Group` = nhl.summary$group[match(nhl.spp$scientific_name,
                                                   nhl.summary$Taxa)],
         `NH Line Season` = nhl.seasons$season[match(nhl.spp$scientific_name,
                                                     nhl.seasons$Taxa)]) %>%
  mutate(`NH Line Group` = ifelse(`NH Line Group` == "earlier",
                                  "Earlier",
                                  ifelse(`NH Line Group` == "later",
                                         "Later", "No Change")),
         `Adult Habitat` = ifelse(`Adult Habitat` == "Demersal",
                                  "D",
                                  ifelse(`Adult Habitat` == "Mesopelagic",
                                         "M", "E")),
         `Cross-shore Distribution` = ifelse(`Cross-shore Distribution` == "Coastal",
                                             "C",
                                             ifelse(`Cross-shore Distribution` == "Oceanic",
                                                    "O", "CO")),
         `Biogeographic Affinity` = ifelse(`Biogeographic Affinity` == "Cool-water",
                                           "C",
                                           ifelse(`Biogeographic Affinity` == "Warm-water",
                                                  "W", "WD")),
         `NH Line Season` = ifelse(`NH Line Season` == "fall_winter",
                                   "FW", "SS"))

all.spp <- merge(calcofi.spp2, nhl.spp2, all=TRUE)
all.spp2 <- all.spp[,c(1,6,2:5,7,14,11,16,8,12,10,15,9,13)]

all.spp2[is.na(all.spp2)] <- "n.p."

all.spp2$Order[all.spp2$Order == "Perciformes/Zoarcoidei"] <- "P.Zoarcoidei"
all.spp2$Order[all.spp2$Order == "Perciformes/Cottoidei"] <- "P.Cottoidei"
all.spp2$Order[all.spp2$Order == "Perciformes/Scorpaenoidei"] <- "P.Scorpaenoidei"

all.spp2$`Common Name`[all.spp2$`Scientific Name` == "Sardinops sagax"] <- "Pacific sardine"
all.spp2$`Common Name`[all.spp2$`Scientific Name` == "Engraulis mordax"] <- "Northern anchovy"
all.spp2$`Common Name`[all.spp2$`Scientific Name` == "Peprilus simillimus"] <- "Butterfish"

# stargazer version
stargazer(all.spp2, summary = F,
          type = "html",
          rownames = F,
          out = "~/Desktop/all.species.htm")

# gt version with references
all.spp2 %>%
  gt() %>%
  tab_options(table.font.names = "Times New Roman",
              table.font.size = 12,
              table.width = 1000) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) %>%
  tab_style(style = cell_text(style = "italic"),
            locations = cells_body(columns = `Scientific Name`)) %>%
  tab_footnote(footnote = "fishbase.se",
               locations = cells_column_labels(columns = c(`Common Name`,
                                                           "Order",
                                                           `Cross-shore Distribution`))) %>%
  tab_footnote(footnote = "fishbase.se",
               locations = cells_body(columns = `Adult Habitat`,
                                      rows = c(1,3:6,9,13,15,17:23,26:28,30,33,35,36,38:42,45:48,51:54,57,62,64,68,69))) %>%
  tab_footnote(footnote = "fishbase.se",
               locations = cells_body(columns = `Biogeographic Affinity`,
                                      rows = c(6,15,17,19,27,35,52))) %>%
  tab_footnote(footnote = "Asch (2015)",
               locations = cells_body(columns = `Adult Habitat`,
                                      rows = c(2,7,8,10:12,14,16,24,25,29,31,32,34,37,43,44,49,50,55,56,58:61,63,65:67,70:73))) %>%
  tab_footnote(footnote = "Hsieh et al. (2005)",
               locations = cells_column_labels(columns = `Cross-shore Distribution`)) %>%
  tab_footnote(footnote = "Hsieh et al. (2005)",
               locations = cells_body(columns = `Biogeographic Affinity`,
                                      rows = c(1,4,5,9,18,20,21,23,26,28,30,33,41,45:48,51,53,54,58:63))) %>%
  tab_footnote(footnote = "Hsieh et al. (2008)",
               locations = cells_body(columns = `Biogeographic Affinity`,
                                      rows = c(2,12,13,16,22,24,29,32,34,36,42:44,49,55:57,64,68,71,72))) %>%
  tab_footnote(footnote = "Hsieh et al. (2009)",
               locations = cells_body(columns = `Biogeographic Affinity`,
                                      rows = c(3,7,8,10,11,14,25,31,37:40,50,65:67,69,70,73))) %>%
  tab_footnote(footnote = html("Columns: Adult Habitat (D, demersal; M, mesopelagic; E, epipelagic). Cross-shore Distribution (C, coastal; O, oceanic; CO, coastal-oceanic). Biogeographic Affinity (C, cool-water; W, warm-water; WD, widely distributed in the Northeast Pacific). MonthMax = calculated month of maximum mean larval abundance. Season = primary spawning season based on mean larval abundance in October-March (FW, fall/winter) compared to April-September (SS, spring/summer). r = Pearson correlation between CTa and time. Group = direction of phenological change based on r ≤ -0.33 or ≥ 0.33. Trend = estimated phenological trend (d y<sup>-1</sup>).")) %>%
  tab_footnote(footnote = "Suborders: P.Suborder = Perciformes/Suborder.") %>%
  tab_footnote(footnote = "n.p. = not present in that analysis.") %>%
  tab_options(footnotes.multiline = F) %>%
  gtsave("all-spp-table.html", "~/Desktop")

# 1.14.25 Asch (2015) Time Period =============================================
comparison = data.frame(names = c("Asch_2008", "New_2008", "New_2022"),
                        Trend = c("-0.12", "-0.14", "-0.18"),
                        `F` = c(6.5, 4.31, 11.77),
                        df = c(288, 322, 379),
                        p = c("< 0.05", "< 0.05", "< 0.001"),
                        Earlier = c("20 (39%)", "25 (44%)", "27 (47%)"),
                        `No Linear Change` = c("22 (43%)", "20 (35%)", "21 (37%)"),
                        Later = c("9 (18%)", "12 (21%)", "9 (16%)"))

comparison %>%
  gt(rowname_col = "names") %>%
  tab_options(table.font.names = "Times New Roman",
              table.font.size = 20,
              table.width = 800) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) %>%
  tab_style(style = cell_text(align = "center"),
            locations = list(cells_body(),
                          cells_column_labels())) %>%
  tab_style(style = cell_text(style = "italic"),
            locations = cells_column_labels(columns = p)) %>%
  tab_style(style = cell_text(style = "italic"),
            locations = cells_body(columns = p)) %>%
  tab_footnote(html("Columns: Trend = estimated assemblage phenological trend (d y<sup>-1</sup>) from slope of linear regression. F = F statistic, df = degrees of freedom, <em>p</em> = p-value from that regression. Earlier, No Linear Change, Later = number (percentage) of species falling into that category of phenological change.")) %>%
  tab_footnote("Rows: Asch_2008 = Asch (2015) data. New_2008 = current study data analyzed over the Asch (2015) time period. New_2022 = current study data analyzed over the extended time period.") %>%
  gtsave("asch-comparison-table.html", "~/Desktop")


# CalCOFI Uniform Months Only =================================================
## Figure description: CalCOFI temperature and mesozooplankton CTa against time when calculated using uniform months only. For SI.

## Read in data needed to make this figure
calcofi.uniform.t <- read.csv(here("data", "calcofi-uniform-temp.csv"))
calcofi.uniform.zdv <- read.csv(here("data", "calcofi-uniform-zdv.csv"))

uniform.t <- calcofi.uniform.t %>%
  ggplot(aes(x = as.numeric(as.character(decade)),
             y = ct_anomaly_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_errorbar(aes(ymin = ct_anomaly_days - se_cta_days,
                    ymax = ct_anomaly_days + se_cta_days),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  scale_y_continuous(limits = c(-45, 45),
                     breaks = c(-40, -20, 0, 20, 40)) +
  scale_x_continuous(limits = c(1950, 2020),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  my.ggplot.theme +
  labs(x = "Decade",
       y = "CT Anomaly (days)",
       title = "CalCOFI Temperature Uniform Months Only")

uniform.zdv <- calcofi.uniform.zdv %>%
  ggplot(aes(x = as.numeric(as.character(decade)),
             y = ct_anomaly_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_errorbar(aes(ymin = ct_anomaly_days - cta_se,
                    ymax = ct_anomaly_days + cta_se),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  scale_y_continuous(breaks = c(-30, -15, 0, 15, 30)) +
  scale_x_continuous(limits = c(1950, 2020),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  labs(x = "Decade",
       y = "CT Anomaly (days)",
       title = "CalCOFI ZDV Uniform Months Only")

uniform <- ggarrange(uniform.t, uniform.zdv,
                     labels = c("(a)", "(b)"),
                     font.label = list(size = 20),
                     # hjust = -1,
                     vjust = 1.25,
                     ncol = 2, nrow = 1, align = "hv")

ggsave(filename = "uniform.png",
       plot = uniform,
       path = here("figures"),
       height = 7, width = 12, units = "in")

# 1.23.25 CalCOFI 3y and 10y CTa x Time =======================================

## Read in data
cta.time.calcofi <- read.csv(here("data", "fig1.calcofi.csv"))
cta.time.calcofi.3y <- read.csv(here("data", "calcofi-3y.csv"))

cta.time.nhl <- read.csv(here("data", "fig1.nhl.csv"))
cta.time.nhl.3y <- read.csv(here("data", "nhl-3y.csv"))

## Convert CalCOFI periods to year labels and combine dataframes
periods <- data.frame(year = c(1950:2023),
                      period = c(rep(1, 2), rep(2:25, each = 3)),
label = c(rep(1950, 2), rep(seq(from = 1953, to = 2023, by = 3), each = 3)))

cta.time.calcofi.3y <- cta.time.calcofi.3y %>%
  mutate(year = periods$label[match(period, periods$period)]) %>%
  select(!period) %>%
  mutate(bins = "Comparison",
         color = "#808080")
cta.time.calcofi <- cta.time.calcofi %>%
  dplyr::rename(year = Decade) %>%
  select(!X) %>%
  mutate(bins = "Full",
         color = "#000000")

together.calcofi <- rbind(cta.time.calcofi, cta.time.calcofi.3y)

## Combine NHL dataframes
cta.time.nhl.3y <- cta.time.nhl.3y %>%
  select(!X) %>%
  mutate(end = "2017",
         color = "#808080")
cta.time.nhl <- cta.time.nhl %>%
  select(!X) %>%
  mutate(end = "2023",
         color = "#000000")

together.nhl <- rbind(cta.time.nhl, cta.time.nhl.3y)

## Combine 3y dataframes
combo.calcofi.3y <- cta.time.calcofi.3y %>%
  select(!bins) %>%
  mutate(loc = "CalCOFI",
         color = "#000000")

periods <- data.frame(period = 1:9,
                      label = seq(from = 1998, to = 2024, by = 3))

combo.nhl.3y <- cta.time.nhl.3y %>%
  select(!end) %>%
  mutate(year = periods$label[match(Period_3y, periods$period)],
         loc = "NHL",
         color = "#808080") %>%
  select(!Period_3y)

combo.3y <- rbind(combo.calcofi.3y, combo.nhl.3y)

## Make the four CalCOFI panels
cta.time.1 <-
  together.calcofi %>%
  filter(group == "assemblage") %>%
  ggplot(aes(x = year,
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color), linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_color_identity(labels = c("Full", "Comparison"),
                       guide = "legend") +
  scale_y_continuous(limits = c(-20, 20),
                     breaks = c(-10, 10)) +
  scale_x_continuous(limits = c(1950, 2020),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  theme(legend.position = c(.38, .8),
        legend.background = element_blank(),
        legend.text = element_text(size = 18, family = "Helvetica"),
        legend.title = element_blank()) +
  labs(title = "CalCOFI Assemblage",
       x = "Year",
       y = "CTa (days)")
cta.time.2 <- together.calcofi %>%
  filter(group == "earlier") %>%
  ggplot(aes(x = year,
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color), linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_color_identity() +
  scale_y_continuous(limits = c(-90, 90),
                     breaks = c(-60, 0, 60)) +
  scale_x_continuous(limits = c(1950, 2020),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  labs(title = "",
       x = "Year",
       y = "CTa (days)")
cta.time.3 <- together.calcofi %>%
  filter(group %in% c("no_change", "noChange")) %>%
  ggplot(aes(x = year,
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color), linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_color_identity() +
  scale_y_continuous(limits = c(-40, 40),
                     breaks = c(-20, 0, 20)) +
  scale_x_continuous(limits = c(1950, 2020),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  labs(title = "",
       x = "Year",
       y = "CTa (days)")
cta.time.4 <-
  together.calcofi %>%
  filter(group == "later") %>%
  ggplot(aes(x = year,
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color), linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_color_identity() +
  scale_y_continuous(limits = c(-60, 60),
                     breaks = c(-30, 0, 30)) +
  scale_x_continuous(limits = c(1950, 2020),
                     breaks = seq(from = 1950, to = 2020, by = 10)) +
  labs(title = "",
       x = "Year",
       y = "CTa (days)")

cta.time.5 <-
  together.nhl %>%
  filter(group == "assemblage") %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color), linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  scale_color_identity(guide = "legend") +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-20, 20),
                     breaks = c(-10, 0, 10)) +
  scale_x_continuous(limits = c(.5, 9.5),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "NH Line Assemblage",
       x = "3-y Period",
       y = "CTa (days)")
cta.time.6 <- together.nhl %>%
  filter(group == "earlier") %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color), linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  scale_color_identity() +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-90, 90),
                     breaks = c(-60, 0, 60)) +
  scale_x_continuous(limits = c(.5, 9.5),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "",
       x = "3-y Period",
       y = "CTa (days)")
cta.time.7 <- together.nhl %>%
  filter(group == "no_change") %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color), linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  scale_color_identity() +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-40, 40),
                     breaks = c(-20, 0, 20)) +
  scale_x_continuous(limits = c(.5, 9.5),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "",
       x = "3-y Period",
       y = "CTa (days)")
cta.time.8 <- together.nhl %>%
  filter(group == "later") %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color), linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_color_identity() +
  scale_y_continuous(limits = c(-60, 60),
                     breaks = c(-30, 0, 30)) +
  scale_x_continuous(limits = c(.5, 9.5),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "",
       x = "3-y Period",
       y = "CTa (days)")

## Put panels together
cta.time <- ggarrange(cta.time.1, cta.time.5, cta.time.2, cta.time.6,
                      cta.time.3, cta.time.7, cta.time.4, cta.time.8,
                      labels = c("(a)", "(b)", "(c)", "(d)",
                                 "(e)", "(f)", "(g)", "(h)"),
                      font.label = list(size = 20),
                      # hjust = -1,
                      vjust = 1.25,
                      ncol = 2, nrow = 4, align = "hv")

ggsave(filename = "cta.time.3y.png",
       plot = cta.time,
       path = ("~/Desktop"), height = 13, width = 12, units = "in")


## Plot the CalCOFI and NHL 3y analyses together
cta.time.9 <- combo.3y %>%
  filter(group == "assemblage") %>%
  ggplot(aes(x = year,
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color), linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  scale_color_identity(labels = c("CalCOFI", "NHL"),
                       guide = "legend") +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-20, 20),
                     breaks = c(-10, 10)) +
  scale_x_continuous(limits = c(1997, 2018),
                     breaks = seq(from = 1997, to = 2018, by = 3)) +
  theme(legend.position = c(.35, .9),
        legend.background = element_blank(),
        legend.text = element_text(size = 18, family = "Helvetica"),
        legend.title = element_blank()) +
  labs(title = "Assemblage",
       x = "3-y Period",
       y = "CTa (days)")
cta.time.10 <- combo.3y %>%
  filter(group == "earlier") %>%
  ggplot(aes(x = year,
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color), linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  scale_color_identity() +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-60, 60),
                     breaks = c(-30, 0, 30)) +
  scale_x_continuous(limits = c(1997, 2018),
                     breaks = seq(from = 1997, to = 2018, by = 3)) +
  labs(title = "Species Shifting Earlier",
       x = "3-y Period",
       y = "CTa (days)")
cta.time.11 <- combo.3y %>%
  filter(group %in% c("no_change", "noChange")) %>%
  ggplot(aes(x = year,
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color),linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  scale_color_identity() +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-45, 45),
                     breaks = c(-20, 0, 20)) +
  scale_x_continuous(limits = c(1997, 2018),
                     breaks = seq(from = 1997, to = 2018, by = 3)) +
  labs(title = "Species With No Linear Change",
       x = "3-y Period",
       y = "CTa (days)")
cta.time.12 <- combo.3y %>%
  filter(group == "later") %>%
  ggplot(aes(x = year,
             y = mean_cta_days,
             group = color)) +
  geom_point(aes(color = color), size = 2.5) +
  geom_path(aes(color = color),linewidth = 1.2) +
  geom_linerange(aes(color = color,
                     ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  scale_color_identity() +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-45, 45),
                     breaks = c(-20, 0, 20)) +
  scale_x_continuous(limits = c(1997, 2018),
                     breaks = seq(from = 1997, to = 2018, by = 3)) +
  labs(title = "Species Shifting Later",
       x = "3-y Period",
       y = "CTa (days)")

cta.time.3yonly <- ggarrange(cta.time.9, cta.time.10, cta.time.11, cta.time.12,
                      labels = c("(a)", "(b)", "(c)", "(d)"),
                      font.label = list(size = 20),
                      # hjust = -1,
                      vjust = 1.25,
                      ncol = 2, nrow = 2, align = "hv")

ggsave(filename = "cta.time.3yonly.png",
       plot = cta.time.3yonly,
       path = ("~/Desktop"), height = 8, width = 12, units = "in")

# (x) Dropped PCA Covariates =================================================
## Figure description: Plots against phenological trend for covariates dropped from the PCA (fishing status, adult trophic level, log10 frequency of abundance, ln(x+1) amplitude of the seasonal cycle. Two columns, left CalCOFI, right NH Line. For SI.

## Read in data needed to make this figure
calcofiNumCovs <- read.csv('/Volumes/petrik-lab/kchen/datasets/calcofi-numreg.csv')
nhlNumCovs <- read.csv('/Volumes/petrik-lab/kchen/datasets/nhl-numreg.csv')

designTrend.calcofi <- read.csv("/Volumes/petrik-lab/kchen/datasets/designTrend.calcofi.csv")
design.spp.calcofi <- read.csv("/Volumes/petrik-lab/kchen/datasets/design.spp.calcofi.csv")

designTrend.nhl <- read.csv("/Volumes/petrik-lab/kchen/datasets/designTrend.nhl.csv")
design.spp.nhl <- read.csv("/Volumes/petrik-lab/kchen/datasets/design.spp.nhl.csv")

## Start with the three numeric covs
calcofiDropTL <- ggplot(calcofiNumCovs,
                        aes(x = adultTL,
                            y = trend_dy)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  my.ggplot.theme +
  labs(x = "Adult Trophic Level",
       y = expression(paste("Trend in Phenology (d ", y^-1, ")")),,
       title = "CalCOFI Assemblage")

calcofiDropFreq <- ggplot(calcofiNumCovs,
                          aes(x = logFreq,
                              y = trend_dy)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  my.ggplot.theme +
  labs(x = "log10 Frequency of Larval Abundance",
       y = expression(paste("Trend in Phenology (d ", y^-1, ")")),,
       title = "")

calcofiDropAmp <- ggplot(calcofiNumCovs,
                         aes(x = logAmplitude,
                             y = trend_dy)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  my.ggplot.theme +
  labs(x = "ln(x+1) Amplitude of Seasonal Cycle",
       y = expression(paste("Trend in Phenology (d ", y^-1, ")")),,
       title = "")

nhlDropTL <- ggplot(nhlNumCovs,
                    aes(x = adultTL,
                        y = trend_dy)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  my.ggplot.theme +
  labs(x = "Adult Trophic Level",
       y = expression(paste("Trend in Phenology (d ", y^-1, ")")),,
       title = "NH Line Assemblage")

nhlDropFreq <- ggplot(nhlNumCovs,
                      aes(x = logFreq,
                          y = trend_dy)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  my.ggplot.theme +
  labs(x = "log10 Frequency of Larval Abundance",
       y = expression(paste("Trend in Phenology (d ", y^-1, ")")),,
       title = "")

nhlDropAmp <- ggplot(nhlNumCovs,
                     aes(x = logAmplitude,
                         y = trend_dy)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  my.ggplot.theme +
  labs(x = "ln(x+1) Amplitude of Seasonal Cycle",
       y = expression(paste("Trend in Phenology (d ", y^-1, ")")),,
       title = "")

## Get fishing status
calcofiFishing <- design.spp.calcofi %>%
  select(fishing_status, trend_in_phenology)

calcofiFishing.means <- calcofiFishing %>%
  group_by(fishing_status) %>%
  reframe(meanTrend = mean(trend_in_phenology))

nhlFishing <- design.spp.nhl %>%
  select(fishing_status, trend_NHL_ddec)

nhlFishing.means <- nhlFishing %>%
  group_by(fishing_status) %>%
  reframe(meanTrend = mean(trend_NHL_ddec))

# set colors. these colors are taken from the grafify colorblind-friendly palette "kelly", which can be viewed with library(grafify) and plot_grafify_palette(palette = "kelly")

calcofiDropFishing <- ggplot(calcofiFishing,
                             aes(y = trend_in_phenology,
                                 x = fishing_status)) +
  geom_boxplot(aes(fill = fishing_status)) +
  scale_fill_grafify(palette = "kelly") +
  geom_point() +
  geom_point(data = calcofiFishing.means,
             aes(x = fishing_status, y = meanTrend),
             size = 4, fill = "red", shape = 23) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  labs(y = expression(paste("Trend in Phenology (d ", y^-1, ")")))

nhlDropFishing <- ggplot(nhlFishing,
                         aes(y = trend_NHL_ddec,
                             x = fishing_status)) +
  geom_boxplot(aes(fill = fishing_status)) +
  scale_fill_grafify(palette = "kelly") +
  geom_point() +
  geom_point(data = nhlFishing.means,
             aes(x = fishing_status, y = meanTrend),
             size = 4, fill = "red", shape = 23) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  scale_y_continuous(limits = c(-6, 6),
                     breaks = seq(from = -6, to = 6, by = 2)) +
  my.ggplot.theme +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  labs(y = expression(paste("Trend in Phenology (d ", y^-1, ")")))

droppedCovs <- ggarrange(calcofiDropTL, nhlDropTL,
                         calcofiDropFreq, nhlDropFreq,
                         calcofiDropAmp, nhlDropAmp,
                         calcofiDropFishing, nhlDropFishing,
                         labels = c("(a)", "(b)", "(c)", "(d)",
                                    "(e)", "(f)", "(g)", "(h)"),
                         font.label = list(size = 20),
                         # hjust = -1,
                         vjust = 1.1,
                         ncol = 2, nrow = 4, align = "hv")

ggsave(filename = "droppedCovs.png",
       plot = droppedCovs,
       path = "~/Desktop",
       height = 13, width = 12, units = "in")


# (x) NH Line MonthMax GAM ========================================================
## Figure description: GAMs of NH Line MonthMax against correlation CTa x time and phenological trend. For SI

## Read in data needed to make this figure
nhlNumCovs <- read.csv('/Volumes/petrik-lab/kchen/datasets/nhl-numreg.csv')
calcofiNumCovs <- read.csv('/Volumes/petrik-lab/kchen/datasets/calcofi-numreg.csv')

## fit a GAM to NHL monthMax x trend. use bs = 'cc' to connect December to January, and set k = 4
G1 <- gam(trend_dy ~ s(monthMax, fx = T, k = 4, bs = "cc"), data = nhlNumCovs)
G1_p <- predict_gam(G1)

nhl.trendGam <-
  ggplot(data = G1_p,
         aes(x = monthMax, y = fit)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  geom_smooth_ci(ci_alpha = 0.4,
                 linewidth = 1.5) +
  geom_point(data = nhlNumCovs,
             aes(x = monthMax, y = trend_dy),
             size = 2.5) +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  scale_y_continuous(limits = c(-7.5, 2.5)) +
  my.ggplot.theme +
  labs(x = "Month of Maximum Larval Abundance",
       y = expression(paste("Trend in Phenology (d ", y^-1, ")")),
       title = "")


## fit a GAM to NHL monthMax x cor. use bs = 'cc' to connect December to January, and set k = 4
G2 <- gam(cor ~ s(monthMax, fx = T, k = 4, bs = "cc"), data = nhlNumCovs)
G2_p <- predict_gam(G2)

nhl.corGam <-
  ggplot(data = G2_p,
         aes(x = monthMax, y = fit)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  geom_smooth_ci(ci_alpha = 0.4,
                 linewidth = 1.5) +
  geom_point(data = nhlNumCovs,
             aes(x = monthMax, y = cor),
             size = 2.5) +
  scale_x_continuous(breaks = seq(1, 12, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  my.ggplot.theme +
  labs(x = "Month of Maximum Larval Abundance",
       y = "Correlation between CTa and time (r)",
       title = "NH Line Assemblage")


## plot together
GAMs <- ggarrange(nhl.corGam, nhl.trendGam,
                  labels = c("(a)", "(b)"),
                  font.label = list(size = 20),
                  # hjust = -1,
                  vjust = 1.25,
                  ncol = 2, nrow = 1, align = "hv")

ggsave(filename = "monthMaxGAMs.png",
       plot = GAMs,
       path = "~/Desktop",
       height = 7, width = 10, units = "in")