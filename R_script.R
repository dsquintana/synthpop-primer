# Synthetic datasets: A non-technical primer for the biobehavioral sciences
# Author: Daniel S. Quintana

# Correspondence to Daniel S. Quintana, NORMENT KG Jebsen Centre for Psychosis Research,
# University of Oslo
# Email: daniel.quintana@medisin.uio.no

# Required datasets to reproduce these results can be downloaded at: https://osf.io/z524n/

# Load required packages

# This is a function that will check to see if packages are installed.
# If they are not, they will be installed.
# After checking, they will be loaded into the R session
# Source: https://gist.github.com/stevenworthington/3178163

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("synthpop", "tidyverse", "cowplot", "car", "ggstatsplot")
ipak(packages)


################################

### Manuscript example 1: Oxytocin and sprituality

ot_dat <- read_csv("ot_dat.csv") # Loads data

ot_dat <- ot_dat %>%
  rename(
    OT_condition = OT_COND,
    rel_affiliation = rel_aff_cat,
    spirituality = spi_1_L
  )  # Renames the variables for easier figure interpretation

## Figure 1a ##

ot_sim <- syn(ot_dat, seed = 1337) # Creates synthetic data

ot_com <- compare(
  ot_sim,
  ot_dat,
  vars = c("OT_condition", "Age",
           "spirituality", "rel_affiliation"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
) # Visual comparison of original and synthetic datasets

fig_1a <- ot_com$plots # Extracts plots from the "ot_com" object

fig_1a <- fig_1a +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_1a <- fig_1a +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") # Renames legend title

fig_1a

#####

## Supplementary Figure 1 ##

ot_dat2 <-
  as.data.frame(ot_dat) # Create new dataframe for multicompare

s_fig_1 <- multi.compare(ot_sim, ot_dat2,
                          var = "OT_condition", by = "spirituality")

s_fig_1 <- s_fig_1 +
  labs(fill = "Nasal spray \n condition") + # Renames legend
  labs(y = "Counts") + # Relabels y-axis 
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

s_fig_1

## Supplementary Figure 2 ##

s_fig_2 <- multi.compare(ot_sim, ot_dat2,
                          var = "OT_condition", 
                          by = "rel_affiliation",
                          cont.type = "box")

s_fig_2 <- s_fig_2 +
  labs(fill = "Nasal spray \n condition") + # Renames legend
  labs(y = "Count") + # Relabels y-axis 
  scale_y_continuous(expand = c(0, 0)) + # Forces y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

s_fig_2

## Supplementary Figure 3 ##

# Create new merged dataset to create figure

dat_m <- read_csv("ot_dat.csv") 
ot_sim_m <- syn(dat_m, seed = 1337)
ot_sim_df <- ot_sim_m$syn
dat_m_df  <- cbind(data = "observed", dat_m)
ot_sim_dff <- cbind(data = "synthetic", ot_sim_df)
merge_dat <- bind_rows(dat_m_df, ot_sim_dff)

s_fig_3 <- ggstatsplot::grouped_ggscatterstats(
  data = merge_dat,
  x = spi_1_L,
  y = Age,
  conf.level = 0.95,
  k = 2, # no. of decimal places in the results
  xlab = "Spirituality",
  bf.message = FALSE,
  grouping.var = data, # grouping variable
  title.prefix = "Dataset",
  marginal.type = "density",
  ggtheme = ggplot2::theme_minimal(),
  messages = FALSE,
  nrow = 2,
  title.text = ""
)
s_fig_3

#######

## Check for replicated unique units ##

ru <- replicated.uniques(ot_sim, ot_dat)

ru

#######

## t-test

a = t.test(spirituality ~ OT_condition,
           data = ot_dat,
           var.equal = FALSE) # Welch's t-test
a 

a1 = lm(ot_dat$spirituality ~ 1 +
          ot_dat$OT_condition) # Linear model equivalent of above t-test

summary(a1) # Confirming results are the same as the t-test
confint(a1)

s_lm <- lm.synds(spirituality ~ 1 +
                   OT_condition,
                 data = ot_sim) # Linear model equivalent in synthetic data

syn <- summary(s_lm) # Synthetic linear model results
syn

t_test_com <- compare(
  s_lm,
  ot_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
) # A comparison of the linear models

t_test_com

fig_1b <- t_test_com$ci.plot


fig_1b <- fig_1b + ggtitle("") +
  theme(axis.text.y = element_blank())

fig_1b <- fig_1b + theme_half_open() +
  background_grid()

fig_1b <- fig_1b +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_1b <- fig_1b +
  annotate("text",
           x = 1,
           y = -1,
           label = "Nasal spray condition")
fig_1b

## Correlation

a_cor = cor.test(ot_dat$Age, ot_dat$spirituality,
                 method = "pearson") # Calculate correlation
a_cor

b_cor = lm(scale(ot_dat$Age) ~ 1 +
             scale(ot_dat$spirituality)) # Linear model equivalent of correlation

summary(b_cor) # Linear model results
confint(b_cor) # Print confidence intervals for linear model coefficients

s_cor <- lm.synds(scale(Age) ~ 1 +
                    scale(spirituality),
                  data = ot_sim) # Linear model equivalent in synthetic data

syn_cor <- summary(s_cor) # Results of linear model equivalent in synthetic data
syn_cor

cor_com <- compare(
  s_cor,
  ot_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

cor_com

fig_1c <- cor_com$ci.plot # Extract plot from "cor_com" object

fig_1c <- fig_1c + ggtitle("") + # Remove title
  theme(axis.text.y = element_blank()) # Remove y-axis text

fig_1c <- fig_1c + theme_half_open() +
  background_grid() # Apply new theme

fig_1c <- fig_1c +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL, name = "") # Remove x-axis text

fig_1c <- fig_1c +
  annotate("text",
           x = 1,
           y = 0.7,
           label = "Spirituality") # Add label to plot
fig_1c

#### Ancova


anc = car::Anova(aov(spirituality ~ OT_condition + rel_affiliation, 
                     data = ot_dat))
anc

anc_lm <- lm(spirituality ~ 1 +
               OT_condition + rel_affiliation,
             data = ot_dat) # Linear model equivalent of ANOVA

summary(anc_lm) # Results from linear model

# Testing for main effect of group to confirm equivalancy

null_rel = lm(spirituality ~ 1 +
                rel_affiliation,
              data = ot_dat) # Null model without OT condition

result_rel = anova(null_rel, anc_lm) # Comparison of null and full model
result_rel # Comparison of null and full model, yielding the same F statistic and p-value


s_ancova <-
  lm.synds(spirituality ~ 1 + OT_condition + rel_affiliation,
           data = ot_sim) # Linear model equivalent in synthetic data

syn_anc <- summary(s_ancova) # Results from synthetic linear model
syn_anc

anc_com <- compare(
  s_ancova,
  ot_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  plot.intercept = FALSE,
  lcol = c("#62B6CB", "#1B4965")
) # Comparison of linear and synthetic model

anc_com

anc_plot <-  anc_com$ci.plot # Extract plot from "anc_com" object

fig_1d <- anc_plot + ggtitle("") +
  theme(axis.text.y = element_blank()) # Remove title and y-axis text

fig_1d <- fig_1d + theme_half_open() +
  background_grid() # Apply theme

fig_1d <- fig_1d +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL, name = "") # Remove x-axis text

fig_1d <- fig_1d +
  annotate("text",
           x = 1.02,
           y = -6.5,
           label = "Nasal spray condition") +
  annotate("text",
           x = 2.02,
           y = -1.8,
           label = "Religious affiliation") # Add labels
fig_1d

# Construct figure 1

p1_top <- plot_grid(fig_1a,
                    labels = c('A'),
                    ncol = 1,
                    label_size = 12) # Create top panel


p1_bottom <- plot_grid(
  fig_1b + theme(legend.position = "none"),
  fig_1c + theme(legend.position = "none"),
  fig_1d + theme(legend.position = "none"),
  labels = c('B', 'C', 'D'),
  ncol = 3,
  rel_widths = c(1, 1, 1),
  label_size = 12
) # Create top panel with stripped legends

legend <- get_legend(fig_1c + theme(legend.box.margin = margin(0, 0, 0, 12))) # Extract legend and create some space

p1_bottom <- plot_grid(p1_bottom,
                       NULL,
                       legend,
                       NULL,
                       ncol = 4,
                       rel_widths = c(3, 0.1, .2, .1)) # Add legend and some more space

p1_bottom <- plot_grid(NULL,
                       p1_bottom,
                       NULL,
                       ncol = 4,
                       rel_widths = c(0.2, 2.2, 0.2)) # Adding a little more space

fig1 <- plot_grid(p1_top, p1_bottom,
                  nrow = 2) # Putting it all together

fig1 # Print at 14 x 6 inches for same dimensions as manuscript

######

### Manuscript example 1: Oxytocin concentrations and theory of mind performance

# Original data source: https://data.mendeley.com/datasets/h3f6ywpd5t/1

b_dat <- read_csv("blood.csv") # Import data

vars_b <- c("EQ", "RMET", "OT", "Sex")
b_dat <- b_dat[, vars_b] # Select variables of interest

b_dat_s <- syn(b_dat, seed = 738) # Create synthetic dataset

fig_2a <- compare(
  b_dat_s,
  b_dat,
  stat = "counts",
  breaks = 12,
  ncol = 2,
  cols = c("#62B6CB", "#1B4965")
) # Compare datasets

fig_2a <- fig_2a$plots # Extract plots from "Fig_2a" object

fig_2a <- fig_2a +
  scale_y_continuous(expand = c(0, 0)) + # Force y-axis to start at zero
  theme_minimal_hgrid(12) # Apply theme

fig_2a <- fig_2a +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank()) +
  labs(fill = "Dataset")

#####
# Check for replicated unique values

ru_b <- replicated.uniques(b_dat_s, b_dat)
ru_b
####

# RMET

rmet_m <- lm(RMET ~ 1 +
               OT + Sex,
             data = b_dat) # linear model

rmet_m_s <- summary(rmet_m) # Result from linear model
rmet_m_s

rmet_m_syn <- lm.synds(RMET ~ 1 + OT + Sex,
                       data = b_dat_s) # Equivalent linear model in synthetic data

rmet_m_syn_s <- summary(rmet_m_syn) # Results from linear model in synthetic data
rmet_m_syn_s

fig_2b <- compare(
  rmet_m_syn,
  b_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
) # Comparison of linear models
fig_2b

fig_2b <-  fig_2b$ci.plot # Extract plot from the "fig_2b" object

fig_2b <- fig_2b + ggtitle("") +
  theme(axis.text.y = element_blank())  # Remove title and y-axis text

fig_2b <- fig_2b + theme_half_open() +
  background_grid() # Add theme

fig_2b <- fig_2b +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL, name = "Coefficient") # Remove x-axis text

fig_2b <- fig_2b +
  annotate("text",
           x = 2,
           y = -1,
           label = "Oxytocin concentration") +
  annotate("text",
           x = 1,
           y = -0.35,
           label = "Sex") # Add labels

# Plot grid

fig2 <- plot_grid(
  fig_2a,
  NULL,
  fig_2b,
  labels = c('A', '', 'B'),
  ncol = 3,
  rel_widths = c(2, 0.15, 1.5)
) # Combine plots (some space was added in between the plots)

fig2 # Print at 14 x 5 inches for same dimensions as manuscript


