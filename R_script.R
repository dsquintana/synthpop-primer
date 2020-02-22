# Synthetic datasets: A non-technical primer for the biobehavioral sciences
# Author: Daniel S. Quintana

# Correspondence to Daniel S. Quintana, NORMENT KG Jebsen Centre for Psychosis Research,
# University of Oslo
# Email: daniel.quintana@medisin.uio.no

# A number of analyses are repeated for different datasets. On first instance of an analysis, the commands are commented

# Occasionally the following error will be encountered: 
# "Error in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : polygon edge not found"
# This is a known intermittent error with Rstudio https://github.com/tidyverse/ggplot2/issues/2252
# If this occurs, confirm you are using the most recent Rstudio version
# Alternatively, simply try re-running the analysis, as it should work again on a 2nd or 3rd attempt

# Load required packages

# This is a function that will check to see if CRAN packages are installed.
# If they are not, they will be installed.
# After checking, they will be loaded into the R session
# Source: https://gist.github.com/stevenworthington/3178163

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("synthpop", "tidyverse", "cowplot", "car",
              "simstudy", "mice", "StatMeasures")
ipak(packages)

# The following is a function from the "faux" package https://debruine.github.io/faux/index.html
# This had to be function had to be entered manually, as faux is not yet available on CRAN
# This function simulates correlated variables
# This function was pasted from version 0.0.0.9018

rnorm_pre <- function (x, mu = 0, sd = 1, r = 0, empirical = FALSE) 
{
  if (!is.vector(x)) 
    stop("x must be a vector")
  if (!is.numeric(x)) 
    stop("x must be numeric")
  if (length(x) < 3) 
    stop("x must have length > 2")
  n <- length(x)
  if (!empirical) {
    sample_params <- sample_from_pop(n, mu, sd, r)
    mu <- sample_params$mu
    sd <- sample_params$sd
    r <- sample_params$r
  }
  y <- stats::rnorm(n)
  z <- r * scale(x)[, 1] + sqrt(1 - r^2) * scale(stats::resid(stats::lm(y ~ 
                                                                          x)))[, 1]
  yresult <- mu + sd * z
  return(yresult)
}

################################

## Manuscript example 1: Oxytocin and sprituality

ot_dat <- read_csv("ot_dat.csv") # Loads data

ot_dat <- ot_dat %>%
  rename(
    OT_condition = OT_COND,
    rel_affiliation = rel_aff_cat,
    spirituality = spi_1_L
  )  # Renames the variables for easier figure interpretation

## Figure 1a 

ot_sim <- syn(ot_dat, seed = 1337) # Creates synthetic data

ot_com <- compare(
  ot_sim, # The synthetic dataset
  ot_dat, # The original dataset
  vars = c("OT_condition", "Age",
           "spirituality", "rel_affiliation"), # The variables for comparison
  print.coef = TRUE, # Print tables of estimates for original and synthetic data
  ncol = 4, # The number of columns in the plot
  breaks = 16, # Gaps between columns 
  stat = "counts", # Present the raw counts for each variable
  cols = c("#62B6CB", "#1B4965") # Setting the colours in the plot
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

fig_1a # Prints figure

#####

## Supplementary Figures 

# To generate figures, see full R script on the project's OSF page https://osf.io/z524n/ 
# This section of the analysis was removed from the Rstudio server instance due to loading constraints

#####

## Check for replicated unique units 

ru <- replicated.uniques(ot_sim, ot_dat)
ru # This prints the number of unique cases

#####

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
  s_lm, # Results from the synthetic linear model
  ot_dat, # The original dataset
  lwd = 1.5, # The type of line in the plot
  lty = 1, # The width of line in the plot
  point.size = 4, # The size of the symbols used in the plot
  lcol = c("#62B6CB", "#1B4965") # Set the colours
) # A comparison of the linear models

t_test_com

fig_1b <- t_test_com$ci.plot # Extract the plot from this object


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

#####

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

#####

## Ancova

anc = car::Anova(aov(spirituality ~ OT_condition + rel_affiliation, 
                     data = ot_dat)) # Perform ANCOVA
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

#####

## Construct figure 1

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

## Prepare data for sharing

ot_synthetic_label <- sdc(ot_sim, ot_dat, 
                          label = "FAKE_DATA") # Adds a "FAKE_DATA" label

ot_synthetic_dat <- ot_synthetic_label$syn # Extracts the synthetic data to a dataframe for sharing

## Manuscript example 2: Sociosexuality and self-rated attractiveness

## Original data source: https://osf.io/6bk3w/

socio_dat <- read_csv("socio.csv") # Import data

socio_dat <- socio_dat %>% drop_na() # Drop NAs

socio_dat <-  socio_dat %>% filter(sex %in% c("male", "female", "intersex")) 

socio_dat_s <- syn(socio_dat, seed = 122) # Create synthetic dataset

fig_2a <- compare(
  socio_dat_s,
  socio_dat,
  breaks = 12,
  ncol = 7,
  nrow = 2,
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

fig_2a

# Models

socio_lm <- lm(behavior2 ~ 1 +
                 sra + age + lab,
               data = socio_dat) # linear model

socio_dat_sum <- summary(socio_lm) # Result from linear model
socio_dat_sum

socio_lm_syn <- lm.synds(behavior2 ~ 1 + sra + age + lab,
                         data = socio_dat_s) # Equivalent linear model in synthetic data


socio_lm_syn_s <- summary(socio_lm_syn) # Results from linear model in synthetic data
socio_lm_syn_s

fig_2b <- compare(
  socio_lm_syn,
  socio_dat,
  breaks = 12,
  ncol = 7,
  nrow = 2,
  cols = c("#62B6CB", "#1B4965")
) # Compare datasets

fig_2b

fig_2b <-  fig_2b$ci.plot # Extract plot from the "fig_3b" object

fig_2b <- fig_2b + ggtitle("") +
  theme(axis.text.y = element_blank())  # Remove title and y-axis text

fig_2b <- fig_2b + theme_half_open() +
  background_grid() # Add theme

fig_2b <- fig_2b +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL, name = "Coefficient") # Remove x-axis text

fig_2b <- fig_2b +
  annotate("text",
           x = 3,
           y = 13.8,
           label = "SRA") +
  annotate("text",
           x = 2,
           y = 27.5,
           label = "Age") +
  annotate("text",
           x = 1,
           y = 1.1,
           label = "Location")
fig_2b

## Detect replicated individuals and prepare synthetic dataset for sharing

dim(socio_dat_s$syn) # Rows and columns before removal of replicated uniques

socio_dat_s_sdc <- sdc(socio_dat_s, socio_dat, 
                       label = "FAKE_DATA", 
                       rm.replicated.uniques = TRUE) # Remove replicated uniques and add FAKE label

dim(socio_dat_s_sdc$syn) # Number of rows and columns AFTER removal of replicated uniques

socio_synthetic_dat <- socio_dat_s_sdc$syn # Extracts the synthetic data (replicated uniques removed) to a dataframe for sharing

# Regression model with uniques excluded

socio_lm_syn_ue <- lm.synds(behavior2 ~ 1 + sra + age + lab,
                            data = socio_dat_s_sdc) # Equivalent linear model in synthetic data

socio_lm_syn_ue_s <- summary(socio_lm_syn_ue) # Results from linear model in synthetic data
socio_lm_syn_ue_s


fig_2c <- compare(
  socio_lm_syn_ue,
  socio_dat,
  breaks = 12,
  ncol = 7,
  nrow = 2,
  cols = c("#62B6CB", "#1B4965")
) # Compare datasets

fig_2c

fig_2c <-  fig_2c$ci.plot # Extract plot from the "fig_3b" object

fig_2c <- fig_2c + ggtitle("") +
  theme(axis.text.y = element_blank())  # Remove title and y-axis text

fig_2c <- fig_2c + theme_half_open() +
  background_grid() # Add theme

fig_2c <- fig_2c +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL, name = "Coefficient") # Remove x-axis text

fig_2c <- fig_2c +
  annotate("text",
           x = 3,
           y = 13.8,
           label = "SRA") +
  annotate("text",
           x = 2,
           y = 27.5,
           label = "Age") +
  annotate("text",
           x = 1,
           y = 1.1,
           label = "Location")
fig_2c

# Create figure 2 plot
# First create regression model panels

fig2bc <- plot_grid(
  fig_2b,
  NULL,
  fig_2c,
  labels = c('B', '', 'C'),
  ncol = 1,
  rel_heights = c(1, 0.01, 1.)
) 

fig2bc

fig2 <- plot_grid(
  fig_2a,
  NULL,
  fig2bc,
  labels = c('A', '', ''),
  ncol = 3,
  rel_widths = c(5, 0.15, 1.5)
) # Combine plots (some space was added in between the plots)

fig2

## Manuscript example 3: Heart rate variability (HRV) and fitness in a series of simulated datasets

## Simulated datasets have been labelled by codes
# A1: 40 cases, normal distribution, none missing
# A2: 40 cases, normal distribution, 5% missing
# A3: 40 cases, normal distribution, 20% missing
# A4: 40 cases, low skew, none missing
# A5: 40 cases, low skew, 5% missing
# A6: 40 cases, low skew, 20% missing
# A7: 40 cases, high skew, none missing
# A8: 40 cases, high skew, 5% missing
# A9: 40 cases, high skew, 20% missing

# B1: 100 cases, normal distribution, none missing
# B2: 100 cases, normal distribution, 5% missing
# B3: 100 cases, normal distribution, 20% missing
# B4: 100 cases, low skew, none missing
# B5: 100 cases, low skew, 5% missing
# B6: 100 cases, low skew, 20% missing
# B7: 100 cases, high skew, none missing
# B8: 100 cases, high skew, 5% missing
# B9: 100 cases, high skew, 20% missing

# C1: 10000 cases, normal distribution, none missing
# C2: 10000 cases, normal distribution, 5% missing
# C3: 10000 cases, normal distribution, 20% missing
# C4: 10000 cases, low skew, none missing
# C5: 10000 cases, low skew, 5% missing
# C6: 10000 cases, low skew, 20% missing
# C7: 10000 cases, high skew, none missing
# C8: 10000 cases, high skew, 5% missing
# C9: 10000 cases, high skew, 20% missing

#### Analysis A1 #####

# A1: 40 cases, normal distribution, none missing

# Create the simulated data using the "simstudy" package

set.seed(1337) # Set seed so that the data is reproducible 

A1_dat <- defData(varname = "nr", dist = "nonrandom", 
                  formula=20, id = "idnum") # Create ID
A1_dat <- defData(A1_dat, varname="heart_rate", dist="normal", 
                  formula = "nr + 55", 
                  variance = 5, link = "identity") # Create normally distributed heart rate data
A1_dat <- defData(A1_dat, varname="weight", dist="normal", 
                  formula = "nr + 65", 
                  variance = 6, link = "identity") # Create normally distributed weight data
A1_dat <- defData(A1_dat, varname="fitness", dist="normal", 
                  formula = "nr + 2", 
                  variance = 2, link = "identity") # Create normally distributed fitness data
A1_dat <- defData(A1_dat, varname = "hrv", dist = "normal" , 
                  formula="nr + 10", 
                  variance = 4, link = "identity") # Create normally distributed HRV data

knitr::kable(A1_dat) # A summary table of the variables

A1_dat <- genData(40, A1_dat) # Generates the dataset with 40 cases

outliers(A1_dat$hrv)  

# All values greater than 75th percentile value + 1.5 times the IQR or lesser than 25th percentile value - 1.5 times the IQR, are tagged as outliers.

A1_dat$fitness <- rnorm_pre(A1_dat$hrv, 
                            mu = 22, sd = 1, r = 0.3, 
                            empirical = TRUE) # Transforms the fitness data so that it's correlated with HRV (r = .3)

cor.test(A1_dat$fitness, A1_dat$hrv) # Confirms the correlation

hist(A1_dat$fitness) # Plots a histogram
mean(A1_dat$fitness) # Prints the mean

A1_dat <- select(A1_dat, heart_rate, weight, 
                 fitness, hrv) # Creates a new dataset with the variables of interest

A1_dat_sim <- syn(A1_dat, seed = 1337) # Creates synthetic data

A1_com <- compare(
  A1_dat_sim,
  A1_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)

A1_com

fig_S1a <- A1_com$plots # Extracts plots from the "A1_com" object

fig_S1a <- fig_S1a +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S1a <- fig_S1a +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("No skew and no missing data") # Names  title

fig_S1a 

A1_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = A1_dat_sim) # Linear model equivalent in synthetic data
summary(A1_com_lm)

A1_comp <- compare(
  A1_com_lm,
  A1_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

A1_comp

fig_a1 <- A1_comp$ci.plot

fig_a1 <- fig_a1 + ggtitle("")

fig_a1 <- fig_a1 + labs(x = "", y = "")

fig_a1 <- fig_a1 + theme_half_open() +
  background_grid()

fig_a1 <- fig_a1 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_a1 <- fig_a1 + ggtitle("No skew and no missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))

fig_a1

#### A2 #####

# A2: 40 cases, normal distribution, 5% missing

set.seed(1337)
A2_dat <- ampute(A1_dat, prop = 0.05, bycases = FALSE)
A2_dat <-A2_dat$amp

outliers(A2_dat$hrv)

A2_dat_sim <- syn(A2_dat, seed = 1337) # Creates synthetic data

A2_com <- compare(
  A2_dat_sim,
  A2_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
A2_com

fig_S2a <- A2_com$plots # Extracts plots from the "ot_com" object

fig_S2a <- fig_S2a +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S2a <- fig_S2a +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("No skew and 5% missing data")

fig_S2a 


A2_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = A2_dat_sim) # Linear model equivalent in synthetic data
summary(A2_com_lm)

A2_comp <- compare(
  A2_com_lm,
  A2_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

A2_comp

fig_a2 <- A2_comp$ci.plot

fig_a2 <- fig_a2 + ggtitle("")


fig_a2 <- fig_a2 + labs(x = "", y = "")

fig_a2 <- fig_a2 + theme_half_open() +
  background_grid()

fig_a2 <- fig_a2 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_a2 <- fig_a2 + ggtitle("No skew and 5% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))

fig_a2

#######

# A3: 40 cases, normal distribution, 20% missing

set.seed(1337)
A3_dat <- ampute(A1_dat, prop = 0.20, bycases = FALSE)
A3_dat <-A3_dat$amp

outliers(A3_dat$hrv)

A3_dat_sim <- syn(A3_dat, seed = 1337) # Creates synthetic data

A3_com <- compare(
  A3_dat_sim,
  A3_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
A3_com

fig_S3a <- A3_com$plots # Extracts plots from the "ot_com" object

fig_S3a <- fig_S3a +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S3a <- fig_S3a +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("No skew and 20% missing data")

fig_S3a 

A3_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = A3_dat_sim) # Linear model equivalent in synthetic data
summary(A3_com_lm)

A3_comp <- compare(
  A3_com_lm,
  A3_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

A3_comp

fig_a3 <- A3_comp$ci.plot

fig_a3 <- fig_a3 + ggtitle("")

fig_a3 <- fig_a3 + labs(x = "", y = "")

fig_a3 <- fig_a3 + theme_half_open() +
  background_grid()

fig_a3 <- fig_a3 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_a3 <- fig_a3 + ggtitle("No skew and 20% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_a3

## A4: 40 cases, low skew, none missing

set.seed(1337)
A4_dat <- defData(varname = "nr", dist = "nonrandom", 
                  formula=20, id = "idnum")
A4_dat <- defData(A4_dat, varname="heart_rate", dist="normal", 
                  formula = "nr + 55", variance = 5, link = "identity")
A4_dat <- defData(A4_dat, varname="weight", dist="normal", 
                  formula = "nr + 65", variance = 6, link = "identity")
A4_dat <- defData(A4_dat, varname="fitness", dist="normal", 
                  formula = "nr + 2", variance = 2, link = "identity")
A4_dat <- defData(A4_dat, varname = "hrv", dist = "gamma" , 
                  formula="nr + 10", variance = 0.1, link = "identity") # Low skew

knitr::kable(A4_dat)

A4_dat <- genData(40, A4_dat)

outliers(A4_dat$hrv)

A4_dat$fitness <- rnorm_pre(A4_dat$hrv, 
                            mu = 22, sd = 1, r = 0.3, empirical = TRUE)

cor.test(A4_dat$fitness, A4_dat$hrv)

A4_dat <- select(A4_dat, heart_rate, weight, 
                 fitness, hrv)

A4_dat_sim <- syn(A4_dat, seed = 1337) # Creates synthetic data

A4_com <- compare(
  A4_dat_sim,
  A4_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)

A4_com

fig_S4a <- A4_com$plots # Extracts plots from the "ot_com" object

fig_S4a <- fig_S4a +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S4a <- fig_S4a +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("Low skew and no missing data")

fig_S4a 

A4_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = A4_dat_sim) # Linear model equivalent in synthetic data
summary(A4_com_lm)

A4_comp <- compare(
  A4_com_lm,
  A4_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

A4_comp

fig_a4 <- A4_comp$ci.plot

fig_a4 <- fig_a4 + ggtitle("")

fig_a4 <- fig_a4 + labs(x = "", y = "")

fig_a4 <- fig_a4 + theme_half_open() +
  background_grid()

fig_a4 <- fig_a4 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_a4 <- fig_a4 + ggtitle("Low skew and no missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_a4

# A5: 40 cases, low skew, 5% missing

set.seed(1337)
A5_dat <- ampute(A4_dat, prop = 0.05, bycases = FALSE)
A5_dat <- A5_dat$amp

outliers(A5_dat$hrv)

A5_dat_sim <- syn(A5_dat, seed = 1337) # Creates synthetic data

A5_com <- compare(
  A5_dat_sim,
  A5_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
A5_com

fig_S5a <- A5_com$plots # Extracts plots from the "ot_com" object

fig_S5a <- fig_S5a +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S5a <- fig_S5a +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("Low skew and 5% missing data")

fig_S5a 

A5_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = A5_dat_sim) # Linear model equivalent in synthetic data
summary(A5_com_lm)

A5_comp <- compare(
  A5_com_lm,
  A5_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

A5_comp

fig_a5 <- A5_comp$ci.plot

fig_a5 <- fig_a5 + ggtitle("")

fig_a5 <- fig_a5 + labs(x = "", y = "")

fig_a5 <- fig_a5 + theme_half_open() +
  background_grid()

fig_a5 <- fig_a5 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_a5 <- fig_a5 + ggtitle("Low skew and 5% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_a5

# A6: 40 cases, low skew, 20% missing

set.seed(1337)
A6_dat <- ampute(A4_dat, prop = 0.20, bycases = FALSE)
A6_dat <-A6_dat$amp

outliers(A6_dat$hrv)

A6_dat_sim <- syn(A6_dat, seed = 1337) # Creates synthetic data

A6_com <- compare(
  A6_dat_sim,
  A6_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
A6_com

fig_S6a <- A6_com$plots # Extracts plots from the "ot_com" object

fig_S6a <- fig_S6a +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S6a <- fig_S6a +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("Low skew and 20% missing data")

fig_S6a 

A6_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = A6_dat_sim) # Linear model equivalent in synthetic data
summary(A6_com_lm)

A6_comp <- compare(
  A6_com_lm,
  A6_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

A6_comp

fig_a6 <- A6_comp$ci.plot

fig_a6 <- fig_a6 + ggtitle("")

fig_a6 <- fig_a6 + labs(x = "", y = "")

fig_a6 <- fig_a6 + theme_half_open() +
  background_grid()

fig_a6 <- fig_a6 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_a6 <- fig_a6 + ggtitle("Low skew and 20% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_a6

# A7: 40 cases, high skew, none missing

set.seed(1337)
A7_dat <- defData(varname = "nr", dist = "nonrandom", 
                  formula=20, id = "idnum")
A7_dat <- defData(A7_dat, varname="heart_rate", dist="normal", 
                  formula = "nr + 55", variance = 5, link = "identity")
A7_dat <- defData(A7_dat, varname="weight", dist="normal", 
                  formula = "nr + 65", variance = 6, link = "identity")
A7_dat <- defData(A7_dat, varname="fitness", dist="normal", 
                  formula = "nr + 2", variance = 2, link = "identity")
A7_dat <- defData(A7_dat, varname = "hrv", dist = "gamma" , 
                  formula="nr + 10", variance = 0.5, link = "identity") # High skew

knitr::kable(A4_dat)

A7_dat <- genData(40, A7_dat)

outliers(A7_dat$hrv)

A7_dat$fitness <- rnorm_pre(A7_dat$hrv, 
                            mu = 22, sd = 1, r = 0.3, empirical = TRUE)

cor.test(A7_dat$fitness, A7_dat$hrv)

A7_dat <- select(A7_dat, heart_rate, weight, 
                 fitness, hrv)

A7_dat_sim <- syn(A7_dat, seed = 1337) # Creates synthetic data

A7_com <- compare(
  A7_dat_sim,
  A7_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
A7_com

fig_S7a <- A7_com$plots # Extracts plots from the "ot_com" object

fig_S7a <- fig_S7a +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S7a <- fig_S7a +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("High skew and no missing data")
fig_S7a

A7_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = A7_dat_sim) # Linear model equivalent in synthetic data
summary(A7_com_lm)

A7_comp <- compare(
  A7_com_lm,
  A7_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

A7_comp

fig_a7 <- A7_comp$ci.plot

fig_a7 <- fig_a7 + ggtitle("")

fig_a7 <- fig_a7 + labs(x = "", y = "")

fig_a7 <- fig_a7 + theme_half_open() +
  background_grid()

fig_a7 <- fig_a7 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_a7 <- fig_a7 + ggtitle("High skew and no missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))

fig_a7

# A8: 40 cases, high skew, 5% missing

set.seed(1337)
A8_dat <- ampute(A7_dat, prop = 0.05, bycases = FALSE)
A8_dat <- A8_dat$amp

outliers(A8_dat$hrv)

A8_dat_sim <- syn(A8_dat, seed = 1337) # Creates synthetic data

A8_com <- compare(
  A8_dat_sim,
  A8_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
A8_com

fig_S8a <- A8_com$plots # Extracts plots from the "ot_com" object

fig_S8a <- fig_S8a +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S8a <- fig_S8a +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("High skew and 5% missing data")
fig_S8a

A8_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = A8_dat_sim) # Linear model equivalent in synthetic data
summary(A8_com_lm)

A8_comp <- compare(
  A8_com_lm,
  A8_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

A8_comp

fig_a8 <- A8_comp$ci.plot

fig_a8 <- fig_a8 + ggtitle("")

fig_a8 <- fig_a8 + labs(x = "", y = "")

fig_a8 <- fig_a8 + theme_half_open() +
  background_grid()

fig_a8 <- fig_a8 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_a8 <- fig_a8 + ggtitle("High skew and 5% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_a8

# A9: 40 cases, high skew, 20% missing

set.seed(1337)
A9_dat <- ampute(A7_dat, prop = 0.2, bycases = FALSE)
A9_dat <- A9_dat$amp

outliers(A9_dat$hrv)

A9_dat_sim <- syn(A9_dat, seed = 1337) # Creates synthetic data

A9_com <- compare(
  A9_dat_sim,
  A9_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
A9_com

fig_S9a <- A9_com$plots # Extracts plots from the "ot_com" object

fig_S9a <- fig_S9a +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S9a <- fig_S9a +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("High skew and 20% missing data")

fig_S9a

A9_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = A9_dat_sim) # Linear model equivalent in synthetic data
summary(A9_com_lm)

A9_comp <- compare(
  A9_com_lm,
  A9_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

A9_comp

fig_a9 <- A9_comp$ci.plot

fig_a9 <- fig_a9 + ggtitle("")

fig_a9 <- fig_a9 + labs(x = "", y = "")

fig_a9 <- fig_a9 + theme_half_open() +
  background_grid()

fig_a9 <- fig_a9 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_a9 <- fig_a9 + ggtitle("High skew and 20% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_a9

figure_S4 <- plot_grid(
  fig_S1a, fig_S2a, fig_S3a, fig_S4a,
    fig_S5a, fig_S6a, fig_S7a, fig_S8a, fig_S9a,
  nrow = 9
) # Create a plot grid visualising general utility


# To save the plot at a readable size, remove the "uncomment" (i.e., remove the leading #) in the next line
# ggsave("figure_S4.pdf", width = 20, height = 22, units = "in")

figure_S7 <- plot_grid(
  fig_a1, fig_a2, fig_a3, 
  fig_a4, fig_a5, fig_a6, 
  fig_a7, fig_a8, fig_a9,
  nrow = 3
)

# To save the plot at a readable size, remove the "uncomment" (i.e., remove the leading #) in the next line
# ggsave("figureS3.pdf", width = 7, height = 4.5, units = "in")

#### B1 #####

# B1: 100 cases, normal distribution, none missing

set.seed(1337)
B1_dat <- defData(varname = "nr", dist = "nonrandom", 
                  formula=20, id = "idnum")
B1_dat <- defData(B1_dat, varname="heart_rate", dist="normal", 
                  formula = "nr + 55", variance = 5, link = "identity")
B1_dat <- defData(B1_dat, varname="weight", dist="normal", 
                  formula = "nr + 65", variance = 6, link = "identity")
B1_dat <- defData(B1_dat, varname="fitness", dist="normal", 
                  formula = "nr + 2", variance = 2, link = "identity")
B1_dat <- defData(B1_dat, varname = "hrv", dist = "normal" , 
                  formula="nr + 10", variance = 4, link = "identity") # Normal dist

knitr::kable(B1_dat)

B1_dat <- genData(100, B1_dat)

outliers(B1_dat$hrv)

B1_dat$fitness <- rnorm_pre(B1_dat$hrv, 
                            mu = 22, sd = 1, r = 0.3, empirical = TRUE)

cor.test(B1_dat$fitness, B1_dat$hrv)

B1_dat <- select(B1_dat, heart_rate, weight, 
                 fitness, hrv)

B1_dat_sim <- syn(B1_dat, seed = 1337) # Creates synthetic data

B1_com <- compare(
  B1_dat_sim,
  B1_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)

B1_com

fig_S1b <- B1_com$plots # Extracts plots from the "ot_com" object

fig_S1b <- fig_S1b +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S1b <- fig_S1b +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("No skew and no missing data")

fig_S1b

B1_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = B1_dat_sim) # Linear model equivalent in synthetic data
summary(B1_com_lm)

B1_comp <- compare(
  B1_com_lm,
  B1_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

B1_comp

fig_b1 <- B1_comp$ci.plot

fig_b1 <- fig_b1 + ggtitle("")

fig_b1 <- fig_b1 + labs(x = "", y = "")

fig_b1 <- fig_b1 + theme_half_open() +
  background_grid()

fig_b1 <- fig_b1 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_b1 <- fig_b1 + ggtitle("No skew and no missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))

fig_b1

#### B2 #####

# B2: 100 cases, normal distribution, 5% missing

set.seed(1337)
B2_dat <- ampute(B1_dat, prop = 0.05, bycases = FALSE)
B2_dat <-B2_dat$amp

outliers(B2_dat$hrv)

B2_dat_sim <- syn(B2_dat, seed = 1337) # Creates synthetic data

B2_com <- compare(
  B2_dat_sim,
  B2_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)

B2_com

fig_S2b <- B2_com$plots # Extracts plots from the "ot_com" object

fig_S2b <- fig_S2b +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S2b <- fig_S2b +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("No skew and no 5% missing data")

fig_S2b

B2_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = B2_dat_sim) # Linear model equivalent in synthetic data
summary(B2_com_lm)

B2_comp <- compare(
  B2_com_lm,
  B2_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

B2_comp

fig_b2 <- B2_comp$ci.plot

fig_b2 <- fig_b2 + ggtitle("")

fig_b2 <- fig_b2 + labs(x = "", y = "")

fig_b2 <- fig_b2 + theme_half_open() +
  background_grid()

fig_b2 <- fig_b2 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_b2 <- fig_b2 + ggtitle("No skew and 5% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))

fig_b2

#######
# B3: 100 cases, normal distribution, 20% missing

set.seed(1337)
B3_dat <- ampute(B1_dat, prop = 0.20, bycases = FALSE)
B3_dat <-B3_dat$amp

outliers(B3_dat$hrv)

B3_dat_sim <- syn(B3_dat, seed = 1337) # Creates synthetic data

B3_com <- compare(
  B3_dat_sim,
  B3_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
B3_com

fig_S3b <- B3_com$plots # Extracts plots from the "ot_com" object

fig_S3b <- fig_S3b +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S3b <- fig_S3b +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("No skew and 20% missing data")

fig_S3b

B3_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = B3_dat_sim) # Linear model equivalent in synthetic data
summary(B3_com_lm)

B3_comp <- compare(
  B3_com_lm,
  B3_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

B3_comp

fig_b3 <- B3_comp$ci.plot

fig_b3 <- fig_b3 + ggtitle("")

fig_b3 <- fig_b3 + labs(x = "", y = "")

fig_b3 <- fig_b3 + theme_half_open() +
  background_grid()

fig_b3 <- fig_b3 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_b3 <- fig_b3 + ggtitle("No skew and 20% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))

fig_b3

## B4: 100 cases, low skew, none missing

set.seed(1337)
B4_dat <- defData(varname = "nr", dist = "nonrandom", 
                  formula=20, id = "idnum")
B4_dat <- defData(B4_dat, varname="heart_rate", dist="normal", 
                  formula = "nr + 55", variance = 5, link = "identity")
B4_dat <- defData(B4_dat, varname="weight", dist="normal", 
                  formula = "nr + 65", variance = 6, link = "identity")
B4_dat <- defData(B4_dat, varname="fitness", dist="normal", 
                  formula = "nr + 2", variance = 2, link = "identity")
B4_dat <- defData(B4_dat, varname = "hrv", dist = "gamma" , 
                  formula="nr + 10", variance = 0.1, link = "identity") # Low skew

knitr::kable(B4_dat)

B4_dat <- genData(100, B4_dat)

outliers(B4_dat$hrv)

B4_dat$fitness <- rnorm_pre(B4_dat$hrv, 
                            mu = 22, sd = 1, r = 0.3, empirical = TRUE)

cor.test(B4_dat$fitness, B4_dat$hrv)

B4_dat <- select(B4_dat, heart_rate, weight, 
                 fitness, hrv)

B4_dat_sim <- syn(B4_dat, seed = 1337) # Creates synthetic data

B4_com <- compare(
  B4_dat_sim,
  B4_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
B4_com

fig_S4b <- B4_com$plots # Extracts plots from the "ot_com" object

fig_S4b <- fig_S4b +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S4b <- fig_S4b +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("Low skew and no missing data")

fig_S4b

B4_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = B4_dat_sim) # Linear model equivalent in synthetic data
summary(B4_com_lm)

B4_comp <- compare(
  B4_com_lm,
  B4_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

B4_comp

fig_b4 <- B4_comp$ci.plot

fig_b4 <- fig_b4 + ggtitle("")

fig_b4 <- fig_b4 + labs(x = "", y = "")

fig_b4 <- fig_b4 + theme_half_open() +
  background_grid()

fig_b4 <- fig_b4 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_b4 <- fig_b4 + ggtitle("Low skew and no missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_b4

# B5: 100 cases, low skew, 5% missing

set.seed(1337)
B5_dat <- ampute(B4_dat, prop = 0.05, bycases = FALSE)
B5_dat <- B5_dat$amp

outliers(B5_dat$hrv)

B5_dat_sim <- syn(B5_dat, seed = 1337) # Creates synthetic data

B5_com <- compare(
  B5_dat_sim,
  B5_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
B5_com

fig_S5b <- B5_com$plots # Extracts plots from the "ot_com" object

fig_S5b <- fig_S5b +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S5b <- fig_S5b +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("Low skew and 5% missing data")

fig_S5b

B5_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = B5_dat_sim) # Linear model equivalent in synthetic data
summary(B5_com_lm)

B5_comp <- compare(
  B5_com_lm,
  B5_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

B5_comp

fig_b5 <- B5_comp$ci.plot

fig_b5 <- fig_b5 + ggtitle("")

fig_b5 <- fig_b5 + labs(x = "", y = "")

fig_b5 <- fig_b5 + theme_half_open() +
  background_grid()

fig_b5 <- fig_b5 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_b5 <- fig_b5 + ggtitle("Low skew and 5% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_b5

# B6: 100 cases, low skew, 20% missing

set.seed(1337)
B6_dat <- ampute(B4_dat, prop = 0.20, bycases = FALSE)
B6_dat <-B6_dat$amp

outliers(B6_dat$hrv)

B6_dat_sim <- syn(B6_dat, seed = 1337) # Creates synthetic data

B6_com <- compare(
  B6_dat_sim,
  B6_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
B6_com

fig_S6b <- B6_com$plots # Extracts plots from the "ot_com" object

fig_S6b <- fig_S6b +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S6b <- fig_S6b +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("Low skew and 20% missing data")

fig_S6b

B6_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = B6_dat_sim) # Linear model equivalent in synthetic data
summary(B6_com_lm)

B6_comp <- compare(
  B6_com_lm,
  B6_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

B6_comp

fig_b6 <- B6_comp$ci.plot

fig_b6 <- fig_b6 + ggtitle("")

fig_b6 <- fig_b6 + labs(x = "", y = "")

fig_b6 <- fig_b6 + theme_half_open() +
  background_grid()

fig_b6 <- fig_b6 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_b6 <- fig_b6 + ggtitle("Low skew and 20% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_b6

# B7: 100 cases, high skew, none missing

set.seed(1337)
B7_dat <- defData(varname = "nr", dist = "nonrandom", 
                  formula=20, id = "idnum")
B7_dat <- defData(B7_dat, varname="heart_rate", dist="normal", 
                  formula = "nr + 55", variance = 5, link = "identity")
B7_dat <- defData(B7_dat, varname="weight", dist="normal", 
                  formula = "nr + 65", variance = 6, link = "identity")
B7_dat <- defData(B7_dat, varname="fitness", dist="normal", 
                  formula = "nr + 2", variance = 2, link = "identity")
B7_dat <- defData(B7_dat, varname = "hrv", dist = "gamma" , 
                  formula="nr + 10", variance = 0.5, link = "identity") # High skew

knitr::kable(B7_dat)

B7_dat <- genData(100, B7_dat)

outliers(B7_dat$hrv)

B7_dat$fitness <- rnorm_pre(B7_dat$hrv, 
                            mu = 22, sd = 1, r = 0.3, empirical = TRUE)

cor.test(B7_dat$fitness, B7_dat$hrv)

B7_dat <- select(B7_dat, heart_rate, weight, 
                 fitness, hrv)

B7_dat_sim <- syn(B7_dat, seed = 1337) # Creates synthetic data

B7_com <- compare(
  B7_dat_sim,
  B7_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
B7_com

fig_S7b <- B7_com$plots # Extracts plots from the "ot_com" object

fig_S7b <- fig_S7b +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S7b <- fig_S7b +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("High skew and no missing data")

fig_S7b

B7_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = B7_dat_sim) # Linear model equivalent in synthetic data
summary(B7_com_lm)

B7_comp <- compare(
  B7_com_lm,
  B7_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

B7_comp

fig_b7 <- B7_comp$ci.plot

fig_b7 <- fig_b7 + ggtitle("")

fig_b7 <- fig_b7 + labs(x = "", y = "")

fig_b7 <- fig_b7 + theme_half_open() +
  background_grid()

fig_b7 <- fig_b7 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_b7 <- fig_b7 + ggtitle("High skew and no missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_b7

# B8: 100 cases, high skew, 5% missing

set.seed(1337)
B8_dat <- ampute(B7_dat, prop = 0.05, bycases = FALSE)
B8_dat <- B8_dat$amp

outliers(B8_dat$hrv)

B8_dat_sim <- syn(B8_dat, seed = 1337) # Creates synthetic data

B8_com <- compare(
  B8_dat_sim,
  B8_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
B8_com

fig_S8b <- B8_com$plots # Extracts plots from the "ot_com" object

fig_S8b <- fig_S8b +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S8b <- fig_S8b +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("High skew and 5% missing data")

fig_S8b

B8_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = B8_dat_sim) # Linear model equivalent in synthetic data
summary(B8_com_lm)

B8_comp <- compare(
  B8_com_lm,
  B8_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

B8_comp

fig_b8 <- B8_comp$ci.plot

fig_b8 <- fig_b8 + ggtitle("")

fig_b8 <- fig_b8 + labs(x = "", y = "")

fig_b8 <- fig_b8 + theme_half_open() +
  background_grid()

fig_b8 <- fig_b8 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_b8 <- fig_b8 + ggtitle("High skew and 5% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_b8

# B9: 100 cases, high skew, 20% missing

set.seed(1337)
B9_dat <- ampute(B7_dat, prop = 0.2, bycases = FALSE)
B9_dat <- B9_dat$amp

outliers(B9_dat$hrv)

B9_dat_sim <- syn(B9_dat, seed = 1337) # Creates synthetic data

B9_com <- compare(
  B9_dat_sim,
  B9_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
B9_com

fig_S9b <- B9_com$plots # Extracts plots from the "ot_com" object

fig_S9b <- fig_S9b +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S9b <- fig_S9b +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("High skew and 20% missing data")

fig_S9b

B9_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = B9_dat_sim) # Linear model equivalent in synthetic data
summary(B9_com_lm)

B9_comp <- compare(
  B9_com_lm,
  B9_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

B9_comp

fig_b9 <- B9_comp$ci.plot

fig_b9 <- fig_b9 + ggtitle("")

fig_b9 <- fig_b9 + labs(x = "", y = "")

fig_b9 <- fig_b9 + theme_half_open() +
  background_grid()

fig_b9 <- fig_b9 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_b9 <- fig_b9 + ggtitle("High skew and 20% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_b9

figure_S5 <- plot_grid(
  fig_S1b, fig_S2b, fig_S3b, fig_S4b,
  fig_S5b, fig_S6b, fig_S7b, fig_S8b, fig_S9b,
  nrow = 9
) # Create a plot grid visualising general utility

# To save the plot at a readable size, remove the "uncomment" (i.e., remove the leading #) in the next line
# ggsave("figure_S5.pdf", width = 20, height = 22, units = "in")


figure_3 <- plot_grid(
  fig_b1, fig_b2, fig_b3, 
  fig_b4, fig_b5, fig_b6, 
  fig_b7, fig_b8, fig_b9,
  nrow = 3
)

# To save the plot at a readable size, remove the "uncomment" (i.e., remove the leading #) in the next line
# ggsave("figure_3.pdf", width = 7, height = 4.5, units = "in")

#######################
#### C1 #####

# C1: 10000 cases, normal distribution, none missing

set.seed(1337)
C1_dat <- defData(varname = "nr", dist = "nonrandom", 
                  formula=20, id = "idnum")
C1_dat <- defData(C1_dat, varname="heart_rate", dist="normal", 
                  formula = "nr + 55", variance = 5, link = "identity")
C1_dat <- defData(C1_dat, varname="weight", dist="normal", 
                  formula = "nr + 65", variance = 6, link = "identity")
C1_dat <- defData(C1_dat, varname="fitness", dist="normal", 
                  formula = "nr + 2", variance = 2, link = "identity")
C1_dat <- defData(C1_dat, varname = "hrv", dist = "normal" , 
                  formula="nr + 10", variance = 4, link = "identity") # Normal dist

knitr::kable(C1_dat)

C1_dat <- genData(10000, C1_dat)

outliers(C1_dat$hrv)

C1_dat$fitness <- rnorm_pre(C1_dat$hrv, 
                            mu = 22, sd = 1, r = 0.3, empirical = TRUE)

cor.test(C1_dat$fitness, C1_dat$hrv)

C1_dat <- select(C1_dat, heart_rate, weight, 
                 fitness, hrv)

C1_dat_sim <- syn(C1_dat, seed = 1337) # Creates synthetic data

C1_com <- compare(
  C1_dat_sim,
  C1_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
C1_com

fig_S1c <- C1_com$plots # Extracts plots from the "ot_com" object

fig_S1c <- fig_S1c +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S1c <- fig_S1c +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("No skew and no missing data")

fig_S1c

C1_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = C1_dat_sim) # Linear model equivalent in synthetic data
summary(C1_com_lm)

C1_comp <- compare(
  C1_com_lm,
  C1_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

C1_comp

fig_c1 <- C1_comp$ci.plot

fig_c1 <- fig_c1 + ggtitle("")

fig_c1 <- fig_c1 + labs(x = "", y = "")

fig_c1 <- fig_c1 + theme_half_open() +
  background_grid()

fig_c1 <- fig_c1 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_c1 <- fig_c1 + ggtitle("No skew and no missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_c1

#### C2 #####

# C2: 10000 cases, normal distribution, 5% missing

set.seed(1337)
C2_dat <- ampute(C1_dat, prop = 0.05, bycases = FALSE)
C2_dat <-C2_dat$amp

outliers(C2_dat$hrv)

C2_dat_sim <- syn(C2_dat, seed = 1337) # Creates synthetic data

C2_com <- compare(
  C2_dat_sim,
  C2_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
C2_com

fig_S2c <- C2_com$plots # Extracts plots from the "ot_com" object

fig_S2c <- fig_S2c +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S2c <- fig_S2c +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("No skew and 5% missing data")

fig_S2c

C2_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = C2_dat_sim) # Linear model equivalent in synthetic data
summary(C2_com_lm)

C2_comp <- compare(
  C2_com_lm,
  C2_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

C2_comp

fig_c2 <- C2_comp$ci.plot

fig_c2 <- fig_c2 + ggtitle("")

fig_c2 <- fig_c2 + labs(x = "", y = "")

fig_c2 <- fig_c2 + theme_half_open() +
  background_grid()

fig_c2 <- fig_c2 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_c2 <- fig_c2 + ggtitle("No skew and no 5% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_c2

#######
# C3: 10000 cases, normal distribution, 20% missing

set.seed(1337)
C3_dat <- ampute(C1_dat, prop = 0.20, bycases = FALSE)
C3_dat <- C3_dat$amp

outliers(C3_dat$hrv)

C3_dat_sim <- syn(C3_dat, seed = 1337) # Creates synthetic data

C3_com <- compare(
  C3_dat_sim,
  C3_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
C3_com

fig_S3c <- C3_com$plots # Extracts plots from the "ot_com" object

fig_S3c <- fig_S3c +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S3c <- fig_S3c +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("No skew and 20% missing data")

fig_S3c

C3_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = C3_dat_sim) # Linear model equivalent in synthetic data
summary(C3_com_lm)

C3_comp <- compare(
  C3_com_lm,
  C3_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

C3_comp

fig_c3 <- C3_comp$ci.plot

fig_c3 <- fig_c3 + ggtitle("")

fig_c3 <- fig_c3 + labs(x = "", y = "")

fig_c3 <- fig_c3 + theme_half_open() +
  background_grid()

fig_c3 <- fig_c3 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_c3 <- fig_c3 + ggtitle("No skew and no 20% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_c3

# C4: 10000 cases, low skew, none missing

set.seed(1337)
C4_dat <- defData(varname = "nr", dist = "nonrandom", 
                  formula=20, id = "idnum")
C4_dat <- defData(C4_dat, varname="heart_rate", dist="normal", 
                  formula = "nr + 55", variance = 5, link = "identity")
C4_dat <- defData(C4_dat, varname="weight", dist="normal", 
                  formula = "nr + 65", variance = 6, link = "identity")
C4_dat <- defData(C4_dat, varname="fitness", dist="normal", 
                  formula = "nr + 2", variance = 2, link = "identity")
C4_dat <- defData(C4_dat, varname = "hrv", dist = "gamma" , 
                  formula="nr + 10", variance = 0.1, link = "identity") # Low skew

knitr::kable(C4_dat)

C4_dat <- genData(10000, C4_dat)

outliers(C4_dat$hrv)

C4_dat$fitness <- rnorm_pre(C4_dat$hrv, 
                            mu = 22, sd = 1, r = 0.3, empirical = TRUE)

cor.test(C4_dat$fitness, C4_dat$hrv)

C4_dat <- select(C4_dat, heart_rate, weight, 
                 fitness, hrv)

C4_dat_sim <- syn(C4_dat, seed = 1337) # Creates synthetic data

C4_com <- compare(
  C4_dat_sim,
  C4_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
C4_com

fig_S4c <- C4_com$plots # Extracts plots from the "ot_com" object

fig_S4c <- fig_S4c +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S4c <- fig_S4c +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("Low skew and no missing data")

fig_S4c

C4_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = C4_dat_sim) # Linear model equivalent in synthetic data
summary(C4_com_lm)

C4_comp <- compare(
  C4_com_lm,
  C4_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

C4_comp

fig_c4 <- C4_comp$ci.plot

fig_c4 <- fig_c4 + ggtitle("")

fig_c4 <- fig_c4 + labs(x = "", y = "")

fig_c4 <- fig_c4 + theme_half_open() +
  background_grid()

fig_c4 <- fig_c4 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_c4 <- fig_c4 + ggtitle("Low skew and no missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_c4

# C5: 10000 cases, low skew, 5% missing

set.seed(1337)
C5_dat <- ampute(C4_dat, prop = 0.05, bycases = FALSE)
C5_dat <- C5_dat$amp

outliers(C5_dat$hrv)

C5_dat_sim <- syn(C5_dat, seed = 1337) # Creates synthetic data

C5_com <- compare(
  C5_dat_sim,
  C5_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
C5_com

fig_S5c <- C5_com$plots # Extracts plots from the "ot_com" object

fig_S5c <- fig_S5c +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S5c <- fig_S5c +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("Low skew and 5% missing data")

fig_S5c

C5_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = C5_dat_sim) # Linear model equivalent in synthetic data
summary(C5_com_lm)

C5_comp <- compare(
  C5_com_lm,
  C5_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

C5_comp

fig_c5 <- C5_comp$ci.plot

fig_c5 <- fig_c5 + ggtitle("")

fig_c5 <- fig_c5 + labs(x = "", y = "")

fig_c5 <- fig_c5 + theme_half_open() +
  background_grid()

fig_c5 <- fig_c5 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_c5 <- fig_c5 + ggtitle("Low skew and 5% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_c5

# C6: 10000 cases, low skew, 20% missing

set.seed(1337)
C6_dat <- ampute(C4_dat, prop = 0.20, bycases = FALSE)
C6_dat <-C6_dat$amp

outliers(C6_dat$hrv)

C6_dat_sim <- syn(C6_dat, seed = 1337) # Creates synthetic data

C6_com <- compare(
  C6_dat_sim,
  C6_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
C6_com

fig_S6c <- C6_com$plots # Extracts plots from the "ot_com" object

fig_S6c <- fig_S6c +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S6c <- fig_S6c +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("Low skew and 20% missing data")

fig_S6c

C6_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = C6_dat_sim) # Linear model equivalent in synthetic data
summary(C6_com_lm)

C6_comp <- compare(
  C6_com_lm,
  C6_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

C6_comp

fig_c6 <- C6_comp$ci.plot

fig_c6 <- fig_c6 + ggtitle("")

fig_c6 <- fig_c6 + labs(x = "", y = "")

fig_c6 <- fig_c6 + theme_half_open() +
  background_grid()

fig_c6 <- fig_c6 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_c6 <- fig_c6 + ggtitle("Low skew and 20% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_c6

# C7: 10000 cases, high skew, none missing

set.seed(1337)
C7_dat <- defData(varname = "nr", dist = "nonrandom", 
                  formula=20, id = "idnum")
C7_dat <- defData(C7_dat, varname="heart_rate", dist="normal", 
                  formula = "nr + 55", variance = 5, link = "identity")
C7_dat <- defData(C7_dat, varname="weight", dist="normal", 
                  formula = "nr + 65", variance = 6, link = "identity")
C7_dat <- defData(C7_dat, varname="fitness", dist="normal", 
                  formula = "nr + 2", variance = 2, link = "identity")
C7_dat <- defData(C7_dat, varname = "hrv", dist = "gamma" , 
                  formula="nr + 10", variance = 0.5, link = "identity") # High skew

knitr::kable(C7_dat)

C7_dat <- genData(10000, C7_dat)

outliers(C7_dat$hrv)

C7_dat$fitness <- rnorm_pre(C7_dat$hrv, 
                            mu = 22, sd = 1, r = 0.3, empirical = TRUE)

cor.test(C7_dat$fitness, C7_dat$hrv)

C7_dat <- select(C7_dat, heart_rate, weight, 
                 fitness, hrv)

C7_dat_sim <- syn(C7_dat, seed = 1337) # Creates synthetic data

C7_com <- compare(
  C7_dat_sim,
  C7_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
C7_com

fig_S7c <- C7_com$plots # Extracts plots from the "ot_com" object

fig_S7c <- fig_S7c +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S7c <- fig_S7c +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("High skew and no missing data")

fig_S7c

C7_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = C7_dat_sim) # Linear model equivalent in synthetic data
summary(C7_com_lm)

C7_comp <- compare(
  C7_com_lm,
  C7_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

C7_comp

fig_c7 <- C7_comp$ci.plot

fig_c7 <- fig_c7 + ggtitle("")

fig_c7 <- fig_c7 + labs(x = "", y = "")

fig_c7 <- fig_c7 + theme_half_open() +
  background_grid()

fig_c7 <- fig_c7 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_c7 <- fig_c7 + ggtitle("High skew and no missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_c7

# C8: 10000 cases, high skew, 5% missing

set.seed(1337)
C8_dat <- ampute(C7_dat, prop = 0.05, bycases = FALSE)
C8_dat <- C8_dat$amp

outliers(C8_dat$hrv)

C8_dat_sim <- syn(C8_dat, seed = 1337) # Creates synthetic data

C8_com <- compare(
  C8_dat_sim,
  C8_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
C8_com

fig_S8c <- C8_com$plots # Extracts plots from the "ot_com" object

fig_S8c <- fig_S8c +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S8c <- fig_S8c +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("High skew and 5% missing data")

fig_S8c

C8_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = C8_dat_sim) # Linear model equivalent in synthetic data
summary(C8_com_lm)

C8_comp <- compare(
  C8_com_lm,
  C8_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

C8_comp

fig_c8 <- C8_comp$ci.plot

fig_c8 <- fig_c8 + ggtitle("")

fig_c8 <- fig_c8 + labs(x = "", y = "")

fig_c8 <- fig_c8 + theme_half_open() +
  background_grid()

fig_c8 <- fig_c8 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_c8 <- fig_c8 + ggtitle("High skew and 5% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_c8

# C9: 10000 cases, high skew, 20% missing

set.seed(1337)
C9_dat <- ampute(C7_dat, prop = 0.2, bycases = FALSE)
C9_dat <- C9_dat$amp

outliers(C9_dat$hrv)

C9_dat_sim <- syn(C9_dat, seed = 1337) # Creates synthetic data

C9_com <- compare(
  C9_dat_sim,
  C9_dat,
  vars = c("heart_rate", "weight", 
           "fitness", "hrv"),
  print.coef = TRUE,
  ncol = 4,
  breaks = 16,
  stat = "counts",
  cols = c("#62B6CB", "#1B4965")
)
C9_com

fig_S9c <- C9_com$plots # Extracts plots from the "ot_com" object

fig_S9c <- fig_S9c +
  scale_y_continuous(expand = c(0, 0)) + # Forces the y-axis to start at zero
  theme_minimal_hgrid(12) # Applies a theme from the 'cowplot' package

fig_S9c <- fig_S9c +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        # Adjusts x-axis tick labels
        axis.title.x = element_blank()) + # Removes x-axis title
  labs(fill = "Dataset") +
  ggtitle("High skew and 20% missing data")

fig_S9c

C9_com_lm <- lm.synds(scale(fitness) ~ 1 +
                        scale(hrv),
                      data = C9_dat_sim) # Linear model equivalent in synthetic data
summary(C9_com_lm)

C9_comp <- compare(
  C9_com_lm,
  C9_dat,
  lwd = 1.5,
  lty = 1,
  point.size = 4,
  lcol = c("#62B6CB", "#1B4965")
)

C9_comp

fig_c9 <- C9_comp$ci.plot

fig_c9 <- fig_c9 + ggtitle("")

fig_c9 <- fig_c9 + labs(x = "", y = "")

fig_c9 <- fig_c9 + theme_half_open() +
  background_grid()

fig_c9 <- fig_c9 +
  theme(axis.text.y = element_blank()) +
  scale_x_discrete(breaks = NULL)

fig_c9 <- fig_c9 + ggtitle("High skew and 5% missing data") +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold"))
fig_c9

figure_S6 <- plot_grid(
  fig_S1c, fig_S2c, fig_S3c, fig_S4c,
  fig_S5c, fig_S6c, fig_S7c, fig_S8c, fig_S9c,
  nrow = 9
) # Create a plot grid visualising general utility

# To save the plot at a readable size, remove the "uncomment" (i.e., remove the leading #) in the next line
# ggsave("figure_S6.pdf", width = 20, height = 22, units = "in")


figure_S8 <- plot_grid(
  fig_c1, fig_c2, fig_c3, 
  fig_c4, fig_c5, fig_c6, 
  fig_c7, fig_c8, fig_c9,
  nrow = 3
)

# To save the plot at a readable size, remove the "uncomment" (i.e., remove the leading #) in the next line
# ggsave("figure_S8.pdf", width = 7, height = 4.5, units = "in")

### END OF SCRIPT ###