---
title: "Code_to_reproduce_Figures"
author: "Cynthia Becker"
date: "2023-03-14"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 3
    number_sections: true
    theme: united
editor_options: 
  markdown: 
    wrap: 73
---



# Setup

## Install necessary packages


```r
# For Maps
library(sf); packageVersion("sf")
library(ggspatial); packageVersion("ggspatial")

# For data wrangling
library(plyr); packageVersion("plyr")
library(dplyr); packageVersion("dplyr")
library(tidyverse); packageVersion("tidyverse")
library(broom); packageVersion("broom")
library(readxl); packageVersion("readxl")

# For microbiome and statistical analyses
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(fantaxtic); packageVersion("fantaxtic")
library(rstatix); packageVersion("rstatix")
library(FactoMineR); packageVersion("FactoMineR")
library(missMDA); packageVersion("missMDA")
library(breakaway); packageVersion("breakaway")
library(corncob); packageVersion("corncob")
library(speedyseq); packageVersion("speedyseq")
library(randomForest); packageVersion("randomForest")
library(lme4); packageVersion("lme4") # linear mixed effects modeling
library(lmerTest); packageVersion("lmerTest") # linear mixed effects modeling with p-value calculation
library(sjPlot) #may be useful for plotting model results from lmer with plot_model(nh4_model, type = "est")
library(corrplot); packageVersion("corrplot")


# For visualization
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw()) # get rid of the gray background
library(knitr); packageVersion("knitr")
library(rcartocolor); packageVersion("rcartocolor")
library(ggpubr); packageVersion("ggpubr")
library(factoextra); packageVersion("factoextra")
library(ggpmisc); packageVersion("ggpmisc")
library(pals); packageVersion("pals")
library(patchwork); packageVersion("patchwork")
```

## Read in prepped data

The shapefiles for making the map can be downloaded at the following
links. To use this code chunk for your own analysis, change the paths
defined in `st_read()`

-   United States and Territories
    <https://earthworks.stanford.edu/catalog/stanford-vt021tk4894>

-   US Virgin Islands and Puerto Rico Habitats
    <https://products.coastalscience.noaa.gov/collections/benthic/e95usvi_pr/>

-   US National Parks
    <https://public-nps.opendata.arcgis.com/datasets/nps::nps-land-resources-division-boundary-and-tract-data-service/explore?layer=2&location=0.239390%2C-12.497900%2C2.00>


```r
# path to the shapefiles. Download the shapefiles linked above and change the path to use this code
usa <- st_read("/Users/cynthiabecker/Documents/Apprill_lab/USVI_Projects/SiteMap/stanford-vt021tk4894-shapefile/", "vt021tk4894")
sttstj <- st_read("/Users/cynthiabecker/Documents/Apprill_lab/USVI_Projects/SiteMap/stsj_fin")
nps <- st_read("/Users/cynthiabecker/Documents/Apprill_lab/USVI_Projects/SiteMap/NPS_-_Land_Resources_Division_Boundary_and_Tract_Data_Service")
```

Environmental Data


```r
envdata <- read_xlsx("data/Env_Data_USVI.xlsx", na = "NA")
reefs <- read.table("data/St.John_sites_latlon.txt", sep = "\t", header = TRUE)
precip <- read.csv("data/EastEndVI_Precip.csv", sep = ",")
temperaturepath <- "data/StJohn_temperature.xlsx"
```

Taxonomic microbiome data (16S rRNA gene abundances). Import the data
that has been filtered for contaminants and low abundance reads. In this
case we defined that as taxa with fewer than 5 counts across the whole
dataset.


```r
## Read in data - use the low abundance filtered dataset
asv <- read.table("data/ASV_lowabund_filt.txt",sep="\t",header=TRUE, row.names=1)
taxa <- as.matrix(read.table("data/taxonomy_lowabund_filt.txt", sep="\t", header=TRUE, row.names=1))
samples <- read.table("data/metadata_lowabund_filt.txt",sep="\t",header=T,row.names=1)

colnames(asv) <- str_replace_all(colnames(asv), pattern = "[X]", "")

# Create phyloseq object
ASV = otu_table(asv, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
META = sample_data(samples)

ps <- phyloseq(ASV, TAX, META)
ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 14009 taxa and 375 samples ]:
## sample_data() Sample Data:        [ 375 samples by 42 sample variables ]:
## tax_table()   Taxonomy Table:     [ 14009 taxa by 7 taxonomic ranks ]:
## taxa are rows
```

```r
# What are the total number of sequences per sample (mean ± sd)
mean(colSums(ASV))
```

```
## [1] 53473.09
```

```r
sd(colSums(ASV))
```

```
## [1] 9272.262
```

```r
# What are the total number of ASVs per sample (mean ± sd)
mean(colSums(ASV != 0))
```

```
## [1] 306.5947
```

```r
sd(colSums(ASV != 0))
```

```
## [1] 59.19492
```

Note that using filtered data for alpha diveristy can lead to incorrect
estimations. Because I am interested in analyzing the richness of the
communities, I will use the full dataset of \~33k ASVs that have not been
filtered to remove low-abundance ASVs. The "filt" extension in this case
only refers to the fact that these tables have been filtered of DNA
extraction contaminatns (via decontam) and filtered of a few samples that
had abnormally low sequencing depth.


```r
asv.a <- read.table("data/ASV_filt.txt",sep="\t",header=TRUE, row.names=1)
taxa.a <- as.matrix(read.table("data/taxonomy_filt.txt", sep="\t", header=TRUE, row.names=1))
samples.a <- read.table("data/metadata_filt.txt",sep="\t",header=T,row.names=1)

colnames(asv.a) <- str_replace_all(colnames(asv.a), pattern = "[X]", "")

# add seasonality information ("wetdryseason" specifically) to the metadata_filt
samples.prev <- samples %>% 
  mutate(sampleID = rownames(.))

seasoninfo <- samples.a %>% 
  mutate(sampleID = rownames(.)) %>% 
  dplyr::select(sampleID) %>% 
  left_join(samples.prev, by = "sampleID")

rownames(seasoninfo) <- seasoninfo$sampleID

ASV.a = otu_table(asv.a, taxa_are_rows = TRUE)
TAX.a = tax_table(taxa.a)
META.a = sample_data(seasoninfo)

ps.a <- phyloseq(ASV.a, TAX.a, META.a)
ps.a
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 33362 taxa and 375 samples ]:
## sample_data() Sample Data:        [ 375 samples by 43 sample variables ]:
## tax_table()   Taxonomy Table:     [ 33362 taxa by 7 taxonomic ranks ]:
## taxa are rows
```

```r
# What are the total number of sequences per sample (mean ± sd)
mean(colSums(ASV.a))
```

```
## [1] 53606.37
```

```r
sd(colSums(ASV.a))
```

```
## [1] 9276.83
```

```r
# What are the total number of ASVs per sample (mean ± sd)
mean(colSums(ASV.a != 0))
```

```
## [1] 346.6347
```

```r
sd(colSums(ASV.a != 0))
```

```
## [1] 72.75927
```

## Generate functions used in analysis


```r
# function for reading in multiple excel sheets
multiplesheets <- function(fname) {
   
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
    
  # assigning names to data frames
  names(data_frame) <- sheets
    
  # print data frame
  print(data_frame)
}


# create a function that will return the values of interest from the Differential abundance test
# THis is for when there are two variables tested or compared to the control
extractmods1 <- function(model) {
  result <- data.frame("Estimate" = model$coefficients[2:3, 1], 
                        "Std.Error" = model$coefficients[2:3, 2], 
                        "p" = model$coefficients[2:3, 4])
  return(result)
}

# create a function that will return the values of interest
# This one is for when there is only one comparison
extractmods2 <- function(model) {
  result <- data.frame("Estimate" = model$coefficients[2, 1], 
                        "Std.Error" = model$coefficients[2, 2], 
                        "p" = model$coefficients[2, 4])
  return(result)
}
```

# Figure 1

Eight reefs with varying benthic cover in St. John, US Virgin Islands,
were sampled from 2016 through 2022 surrounding two major disturbances.


```r
# site metadata
reefs$Site <- factor(reefs$Site, levels = c("Dittlif", "Cocoloba", "Joels Shoal", "Europa", "Yawzi", "Tektite", "Booby Rock", "Ram Head"))

# Make color palette from colorblind friendly palettes in rcartocolor Safe palette
reefpalette <- c("#88CCEE", "gray60", "#DDCC77", "#117733", 
                 "#332288", "#AA4499", "#44AA99", "tan3")

# Filter to only layers you want
usvi <- usa %>% filter(state == "United States Virgin Islands")
reefZONE <- sttstj %>%
  filter(ZONE %in% c("Forereef", "Reef Crest", "Backreef"))
vicoral <- nps %>% filter(PARKNAME %in% c("Virgin Islands", "Virgin Islands Coral Reef"))

map <- ggplot() +
  geom_sf(data = usvi, fill = "darkgray", color = "darkgray") +
  geom_sf(data = reefZONE, fill = "pink", color = "coral") +
  geom_sf(data = vicoral, fill = NA, color = "tan4", linewidth = 1) +
  coord_sf(xlim = c(-64.8038, -64.685), ylim = c(18.2895, 18.375364), expand = FALSE) +
  geom_point(data = reefs, mapping = aes(x = Lon, y = Lat, fill = Site), pch = 21, size = 4.5) +
  scale_fill_carto_d(palette = "Safe") +
  annotation_scale(location = "bl", width_hint = 0.2) +
  theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = "white"))

coralalgae <- envdata %>%
  filter(site != "USVI Blue Water") %>%
  filter(depthtype == "benthic") %>%
  dplyr::select(Date, site, siteacronym, yearmonth, disturbance, BleachCoral:TurfAlgae) %>%
  mutate(site = factor(site, levels = c("Dittlif", "Cocoloba", "Joels Shoal", "Europa", "Yawzi", "Tektite", "Booby Rock", "Ram Head", "Sand patch"))) %>%
  filter(site != "Sand patch") %>% #no need to calculate sand patch stuff
  mutate(Algae = Ramicrusta + Macroalgae + TurfAlgae) %>% # only interested in total algae
  mutate(Coral = HardCoral + SoftCoral) %>%               # only interested in total coral
  dplyr::select(Date, site, disturbance, siteacronym, yearmonth, Algae, Coral) %>%
  mutate(disturbance = factor(disturbance, levels = c("historic", "hurricane", "disease"))) %>%
  drop_na(Algae) %>%
  gather(key = "organism", value = "percentcover", Algae:Coral)

benthos <- ggplot(coralalgae, aes(x = Date, y = percentcover)) +
  geom_vline(xintercept = as.POSIXct(as.Date("2017-09-06")), colour = "gray", linewidth = 3) + # Irma and Maria
  geom_vline(xintercept = as.POSIXct(as.Date("2020-06-01")), colour = "gray", linewidth = 3) + # Bleaching event
  geom_point(aes(shape = disturbance, color = site), data = coralalgae, size = 4) +
  geom_line(aes(color = site), linewidth = 1.2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_carto_d(palette = "Safe") +
  facet_wrap( ~ organism, scales = "free_y", ncol = 1) +
  labs(y = "Relative cover", x = "Date", color = "Reef", shape = "Disturbance")

ggarrange(map, benthos, labels = c("a.", "b."), widths = c(2,1), common.legend = TRUE, legend = "bottom")
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

<img src="figures/fig-rev-Map-1.png" width="672" />

```r
# Retrieve the linear model results for the trend line so you can report the results
HCtrend <- lm(percentcover ~ Date, data = filter(coralalgae, organism == "Coral"))
summary(HCtrend)
```

```
## 
## Call:
## lm(formula = percentcover ~ Date, data = filter(coralalgae, organism == 
##     "Coral"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.16204 -0.06786 -0.01663  0.05511  0.24768 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)   
## (Intercept)  7.939e-01  2.786e-01   2.849  0.00575 **
## Date        -3.770e-10  1.784e-10  -2.113  0.03819 * 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.09468 on 70 degrees of freedom
## Multiple R-squared:  0.05994,	Adjusted R-squared:  0.04651 
## F-statistic: 4.464 on 1 and 70 DF,  p-value: 0.03819
```

```r
HCtrend$coefficients[2]*31556926 # get the estimate in years (need to multiply the # of seconds in a year, 31556926)
```

```
##        Date 
## -0.01189656
```

```r
Algaetrend <- lm(percentcover ~ Date, data = filter(coralalgae, organism == "Algae"))
summary(Algaetrend)
```

```
## 
## Call:
## lm(formula = percentcover ~ Date, data = filter(coralalgae, organism == 
##     "Algae"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.34357 -0.06658  0.00335  0.07283  0.34459 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -1.540e+00  3.832e-01  -4.019 0.000145 ***
## Date         1.301e-09  2.454e-10   5.304 1.26e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1302 on 70 degrees of freedom
## Multiple R-squared:  0.2867,	Adjusted R-squared:  0.2765 
## F-statistic: 28.13 on 1 and 70 DF,  p-value: 1.26e-06
```

```r
Algaetrend$coefficients[2]*31556926 # get the estimate in years (need to multiply the # of seconds in a year, 31556926)
```

```
##       Date 
## 0.04107202
```

# Supplementary Figure S1 - SCTLD 


```r
sctld <- envdata %>%
  select(Date:month, SCTLDprevalence) %>%
  filter(!is.na(SCTLDprevalence)) %>%
  mutate(site = factor(site, levels = c("Dittlif", "Cocoloba", "Joels Shoal", "Europa", "Yawzi", "Tektite", "Booby Rock", "Ram Head")))

ggplot(sctld, aes(x = site, y = SCTLDprevalence*100, fill = yearmonth)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("black", "gray55")) +
  labs(x = "Reef", y = "Percent of corals with SCTLD", fill = "Year-month", title = "Figure S1. SCTLD Prevalence") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

<img src="figures/fig-rev-unnamed-chunk-2-1.png" width="672" />


# Figures 2 & 3 - Data Prep and Stats Analysis

## Prepare benthic and environmental data for analysis


```r
## Prepare data for manipulation
envdata2 <- envdata %>%       
  filter(site != "USVI Blue Water") %>%
  filter(site != "Sand patch") %>%
  mutate(TON = tn_um - (no2no3_um + nh4_um)) %>%
  mutate(no3_um = no2no3_um - no2_um) %>%
  dplyr::select(-Phaeo_ug_per_l, -no2no3_um, -tn_um, -temp_c, -salinity) %>%
  mutate(cellRatio = hbact_per_ml / (pro_per_ml + syn_per_ml + peuk_per_ml)) %>%
  mutate(Pro_fgC = pro_per_ml * 30) %>%
  mutate(Syn_fgC = syn_per_ml * 100) %>%
  mutate(Hbact_fgC = hbact_per_ml * 10) %>%
  mutate(biomassRatio = Hbact_fgC / (Pro_fgC + Syn_fgC)) %>%
  dplyr::select(-Pro_fgC, -Syn_fgC, -Hbact_fgC) %>%
  mutate(Chl_ug_per_l = ifelse(depthtype == "surface", NA, Chl_ug_per_l)) %>% #remove surface chlorophyll values
  mutate(BenthicAlgae = Ramicrusta + Macroalgae + TurfAlgae) %>% # choose true algae only, since CYAN aren't and ignore CCA b/c it is beneficial
  mutate(coral = HardCoral + SoftCoral) %>%
  mutate(AlgaeToCoral = BenthicAlgae / (HardCoral + SoftCoral)) %>%
  dplyr::select(-BleachCoral, -Other, -CYAN, -CCA, -DiseasedCoral, -Sponge, -Substrate) %>% # no need for non-coral or non-algae groups
  mutate(disturbance = factor(disturbance, levels = c("historic", "hurricane", "disease"))) %>%
  mutate(site = factor(site, levels = c("Dittlif", "Cocoloba", "Joels Shoal", "Europa", "Yawzi", "Tektite", "Booby Rock", "Ram Head"))) %>%
  mutate(disturbance = factor(disturbance, levels = c("historic", "hurricane", "disease"))) 

# Put data into long format
env.long <- envdata2 %>%
  gather(key = "datatype", value = "value", npoc_um:AlgaeToCoral) %>%
  drop_na(value)

# Remove outliers for the likely erroneously high Chl, ammonium, nitrate, and npoc, since this could cause significant relationships when there potentially are none. This way I am being more cautious. I define outliers as ±3 std deviations away from the mean. Remove them.
envdata_olrm <- env.long %>%
  filter(datatype %in% c("Chl_ug_per_l", "nh4_um", "no3_um", "npoc_um")) %>%
  group_by(datatype) %>%
  filter(!(abs(value - median(value)) > 3*sd(value)))

# Add the envdata_olrm back to the env.long dataframe
env.long.olrm <- env.long %>%
  filter(datatype != "Chl_ug_per_l") %>%
  filter(datatype != "nh4_um") %>%
  filter(datatype != "no3_um") %>%
  filter(datatype != "npoc_um") %>%
  bind_rows(envdata_olrm)
```

## Are surface and benthic values distinct for environmental (nutrient and microbial) parameters?

This analysis is important to justify whether or not to combine surface
and reef-depth (also called benthic) data in the analysis of how time,
reef, and disturbances influenced nutrient and microbial parameters.


```r
# For the benthic data. Do a shapiro-wilk test to test for normality. If normal, do an ANOVA followed by a Tukey HSD post-hoc test. 
surfacebenthic.long <- env.long.olrm %>%
  filter(datatype %in% c("po4_um",
                         "silicate_um",
                         "no2_um",
                         "pro_per_ml",
                         "syn_per_ml",
                         "peuk_per_ml",
                         "hbact_per_ml",
                         "TON",
                         "cellRatio",
                         "biomassRatio",
                         "npoc_um",
                         "nh4_um",
                         "no3_um")) %>%
  filter(value != "Inf")

#check data are normally distributed
normal <- surfacebenthic.long %>%
  group_by(datatype, depthtype) %>%
  do(tidy(shapiro.test(.$value))) %>% 
  ungroup() %>%
  mutate(normal = ifelse(p.value < 0.05, "not_normal", "yes_normal")) %>%
  dplyr::select(datatype, depthtype, normal) %>%
  spread(key = "depthtype", value = "normal")

# Check for homogeneity of variances. This is an assumption of the ANOVA test
variance <- surfacebenthic.long %>%
  group_by(datatype) %>%
  levene_test(value ~ disturbance, center = median) %>% 
  mutate(equalvariance = ifelse(p < 0.05, "no", "yes")) %>%
  dplyr::select(datatype, equalvariance)

# t-test results
ttest <- surfacebenthic.long %>%
  group_by(datatype) %>%
  t_test(value ~ depthtype) %>%
  mutate(t.test.sig = ifelse(p < 0.05, "yes", "no")) %>%
  dplyr::select(datatype, p, t.test.sig)

# Do wilcox test, which is a non-parametric alternative to a t-test
wilcox <- surfacebenthic.long %>%
  group_by(datatype) %>%
  wilcox_test(value ~ depthtype) %>%
  mutate(Wilcox.Sig = ifelse(p < 0.05, "yes", "no")) %>%
  dplyr::select(datatype, p, Wilcox.Sig)

surf_benth_envdata_results <- normal %>%
  left_join(variance, by = "datatype") %>%
  left_join(ttest, by = "datatype") %>%
  left_join(wilcox, by = "datatype") %>%
  mutate(test2use = ifelse(surface == "yes_normal" & benthic == "yes_normal" & equalvariance == "yes", "ttest", "wicox")) %>%
  mutate(Significant = ifelse(test2use == "ttest", t.test.sig, Wilcox.Sig)) %>%
  mutate(p = ifelse(test2use == "ttest", p.x, p.y)) %>%
  dplyr::select(datatype, benthic, surface, test2use, Significant, p)

kable(surf_benth_envdata_results)
```



|datatype     |benthic    |surface    |test2use |Significant |       p|
|:------------|:----------|:----------|:--------|:-----------|-------:|
|biomassRatio |not_normal |not_normal |wicox    |no          | 0.99600|
|cellRatio    |not_normal |not_normal |wicox    |no          | 0.98800|
|hbact_per_ml |yes_normal |not_normal |wicox    |yes         | 0.00747|
|nh4_um       |not_normal |not_normal |wicox    |no          | 0.79600|
|no2_um       |yes_normal |yes_normal |ttest    |no          | 0.18800|
|no3_um       |not_normal |not_normal |wicox    |yes         | 0.00478|
|npoc_um      |not_normal |not_normal |wicox    |no          | 0.21400|
|peuk_per_ml  |not_normal |not_normal |wicox    |no          | 0.42400|
|po4_um       |not_normal |not_normal |wicox    |no          | 0.33000|
|pro_per_ml   |not_normal |not_normal |wicox    |no          | 0.33100|
|silicate_um  |yes_normal |not_normal |wicox    |no          | 0.89600|
|syn_per_ml   |yes_normal |not_normal |wicox    |no          | 0.35200|
|TON          |yes_normal |yes_normal |wicox    |no          | 0.13000|

From this, we determine that most parameters are not significantly
different between surface and benthic waters. Since we were most
interested in how time, reef, and disturbances influenced nutrient and
microbial parameters, we chose to combine benthic (reef-depth) and
surface measurements for statistical analyses.

## Linear Mixed Effects Modeling (LMM) test on environmental data (benthic composition, nutrients, microbial abundances)

We hypothesized that disturbances (hurricanes, disease) impacted reef
ecosystems in the USVI, which would cause changes in cover of major coral
and algal groups comprising the benthic habitat, as well as reef water
inorganic and organic nutrients, chlorophyll, and microbial cell
abundances from both reef depth and surface waters. We also hypothesized
that sampling timepoint and reef influenced these metrics, and therefore
controlled for these effects. Toward this end, we used a linear mixed
effects modeling (LMM) approach that tested a fixed effect of disturbance
on each environmental variable while controlling for the random effects
of reef site and sampling time. Controlling for the eight reef sites was
important to account for variability in each reef as they are visually
distinct and may be impacted by different hydrological regimes and
location, which would influence the environmental variables we measured.
Controlling for time was important because the samplings were irregularly
timed year to year and there was a longer stretch of time between the
post-hurricane samplings (late 2017-2018) and the disease samplings
(2020-2022). Due to the irregular samplings, we assumed that impacts to
environmental variables may not be linear over time. To serve as a
reference, we graphed each environmental parameter over time to highlight
potential changes with three main covariates: time, reefs, and
disturbances, and placed these in the supplement.

To do the linear mixed effects modeling, I will write a for loop that
creates a model, investigates the assumptions of the model, and saves the
output of the model for each environmental variable. I commented out all
the code I used for saving the figures, but provided it in case someone
would like to re-do the analysis and investigate the assumptions and
output. To save space, I will hide the code output, which includes
numerous graphs and statistical output. Below, I will provide an example
model summary output and one set of graphs used to assess model
assumptions.


```r
# prepare a dataframe for LMM analysis. I will need to remove the disease prevalence, BenthicAlgae (total algae), coral (total coral) and cellRatio varaibles. Also remove any missing or "Infinite" values. 
env.data.foranalysis <- env.long.olrm %>%
  filter(datatype != "SCTLDprevalence") %>%
  filter(datatype != "BenthicAlgae") %>%
  filter(datatype != "coral") %>%
  filter(datatype != "cellRatio") %>%
  drop_na(value) %>%
  filter(value != "Inf")

variables <- as.vector(unique(env.data.foranalysis$datatype))

# Initialize an empty list to store results
results_list <- list()
CI_list <- list()
  
for (i in 1:length(variables)) {
  cat("Variable:", variables[i], "\n")
  
  vardata <- filter(env.data.foranalysis, datatype == variables[i])
  
  model <- lmerTest::lmer(value ~ disturbance + (1|site) + (1|yearmonth), data = vardata)
  
  cat("Summary of Model for", variables[i], "\n")
  print(summary(model))
  
  # save the predicted and raw values
  pred <- predict(model)
  res <- residuals(model)
  res.pred.data <- data.frame(pred, res, vardata)
  
  var <- variables[i]
  
  # Open a PDF file for saving plots
  #pdf(paste("figures/lmer_figs/", var, "_assumptions.pdf", sep = ""))
  
  # Investigate assumptions of equal variance of residuals along predicted values and normally distributed residuals. I will 
  par(mfrow=c(2,2))
  plot(x = pred, y = res)
  abline(h = 0)
  qqnorm(res) 
  hist(res)
  plot.new()
  text(0.5,0.5, paste(var), cex=2)
  
  # Close the PDF file
  #dev.off()
  
  # Investigate assumptions
  p1 <- ggplot(res.pred.data, aes(x = pred, y = res, color = disturbance)) +
    geom_point() +
    geom_hline(yintercept = 0) +
    facet_wrap(. ~ disturbance) +
    ggtitle(paste(var))
  
  p2 <- ggplot(res.pred.data, aes(x = pred, y = res, color = site)) +
    geom_point() +
    geom_hline(yintercept = 0) +
    facet_wrap(. ~ site) +
    ggtitle(paste(var))
  
  p3 <- ggplot(res.pred.data, aes(x = pred, y = res, color = yearmonth)) +
    geom_point() +
    geom_hline(yintercept = 0) +
    facet_wrap(. ~ yearmonth) +
    ggtitle(paste(var))
  
  #ggsave(paste("figures/lmer_figs/", var, "_assumptions_disturb.pdf", sep = ""), plot = p1)
  #ggsave(paste("figures/lmer_figs/", var, "_assumptions_site.pdf", sep = ""), plot = p2)
  #ggsave(paste("figures/lmer_figs/", var, "_assumptions_yearmonth.pdf", sep = ""), plot = p3)
  
  # Make models with and without the random effects and a null model without disturbance as fixed effect
  notime_model <- lmerTest::lmer(value ~ disturbance + (1|site), data = vardata)
  nosite_model <- lmerTest::lmer(value ~ disturbance + (1|yearmonth), data = vardata)
  null_model <- lmerTest::lmer(value ~ 1 + (1|site) + (1|yearmonth), data = vardata)
  
  # run a likelihood ratio test to evaluate if additional random/fixed effects leads to a better fit of the model
  time_LRT <- anova(model, notime_model)
  site_LRT <- anova(model, nosite_model)
  disturbance_LRT <- anova(model, null_model)
  
  # Save graphs of both the predicted and actual values just out of interest and verification model is appropriate
  predicted <- ggplot(res.pred.data, aes(x = disturbance, y = pred)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(width=.1, height=0), aes(color = disturbance), size = 3, alpha = 0.7) +
    scale_color_brewer(palette = "Dark2") 
  actual <- ggplot(res.pred.data, aes(x = disturbance, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(width=.1, height=0), aes(color = disturbance), size = 3, alpha = 0.7) +
    scale_color_brewer(palette = "Dark2") 
  
  p4 <- ggarrange(predicted, actual, ncol = 2, 
                  labels = c(paste("Predicted", var), paste("Actual", var)), 
                  common.legend = TRUE, widths = c(1,1))
  
  #ggsave(paste("figures/lmer_figs/", var, "_pred_actual.pdf", sep = ""), plot = p4)

  # Plot the estimates with confidence intervals
  fe <- fixef(model)
  ci <- confint(model)

  cife <- data.frame(ci[4:6,], fe, "disturbance" = c("historic", "hurricane", "disease"))
  colnames(cife) <- c("lower", "upper", "fixedeff", "disturbance")
  cife$disturbance <- factor(cife$disturbance, levels = c("historic", "hurricane", "disease"))

  cife_noIntercept <- cife[2:3,]
  
  # Plot the coefficients with the confidence intervals (95%). If they go through zero that indicates less confidence in the change
  p5 <- ggplot(cife_noIntercept, aes(x = disturbance, y = fixedeff, color = disturbance)) +
    geom_hline(yintercept = 0, color = "#1B9E77", size = 1) +
    geom_point(size = 5) +
    geom_errorbar(data = cife_noIntercept, aes(ymin = lower, ymax = upper), width = 0.2, size = 1.5) +
    scale_color_manual(values = c("#D95F02", "#7570B3")) +
    labs(x = "Disturbance", y = "Estimate with 95% CI relative to Historic", title = paste(var))
  
  #ggsave(paste("figures/lmer_figs/", var, "_estimate_CI.pdf", sep = ""), plot = p5, height = 3, width = 3)

  # Add results to a list
  results_list[[i]] <- list(
    variable = variables[i],
    model_summary = summary(model), 
    time_LRT_results = time_LRT, 
    site_LRT_results = site_LRT, 
    disturbance_LRT_results = disturbance_LRT, 
    CI_95 = ci
  )
  
  # Add confidence interval results to a list so I can easily make a lot of graphs after this
  CI_list[[i]] <- list(
    variable = variables[i],
    conf_int = cife_noIntercept
  )
}
```

Example of a model summary and summary figure investigating the model
assumptions of normal distribution of residuals and homogeneity of
variances of the residuals. I will show the model output for turf algae
which was significantly different. The statistical results from the
likelihood ratio test are stored in `disturbance_LRT_results` and the
model summary which includes t-tests between groups is saved in the
`model_summary`. For the graph I show, I am just showing the graphs
assessing model assumptions for the last variable (Nitrate) run in the
above for loop.


```r
# I will print the results for turf algae, which was significantly different due to disturbance (see LRT_results - likelihood ratio test results) and specificaly was different between 
results_list[[12]]$variable
```

```
## [1] "TurfAlgae"
```

```r
results_list[[12]]$disturbance_LRT_results
```

```
## Data: vardata
## Models:
## null_model: value ~ 1 + (1 | site) + (1 | yearmonth)
## model: value ~ disturbance + (1 | site) + (1 | yearmonth)
##            npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)  
## null_model    4 -126.22 -117.12 67.112  -134.22                       
## model         6 -129.53 -115.87 70.765  -141.53 7.3064  2    0.02591 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
results_list[[12]]$model_summary
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: value ~ disturbance + (1 | site) + (1 | yearmonth)
##    Data: vardata
## 
## REML criterion at convergence: -125.7
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.7232 -0.5136 -0.1426  0.4255  4.0372 
## 
## Random effects:
##  Groups    Name        Variance  Std.Dev.
##  yearmonth (Intercept) 0.0021477 0.04634 
##  site      (Intercept) 0.0002534 0.01592 
##  Residual              0.0072020 0.08486 
## Number of obs: 72, groups:  yearmonth, 9; site, 8
## 
## Fixed effects:
##                      Estimate Std. Error      df t value Pr(>|t|)  
## (Intercept)           0.06219    0.03944 4.95951   1.577   0.1762  
## disturbancehurricane  0.05918    0.04860 5.29845   1.218   0.2748  
## disturbancedisease    0.13994    0.04756 5.18140   2.942   0.0308 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) dstrbnch
## dstrbnchrrc -0.795         
## distrbncdss -0.812  0.697
```

```r
# This will print the results for nitrate, which was the last variable tested in the above for loop
par(mfrow=c(2,2))
plot(x = pred, y = res)
abline(h = 0)
qqnorm(res) 
hist(res)
plot.new()
text(0.5,0.5, paste(var), cex=2)
```

<img src="figures/fig-rev-unnamed-chunk-4-1.png" width="672" />

## LMM Model assumption testing

I conducted all the Linear mixed effects models and produced the figures
where I assessed the model assumptions of normal distribution of
residuals and homogeneity of variances of the residuals, but I didn't
adjust the data to fit the model assumptions. In some cases, where the
assumptions were not met, or suspected to not be met, I went back through
each variable and removed outliers that were the likely cause of poor
model fit. After removing the outliers, I re-ran the LMM to see if
removing the outliers led to better model fit. After conducting this
analysis on all tests with slight concerns for the model meeting
assumptions, we chose to only remove outliers on *Ramicrusta* and
picoeukaryotes. This is because removing the outliers led to a better
model fit and matching of assumptions. For all other taxa, removing
outliers did not change model assumptions or statistical results. I will
therefore only show the code for *Ramicrusta* and picoeukaryotes below.

### Picoeukaryotes


```r
# analyze picoeukaryotes. In the first model, there were some clear outliers in the residuals that was leading to non-normal data and right-skewed data. 
# Re-run and check the outliers and results after removing the outliers. 
vardata <- filter(env.data.foranalysis, datatype == "peuk_per_ml")
  
model <- lmerTest::lmer(value ~ disturbance + (1|site) + (1|yearmonth), data = vardata)
# Investigate the model output BEFORE removing outliers
summary(model)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: value ~ disturbance + (1 | site) + (1 | yearmonth)
##    Data: vardata
## 
## REML criterion at convergence: 2362.6
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.9688 -0.6294 -0.1198  0.5223  4.2954 
## 
## Random effects:
##  Groups    Name        Variance Std.Dev.
##  yearmonth (Intercept) 166964   408.6   
##  site      (Intercept)  14395   120.0   
##  Residual              190841   436.9   
## Number of obs: 158, groups:  yearmonth, 10; site, 8
## 
## Fixed effects:
##                      Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)          1610.362    215.809    8.309   7.462 5.89e-05 ***
## disturbancehurricane -396.663    290.176    9.066  -1.367    0.205    
## disturbancedisease   -323.090    296.424    9.562  -1.090    0.302    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) dstrbnch
## dstrbnchrrc -0.715         
## distrbncdss -0.700  0.736
```

```r
# save the predicted and raw values
pred <- predict(model)
res <- residuals(model)
res.pred.data <- data.frame(pred, res, vardata)
  
var <- "peuk_per_ml"
  
# Investigate the residuals to see if the model meets assumptions
par(mfrow=c(2,2))
plot(x = pred, y = res)
abline(h = 0)
qqnorm(res) 
qqline(res)
hist(res)
plot.new()
text(0.5,0.5, paste(var), cex=2)
```

<img src="figures/fig-rev-picoeukaryotes-1.png" width="672" />

```r
# Which samples are the high residuals so I can remove them?
# set the residual number to remove
cutoff = 1500
high.residuals <- res > cutoff
table(high.residuals) # found the top 3, need to remove
```

```
## high.residuals
## FALSE  TRUE 
##   155     3
```

```r
envdata.filt <- res.pred.data[res.pred.data$res < cutoff,]

# Try running the data again with the filtered data and see if it impacts the broad findings:
model.filt <- lmerTest::lmer(value ~ disturbance + (1|site) + (1|yearmonth), data = envdata.filt)
summary(model.filt)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: value ~ disturbance + (1 | site) + (1 | yearmonth)
##    Data: envdata.filt
## 
## REML criterion at convergence: 2253.2
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.22767 -0.73520  0.00262  0.68788  2.84349 
## 
## Random effects:
##  Groups    Name        Variance Std.Dev.
##  yearmonth (Intercept) 161777   402.2   
##  site      (Intercept)  15240   123.4   
##  Residual              120701   347.4   
## Number of obs: 155, groups:  yearmonth, 10; site, 8
## 
## Fixed effects:
##                      Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)          1609.903    210.399    8.600   7.652 4.06e-05 ***
## disturbancehurricane -469.968    278.423    9.091  -1.688    0.125    
## disturbancedisease   -372.218    283.132    9.527  -1.315    0.219    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) dstrbnch
## dstrbnchrrc -0.723         
## distrbncdss -0.711  0.793
```

```r
# save the predicted and raw values of the model WITH OUTLIERS REMOVED
pred.filt <- predict(model.filt)
res.filt <- residuals(model.filt)
res.pred.data <- data.frame(pred.filt, res.filt, envdata.filt)
  
# Investigate the residuals to see if the model meets assumptions WITH OUTLIERS REMOVED
par(mfrow=c(2,2))
plot(x = pred.filt, y = res.filt)
abline(h = 0)
qqnorm(res.filt) 
qqline(res.filt)
hist(res.filt)
plot.new()
text(0.5,0.5, paste(var), cex=2)
```

<img src="figures/fig-rev-picoeukaryotes-2.png" width="672" />

```r
# Happy with new model with outliers removed because the model assumptions (normal residucals and homogeneity of variances of residuals) is met bettter now run the rest of it to find out answers. 

# Make models with and without the random effects and a null model without disturbance as fixed effect
notime_model <- lmerTest::lmer(value ~ disturbance + (1|site), data = envdata.filt)
nosite_model <- lmerTest::lmer(value ~ disturbance + (1|yearmonth), data = envdata.filt)
null_model <- lmerTest::lmer(value ~ 1 + (1|site) + (1|yearmonth), data = envdata.filt)
  
# run a likelihood ratio test to evaluate if additional random/fixed effects leads to a better fit of the model
time_LRT <- anova(model.filt, notime_model)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
site_LRT <- anova(model.filt, nosite_model)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
disturbance_LRT <- anova(model.filt, null_model)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
# Results are different, so use the new set of data
# Plot the estimates with confidence intervals
fe <- fixef(model.filt)
ci <- confint(model.filt)
```

```
## Computing profile confidence intervals ...
```

```r
cife <- data.frame(ci[4:6,], fe, "disturbance" = c("historic", "hurricane", "disease"))
colnames(cife) <- c("lower", "upper", "fixedeff", "disturbance")
cife$disturbance <- factor(cife$disturbance, levels = c("historic", "hurricane", "disease"))

cife_noIntercept <- cife[2:3,]

# Save the picoeukaryote data to add it back with the other model information
cife_Peuk_per_ml <- tibble("variable" = c("peuk_per_ml", "peuk_per_ml"), cife_noIntercept)
```

When comparing the results of the picoeukaryote model before and after
removing the outliers to satisfy the assumptions of normality and
homogeneity of variance of residuals, I found some change in the overall
trends of the results. Specifically, site as a random effect
significantly improved model fit in the new results, while before it
didn't. As such, I will proceed reporting the results of the NEW model as
it removed outliers that were leading to the failure of the model to meet
assumptions.

### Ramicrusta


```r
# analyze Ramicrusta.  In the first model, there were some outliers in the residuals that was leading to non-normal data and right-skewed data. 
# Re-run and check the outliers and results after removing the outliers. 
vardata <- filter(env.data.foranalysis, datatype == "Ramicrusta")
  
model <- lmerTest::lmer(value ~ disturbance + (1|site) + (1|yearmonth), data = vardata)
# Investigate the model output BEFORE removing outliers
summary(model)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: value ~ disturbance + (1 | site) + (1 | yearmonth)
##    Data: vardata
## 
## REML criterion at convergence: -167.1
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.5426 -0.5309 -0.0507  0.3751  3.4239 
## 
## Random effects:
##  Groups    Name        Variance Std.Dev.
##  yearmonth (Intercept) 0.002096 0.04578 
##  site      (Intercept) 0.001789 0.04230 
##  Residual              0.003202 0.05658 
## Number of obs: 72, groups:  yearmonth, 9; site, 8
## 
## Fixed effects:
##                      Estimate Std. Error       df t value Pr(>|t|)
## (Intercept)          0.041917   0.038364 8.657001   1.093    0.304
## disturbancehurricane 0.001466   0.043345 7.304548   0.034    0.974
## disturbancedisease   0.066196   0.042551 7.107646   1.556    0.163
## 
## Correlation of Fixed Effects:
##             (Intr) dstrbnch
## dstrbnchrrc -0.751         
## distrbncdss -0.765  0.741
```

```r
# save the predicted and raw values
pred <- predict(model)
res <- residuals(model)
res.pred.data <- data.frame(pred, res, vardata)
  
var <- "Ramicrusta"

# Investigate the residuals to see if the model meets assumptions
par(mfrow=c(2,2))
plot(x = pred, y = res)
abline(h = 0)
qqnorm(res) 
qqline(res)
hist(res)
plot.new()
text(0.5,0.5, paste(var), cex=2)
```

<img src="figures/fig-rev-Ramicrusta-1.png" width="672" />

```r
# Which samples are the high residuals so I can remove them?
# set the residual number to remove
cutoff = 0.1
high.residuals <- res > cutoff
table(high.residuals) # found the top 3, need to remove
```

```
## high.residuals
## FALSE  TRUE 
##    69     3
```

```r
envdata.filt <- res.pred.data[res.pred.data$res < cutoff,]

# Try running the data again with the filtered data and see if it impacts the broad findings:
model.filt <- lmerTest::lmer(value ~ disturbance + (1|site) + (1|yearmonth), data = envdata.filt)
summary(model.filt)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: value ~ disturbance + (1 | site) + (1 | yearmonth)
##    Data: envdata.filt
## 
## REML criterion at convergence: -191
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -1.65536 -0.62940 -0.00202  0.39117  2.47757 
## 
## Random effects:
##  Groups    Name        Variance Std.Dev.
##  yearmonth (Intercept) 0.001339 0.03659 
##  site      (Intercept) 0.001264 0.03555 
##  Residual              0.001938 0.04402 
## Number of obs: 69, groups:  yearmonth, 9; site, 8
## 
## Fixed effects:
##                      Estimate Std. Error      df t value Pr(>|t|)
## (Intercept)           0.03182    0.03095 9.09999   1.028    0.330
## disturbancehurricane  0.01116    0.03464 7.54591   0.322    0.756
## disturbancedisease    0.06334    0.03410 7.35396   1.858    0.104
## 
## Correlation of Fixed Effects:
##             (Intr) dstrbnch
## dstrbnchrrc -0.746         
## distrbncdss -0.759  0.741
```

```r
# save the predicted and raw values WITH OUTLIERS REMOVED
pred.filt <- predict(model.filt)
res.filt <- residuals(model.filt)
res.pred.data <- data.frame(pred.filt, res.filt, envdata.filt)
  
# Investigate the residuals to see if the model meets assumptions WITH OUTLIERS REMOVED
par(mfrow=c(2,2))
plot(x = pred.filt, y = res.filt)
abline(h = 0)
qqnorm(res.filt) 
qqline(res.filt)
hist(res.filt)
plot.new()
text(0.5,0.5, paste(var), cex=2)
```

<img src="figures/fig-rev-Ramicrusta-2.png" width="672" />

```r
# Happy with new model with outliers removed because the model assumptions (normal residucals and homogeneity of variances of residuals) is met better so now run the rest of it to find out answers. 

# Make models with and without the random effects and a null model without disturbance as fixed effect
notime_model <- lmerTest::lmer(value ~ disturbance + (1|site), data = envdata.filt)
nosite_model <- lmerTest::lmer(value ~ disturbance + (1|yearmonth), data = envdata.filt)
null_model <- lmerTest::lmer(value ~ 1 + (1|site) + (1|yearmonth), data = envdata.filt)
  
  # run a likelihood ratio test to evaluate if additional random/fixed effects leads to a better fit of the model
time_LRT <- anova(model.filt, notime_model)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
site_LRT <- anova(model.filt, nosite_model)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
disturbance_LRT <- anova(model.filt, null_model)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
# Results are different, so use the new set of data
# Plot the estimates with confidence intervals
fe <- fixef(model.filt)
ci <- confint(model.filt)
```

```
## Computing profile confidence intervals ...
```

```r
cife <- data.frame(ci[4:6,], fe, "disturbance" = c("historic", "hurricane", "disease"))
colnames(cife) <- c("lower", "upper", "fixedeff", "disturbance")
cife$disturbance <- factor(cife$disturbance, levels = c("historic", "hurricane", "disease"))

cife_noIntercept <- cife[2:3,]

cife_Ramicrusta <- tibble("variable" = c("Ramicrusta", "Ramicrusta"), cife_noIntercept)
```

When comparing the results of the *Ramicrusta* model before and after
removing the outliers to satisfy the assumptions of normality and
homogeneity of variance of residuals, I found some change in the overall
trends of the results. Specifically, disturbance was significant at
p\<0.05 at impacting *Ramicrusta*. But before outlier removal, the model
assumptions were not as well met and the disturbance was significant at
p\<0.10. Ultimately, I conclude that I will move forward with removing
the outliers as it improved the assumptions of the model, which also
changed the conclusions.

## Figure 2

Benthic composition at eight reefs significantly changed with disturbances in St. John, US Virgin Islands. 


```r
# need confidence interval data from LMM for loop
CI_df <- ldply(CI_list, data.frame)

# remove the ramicrusta and picoeukaryote data from this since you had to adjust those models to fit the assumptionsof a linear mixed effects model
# Then add the ramicrusta and picoeukaryote data back in fromm the updated model
CI_df_all <- CI_df %>%
  filter(variable != "Ramicrusta") %>%
  filter(variable != "peuk_per_ml") %>% # remove the old ramicrusta and picoeukaryote data
  rename(lower = conf_int.lower, 
         upper = conf_int.upper, 
         fixedeff = conf_int.fixedeff, 
         disturbance = conf_int.disturbance) %>% # rename columns so you can bind rows
  bind_rows(cife_Peuk_per_ml) %>%
  bind_rows(cife_Ramicrusta) # add new data

# Figure combos to create
# 1. Ramicrusta, Turf algae, and macroalgae for Figure 2
CI_df_fig2 <- CI_df_all %>%
  filter(variable %in% c("Ramicrusta", "TurfAlgae", "Macroalgae")) %>%
  mutate(variable = factor(variable, levels = c("Ramicrusta", "TurfAlgae", "Macroalgae")))

a <- ggplot(CI_df_fig2, aes(x = disturbance, y = fixedeff, color = disturbance)) +
    geom_hline(yintercept = 0, color = "#1B9E77", size = 1) +
    geom_point(size = 5) +
    geom_errorbar(data = CI_df_fig2, aes(ymin = lower, ymax = upper), width = 0.2, size = 1.5) +
    scale_color_manual(values = c("#D95F02", "#7570B3")) +
    facet_wrap(. ~ variable, scales = "free", ncol = 2) +
    labs(x = "Disturbance", y = "Estimate with 95% CI relative to Historic")

## Conduct a PCA
benthic <- envdata %>%
  filter(site != "USVI Blue Water") %>%
  dplyr::select(Date, site, siteacronym, depthtype, yearmonth, disturbance, BleachCoral, CCA, CYAN, DiseasedCoral, HardCoral, Macroalgae, Other, Ramicrusta, SoftCoral, Sponge, Substrate, TurfAlgae) %>%
  mutate(site = factor(site, levels = c("Dittlif", "Cocoloba", "Joels Shoal", "Europa", "Yawzi", "Tektite", "Booby Rock", "Ram Head", "Sand patch"))) %>%
  drop_na(HardCoral) %>%
  filter(site != "Sand patch") %>% #no need to calculate sand patch stuff
  mutate(Other = Other + BleachCoral) %>%
  dplyr::select(-BleachCoral) %>% #no need for this section
  mutate(disturbance = factor(disturbance, levels = c("historic", "hurricane", "disease")))

# Prep the benthic data 
pca_data_short_final <- as.matrix(benthic) #change from tibble to matrix
rownames(pca_data_short_final) <- pca_data_short_final[,3] #make rownames the unique sites
pca_data_short_final <- pca_data_short_final[,7:17] #select numeric variables
class(pca_data_short_final) <- "numeric" #change from character to numeric values for pca

#Do the PCA
pca <- PCA(pca_data_short_final, scale.unit=TRUE, graph = FALSE) #performs the Principal component analysis #scale.unit=TRUE then data are scaled to unit variance

#All pca values for density plots to acompany the ggplot summary figure
pcaAll <- cbind(benthic, pca$ind$coord[,1:2]) %>%
  mutate(disturbance = factor(disturbance, levels = c("historic", "hurricane", "disease")))

pcaVars <- data.frame(pca$var$coord[,1:2], "variable" = rownames(pca$var$coord))

pcagraph <- ggscatter(pcaAll, x = "Dim.1", y = "Dim.2", color = "disturbance", size = 4, alpha = 0.7, palette = "Dark2") +
  labs(x = "PC1 (27.2%)", y = "PC2 (15.2%)") +
  border()+
  theme(legend.position = "none") +
  geom_segment(data = pcaVars, aes(x = 0, xend = Dim.1*4, y = 0, yend = Dim.2*4), arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = pcaVars, aes(x = Dim.1*4, y = Dim.2*4, label = variable)) +
  coord_fixed()

cx <- ggplot(pcaAll, aes(x = disturbance, y = Dim.1, fill = disturbance)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Dark2") +
  rotate() +
  clean_theme()+
  theme(legend.position = "none")

dy <- ggplot(pcaAll, aes(x = disturbance, y = Dim.2, fill = disturbance)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Dark2") +
  clean_theme() +
  theme(legend.position = "none")

# Arranging the plot
fig2 <- ggarrange(ggarrange(cx, NULL, pcagraph, dy, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(4,1), heights = c(1, 4)), 
          a, labels = c("a.", "b."))
annotate_figure(fig2, top = "Figure 2. Benthic composition at eight reefs significantly changed with disturbances")
```

<img src="figures/fig-rev-unnamed-chunk-5-1.png" width="672" />

```r
#Use vegan for the benthic pca stuff - PCAs are made with euclidian distance so use adonis function but make sure method = "eu"
# Interested in how environmental variables influence benthic composition
# Does site (reef) influence benthic composition
site <- adonis2(pca_data_short_final ~ site, data = benthic, method = "eu")
site[1,]
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = pca_data_short_final ~ site, data = benthic, method = "eu")
##      Df SumOfSqs      R2      F Pr(>F)    
## site  7   1.7506 0.34489 4.8135  0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Does disturbance influence benthic composition
dist <- adonis2(pca_data_short_final ~ disturbance, data = benthic, method = "eu")
dist[1,]
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = pca_data_short_final ~ disturbance, data = benthic, method = "eu")
##             Df SumOfSqs      R2      F Pr(>F)    
## disturbance  2   1.0434 0.20556 8.9268  0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Does time (date, continuous) influence benthic composition
date <- adonis2(pca_data_short_final ~ Date, data = benthic, method = "eu")
date[1,]
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = pca_data_short_final ~ Date, data = benthic, method = "eu")
##      Df SumOfSqs      R2      F Pr(>F)    
## Date  1  0.99073 0.19518 16.977  0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Does sampling timepoint (categorical) influence benthic composition
timepoint <- adonis2(pca_data_short_final ~ yearmonth, data = benthic, method = "eu")
timepoint[1,]
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = pca_data_short_final ~ yearmonth, data = benthic, method = "eu")
##           Df SumOfSqs      R2      F Pr(>F)    
## yearmonth  8   1.7017 0.33525 3.9716  0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Supplementary Figure S3

```r
# 2. Supplementary figure of all the benthic components that werent significantly distinct
# Uses code from a previous chunk and results from the linear mixed effects modeling approach
CI_df_suppfig1 <- CI_df_all %>%
  filter(variable %in% c("HardCoral", "SoftCoral", "AlgaeToCoral")) %>%
  mutate(variable = factor(variable, levels = c("HardCoral", "SoftCoral", "AlgaeToCoral")))

b <- ggplot(CI_df_suppfig1, aes(x = disturbance, y = fixedeff, color = disturbance)) +
    geom_hline(yintercept = 0, color = "#1B9E77", size = 1) +
    geom_point(size = 5) +
    geom_errorbar(data = CI_df_suppfig1, aes(ymin = lower, ymax = upper), width = 0.2, size = 1.5) +
    scale_color_manual(values = c("#D95F02", "#7570B3")) +
    facet_wrap(. ~ variable, scales = "free", ncol = 3) +
    labs(x = "Disturbance", y = "Estimate with 95% CI relative to Historic", title = "Figure S3. Benthic components not significantly different with disturbance")

b
```

<img src="figures/fig-rev-unnamed-chunk-6-1.png" width="672" />

##  Supplementary Figure S5

```r
# 3. Supplementary figure of all the nutrients and cell abundances that weren't significantly different (all but ammonium)
# Uses code from a previous chunk and results from the linear mixed effects modeling approach

CI_df_suppfig2 <- CI_df_all %>%
  filter(variable %in% c("po4_um", "silicate_um", "no2_um", "pro_per_ml", "syn_per_ml", "hbact_per_ml", "TON", "npoc_um", "Chl_ug_per_l", "no3_um", "peuk_per_ml")) %>%
  mutate(variable = factor(variable, levels = c("pro_per_ml", "syn_per_ml", "peuk_per_ml", "hbact_per_ml", "silicate_um", "npoc_um", "TON",  "Chl_ug_per_l", "no2_um", "no3_um", "po4_um")))

c <- ggplot(CI_df_suppfig2, aes(x = disturbance, y = fixedeff, color = disturbance)) +
    geom_hline(yintercept = 0, color = "#1B9E77", size = 1) +
    geom_point(size = 5) +
    geom_errorbar(data = CI_df_suppfig2, aes(ymin = lower, ymax = upper), width = 0.2, size = 1.5) +
    scale_color_manual(values = c("#D95F02", "#7570B3")) +
    facet_wrap(. ~ variable, scales = "free", ncol = 4) +
    labs(x = "Disturbance", y = "Estimate with 95% CI relative to Historic", title = "Figure S5. Nutrients and microbes not significantly different with disturbance")

c
```

<img src="figures/fig-rev-unnamed-chunk-7-1.png" width="672" />




## Figure 3 Code

First investigate the collinearity of the environmental variables, and remove variables that are highly collinear before moving on to the PCA analysis. 

```r
## Look at pairwise correlation between the variables in all the environmental data
data_cor_all <- envdata %>%                        
  dplyr::select(Date, site, siteacronym, depthtype, yearmonth, disturbance, month, temp_c:hbact_per_ml) %>% 
  mutate(TON = tn_um - (no2no3_um + nh4_um)) %>%
  mutate(no3 = no2no3_um - no2_um) %>%
  mutate(site = factor(site, levels = c("Dittlif", "Cocoloba", "Joels Shoal", "Europa", "Yawzi", "Tektite", "Booby Rock", "Ram Head", "Sand patch"))) %>%
  dplyr::select(-Phaeo_ug_per_l, -no2no3_um, -tn_um) %>%
  filter(site != "Sand patch") %>%
  mutate(disturbance = factor(disturbance, levels = c("historic", "hurricane", "disease")))

data_cor_all2 <- data_cor_all[,8:21]                                # select indep. variables
var <- cor(data_cor_all2, use = "complete.obs", method = "pearson")  # correlation matrix  with casewise deletion of missing values

colnames(var) <- colnames(data_cor_all2)                    # rename the row names and column names
rownames(var) <- colnames(data_cor_all2)
corrplot(var, method='number')                              # visualize any multicollinearity
```

<img src="figures/fig-rev-unnamed-chunk-8-1.png" width="672" />

```r
# Remove collinear variables with correlation greater than abs(0.7)
toremove <- var < (-0.7)

# Variables with TRUE are Salinity, silicate, and nitrate

## Based on this, I will remove: Silicate, Salinity, and Nitrate from this analysis (anything with correlation greater than 0.7)

# Also keep all data for analysis
collinear_remove <- data_cor_all %>%
  filter(yearmonth != "2022-06") %>% # remove june 2022 since no cell count data
  dplyr::select(!c(silicate_um, salinity, no3))

## BENTHIC and SURFACE TOGETHER ANALYSIS ##

# I did a test comparing depth of all variables and most were not significantly different in value between the surface and benthic values, so I am going to combine both for an analysis

allmatrix <- as.matrix(collinear_remove) #change from tibble to matrix
rownames(allmatrix) <- allmatrix[,3] #make rownames the unique sites
allmatrix_vars <- allmatrix[,8:18] #select numeric variables
class(allmatrix_vars) <- "numeric" #change from character to numeric values for pca

## impute missing values
replaceMissing_all <- imputePCA(allmatrix_vars)
allmatrix_vars_noNA <- replaceMissing_all$completeObs

#Do the PCA
pca_all<- PCA(allmatrix_vars_noNA, scale.unit=TRUE, graph = FALSE) #performs the Principal component analysis, and replaces missing values. #scale.unit=TRUE then data are scaled to unit variance

# Investigate eigenvalues and make a scree plot to see how many dimensions to plot
# 2 appears sufficient and i will comment out this code for now so it doesn't print in github
# get_eigenvalue(pca_all)
# fviz_eig(pca_all) #1st dimension is most. First 2 should be sufficient

all_pca_viz <- fviz_pca(pca_all,
         geom=c("point"),
         pointshape=21,
         pointsize=3,
         fill.ind=collinear_remove$yearmonth,
         invisible="quali", # gets rid of the unequally sized dots, so all of them are equal in size
         geom.var =c("arrow", "text"),
         col.var = "black",
         repel=TRUE,
         title="PCA of environmental variables",
         legend.title=list(fill="Time")) +
  scale_fill_viridis_d() +
  coord_fixed() +
  theme_bw()

# Save the graph of only ammonium from the LMM analysis
# Need the CI_df_all data frame
nh4_ci <- CI_df_all %>%
  filter(variable == "nh4_um")

nh4 <- ggplot(nh4_ci, aes(x = disturbance, y = fixedeff, color = disturbance)) +
    geom_hline(yintercept = 0, color = "#1B9E77", size = 1) +
    geom_point(size = 5) +
    geom_errorbar(data = nh4_ci, aes(ymin = lower, ymax = upper), width = 0.2, size = 1.5) +
    scale_color_manual(values = c("#D95F02", "#7570B3")) +
    labs(x = "Disturbance", y = "Estimate with 95% CI relative to Historic", title = "b. Ammonium")

# Put the 2 graphs together using ggpubr
ggarrange(all_pca_viz, ggarrange(nh4, NULL, ncol = 1), 
                    labels = c("a.", "b."), 
          widths = c(2,1)) 
```

<img src="figures/fig-rev-unnamed-chunk-8-2.png" width="672" />

```r
#Use vegan for the PERMANOVA tests on benthic pca stuff - PCAs are made with euclidian distance so use adonis function but make sure method = "eu
# Report statistics to go with graph
# Does site (reef) influence envionrmental variables?
site <- adonis2(allmatrix_vars_noNA ~ site, data = collinear_remove, method = "eu")
site[1,]
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = allmatrix_vars_noNA ~ site, data = collinear_remove, method = "eu")
##      Df   SumOfSqs      R2      F Pr(>F)  
## site  7 3.7868e+11 0.10618 2.5796  0.011 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Does sampling timepoint (yearmonth) influence environmental variables?
yearmonth <- adonis2(allmatrix_vars_noNA ~ yearmonth, data = collinear_remove, method = "eu")
yearmonth[1,]
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = allmatrix_vars_noNA ~ yearmonth, data = collinear_remove, method = "eu")
##           Df   SumOfSqs      R2      F Pr(>F)    
## yearmonth  9 1.5351e+12 0.43046 12.597  0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Does disturbance influence environmental variables?
dist <- adonis2(allmatrix_vars_noNA ~ disturbance, data = collinear_remove, method = "eu")
dist[1,]
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = allmatrix_vars_noNA ~ disturbance, data = collinear_remove, method = "eu")
##             Df   SumOfSqs       R2      F Pr(>F)  
## disturbance  2 1.0804e+11 0.030294 2.4524  0.075 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

# Figure 4 - Data Prep and Stats for Alpha & Beta Diversity

## Alpha diversity

Make sure to use the `ps.a` physloseq object since that one has all the
low abundaance ASVs RETAINED. This is important since I am trying to
estimate the true richness, so I should not use the low-abunndance
filtered data.

### Does depth influence alpha diversity?

Use this analysis to justify grouping surface and benthic (reef-depth)
samples together


```r
#set levels
sample_data(ps.a)$date <- as.Date(sample_data(ps.a)$date, tryFormats = c("%m/%d/%y"))
sample_data(ps.a)$site <- factor(sample_data(ps.a)$site, levels = c("Dittlif", "Cocoloba", "Joels Shoal", "Europa", "Yawzi", "Tektite", "Booby Rock", "Ram Head", "Sand patch"))
sample_data(ps.a)$disturbance <- factor(sample_data(ps.a)$disturbance, levels = c("historic", "hurricane", "disease"))
sample_data(ps.a)$depthtype <- factor(sample_data(ps.a)$depthtype, levels = c("benthic", "surface"))

### Does Alpha diversity change with depth?? ###
# estimate richness
richness <- ps.a %>% breakaway        # estimate true richness

# Do inference on richness - I hypothesize that richness changes with depth
meta_richness <- ps.a %>%
  sample_data %>%
  as_tibble %>%
  mutate("sample_names" = ps.a %>% sample_names) %>%
  left_join(summary(richness),
            by = "sample_names") 

# Test the effect of depth on estimated richness
bt_depth_fixed <- betta(formula = estimate ~ depthtype, 
                      ses = error, data = meta_richness)

alphadepth <- ggplot(meta_richness, aes(x = depthtype, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width=.1, height=0), aes(colour = depthtype, shape = disturbance), size = 3, alpha = 0.7) +
  labs(y = "Breakaway estimate of richness", x = "Disturbance", color = "Depth", shape = "Disturbance", 
       title = "Figure S8") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Print the statistical result answering "Does depth influence alpha diversity?"
bt_depth_fixed$table
```

```
##                  Estimates Standard Errors p-values
## (Intercept)      360.78226        3.765374        0
## depthtypesurface -27.56523        5.390113        0
```

```r
bt_depth_fixed$global[2]
```

```
## [1] 0
```

Based on this, alpha diversity significantly changes with depth,
therefore, I will split up the alpha diversity analyses and proceed first
to test changes on benthic depth seawater alpha diversity. Then I will
analyze surface depth alpha diversity.

### Benthic - Alpha diversity

Does disturbance, reef, and sampling timepoint influence alpha diversity?


```r
## BENTHIC ANALYSIS FIRST #####
ps.a_benthic <- ps.a %>%
  subset_samples(site != "Sand patch") %>%  # get rid of sandpatch
  subset_samples(depth == "Benthic")        # only interested in benthic data

benthic_richness <- ps.a_benthic %>% breakaway        # estimate true richness

# Do inference on richness - I hypothesize that richness changes by reef and disturbance
benthic_meta_richness <- ps.a_benthic %>%
  sample_data %>%
  as_tibble %>%
  mutate("sample_names" = ps.a_benthic %>% sample_names) %>%
  left_join(summary(benthic_richness),
            by = "sample_names")

# I hypothesize that alpha diversity cahnges with disturbance - IT DOES - This is for a main figure
set.seed(100) #set seed for reproducibility
bt_disturbance_fixed <- betta(formula = estimate ~ disturbance, 
                      ses = error, data = benthic_meta_richness)
bt_disturbance_fixed$table
```

```
##                      Estimates Standard Errors p-values
## (Intercept)          365.34338         6.61298    0.000
## disturbancehurricane  28.16266        13.15232    0.032
## disturbancedisease   -29.51260        10.69924    0.006
```

```r
bt_disturbance_fixed$global[2] # get global p-value to see if there are overall changes across disturbance
```

```
## [1] 0
```

```r
# I hypothesize that reef (site) will have a significant impact on alpha diversity - IT DOES
bt_site_fixed <- betta(formula = estimate ~ site, 
                      ses = error, data = benthic_meta_richness)
bt_site_fixed$global[2] # get global p-value to see if there are changes across all reefs
```

```
## [1] 0
```

```r
# I hypothesize that sampling timepoint (yearmonth, coded as year and month of sampling event) will significanlty influence alpha diversity because of seasonal fluctuations in temperature and rainfall and local conditions and other environmental features that may not be accounted for in the sampling regime. 
# Overall significant but not different between events. 
bt_yearmonth_fixed <- betta(formula = estimate ~ yearmonth, 
                      ses = error, data = benthic_meta_richness)
bt_yearmonth_fixed$global[2] # get global p-value to see if there are changes across all reefs
```

```
## [1] 0
```

```r
# check if the richness changes with disturbance within each reef site (adding in random effects)
bt_disturbance_fixed_site_random <- betta_random(chats = benthic_meta_richness$estimate,
                                       ses = benthic_meta_richness$error,
                                       X = model.matrix(~ disturbance, data = benthic_meta_richness),
                                       groups = benthic_meta_richness$site)
bt_disturbance_fixed_site_random$global[2]
```

```
## [1] 0
```

```r
## Make graphs that display the different alpha diversity metrics:

benthicdisturbance <- ggplot(benthic_meta_richness, aes(x = disturbance, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width=.1, height=0), aes(colour = disturbance, shape = disturbance), size = 3, alpha = 0.7) +
  labs(y = "Breakaway estimate of richness", x = "Disturbance", color = "Disturbance", 
       title = "Fig 4. Benthic richness") +
  scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

benthicreef <- ggplot(benthic_meta_richness, aes(x = site, y = estimate, fill = disturbance)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Breakaway estimate of richness", x = "Reef", fill = "Disturbance") +
  scale_fill_brewer(palette = "Dark2") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

benthictime <- ggplot(benthic_meta_richness, aes(x = yearmonth, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width=.1, height=0), aes(color = site, shape = disturbance), size = 3) +
  labs(y = "Breakaway estimate of richness", x = "Sampling Timepoint", color = "Reef", shape = "Disturbance") +
  scale_color_carto_d(palette = "Safe") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
```

The output above includes the statistical results that will be
incorporated into the main alpha diversity figure 4 and the supplementary
alpha diversity figure S9.

### Surface - Alpha diversity

Does disturbance, reef, and sampling timepoint influence alpha diversity?


```r
## SURFACE DATA - USE Alpha diversity phyloseq object from last chunk
ps.a_surface <- ps.a %>%
  subset_samples(site != "Sand patch") %>%  # get rid of sandpatch
  subset_samples(depth == "Surface")         # only interested in benthic data

surface_richness <- ps.a_surface %>% breakaway        # estimate true richness

# Do inference on richness - I hypothesize that richness changes by reef and disturbance
surface_meta_richness <- ps.a_surface %>%
  sample_data %>%
  as_tibble %>%
  mutate("sample_names" = ps.a_surface %>% sample_names) %>%
  left_join(summary(surface_richness),
            by = "sample_names")

# I hypothesize that alpha diversity cahnges with disturbance - IT DOES - This is for a main figure
set.seed(100) #set seed for reproducibility
bt_disturbance_fixed <- betta(formula = estimate ~ disturbance, 
                      ses = error, data = surface_meta_richness)
bt_disturbance_fixed$table
```

```
##                       Estimates Standard Errors p-values
## (Intercept)          332.002614        4.039863    0.000
## disturbancehurricane  18.616862        7.625434    0.015
## disturbancedisease    -9.512116        6.709065    0.156
```

```r
bt_disturbance_fixed$global[2] # get global p-value to see if there are overall changes across disturbance
```

```
## [1] 0
```

```r
# I hypothesize that reef (site) will have a significant impact on alpha diversity - IT DOES
bt_site_fixed <- betta(formula = estimate ~ site, 
                      ses = error, data = surface_meta_richness)
bt_site_fixed$global[2] # get global p-value to see if there are changes across all reefs
```

```
## [1] 0
```

```r
# I hypothesize that sampling timepoint (yearmonth, coded as year and month of sampling event) will significanlty influence alpha diversity because of seasonal fluctuations in temperature and rainfall and local conditions and other environmental features that may not be accounted for in the sampling regime. 
# Overall significant but not different between events. 
bt_yearmonth_fixed <- betta(formula = estimate ~ yearmonth, 
                      ses = error, data = surface_meta_richness)
bt_yearmonth_fixed$global[2] # get global p-value to see if there are changes across all reefs
```

```
## [1] 0
```

```r
# check if the richness changes with disturbance within each reef site (adding in random effects)
bt_disturbance_fixed_site_random <- betta_random(chats = surface_meta_richness$estimate,
                                       ses = surface_meta_richness$error,
                                       X = model.matrix(~ disturbance, data = surface_meta_richness),
                                       groups = surface_meta_richness$site)
bt_disturbance_fixed_site_random$global[2]
```

```
## [1] 0
```

```r
## Make graphs that display the different alpha diversity metrics:

surfacedisturbance <- ggplot(surface_meta_richness, aes(x = disturbance, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width=.1, height=0), aes(colour = disturbance, shape = disturbance), size = 3, alpha = 0.7) +
  labs(y = "Breakaway estimate of richness", x = "Disturbance", color = "Disturbance", 
       title = "Fig 4. Surface richness") +
  scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

surfacereef <- ggplot(surface_meta_richness, aes(x = site, y = estimate, fill = disturbance)) +
  geom_boxplot(outlier.shape = NA) +
  labs(y = "Breakaway estimate of richness", x = "Reef", fill = "Disturbance") +
  scale_fill_brewer(palette = "Dark2") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

surfacetime <- ggplot(surface_meta_richness, aes(x = yearmonth, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width=.1, height=0), aes(color = site, shape = disturbance), size = 3) +
  labs(y = "Breakaway estimate of richness", x = "Sample Timepoint", color = "Reef", shape = "Disturbance") +
  scale_color_carto_d(palette = "Safe") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

surfaceszn <- ggplot(surface_meta_richness, aes(x = wetdryseason, y = estimate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width=.1, height=0), aes(color = wetdryseason, shape = disturbance), size = 3, alpha = 0.7) +
  labs(y = "Breakaway estimate of richness", x = "Season", color = "Season", shape = "Disturbance", title = "Surface richness") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank()) 
```

The output above includes the statistical results that will be
incorporated into the main alpha diversity figure 4 and the supplementary
alpha diversity figure S9.

## Beta diversity

### Does depth influence beta diversity?

Make sure to use the `ps` phyloseq object that has been low abundance
filtered so it has about 14,000 taxa. This will make the analyses run
more quickly because about half of the taxa that only show up for a total
of 5 or fewer counts are incredibly sparse and less influential on
overall beta diversity.


```r
sample_data(ps)$date <- as.Date(sample_data(ps)$date, tryFormats = c("%m/%d/%y"))
sample_data(ps)$site <- factor(sample_data(ps)$site, levels = c("Dittlif", "Cocoloba", "Joels Shoal", "Europa", "Yawzi", "Tektite", "Booby Rock", "Ram Head", "Sand patch"))
sample_data(ps)$disturbance <- factor(sample_data(ps)$disturbance, levels = c("historic", "hurricane", "disease"))
sample_data(ps)$depthtype <- factor(sample_data(ps)$depthtype, levels = c("benthic", "surface"))

# remove sand patch from dataset
ps.beta <- ps %>% subset_samples(site != "Sand patch")

# transform the data
ps.beta_clr <- microbiome::transform(ps.beta, 'clr')
ps.beta_clr_euc <- ordinate(ps.beta_clr, "RDA", "euclidean")

# plot PCA
betadepth <- plot_ordination(ps.beta_clr, ps.beta_clr_euc, 
                             type="samples", color="depthtype", shape = "disturbance") +
  coord_fixed() +
  geom_point(size = 3, alpha = 0.7) +
  labs(color = "Depth", shape = "Disturbance", title = "Figure S8. Beta diversity") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank())

# Do permanova (adonis test) to evaluate effect of reef and disturbance on composition.
clr_meta <- as(sample_data(ps.beta_clr), "data.frame")
clr_otu <- otu_table(ps.beta_clr) %>% as.matrix %>% t() #samples must be rows, so transpose

#get dissimilarity matrix 
euc_diss <- vegdist(clr_otu, method = "euclidean")

# permanova test to evaluate if depth is significantly structures microbial community beta diversity
set.seed(10)
adonis2(formula = euc_diss ~ depthtype, data = clr_meta, permutations = 999) 
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = euc_diss ~ depthtype, data = clr_meta, permutations = 999)
##            Df SumOfSqs     R2      F Pr(>F)
## depthtype   1     2861 0.0035 1.2192  0.113
## Residual  347   814187 0.9965              
## Total     348   817048 1.0000
```

From this, I conclude that depth is not a strong driver of changes in
microbial community composition. I will proceed with beta diversity
analyses that combine both surface and benthic depth reef water
microbiome samples.

### Disturbance, time, reef on beta diversity


```r
# Investigate changes and drivers to beta diversity
betareef <- plot_ordination(ps.beta_clr, ps.beta_clr_euc, type="samples", color = "disturbance") +
  coord_fixed() +
  facet_wrap(~ site, ncol = 4) +
  geom_point(size = 3, alpha =  0.7) +
  scale_color_brewer(palette = "Dark2") +
  labs(color = "Disturbance", title = "Figure S7. PCA on seawater microbial community faceted by reef") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank())

betamain <- plot_ordination(ps.beta_clr, ps.beta_clr_euc, type="samples", color = "yearmonth", shape = "disturbance") +
  coord_fixed() +
  geom_point(size = 3, alpha =  0.7) +
  scale_color_viridis_d() + 
  labs(color = "Sampling Timepoin", title = "Beta diversity in benthic & surface sw") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        legend.position = "bottom")

# Does sampling timepoint influence microbiome composition?
timepoint <- adonis2(formula = euc_diss ~ yearmonth, data = clr_meta, permutations = 999)
timepoint[1,]
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = euc_diss ~ yearmonth, data = clr_meta, permutations = 999)
##           Df SumOfSqs      R2      F Pr(>F)    
## yearmonth 10   216935 0.26551 12.218  0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Does disturbance influence microbiome composition?
dist <- adonis2(formula = euc_diss ~ disturbance, data = clr_meta, permutations = 999) 
dist[1,]
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = euc_diss ~ disturbance, data = clr_meta, permutations = 999)
##             Df SumOfSqs       R2      F Pr(>F)    
## disturbance  2    49912 0.061089 11.256  0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Does reef site influence microbiome composition?
reef <- adonis2(formula = euc_diss ~ site, data = clr_meta, permutations = 999) 
reef[1,]
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = euc_diss ~ site, data = clr_meta, permutations = 999)
##      Df SumOfSqs       R2      F Pr(>F)    
## site  7    38855 0.047556 2.4323  0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# Does disturbance nested within site influence microbiome composition?
nested <- adonis2(formula = euc_diss ~ site/disturbance, data = clr_meta, permutations = 999) #disturbance nested within site
nested[1,]
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = euc_diss ~ site/disturbance, data = clr_meta, permutations = 999)
##      Df SumOfSqs       R2      F Pr(>F)    
## site  7    38855 0.047556 2.5981  0.001 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Figure 4


```r
# Figure 4
# Microbial beta diversity, alpha diversity benthic, alpha diversity surface
ggarrange(betamain, 
          ggarrange(benthicdisturbance, surfacedisturbance, labels = c("b.", "c."), ncol = 2, common.legend = TRUE, legend = "bottom"),
          ncol = 2, labels = c("a.", "b."))
```

<img src="figures/fig-rev-diversity-1.png" width="672" />

## Supplementary Figure S7


```r
# Supplementary Figure S7
# Microbial alpha and beta diversity across disturbbances within individual reefs
ggarrange(betareef, 
          common.legend = TRUE, legend = "bottom")
```

<img src="figures/fig-rev-unnamed-chunk-9-1.png" width="672" />

## Supplementary Figure S8


```r
# Supplementary Figure S8
# Microbial beta diversity, alpha diversity depth tests and results
ggarrange(betadepth, alphadepth, ncol = 2, labels = c("a.", "b."), common.legend = TRUE, widths = c(2,1))
```

<img src="figures/fig-rev-unnamed-chunk-10-1.png" width="672" />

## Supplementary Figure S9


```r
# put the supplementary alpha diversity graphs together 
ggarrange(ggarrange(benthicreef, surfacereef, common.legend = TRUE, 
                    labels = c("a. Benthic", "b. Surface")), 
          ggarrange(benthictime, surfacetime, common.legend = TRUE, 
                    labels = c("c. Benthic", "d. Surface"), legend = "bottom"), 
          ncol = 1)
```

<img src="figures/fig-rev-unnamed-chunk-11-1.png" width="672" />

# Figures 5 & 6 - Data Prep and Stats - Sensitive and Predictive Analysis

One of the main findings was the enrichment of ammonium following
disturbances. Ammonium is readily assimilated by microorganisms and may
be a key mechanism leading to shifts in microbial community composition
in the water column overlying reefs. Toward this end, I will conduct a
differential abundance test of taxa that are differentialy abundant in
response to ammonium, with the idea that this will identify taxa that are
sensitive of ammonium. I will conduct a second test using random forest
machine learning approaches to identify taxa that are predictive of
ammonium. I will focus the whole analysis on benthic-depth reef water
microbial communities as I predict that they are more readily influenced
by and influence the reef benthic habitat.

## Prep data for differential abundance test - corncob


```r
# Prep the phyloseq object - represents low-abundance filtered data
sample_data(ps)$site <- factor(sample_data(ps)$site, levels = c("Dittlif", "Cocoloba", "Joels Shoal", "Europa", "Yawzi", "Tektite", "Booby Rock", "Ram Head", "Sand patch"))
sample_data(ps)$disturbance <- factor(sample_data(ps)$disturbance, levels = c("historic", "hurricane", "disease"))
sample_data(ps)$yearmonth <- factor(sample_data(ps)$yearmonth, levels = c("2016-06", "2016-10", "2017-03", "2017-07", "2017-11", "2018-04", "2018-11", "2020-08", "2021-01", "2021-10", "2022-06"))

# Remove Ammonium values that were outliers and weren't in the environmental analysis. Also remove NAs
sample_data(ps)$missingnh4 <- is.na(sample_data(ps)$nh4)
sample_data(ps)$outliernh4 <- sample_data(ps)$nh4 > 0.8
ps.nh4.narm <- ps %>%
  subset_samples(missingnh4 == "FALSE") %>%
  subset_samples(outliernh4 == "FALSE")

# BENTHIC DATA - You are only analyzing the benthic data for the differential abundance and random forest stuff. It was not significantly different from the surface, and I would like to use the surface waters to test the RF model. 
ps.benthic.nh4.narm <- ps.nh4.narm %>%  # name with nh4 so you know it is for the ammonium test
  subset_samples(site != "Sand patch") %>%
  subset_samples(depth == "Benthic")

# Try renaming all the NA taxa using the fantaxtic package Natasha sent me
ps.benthic.nh4.narm <- name_na_taxa(ps.benthic.nh4.narm, 
                                        na_label = "Unclassified <tax> (<rank>)")
```

## Differential abundance test in corncob


```r
# What microorganisms are significantly changing with changes in ammonium? 

set.seed(1)
nh4_narm.da <- differentialTest(formula = ~ nh4, 
                             phi.formula = ~ nh4,
                             formula_null = ~ 1,
                             phi.formula_null = ~ nh4,
                             test = "LRT", boot = FALSE,
                             data = ps.benthic.nh4.narm,
                             fdr_cutoff = 0.05)

# save the nh4_narm.da output and add asv, p.adj values to it, and ultimately taxonomny info
nh4.da.models <- lapply(nh4_narm.da$significant_models, extractmods2)
names(nh4.da.models) <- nh4_narm.da$significant_taxa

# Add ASVs to the taxonomy table and save the significant asvs
tax_table(ps.benthic.nh4.narm)[,7] <- rownames(tax_table(ps.benthic.nh4.narm))
sig.taxonomy.nh4 <- as.data.frame(tax_table(ps.benthic.nh4.narm)[nh4_narm.da$significant_taxa,]) 

# Move the data from a list to a dataframe and add taxonomy info
nh4.da.models.df <- ldply(nh4.da.models, data.frame) %>% 
  left_join(sig.taxonomy.nh4, by = c(".id" = "Species")) %>%
  mutate(genusasv = paste0(Genus, "_(", .id,")"))

### Save these data so you don't have to re-run the model ###
write.table(nh4.da.models.df, "data/Sig_ASVs_NH4_Oct23.txt", sep="\t",row.names = FALSE)
```

I saved the output of all the ASVs that were differentially enriched or
depleted with ammonium. There were 22 of them.

## Supplementary Figure S10 - Differential Abundance output


```r
DAnh4 <- read_delim("data/Sig_ASVs_NH4_Oct23.txt", delim = "\t", col_names = TRUE)

DAnh4_desc <- DAnh4 %>% 
  arrange(desc(Class), Estimate)

# plot the differential abundance results

genusASVorder <- DAnh4_desc %>%
  arrange(desc(Class), Estimate) %>%
  dplyr::select(genusasv) %>%
  unique() 

DAnh4_desc$genusasv <- factor(DAnh4_desc$genusasv, levels = as.vector(genusASVorder$genusasv))
DAnh4_desc$Class <- factor(DAnh4_desc$Class, levels = c("Acidimicrobiia", 
                                                        "Alphaproteobacteria", 
                                                "Bacteroidia",
                                                "Cyanobacteriia",
                                                "Gammaproteobacteria", 
                                                "Unclassified SAR324 clade(Marine group B) (Phylum)", 
                                                "Verrucomicrobiae"))

palcolors <- c("#990F26", "#CC7A88", "#99600F", "#54990F", "#967ACC", "#7ABECC", "#999999")

ggplot(DAnh4_desc, aes(x = genusasv, y = Estimate, color = Class)) +
  geom_errorbar(aes(ymin = Estimate-Std.Error, ymax = Estimate+Std.Error), color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4) +
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Differential Abundance") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(plot.background = element_blank(), 
        legend.position = "bottom") +
  scale_color_manual(values = palcolors)
```

<img src="figures/fig-rev-Corncob_results-1.png" width="672" />

## Prep data for Random Forest Analysis


```r
# Use the same phyloseq object where I removed the samples with NAs for the ammonium values
# BENTHIC FILTERED data (~1000 taxa)
ps.benthic.nh4.narm.filt <- ps.benthic.nh4.narm %>% 
  filter_taxa(function(x) mean(x) > 0.5, TRUE) # keep if avg abundance 0.5 counts

taxa <- ps.benthic.nh4.narm.filt %>% tax_table() %>% rownames

# Get surface taxa based on the asvs in the benthic seawater - This will be for validation of the random forest model
ps.surface.nh4.narm.filt <- ps.nh4.narm %>% 
  subset_samples(site != "Sand patch") %>%
  subset_samples(depth == "Surface") %>%
  mutate_tax_table(ASV = .otu) %>% 
  subset_taxa(ASV %in% taxa)
```

## Random forest model and testing


```r
## RANDOM FOREST ANALYSIS ##
# get the ASVs as columns and samples as rows, which are the predictors
predictors <- t(otu_table(ps.benthic.nh4.narm.filt))

# response variable is ammonium
response <- sample_data(ps.benthic.nh4.narm.filt)$nh4

# combine with the ASV table
resp.pred <- data.frame(response, predictors)

# set seed for reproducibility and this is the main random forest model
set.seed(10)
nh4.RF <- randomForest(response ~ ., data = resp.pred, ntree = 1000, importance = TRUE)

# check output of model
# Explains 66% of variance?!?!?!?!?!? That is way higher than anything I have seen in tutorials
print(nh4.RF)
```

```
## 
## Call:
##  randomForest(formula = response ~ ., data = resp.pred, ntree = 1000,      importance = TRUE) 
##                Type of random forest: regression
##                      Number of trees: 1000
## No. of variables tried at each split: 345
## 
##           Mean of squared residuals: 0.009892131
##                     % Var explained: 66.48
```

```r
# Look at the predicted ammonium compared to actual ammonium in BENTHIC data
predicted.actual.benthic <- data.frame("predicted" = nh4.RF$predicted, 
                                       "actual" = resp.pred$response,
                                       sample_data(ps.benthic.nh4.narm.filt))

benthRFplot <- ggplot(predicted.actual.benthic, aes(x = actual, y = predicted)) +
  geom_point(aes(color = disturbance), size = 3, alpha = 0.8) +
  geom_smooth(method = 'lm', se = FALSE, color = "black") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Actual ammonium", y = "Predicted ammonium", 
       title = "Prediction of ammonium from Random Forest model on benthic data")

# Validate the model on surface data
validation.pred <- t(otu_table(ps.surface.nh4.narm.filt))

validation.resp <- sample_data(ps.surface.nh4.narm.filt)$nh4

valid.pred.resp <- data.frame(validation.resp, validation.pred)

prediction_surface <- predict(nh4.RF, valid.pred.resp[,-1])

predicted.actual.surface <- data.frame("predicted" = prediction_surface, 
                                       "actual" = sample_data(ps.surface.nh4.narm.filt)$nh4,
                                       sample_data(ps.surface.nh4.narm.filt))

surfRFplot <- ggplot(predicted.actual.surface, aes(x = actual, y = predicted)) +
  geom_point(aes(color = disturbance), size = 3, alpha = 0.8) +
  geom_smooth(method = 'lm', se = FALSE, color = "black") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Actual ammonium", y = "Predicted ammonium", 
       title = "Prediction of ammonium in surface water, model trained on benthic data")

# I need to calculate the mean squared error between the actual and predicted values using the model on the test set (surface seawater data)
mse_predict <- mean((prediction_surface-validation.resp)^2) #on the surface data using model trained on benthic data
mse_predict
```

```
## [1] 0.009975977
```

```r
mse_model <- mean((nh4.RF$predicted-resp.pred$response)^2) #on the benthic data - direct from the model with data used to train model
# THIS IS PRETTY LOW COMPARED TO THE ONLINE TUTORIALS I AM FINDING
mse_model
```

```
## [1] 0.009892131
```

## Save the top 50 predictive taxa


```r
# make a data frame with predictor names and the importance value
imp <- importance(nh4.RF)
imp <- data.frame(predictors = rownames(imp), imp)

# order the predictor levels by importance (higher number is more important)
imp.sort.Accuracy <- arrange(imp, desc(X.IncMSE))

# Select the top 50 important taxa based on the Mean Decrease Accurracy metric
imp.50.acc <- imp.sort.Accuracy[1:50,]

# Save a table of the top 50 ASVs
ASV.50.acc <- imp.50.acc$predictors
r <- rownames(tax_table(ps.benthic.nh4.narm.filt)) %in% ASV.50.acc
table <- tax_table(ps.benthic.nh4.narm.filt)[r,]

table <- as.data.frame(table)
table$ASV <- rownames(table)

tableAccuracy <- left_join(table, imp.50.acc, by = c("ASV" = "predictors"))

# Write the table in case someone else wants to generate the figure
write.table(tableAccuracy, "data/RandomForest_NH4_top50_accuracy_benthic.txt", sep = "\t", row.names = FALSE)

# Graph the top 50 taxa
tableAccuracy.genusasv <- tableAccuracy %>%
  mutate(genusasv = paste0(Genus, "_(", ASV,")"))

genusASVorder <- tableAccuracy.genusasv %>%
  arrange(desc(Class), X.IncMSE) %>%
  dplyr::select(genusasv) %>%
  unique() 

tableAccuracy.genusasv$genusasv <- factor(tableAccuracy.genusasv$genusasv, levels = as.vector(genusASVorder$genusasv))
tableAccuracy.genusasv$Class <- factor(tableAccuracy.genusasv$Class, levels = c("Alphaproteobacteria",
                                                                                "Bacteroidia",
                                                                                "Bdellovibrionia",
                                                                                "Cyanobacteriia",
                                                                                "Dadabacteriia",
                                                                                "Gammaproteobacteria",
                                                                                "Unclassified Marinimicrobia (SAR406 clade) (Phylum)",
                                                                                "Unclassified SAR324 clade(Marine group B) (Phylum)",
                                                                                "Thermoplasmata",
                                                                                "Verrucomicrobiae"))


palcolorsRF <- c("#CC7A88", "#99600F", "#CCAA7A", "#54990F", "#3D0F99", "#967ACC", "#3E9FB3", "#7ABECC", "#333333", "#999999")


# Graph all 50 taxa and their percent increase in MSE
top50 <- ggplot(tableAccuracy.genusasv, aes(x = genusasv, y = X.IncMSE, fill = Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = palcolorsRF) +
  theme(plot.background = element_blank()) +
  labs(x = "Taxa", y = "% Increase MSE") +
  theme(legend.position = "bottom")
```

Explanation of Random Forest 50 taxa. Ultimately we want to have as small
a mean squared error as possible in the model. The most important ASVs
are those that reduce the variance in the model. In other words, when the
ASV is removed, it leads to the highest INCREASE in MSE (mean square
error) in the regression-based model. More error is generally bad. This
is the metric that the random forest algorithm uses to identify an
"IMPORTANT" ASV. I just chose to select the top 50 ASVs that are
important based on their % increase in mean square error of the model
when they are removed.

## Supplementary Figure S11


```r
# Combine the top 50 plot with the regression model plots
supp <- ggarrange(ggarrange(benthRFplot, surfRFplot, 
          ncol = 1, nrow = 2, labels = c("a.", "b.")), 
          top50, widths = c(2,3), labels = c("a.", "c."))
```

```
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
```

```r
supp
```

<img src="figures/fig-rev-unnamed-chunk-15-1.png" width="672" />

## Figure 5 - Sensitive AND Predictive taxa


```r
DAtaxa <- read_delim("data/Sig_ASVs_NH4_Oct23.txt", delim = "\t", col_names = TRUE)

RFtaxa <- read_delim("data/RandomForest_NH4_top50_accuracy_benthic.txt", delim = "\t", col_names = TRUE)

RFasvs <- RFtaxa %>% 
  dplyr::select(ASV, X.IncMSE)

shared <- DAtaxa %>%
  left_join(RFasvs, by = c(".id" = "ASV")) %>%
  drop_na(X.IncMSE) # remove NAs, which indicates the microbe was NOT a predictor in the random forest analysis

importantASVs <- shared$.id %>% unique

r <- rownames(tax_table(ps.benthic.nh4.narm.filt)) %in% importantASVs

sharedtaxa <- as.data.frame(tax_table(ps.benthic.nh4.narm.filt)[r,])

## Make the shared graph
# corncob results
# order the taxa how i want them in the plot, with the archaea at the bottom
genusASVorder <- shared %>%
  arrange(desc(Class), Estimate) %>%
  dplyr::select(genusasv) %>%
  unique() 

shared$genusasv <- factor(shared$genusasv, levels = as.vector(genusASVorder$genusasv))
shared$Class <- factor(shared$Class, levels = c("Alphaproteobacteria", 
                                                "Bacteroidia",
                                                "Cyanobacteriia",
                                                "Gammaproteobacteria", 
                                                "Unclassified SAR324 clade(Marine group B) (Phylum)"))


palcolor3 <- c("#CC7A88", "#99600F", "#54990F", "#967ACC","#7ABECC")

a <- ggplot(shared, aes(x = genusasv, y = Estimate, color = Class)) +
  geom_errorbar(aes(ymin = Estimate-Std.Error, ymax = Estimate+Std.Error), color = "black", width = .3, position=position_dodge(.9)) +
  geom_point(size = 4) +
  coord_flip() +
  theme_bw() +
  labs(x = "Taxa", y = "Differential Abundance") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(plot.background = element_blank(), 
        legend.position = "bottom") +
  scale_color_manual(values = palcolor3)

b <- ggplot(shared, aes(x = genusasv, y = X.IncMSE, fill = Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = palcolor3) +
  theme(plot.background = element_blank()) +
  labs(x = "Taxa", y = "% Increase MSE") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

wrap_plots(a, b, widths = c(1,1))
```

<img src="figures/fig-rev-unnamed-chunk-16-1.png" width="672" />

## Figure 6 - Relative abundance of sensitive and predictive taxa

Relative abundance of shared sensitive (differentially abundant) and
predictive (random forests) taxa in relation to ammonium


```r
# Shared is from the section on "How many ASVs were shared?"
# prep the data
ASVshared <- shared %>%
  dplyr::select(.id) %>%
  unique()

ASVshared <- ASVshared$.id # must be a vector

# Take data from only benthic samples, since the main models were tested on benthic data. Could also show benthic AND surface if you want.
ps.shared.ra <- ps.benthic.nh4.narm %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  mutate_tax_table(ASV = .otu) %>%
  subset_taxa(ASV %in% ASVshared) 

ASVorder <- shared %>%
  arrange(Kingdom, desc(Class), Estimate) %>%
  dplyr::select(.id) %>%
  unique()

# plot relative abundance data
shared.ra.long <- ps.shared.ra %>%
  otu_table() %>%
  as_tibble() %>%
  left_join(as_tibble(tax_table(ps.shared.ra)), by = ".otu") %>%
  left_join(as_tibble(sample_data(ps.shared.ra)), by = ".sample") %>%
  mutate(genusasv = paste0(Genus, "_(", .otu, ")")) %>%
  mutate(.otu = factor(.otu, levels = as.vector(ASVorder$.id))) %>%
  mutate(genusasv = factor(genusasv, levels = as.vector(genusASVorder$genusasv))) %>% 
  mutate(Class = factor(Class, levels = c("Alphaproteobacteria", 
                                                "Bacteroidia",
                                                "Cyanobacteriia",
                                                "Gammaproteobacteria", 
                                                "Unclassified SAR324 clade(Marine group B) (Phylum)")))

ggplot(shared.ra.long, aes(x = nh4, y = .abundance)) +
  geom_point(aes(color = disturbance), size = 2, alpha = 0.7) +
  geom_smooth(method = 'lm', se = FALSE, color = "black") +
  facet_wrap(~ Class + genusasv, ncol = 5, scales = "free_y") +
  labs(y = "Relative abundance", x = "Ammonium (uM)", color = "Disturbance") +
  scale_color_brewer(palette = "Dark2") +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position = "bottom")
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

<img src="figures/fig-rev-unnamed-chunk-17-1.png" width="672" />

# Supplementary Figure S12 - Pro vs ammonium


```r
pronh4 <- envdata %>%
  filter(site != "USVI Blue Water") %>%
  dplyr::select(Date, site, siteacronym, yearmonth, disturbance, pro_per_ml, nh4_um) %>%
  mutate(site = factor(site, levels = c("Dittlif", "Cocoloba", "Joels Shoal", "Europa", "Yawzi", "Tektite", "Booby Rock", "Ram Head", "Sand patch"))) %>%
  filter(site != "Sand patch") %>% #no need to calculate sand patch stuff
  mutate(disturbance = factor(disturbance, levels = c("historic", "hurricane", "disease"))) %>%
  filter(nh4_um < 0.8) # remove outlier/high ammonium values removed from environmental data analysis 

ggplot(pronh4, aes(x = nh4_um, y = pro_per_ml)) +
  geom_point(aes(color = disturbance), size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  ylim(0,150000) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Ammonium (uM)", y = "Prochlorococcus abundance (cells ml-1)", color = "Disturbance", title = "Figure S12. Prochlorococcus vs. Ammonium")
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

<img src="figures/fig-rev-pro vs nh4-1.png" width="672" />

```r
lm <- lm(pro_per_ml ~ nh4_um, data = pronh4)
summary(lm)
```

```
## 
## Call:
## lm(formula = pro_per_ml ~ nh4_um, data = pronh4)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -61646 -23576  -3465  14862  79872 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    63548       3650  17.412  < 2e-16 ***
## nh4_um        -79558      14730  -5.401 2.59e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 27710 on 148 degrees of freedom
##   (15 observations deleted due to missingness)
## Multiple R-squared:  0.1646,	Adjusted R-squared:  0.159 
## F-statistic: 29.17 on 1 and 148 DF,  p-value: 2.585e-07
```

# Supplementary Figure S13 - temperature


```r
# hobo and lameshur bay data from NOAA
tempwlam <- multiplesheets(temperaturepath)

tempwlamdf <- ldply(tempwlam, data.frame, .id = "reef") %>%
  filter(reef != "buoy41052") %>%
 mutate(reefname = case_when(reef == "lamv3" ~ "LameshurBay_NOAA",
                             reef == "BR" ~ "Booby Rock", 
                             reef == "CO" ~ "Cocoloba",
                             reef == "DL" ~ "Dittlif",
                             reef == "DP" ~ "Dittlif",
                             reef == "EU" ~ "Europa",
                             reef == "JS" ~ "Joels Shoal",
                             reef == "RB" ~ "Sand patch",
                             reef == "RH" ~ "Ram Head",
                             reef == "TK" ~ "Tektite",
                             reef == "YZ" ~ "Yawzi",
                             ))

tempwlamdf$Time <- as.POSIXct(tempwlamdf$Time, tz = "", format = "%Y/%m/%d %H:%M")
```


```r
trendline <- lm(temperature_daily_mean_C ~ Time, data = tempwlamdf)

ggplot(tempwlamdf, aes(x = Time, y = temperature_daily_mean_C)) +
  geom_vline(xintercept = as.POSIXct(as.Date("2017-09-06")), colour = "darkred", linewidth = 4) + # hurricane
  geom_vline(xintercept = as.POSIXct(as.Date("2020-06-01")), colour = "darkred", linewidth = 4) + # disease emergence
  geom_vline(xintercept = as.POSIXct(as.Date("2016-06-01")), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.POSIXct(as.Date("2016-10-01")), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.POSIXct(as.Date("2017-03-01")), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.POSIXct(as.Date("2017-07-01")), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.POSIXct(as.Date("2017-11-01")), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.POSIXct(as.Date("2018-04-01")), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.POSIXct(as.Date("2018-11-01")), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.POSIXct(as.Date("2020-08-01")), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.POSIXct(as.Date("2021-01-01")), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.POSIXct(as.Date("2021-10-01")), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.POSIXct(as.Date("2022-06-01")), colour = "gray", linewidth = 2) +
  geom_line(data = tempwlamdf, aes(color = reefname, group = format(Time, "%Y-%m"))) +
  stat_poly_line(method = "lm") +
  stat_poly_eq(aes(label = after_stat(eq.label)), label.x = as.Date("2022-01-01")) +
  scale_color_carto_d(palette = "Bold") +
  labs(y = "Mean daily temperature (Celsius)", title = "Figure S13. Temperature at reef depth over seven years", color = "Site")
```

<img src="figures/fig-rev-temp2-1.png" width="672" />

# Supplementary Figure S14 - Precipitation


```r
precip2016to2022 <- precip %>% 
  mutate(DATE =  as.Date(DATE, "%Y-%m-%d")) %>%
  filter(between(DATE, as.Date("2016-01-15"), as.Date("2022-12-15")))

ggplot(precip2016to2022, aes(x = DATE, y = PRCP)) +
  geom_vline(xintercept = as.Date("2017-09-06"), colour = "darkred", linewidth = 4, alpha = 0.8) + # hurricane
  geom_vline(xintercept = as.Date("2020-06-01"), colour = "darkred", linewidth = 4, alpha = 0.8) + # disease emergence
  geom_vline(xintercept = as.Date("2016-06-01"), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.Date("2016-10-01"), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.Date("2017-03-01"), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.Date("2017-07-01"), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.Date("2017-11-01"), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.Date("2018-04-01"), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.Date("2018-11-01"), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.Date("2020-08-01"), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.Date("2021-01-01"), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.Date("2021-10-01"), colour = "gray", linewidth = 2) +
  geom_vline(xintercept = as.Date("2022-06-01"), colour = "gray", linewidth = 2) +
  geom_line(color = "darkblue", linewidth = 1.5) +
  labs(x = "Date", y = "Precipitation", title = "Figure S14. Monthly precipitation at East End, US Virgin Islands")
```

<img src="figures/fig-rev-unnamed-chunk-18-1.png" width="672" />

# Supplementary Figure S2 - benthic cover over time

```r
benthic_vars <- c("Ramicrusta","TurfAlgae","Macroalgae","HardCoral","SoftCoral")
benthic.forFigures <- env.long.olrm %>%
  filter(datatype %in% benthic_vars) %>%
  drop_na(value) %>%
  mutate(datatype = factor(datatype, levels = 
                             c("Ramicrusta","TurfAlgae","Macroalgae","HardCoral","SoftCoral")))

# all values over time (color by reef)
ggplot(benthic.forFigures, aes(x = Date, y = value, color = site)) +
  geom_vline(xintercept = as.POSIXct(as.Date("2017-09-06")), colour = "gray", linewidth = 3) + # Irma and Maria
  geom_vline(xintercept = as.POSIXct(as.Date("2020-06-01")), colour = "gray", linewidth = 3) + # Bleaching event
  geom_point(aes(shape = disturbance), data = benthic.forFigures, size = 4, alpha = 0.8) +
  geom_line(linewidth = 1.2) +
  scale_color_carto_d(palette = "Safe") +
  facet_wrap( ~ datatype, scales = "free_y", ncol = 3) +
  labs(y = "Relative cover", x = "Date", color = "Reef", shape = "Disturbance", title = "Figure S2. Benthic organism cover over time")
```

<img src="figures/fig-rev-extra_supplementary_figs-1.png" width="672" />

# Supplementary Figure S4 - environmental nutrient data over time

```r
nuts_vars <- c("nh4_um","TON","Chl_ug_per_l","po4_um","silicate_um", "npoc_um","no2_um", "no3_um")
nuts.forFigures <- env.long.olrm %>%
  filter(datatype %in% nuts_vars) %>%
  drop_na(value) %>%
  mutate(datatype = factor(datatype, levels = c("nh4_um","TON","Chl_ug_per_l","po4_um","silicate_um", "npoc_um","no2_um", "no3_um")))

# all values over time (color by reef)
ggplot(nuts.forFigures, aes(x = Date, y = value, color = site)) +
  geom_vline(xintercept = as.POSIXct(as.Date("2017-09-06")), colour = "gray", linewidth = 3) + # Irma and Maria
  geom_vline(xintercept = as.POSIXct(as.Date("2020-06-01")), colour = "gray", linewidth = 3) + # Bleaching event
  geom_point(aes(shape = depthtype), data = nuts.forFigures, size = 4, alpha = 0.8) +
  geom_line(linewidth = 1.2) +
  scale_color_carto_d(palette = "Safe") +
  facet_wrap( ~ datatype, scales = "free_y", ncol = 4) +
  ggtitle("Figure S4. Nutrients over time")
```

<img src="figures/fig-rev-unnamed-chunk-19-1.png" width="672" />

# Supplementary Figure S6 - microbial cell abundances over time

```r
cells_vars <- c("pro_per_ml","peuk_per_ml","syn_per_ml","hbact_per_ml")
cells.forFigures <- env.long.olrm %>%
  filter(datatype %in% cells_vars) %>%
  drop_na(value) %>%
  filter(datatype != "biomassRatio") %>%
  mutate(datatype = factor(datatype, levels = c("pro_per_ml","peuk_per_ml","syn_per_ml","hbact_per_ml")))

# all values over time (color by reef)
ggplot(cells.forFigures, aes(x = Date, y = value, color = site)) +
  geom_vline(xintercept = as.POSIXct(as.Date("2017-09-06")), colour = "gray", linewidth = 3) + # Irma and Maria
  geom_vline(xintercept = as.POSIXct(as.Date("2020-06-01")), colour = "gray", linewidth = 3) + # Bleaching event
  geom_point(aes(shape = depthtype), data = cells.forFigures, size = 4, alpha = 0.8) +
  geom_line(linewidth = 1.2) +
  scale_color_carto_d(palette = "Safe") +
  facet_wrap( ~ datatype, scales = "free_y", ncol = 2) +
  labs(y = "Abundance (cells per milliliter)", shape = "Depth", color = "Reef", title = "Figure S6. Cell abundances over time") 
```

<img src="figures/fig-rev-unnamed-chunk-20-1.png" width="672" />



