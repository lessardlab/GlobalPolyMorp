This directory contains a series of R scripts that handle various aspects of data preparation, model fitting, and ecological analysis. Below is an outline of the main files and their functions:

## 00_Functions
Description: Contains all reusable functions utilized throughout the analysis. Make sure to source this script first to ensure that the necessary functions are available for other scripts.

## 01_dataPrep

> This script processes a dataset of species polymorphism occurrences, matching species names and adding a polymorphism identifier (`poly_id`). It first reads species and occurrence data, cleans the data by removing NAs and matching coordinates to continents, and excludes occurrences in oceans and Hawaii. The script then calculates the density of species occurrences at different site coordinates. Climate data (mean annual temperature and precipitation) is extracted for each occurrence point. The data is further refined by removing invalid continent entries, and the final dataset is saved. Additionally, polymorphism data from external files is cleaned and matched with species subfamilies, and the results are exported for further analysis.

## 02_EqualAreaAntAssemblages

> This script processes ant occurrence data to create 1x1-degree spatial assemblages using an equal-area projection. It begins by loading and cleaning the occurrence dataset, removing duplicates, and setting up the spatial projection of the data. The script reprojects the spatial points to an equal-area projection (LAEA) and overlays them on a global 1x1-degree grid. Species presence-absence data is tabulated for each grid cell. The centroids of these grid cells are calculated and assigned to the corresponding ant occurrences. Climate data (mean annual temperature and precipitation) is then extracted for each occurrence point and added to the final dataset, which is prepared for further analysis.


## 03_FitModels

> This script fits linear mixed-effect models to analyze worker caste polymorphism in ants using occurrence and environmental data. After loading dependencies and ensuring data completeness, the script transforms climate variables (temperature and precipitation) to a standardized range. Several generalized linear mixed models (GLMMs) are then fitted using `glmer`, with various combinations of climate variables and random effects based on grid IDs. These models evaluate how climate affects polymorphism across different species, genera, and subfamilies. The script compares model performance using AIC, and the best-fitting models are identified and summarized. Additionally, models excluding specific genera (like *Pheidole* and *Camponotus*) are tested to assess their impact on the results. The script further explores subfamily-specific and genus-specific models, calculates partial AIC values, and stores the outputs for comparison. Lastly, it generates plots showing the relationships between polymorphism, temperature, and precipitation.

## 03_ModelswtRichness

> This script fits linear mixed-effect models to study worker caste polymorphism in ants, incorporating species richness per grid as a covariate. After adding richness and its logarithmic transformation to the dataset, the script fits several generalized linear models (GLMs) to explore the effects of climate variables (temperature and precipitation) and species richness on polymorphism. Different models with combinations of these factors are tested, and the best models are identified using AIC values. Additional models test the exclusion of certain ant genera (*Pheidole* and *Camponotus*) to assess their influence on the results. The script also explores the effects of subfamily-level variation and calculates partial AIC values for comparisons. Finally, the results are visualized using plots to show the relationship between temperature, precipitation, and polymorphism, and the results are saved for further analysis.


## 04_BiomeTaxa

> This script loads ecoregion shapes and transforms them into an equal-area projection to analyze the distribution of ant species across biomes. It creates spatial objects from species occurrence data and overlays them with ecoregion data to associate species with specific biomes. The script defines functions to calculate proportions of worker caste polymorphism and species occurrence counts by biome. These results are used to generate heatmaps that show species richness, polymorphism, and observation counts per genus across different biomes. The script concludes by saving these heatmaps as PNG files and generating a custom legend for the biomes. Finally, it combines the heatmaps and the biome legend into a single image for visualization.

## 05_EcoRegions

> This script focuses on visualizing the interaction between temperature and precipitation globally, analyzing biomes and mapping species richness and distribution using various raster and spatial methods. It uses a PCA (Principal Component Analysis) on climate variables (precipitation and temperature) and creates global maps to show the interaction between these variables. The script also plots and aggregates species distribution data by biome, extracting values for each biome and displaying them in climatic space. Additionally, it generates 3D scatter plots to depict species polymorphism probability in relation to climatic variables. The script further creates boxplots to show the predicted probabilities of polymorphism across different biomes and composes multiple figures, including maps of biomes, heatmaps of biome probabilities, and richness distribution maps for different ant subfamilies (*Formicinae*, *Myrmicinae*, *Dolichoderinae*). These results are visualized and saved as PNG images, with composite images created to combine individual figures into cohesive final outputs for reporting and analysis.
