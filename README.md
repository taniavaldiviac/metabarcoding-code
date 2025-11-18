# Metabarcoding de comunidades de eucariontes, ICMyL-UNAM

### Guía para el uso del repositorio metabarcoding-code

Este instructivo te permite reproducir el análisis y el reporte de Quarto en cualquier servidor usando la carpeta `metabarcoding-code`.  

## Pasos para reproducir el análisis

### 1. Instala Miniconda/Conda en el servidor (si no está instalado)

### 2. Crea y activa el entorno conda

```sh
conda activate r-bio
```

### 3. Instala los paquetes base de R y dependencias del sistema

```sh
conda install -c conda-forge libcurl
```

### 4. Clona tu repositorio

```sh
git clone https://github.com/taniavaldiviac/metabarcoding-code/
cd metabarcoding-code
```

### 5. Abre R en el entorno conda

```sh
R
```

### 6. Instala los paquetes CRAN y Bioconductor desde R

```r
install.packages(c(
  "RColorBrewer", "ape", "colorspace", "data.table", "dplyr", "filesstrings",
  "ggforce", "ggplot2", "ggsci", "here", "insect", "kableExtra", "pander",
  "pandoc", "pheatmap", "purrr", "readr", "rentrez", "scales", "seqinr",
  "sqldf", "strex", "stringr", "svglite", "taxize", "tibble", "tidyr",
  "tidyverse", "vegan", "viridis"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ShortRead", "phyloseq", "dada2", "Biostrings"))
```

### 7. Si usas renv, restaura el entorno

```r
renv::restore()
```

### 8. Si los paquetes base no se ven en renv, agrega la ruta de conda a `.libPaths()`

```r
.libPaths(c(.libPaths(), "/home/bio/envs/r-bio/lib/R/library"))
```

### 9. Renderiza tu reporte Quarto
git push

```sh
quarto render scripts/Report.qmd
```

---

Usa la carpeta `metabarcoding-code` como tu proyecto.  
Sigue estos pasos para reproducir el análisis y reporte en cualquier servidor.
