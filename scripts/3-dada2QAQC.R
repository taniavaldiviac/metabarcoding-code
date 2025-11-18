# dada2 QAQC 
# Metabarcoding en comunidades de eucariontes

## Set up reproducible environment -------------------------------------------------------
# Auto-activar renv (útil cuando se ejecuta con Rscript en servidores sin RStudio)
if (file.exists(file.path("renv","activate.R"))) {
  source(file.path("renv","activate.R"))
}

# Cargamos librerías necesarias. Cada una se usa en distintas partes del pipeline.
library(dada2)     # pipeline principal para denoising
library(tidyverse) # manipulación de datos (incluye readr)
library(seqinr)    # utilidades FASTA/FASTA
library(ShortRead) # lectura/QA de FASTQ
library(digest)    # generación de hashes SHA1 para ASVs

# Helpers reproducibilidad y manejo de rutas
library(here)      # rutas relativas al RProject (más robusto que setwd())
library(fs)        # operaciones con archivos/carpeta (cross-platform)
library(readr)     # lectura/escritura de CSVs (mejor manejo de encoding)

# Fijar semilla si hay pasos aleatorios y opciones generales
set.seed(42)                        
#options(stringsAsFactors = FALSE)   # evitar factores implícitos

# Definición de rutas usando el root del proyecto
# here::here() detecta automáticamente la raíz (donde está el .Rproj)
proj_root <- here::here()
message("Proyecto detectado en: ", proj_root)

# Fijar working directory a la raíz del proyecto
# Esto es útil para código que no usa here() internamente y para Rscript
setwd(proj_root)
message("Working directory fijado a: ", proj_root)

# Definición de rutas de entrada/salida relativas a la raíz del proyecto
fastq_location    <- file.path(proj_root, "for_dada2")   # entrada: FASTQ sin procesar
output_location   <- file.path(proj_root, "final_data")  # salida: logs, csv, rdata
metadata_location <- file.path(proj_root, "metadata")    # metadatos (primer_data.csv, dbs, known_hashes)
primer_csv       <- file.path(metadata_location, "primer_data.csv") # tabla de parámetros por locus/primer  
i <- 1  # índice de fila en primer_data.csv para procesar (cambiar si se desea)
run_name <- "13112025" # identificador de la corrida (cambiar por fecha/ID si se desea)

# Leer tabla de parámetros por locus/primer de forma robusta y validar columnas
# Leer primer_data asegurando que las dos primeras columnas son character
# Obtén nombres de columnas

primer_cols <- readr::read_csv(primer_csv, n_max = 0, show_col_types = FALSE) %>% names()
if(!file.exists(primer_csv)) stop("No se encuentra primer_data.csv en: ", primer_csv)

(primer.data <- readr::read_csv(
  primer_csv,
  col_types = readr::cols(
    !!primer_cols[1] := readr::col_character(),
    !!primer_cols[2] := readr::col_character(),
    .default = readr::col_guess())))

# Lista de columnas que el script requiere; detener con mensaje claro si faltan
required_cols <- c("locus_shorthand","F_qual","R_qual","db_name","max_trim")
missing_cols <- setdiff(required_cols, colnames(primer.data))
if(length(missing_cols) > 0) stop("Faltan columnas en primer_data.csv: ", paste(missing_cols, collapse=", "))

# Guardar información de sesión para que los alumnos puedan reproducir el entorno
writeLines(capture.output(sessionInfo()), file.path(output_location, "logs", "sessionInfo.txt"))

# Loop principal: procesar cada fila de primer.data (cada locus/primer)
#for (i in 1:nrow(primer.data)) {
  # Buscar archivos forward (R1) que contengan el shorthand del locus.
  # Esto permite tener múltiples loci en la carpeta for_dada2 y procesarlos por separado.
  check <- grep(primer.data$locus_shorthand[i],
                sort(list.files(fastq_location, pattern="R1_001.fastq", full.names = TRUE)),
                value = TRUE)
  # Si no hay archivos para ese locus, se salta la iteración (buena práctica en pipelines)
  #if(length(check) > 0) {
    
    ### Lectura de archivos FASTQ ---------------------------------------------------
    # fnFs y fnRs son vectores con rutas a archivos forward y reverse respectivamente
    fnFs <- grep(primer.data$locus_shorthand[i],
                 sort(list.files(fastq_location, pattern="R1_001.fastq", full.names = TRUE)),
                 value = TRUE)
    fnRs <- grep(primer.data$locus_shorthand[i],
                 sort(list.files(fastq_location, pattern="R2_001.fastq", full.names = TRUE)),
                 value = TRUE)
    
    # Construcción de nombres de muestra a partir del nombre de archivo (ajustar si su formato difiere)
    sample.names1 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    sample.names2 <- sapply(strsplit(basename(fnFs), "_"), `[`, 3)
    sample.names <- paste(sample.names1, sample.names2, sep = "_")
    
    # Archivos pareados y muestras únicas
    stopifnot(length(fnFs) == length(fnRs))
    stopifnot(!any(duplicated(sample.names)))
    
    ### Localización de bases de datos y hashes en metadata -------------------------
    # tax_location: ruta al fichero de referencia para assignTaxonomy
    tax_location <- file.path(metadata_location, primer.data$db_name[i])
    
    ### Preparar nombres para archivos filtrados ------------------------------------
    # Los archivos filtrados se escriben en subdirectorio "filtered" dentro de fastq_location
    filtFs <- file.path(fastq_location, "filtered", paste0(sample.names, "_F_filt.fastq"))
    filtRs <- file.path(fastq_location, "filtered", paste0(sample.names, "_R_filt.fastq"))
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names
    
    ### Filtrar muestras vacías -----------------------------------------------------
    # Algunos FASTQ pueden estar vacíos / corruptos; evitamos que rompan pasos posteriores.
    file.empty <- function(filenames) file.info(filenames)$size == 20
    empty_files <- file.empty(fnFs) | file.empty(fnRs)
    fnFs <- fnFs[!empty_files]; fnRs <- fnRs[!empty_files]
    filtFs <- filtFs[!empty_files]; filtRs <- filtRs[!empty_files]
    sample.names <- sample.names[!empty_files]
    
    ### Guardar perfiles de calidad (para inspección manual) ------------------------
    # Guardamos PNGs en logs para que los alumnos vean las curvas de calidad por run/locus
    png(filename=file.path(output_location, "logs", paste0(run_name,"_",primer.data$locus_shorthand[i],"_forward.png")))
    plotQualityProfile(fnFs[seq_len(min(10, length(fnFs)))]) # graficar hasta 10 archivos si existen
    dev.off()
    png(filename=file.path(output_location, "logs", paste0(run_name,"_",primer.data$locus_shorthand[i],"_reverse.png")))
    plotQualityProfile(fnRs[seq_len(min(10, length(fnRs)))])
    dev.off()
    
  ### Calcular longitud de truncamiento basada en calidad --------------------------
    # Se analiza calidad por ciclo y se aplica una ventana móvil; se usa la mediana entre muestras.
    message("Calculating quality trimming length... ", Sys.time())
    n <- 500000
    window_size <- 10
    trimsF <- numeric(length(fnFs))
    names(trimsF) <- fnFs
    for(f in seq_along(fnFs)) {
      srqa <- qa(fnFs[f], n = n)
      df <- srqa[["perCycle"]]$quality
      means <- rowsum(df$Score * df$Count, df$Cycle) / rowsum(df$Count, df$Cycle)
      window_values <- vapply(1:(length(means)-window_size), function(j) mean(means[j:(j+window_size)]), numeric(1))
      where_to_cut_F <- min(which(window_values < 30))
      trimsF[f] <- where_to_cut_F
    }
    where_trim_all_Fs <- median(trimsF, na.rm = TRUE) # truncLen forward global (mediana)
    
    trimsR <- numeric(length(fnRs))
    names(trimsR) <- fnRs
    for(r in seq_along(fnRs)) {
      srqa <- qa(fnRs[r], n = n)
      df <- srqa[["perCycle"]]$quality
      means <- rowsum(df$Score * df$Count, df$Cycle) / rowsum(df$Count, df$Cycle)
      window_values <- vapply(1:(length(means)-window_size), function(j) mean(means[j:(j+window_size)]), numeric(1))
      where_to_cut_R <- min(which(window_values < 30))
      trimsR[r] <- where_to_cut_R
    }
    where_trim_all_Rs <- median(trimsR, na.rm = TRUE)
    
    # Limitar truncamiento al máximo permitido (evita cortar demasiado y perder overlap)
    if(!is.na(where_trim_all_Fs) && where_trim_all_Fs > primer.data$max_trim[i]) {
      message("ADVERTENCIA: truncLen_F calculado (", where_trim_all_Fs, 
              ") excede max_trim (", primer.data$max_trim[i], "). Usando max_trim.")
      where_trim_all_Fs <- primer.data$max_trim[i]
    }
    if(!is.na(where_trim_all_Rs) && where_trim_all_Rs > primer.data$max_trim[i]) {
      message("ADVERTENCIA: truncLen_R calculado (", where_trim_all_Rs, 
              ") excede max_trim (", primer.data$max_trim[i], "). Usando max_trim.")
      where_trim_all_Rs <- primer.data$max_trim[i]
    }
    
    message("Puntos de truncamiento finales: F=", where_trim_all_Fs, " | R=", where_trim_all_Rs)
    message("Finished calculating quality trimming length... ", Sys.time())
    
    ### Distribución de longitudes de lecturas (OPCIONAL - diagnóstico) -------------
    # Este paso es útil para detectar problemas de secuenciación, pero no afecta el pipeline.
    # Puedes comentarlo si tienes muchas muestras y quieres acelerar el análisis.
    
    PLOT_RAW_LENGTHS <- TRUE  # Cambiar a FALSE para saltar este paso
    
    if (PLOT_RAW_LENGTHS) {
      message("Calculando distribución de longitudes crudas... ", Sys.time())
      
      sample_n_reads <- 100000
      get_lengths <- function(fq_file, n = sample_n_reads) {
        streamer <- ShortRead::FastqStreamer(fq_file, n = n)
        on.exit(close(streamer))
        chunk <- ShortRead::yield(streamer)
        if (length(chunk) == 0) return(integer(0))
        ShortRead::width(ShortRead::sread(chunk))
      }
      
      lengths_F <- purrr::map_df(fnFs, ~ {
        L <- get_lengths(.x)
        tibble::tibble(ReadType = "Forward",
                       File = basename(.x),
                       Length = L)
      })
      lengths_R <- purrr::map_df(fnRs, ~ {
        L <- get_lengths(.x)
        tibble::tibble(ReadType = "Reverse",
                       File = basename(.x),
                       Length = L)
      })
      lengths_all <- dplyr::bind_rows(lengths_F, lengths_R)
      
      png(file.path(output_location, "logs",
                    paste0(run_name,"_",primer.data$locus_shorthand[i],"_raw_read_lengths.png")),
          width = 2000, height = 1200, res=300)
      print(
        ggplot2::ggplot(lengths_all, ggplot2::aes(x = Length, fill = ReadType)) +
          ggplot2::geom_histogram(alpha = 0.6, position = "identity", binwidth = 1, color = "black") +
          ggplot2::facet_wrap(~ ReadType, ncol = 1, scales = "free_y") +
          ggplot2::labs(title = "Distribución de longitudes de lecturas crudas (primeras 100k por archivo)",
                        x = "Longitud (bp)", y = "Frecuencia") +
          ggplot2::theme_minimal()
      )
      dev.off()
      message("Distribución de longitudes crudas guardada. ", Sys.time())
    } else {
      message("Saltando análisis de longitudes crudas (PLOT_RAW_LENGTHS = FALSE)")
    }
    
    ### Filtrado y trimming con dada2 ------------------------------------------------
    message("Starting filter and trim... ", Sys.time())
    message("Parámetros de filtrado:")
    message("  - truncLen: F=", where_trim_all_Fs, ", R=", where_trim_all_Rs)
    message("  - maxEE: 2,2 (máximo errores esperados)")
    message("  - maxN: 0 (no permite bases ambiguas)")
    message("  - truncQ: 2 (trunca en calidad < 2)")
    message("  - rm.phix: TRUE (remueve contaminación PhiX)")
    
    # Validar que los puntos de truncamiento son razonables
    if (is.na(where_trim_all_Fs) || where_trim_all_Fs < 50) {
      stop("truncLen_F muy bajo o NA (", where_trim_all_Fs, "). Revisa calidad de Forward.")
    }
    if (is.na(where_trim_all_Rs) || where_trim_all_Rs < 50) {
      stop("truncLen_R muy bajo o NA (", where_trim_all_Rs, "). Revisa calidad de Reverse.")
    }
    
    # Ejecutar filterAndTrim con manejo de errores
    out <- tryCatch({
      filterAndTrim(
        fnFs, filtFs, fnRs, filtRs, 
        truncLen = c(where_trim_all_Fs, where_trim_all_Rs),
        maxN = 0,           # elimina lecturas con Ns
        maxEE = c(2, 2),    # máximo 2 errores esperados por lectura
        truncQ = 2,         # trunca lecturas con Q < 2
        rm.phix = TRUE,     # remueve PhiX (control Illumina)
        compress = TRUE,    # comprime FASTQs filtrados (ahorra espacio)
        multithread = TRUE, # paraleliza por núcleo
        matchIDs = TRUE     # asegura pareado F/R por ID
      )
    }, error = function(e) {
      stop("Error en filterAndTrim: ", conditionMessage(e), 
           "\nRevisa que los archivos FASTQ existan y no estén corruptos.")
    })
    
    message("Finished filter and trim. ", Sys.time())
    
    # Estadísticas de filtrado
    total_input <- sum(out[, "reads.in"])
    total_output <- sum(out[, "reads.out"])
    pct_retained <- round(total_output / total_input * 100, 2)
    
    message("\n=== Resumen de filtrado ===")
    message("Total lecturas input:  ", format(total_input, big.mark = ","))
    message("Total lecturas output: ", format(total_output, big.mark = ","))
    message("Retención global:      ", pct_retained, "%")
    message("Muestras procesadas:   ", nrow(out))
    
    # Advertir si la retención es muy baja
    if (pct_retained < 50) {
      warning("Retención < 50%! Considera relajar parámetros (maxEE, truncLen) o revisar calidad.")
    }
    
    # Identificar muestras con retención muy baja (< 30%)
    low_retention <- out[, "reads.out"] / out[, "reads.in"] < 0.3
    if (any(low_retention)) {
      problem_samples <- rownames(out)[low_retention]
      warning("Muestras con retención < 30%: ", paste(problem_samples, collapse = ", "))
    }
    
    # Guardar tabla de filtrado para inspección
    write.csv(out, 
              file.path(output_location, "logs", 
                        paste0(run_name, "_", primer.data$locus_shorthand[i], "_filterAndTrim_stats.csv")),
              row.names = TRUE)
    message("Estadísticas de filtrado guardadas en logs/")
    
    ### Dereplicación ----------------------------------------------------------------
    # La dereplicación identifica secuencias únicas en cada muestra y mantiene un
    # registro de su abundancia. Esto reduce redundancia y optimiza pasos posteriores.
    # Es un paso crítico que NO elimina información biológica.
    
    message("Iniciando dereplicación... ", Sys.time())
    
    # Conservamos solo archivos filtrados que realmente existen
    # (protege contra fallos en filterAndTrim que no hayan creado ciertos archivos)
    exists <- file.exists(filtFs) & file.exists(filtRs)
    
    if (!any(exists)) {
      stop("Ningún archivo filtrado encontrado. Revisa que filterAndTrim haya funcionado correctamente.")
    }
    
    if (sum(!exists) > 0) {
      message("ADVERTENCIA: ", sum(!exists), " muestras sin archivos filtrados (saltadas).")
      dropped_samples <- sample.names[!exists]
      message("Muestras eliminadas: ", paste(dropped_samples, collapse = ", "))
    }
    
    # Filtrar vectores para mantener solo archivos existentes
    filtFs <- filtFs[exists]
    filtRs <- filtRs[exists]
    sample.names <- sample.names[exists]
    
    # Ejecutar dereplicación con verbose para seguimiento
    message("Dereplicando ", length(filtFs), " muestras Forward...")
    derepFs <- derepFastq(filtFs, verbose = TRUE)
    
    message("Dereplicando ", length(filtRs), " muestras Reverse...")
    derepRs <- derepFastq(filtRs, verbose = TRUE)
    
    # Asignar nombres de muestra a los objetos dereplicados
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names
    
    # Diagnóstico rápido: inspeccionar primera muestra
    if (length(derepFs) > 0) {
      n_unique_F <- length(derepFs[[1]]$uniques)
      top5_F <- head(sort(derepFs[[1]]$uniques, decreasing = TRUE), 5)
      lengths_F <- unique(nchar(names(derepFs[[1]]$uniques)))
      
      message("\n=== Diagnóstico de dereplicación (muestra: ", sample.names[1], ") ===")
      message("Secuencias únicas Forward: ", format(n_unique_F, big.mark = ","))
      message("Top 5 abundancias: ", paste(top5_F, collapse = ", "))
      message("Longitudes observadas: ", paste(lengths_F, collapse = ", "), " bp")
      
      # Calcular ratio de compresión (útil para entender diversidad)
      total_reads_F <- sum(derepFs[[1]]$uniques)
      compression_ratio <- round(n_unique_F / total_reads_F * 100, 2)
      message("Ratio de compresión: ", compression_ratio, "% (", 
              n_unique_F, " únicas / ", format(total_reads_F, big.mark = ","), " lecturas)")
      
      # Advertir si hay poca diversidad (pocas secuencias únicas)
      if (compression_ratio < 1) {
        message("INFO: Baja diversidad detectada (< 1% secuencias únicas). Esto es normal en muestras con especies dominantes.")
      }
      if (compression_ratio > 50) {
        warning("Alta diversidad inusual (> 50% secuencias únicas). Revisa calidad de filtrado o contaminación.")
      }
    } else {
      stop("No hay muestras dereplicadas. El pipeline no puede continuar.")
    }
    
    # Guardar estadísticas de dereplicación
    derep_stats <- data.frame(
      Sample = sample.names,
      Unique_seqs_F = sapply(derepFs, function(x) length(x$uniques)),
      Total_reads_F = sapply(derepFs, function(x) sum(x$uniques)),
      Unique_seqs_R = sapply(derepRs, function(x) length(x$uniques)),
      Total_reads_R = sapply(derepRs, function(x) sum(x$uniques)),
      stringsAsFactors = FALSE
    )
    
    derep_stats$Compression_ratio_F <- with(derep_stats, 
                                            round(Unique_seqs_F / Total_reads_F * 100, 2))
    derep_stats$Compression_ratio_R <- with(derep_stats, 
                                            round(Unique_seqs_R / Total_reads_R * 100, 2))
    
    write.csv(derep_stats,
              file.path(output_location, "logs",
                        paste0(run_name, "_", primer.data$locus_shorthand[i], "_derep_stats.csv")),
              row.names = FALSE)
    
    message("Estadísticas de dereplicación guardadas en logs/")
    message("Dereplicación completada. ", Sys.time())
    
    ### Aprendizaje de modelos de error (forzar convergencia) -----------------------
    # El modelo de error es crítico en dada2: aprende las tasas reales de error de
    # secuenciación (A->C, A->G, etc.) para distinguir variantes biológicas de errores técnicos.
    # Usamos parámetros robustos para forzar convergencia en datasets pequeños/medianos.
    
    message("\n=== Aprendiendo modelos de error ===")
    message("Inicio: ", Sys.time())
    
    # Parámetros para forzar convergencia
    nbases_target <- 1e8   # 100 millones de bases (aumentar si tienes > 50 muestras)
    max_consist   <- 25    # más iteraciones que default (10) para convergencia
    
    # Validación previa: asegurar que hay suficientes datos
    total_bases_F <- sum(sapply(filtFs, function(f) {
      if (file.exists(f)) file.info(f)$size else 0
    }))
    total_bases_R <- sum(sapply(filtRs, function(r) {
      if (file.exists(r)) file.info(r)$size else 0
    }))
    
    message("Datos disponibles para aprendizaje:")
    message("  - Forward: ~", format(total_bases_F / 1e6, digits = 2), " MB")
    message("  - Reverse: ~", format(total_bases_R / 1e6, digits = 2), " MB")
    
    if (total_bases_F < 1e6 || total_bases_R < 1e6) {
      warning("Datos muy limitados para aprender errores (< 1 MB). El modelo puede ser inestable.")
    }
    
    # Aprender errores Forward
    message("\nAprendiendo errores Forward...")
    message("  - Objetivo: ", format(nbases_target, scientific = FALSE), " bases")
    message("  - MAX_CONSIST: ", max_consist, " iteraciones")
    
    errF <- tryCatch({
      learnErrors(
        filtFs,
        nbases      = nbases_target,
        randomize   = TRUE,      # muestrea aleatoriamente (evita sesgos por orden)
        multithread = TRUE,      # paraleliza por núcleo
        MAX_CONSIST = max_consist,
        verbose     = TRUE
      )
    }, error = function(e) {
      stop("Error al aprender modelo Forward: ", conditionMessage(e),
           "\nPosibles causas: datos insuficientes, archivos corruptos, o calidad muy baja.")
    })
    
    # Verificar que el modelo se creó correctamente
    if (!is.null(errF)) {
      message("Modelo de error Forward: ✓ Creado exitosamente")
      
      # Verificar estructura básica del modelo
      if (!is.null(errF$trans) && !is.null(errF$err_out)) {
        n_transitions <- length(unique(errF$trans))
        message("  - Transiciones aprendidas: ", n_transitions)
        message("  - Rango de tasas de error: ", 
                round(min(errF$err_out, na.rm = TRUE), 5), " - ", 
                round(max(errF$err_out, na.rm = TRUE), 5))
        
        # Validación: tasas de error razonables (0.0001 - 0.1)
        if (any(errF$err_out > 0.1, na.rm = TRUE)) {
          warning("Tasas de error inusualmente altas (> 10%). Revisa calidad de datos.")
        }
      } else {
        warning("Modelo Forward tiene estructura incompleta. Revisa mensajes de learnErrors.")
      }
    } else {
      stop("Modelo de error Forward es NULL. Revisa que haya datos suficientes.")
    }
    
    # Aprender errores Reverse
    message("\nAprendiendo errores Reverse...")
    message("  - Objetivo: ", format(nbases_target, scientific = FALSE), " bases")
    message("  - MAX_CONSIST: ", max_consist, " iteraciones")
    
    errR <- tryCatch({
      learnErrors(
        filtRs,
        nbases      = nbases_target,
        randomize   = TRUE,
        multithread = TRUE,
        MAX_CONSIST = max_consist,
        verbose     = TRUE
      )
    }, error = function(e) {
      stop("Error al aprender modelo Reverse: ", conditionMessage(e),
           "\nPosibles causas: datos insuficientes, archivos corruptos, o calidad muy baja.")
    })
    
    # Verificar modelo Reverse
    if (!is.null(errR)) {
      message("Modelo de error Reverse: ✓ Creado exitosamente")
      
      if (!is.null(errR$trans) && !is.null(errR$err_out)) {
        n_transitions <- length(unique(errR$trans))
        message("  - Transiciones aprendidas: ", n_transitions)
        message("  - Rango de tasas de error: ", 
                round(min(errR$err_out, na.rm = TRUE), 5), " - ", 
                round(max(errR$err_out, na.rm = TRUE), 5))
        
        if (any(errR$err_out > 0.1, na.rm = TRUE)) {
          warning("Tasas de error inusualmente altas (> 10%). Revisa calidad de datos.")
        }
      } else {
        warning("Modelo Reverse tiene estructura incompleta. Revisa mensajes de learnErrors.")
      }
    } else {
      stop("Modelo de error Reverse es NULL. Revisa que haya datos suficientes.")
    }
    
    message("\nAprendizaje de errores completado. ", Sys.time())
    
    # Guardar plots de diagnóstico (líneas negras = tasas esperadas, puntos rojos = observadas)
    message("Generando gráficos de modelos de error...")
    
    dir.create(file.path(output_location, "logs"), recursive = TRUE, showWarnings = FALSE)
    
    png(file.path(output_location, "logs", 
                  paste0(run_name, "_", primer.data$locus_shorthand[i], "_errF.png")), 
        width = 1200, height = 1200, res = 150)
    tryCatch({
      print(plotErrors(errF, nominalQ = TRUE))
    }, error = function(e) {
      warning("No se pudo generar gráfico de errores Forward: ", conditionMessage(e))
    })
    dev.off()
    
    png(file.path(output_location, "logs", 
                  paste0(run_name, "_", primer.data$locus_shorthand[i], "_errR.png")), 
        width = 1200, height = 1200, res = 150)
    tryCatch({
      print(plotErrors(errR, nominalQ = TRUE))
    }, error = function(e) {
      warning("No se pudo generar gráfico de errores Reverse: ", conditionMessage(e))
    })
    dev.off()
    
    message("Gráficos de modelos de error guardados en logs/")
    
    # Diagnóstico rápido: tasas de error promedio por transición
    # errF$err_out es una matriz [16 transiciones × 41 Q-scores]
    # Los nombres de transiciones están en rownames, no en errF$trans
    
    if (!is.null(errF) && is.matrix(errF$err_out)) {
      # Obtener nombres de transiciones (filas de la matriz)
      transition_names <- rownames(errF$err_out)
      
      if (is.null(transition_names)) {
        warning("errF$err_out no tiene nombres de fila. No se puede calcular resumen de errores.")
        err_summary_F <- data.frame(
          Direction = "Forward",
          A2C = NA, A2G = NA, A2T = NA, C2A = NA,
          stringsAsFactors = FALSE
        )
      } else {
        # Calcular promedio por fila (transición) excluyendo transiciones correctas (A2A, C2C, etc.)
        # Promediamos a través de todos los Q-scores (columnas)
        err_summary_F <- data.frame(
          Direction = "Forward",
          A2C = if("A2C" %in% transition_names) mean(errF$err_out["A2C", ], na.rm = TRUE) else NA,
          A2G = if("A2G" %in% transition_names) mean(errF$err_out["A2G", ], na.rm = TRUE) else NA,
          A2T = if("A2T" %in% transition_names) mean(errF$err_out["A2T", ], na.rm = TRUE) else NA,
          C2A = if("C2A" %in% transition_names) mean(errF$err_out["C2A", ], na.rm = TRUE) else NA,
          C2G = if("C2G" %in% transition_names) mean(errF$err_out["C2G", ], na.rm = TRUE) else NA,
          C2T = if("C2T" %in% transition_names) mean(errF$err_out["C2T", ], na.rm = TRUE) else NA,
          G2A = if("G2A" %in% transition_names) mean(errF$err_out["G2A", ], na.rm = TRUE) else NA,
          G2C = if("G2C" %in% transition_names) mean(errF$err_out["G2C", ], na.rm = TRUE) else NA,
          G2T = if("G2T" %in% transition_names) mean(errF$err_out["G2T", ], na.rm = TRUE) else NA,
          T2A = if("T2A" %in% transition_names) mean(errF$err_out["T2A", ], na.rm = TRUE) else NA,
          T2C = if("T2C" %in% transition_names) mean(errF$err_out["T2C", ], na.rm = TRUE) else NA,
          T2G = if("T2G" %in% transition_names) mean(errF$err_out["T2G", ], na.rm = TRUE) else NA,
          stringsAsFactors = FALSE
        )
        
        message("Transiciones Forward disponibles: ", paste(transition_names, collapse = ", "))
      }
    } else {
      err_summary_F <- data.frame(
        Direction = "Forward",
        A2C = NA, A2G = NA, A2T = NA, C2A = NA,
        stringsAsFactors = FALSE
      )
      warning("errF$err_out no es una matriz. No se puede calcular resumen de errores Forward.")
    }
    
    if (!is.null(errR) && is.matrix(errR$err_out)) {
      transition_names <- rownames(errR$err_out)
      
      if (is.null(transition_names)) {
        warning("errR$err_out no tiene nombres de fila. No se puede calcular resumen de errores.")
        err_summary_R <- data.frame(
          Direction = "Reverse",
          A2C = NA, A2G = NA, A2T = NA, C2A = NA,
          stringsAsFactors = FALSE
        )
      } else {
        err_summary_R <- data.frame(
          Direction = "Reverse",
          A2C = if("A2C" %in% transition_names) mean(errR$err_out["A2C", ], na.rm = TRUE) else NA,
          A2G = if("A2G" %in% transition_names) mean(errR$err_out["A2G", ], na.rm = TRUE) else NA,
          A2T = if("A2T" %in% transition_names) mean(errR$err_out["A2T", ], na.rm = TRUE) else NA,
          C2A = if("C2A" %in% transition_names) mean(errR$err_out["C2A", ], na.rm = TRUE) else NA,
          C2G = if("C2G" %in% transition_names) mean(errR$err_out["C2G", ], na.rm = TRUE) else NA,
          C2T = if("C2T" %in% transition_names) mean(errR$err_out["C2T", ], na.rm = TRUE) else NA,
          G2A = if("G2A" %in% transition_names) mean(errR$err_out["G2A", ], na.rm = TRUE) else NA,
          G2C = if("G2C" %in% transition_names) mean(errR$err_out["G2C", ], na.rm = TRUE) else NA,
          G2T = if("G2T" %in% transition_names) mean(errR$err_out["G2T", ], na.rm = TRUE) else NA,
          T2A = if("T2A" %in% transition_names) mean(errR$err_out["T2A", ], na.rm = TRUE) else NA,
          T2C = if("T2C" %in% transition_names) mean(errR$err_out["T2C", ], na.rm = TRUE) else NA,
          T2G = if("T2G" %in% transition_names) mean(errR$err_out["T2G", ], na.rm = TRUE) else NA,
          stringsAsFactors = FALSE
        )
        
        message("Transiciones Reverse disponibles: ", paste(transition_names, collapse = ", "))
      }
    } else {
      err_summary_R <- data.frame(
        Direction = "Reverse",
        A2C = NA, A2G = NA, A2T = NA, C2A = NA,
        stringsAsFactors = FALSE
      )
      warning("errR$err_out no es una matriz. No se puede calcular resumen de errores Reverse.")
    }
    
    err_summary <- rbind(err_summary_F, err_summary_R)
    write.csv(err_summary,
              file.path(output_location, "logs",
                        paste0(run_name, "_", primer.data$locus_shorthand[i], "_error_rates_summary.csv")),
              row.names = FALSE)
    
    message("\n=== Resumen de tasas de error (promedio) ===")
    print(err_summary)
    message("CSV guardado en logs/")
    
    
    ### Inferencia de muestras -------------------------------------------------------
    # Inferimos ASVs usando las tasas de error aprendidas
    dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
    dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
    
    ### Unión de pares (merge) ------------------------------------------------------
    
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 12, verbose = TRUE)
    
    ### Construir tabla de secuencias (ASV table) -----------------------------------
    seqtab <- makeSequenceTable(mergers)
    
    ### Eliminación de quimeras ------------------------------------------------------
    seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
    freq.nochim <- sum(seqtab.nochim) / sum(seqtab) # proporción de lecturas no quiméricas
    
    message("Frecuencia de lecturas no quiméricas: ", round(freq.nochim * 100, 2), "%")
    message("ASVs totales después de remover quimeras: ", ncol(seqtab.nochim))
    
    ### Visualizar distribución de longitudes ANTES del filtro ------------------------
    seq_lengths <- nchar(colnames(seqtab.nochim))
    median_len <- median(seq_lengths)
    sd_len <- sd(seq_lengths)
    
    message("Longitud mediana de ASVs: ", median_len, " bp")
    message("Desviación estándar: ", round(sd_len, 2), " bp")
    
    png(file.path(output_location, "logs",
                  paste0(run_name, "_", primer.data$locus_shorthand[i], "_seq_lengths_distribution.png")),
        width = 1800, height = 1200, res = 300)
    print(
      ggplot2::ggplot(data = data.frame(Length = seq_lengths), ggplot2::aes(x = Length)) +
        ggplot2::geom_histogram(binwidth = 1, fill = "#69b3a2", color = "#2d5f4f", alpha = 0.8) +
        ggplot2::geom_vline(aes(xintercept = median_len), 
                            color = "#e63946", linetype = "dashed", linewidth = 1) +
        ggplot2::geom_vline(aes(xintercept = median_len - sd_len), 
                            color = "#f77f00", linetype = "dotted", linewidth = 0.8) +
        ggplot2::geom_vline(aes(xintercept = median_len + sd_len), 
                            color = "#f77f00", linetype = "dotted", linewidth = 0.8) +
        ggplot2::annotate("text", x = median_len, y = Inf, 
                          label = paste0("Mediana: ", round(median_len, 1), " bp"), 
                          vjust = 1.5, hjust = -0.1, color = "#e63946", fontface = "bold") +
        ggplot2::annotate("text", x = median_len + sd_len, y = Inf, 
                          label = paste0("±SD: ", round(sd_len, 1), " bp"), 
                          vjust = 3, hjust = -0.1, color = "#f77f00", size = 3.5) +
        ggplot2::labs(
          title = "Distribución de Longitudes de Secuencias ASV (Sin Filtrar)",
          subtitle = paste0("n = ", length(seq_lengths), " secuencias"),
          x = "Longitud de Secuencia (bp)", 
          y = "Frecuencia"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()
        )
    )
    dev.off()
    message("Gráfico de distribución de longitudes (sin filtrar) guardado. ", Sys.time())
    
    ### Filtrado por longitud --------------------------------------------------------
    min_len <- median_len - sd_len 
    max_len <- median_len + sd_len
    
    message("Rango de longitudes aceptado: ", round(min_len, 1), " - ", round(max_len, 1), " bp")
    
    indexes.to.keep <- which((nchar(colnames(seqtab.nochim)) <= max_len) & 
                               (nchar(colnames(seqtab.nochim)) >= min_len))
    cleaned.seqtab.nochim <- seqtab.nochim[, indexes.to.keep, drop = FALSE]
    filteredout.seqtab.nochim <- seqtab.nochim[, -indexes.to.keep, drop = FALSE]
    
    message("ASVs antes del filtro: ", ncol(seqtab.nochim))
    message("ASVs después del filtro: ", ncol(cleaned.seqtab.nochim))
    message("ASVs filtradas por longitud: ", ncol(filteredout.seqtab.nochim))
    
    # Guardar ASVs filtradas para inspección
    if (ncol(filteredout.seqtab.nochim) > 0) {
      write.csv(filteredout.seqtab.nochim, 
                file.path(output_location, "logs", 
                          paste0(run_name, "_", primer.data$locus_shorthand[i], "_filtered_out_asv.csv")))
      message("ASVs filtradas guardadas en logs/")
    }
    
    ### Visualizar distribución de longitudes DESPUÉS del filtro ---------------------
    if (ncol(cleaned.seqtab.nochim) > 0) {
      seq_lengths_clean <- nchar(colnames(cleaned.seqtab.nochim))
      median_len_clean <- median(seq_lengths_clean)
      sd_len_clean <- sd(seq_lengths_clean)
      
      png(file.path(output_location, "logs",
                    paste0(run_name, "_", primer.data$locus_shorthand[i], "_cleaned_seq_lengths_distribution.png")),
          width = 1800, height = 1200, res = 300)
      print(
        ggplot2::ggplot(data = data.frame(Length = seq_lengths_clean), ggplot2::aes(x = Length)) +
          ggplot2::geom_histogram(binwidth = 1, fill = "#06d6a0", color = "#118ab2", alpha = 0.8) +
          ggplot2::geom_vline(aes(xintercept = median_len_clean), 
                              color = "#ef476f", linetype = "dashed", linewidth = 1) +
          ggplot2::geom_vline(aes(xintercept = median_len_clean - sd_len_clean), 
                              color = "#ffd166", linetype = "dotted", linewidth = 0.8) +
          ggplot2::geom_vline(aes(xintercept = median_len_clean + sd_len_clean), 
                              color = "#ffd166", linetype = "dotted", linewidth = 0.8) +
          ggplot2::annotate("text", x = median_len_clean, y = Inf, 
                            label = paste0("Mediana: ", round(median_len_clean, 1), " bp"), 
                            vjust = 1.5, hjust = -0.1, color = "#ef476f", fontface = "bold") +
          ggplot2::annotate("text", x = median_len_clean + sd_len_clean, y = Inf, 
                            label = paste0("±SD: ", round(sd_len_clean, 1), " bp"), 
                            vjust = 3, hjust = -0.1, color = "#ffd166", size = 3.5) +
          ggplot2::labs(
            title = "Distribución de Longitudes de Secuencias ASV (Filtradas)",
            subtitle = paste0("n = ", length(seq_lengths_clean), " secuencias | ",
                              ncol(filteredout.seqtab.nochim), " ASVs removidas"),
            x = "Longitud de Secuencia (bp)", 
            y = "Frecuencia"
          ) +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank()
          )
      )
      dev.off()
      message("Gráfico de distribución de longitudes (filtradas) guardado. ", Sys.time())
    } else {
      warning("No quedaron ASVs después del filtro por longitud!")
    }
    
    ### Seguimiento de reads a través del pipeline -----------------------------------
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out[exists, ],
                   sapply(dadaFs, getN),
                   sapply(dadaRs, getN),
                   sapply(mergers, getN),
                   rowSums(seqtab.nochim),
                   rowSums(cleaned.seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "length_filter")
    rownames(track) <- sample.names
    
    # Guardar tabla de tracking
    write.csv(track, 
              file.path(output_location, "logs", 
                        paste0(run_name, "_", primer.data$locus_shorthand[i], "_read_tracking.csv")),
              row.names = TRUE)
    
    message("\n=== Resumen de lecturas por paso ===")
    print(head(track))
    message("\nEstadísticas de retención:")
    message("- Promedio input → filtered: ", round(mean(track[,"filtered"] / track[,"input"] * 100), 1), "%")
    message("- Promedio filtered → merged: ", round(mean(track[,"merged"] / track[,"filtered"] * 100), 1), "%")
    message("- Promedio merged → nonchim: ", round(mean(track[,"nonchim"] / track[,"merged"] * 100), 1), "%")
    message("- Promedio nonchim → length_filter: ", round(mean(track[,"length_filter"] / track[,"nonchim"] * 100), 1), "%")
    message("Total ASVs eliminadas por longitud: ", ncol(seqtab.nochim) - ncol(cleaned.seqtab.nochim))
    
    ### Definir rutas de salida ------------------------------------------------------
    ASV_file <- file.path(output_location, "csv_output", 
                          paste0(run_name, "_", primer.data$locus_shorthand[i], "_ASV_table.csv"))
    taxonomy_file <- file.path(output_location, "csv_output", 
                               paste0(run_name, "_", primer.data$locus_shorthand[i], "_taxonomy_output.csv"))
    bootstrap_file <- file.path(output_location, "csv_output", 
                                paste0(run_name, "_", primer.data$locus_shorthand[i], "_tax_bootstrap.csv"))
    
    dir.create(dirname(ASV_file), recursive = TRUE, showWarnings = FALSE)
    
    ### Generar tabla ASV en formato largo (usando cleaned) --------------------------
    current_asv <- as.data.frame(cleaned.seqtab.nochim) %>%
      tibble::rownames_to_column("Sample_name") %>%
      tidyr::pivot_longer(cols = -Sample_name, names_to = "Sequence", values_to = "nReads") %>%
      dplyr::filter(nReads > 0) %>%
      dplyr::mutate(Label = primer.data$locus_shorthand[i]) %>%
      dplyr::relocate(Label, .after = Sample_name)
    
    readr::write_csv(current_asv, ASV_file)
    message("Tabla ASV (filtrada por longitud) guardada en: ", ASV_file)
    
    ### Asignación de taxonomía a ASVs ---------------------------------------
    dir.create(dirname(taxonomy_file), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(bootstrap_file), recursive = TRUE, showWarnings = FALSE)
    
    if (!file.exists(tax_location)) {
      stop("No se encuentra la base de taxonomía: ", tax_location)
    }
    
    if (ncol(seqtab.nochim) == 0) {
      warning("No hay ASVs para clasificar; se escribirán CSV vacíos.")
      joined_old_new_taxa <- tibble::tibble(
        Sequence = character(),
        Kingdom = character(), Phylum = character(), Class = character(),
        Order = character(), Family = character(), Genus = character(),
        Species = character()
      )
      ready_new_taxa <- tibble::tibble(Sequence = character())
      readr::write_csv(joined_old_new_taxa, taxonomy_file)
      readr::write_csv(ready_new_taxa, bootstrap_file)
    } else {
      seqs <- colnames(cleaned.seqtab.nochim)
      
      message("Asignando taxonomía a ", length(seqs), " ASVs ...")
      at <- tryCatch(
        assignTaxonomy(seqs, tax_location,
                       tryRC = TRUE, multithread = TRUE,
                       outputBootstraps = TRUE, verbose = TRUE),
        error = function(e) stop("Fallo assignTaxonomy: ", conditionMessage(e))
      )
      
      # Normalizar salida (lista con $tax/$boot o matriz)
      if (is.list(at) && !is.null(at$tax)) {
        tax_mat  <- at$tax
        boot_mat <- at$boot
      } else {
        tax_mat  <- at
        boot_mat <- NULL
      }
      
      tax_df <- as.data.frame(tax_mat) %>%
        tibble::rownames_to_column("Sequence")
      # Quitar prefijo "tax." si viniera
      nc <- setdiff(names(tax_df), "Sequence")
      names(tax_df)[match(nc, names(tax_df))] <- sub("^tax\\.", "", nc)
      
      # Reordenar/asegurar columnas estándar
      ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
      for (rk in ranks) if (!rk %in% names(tax_df)) tax_df[[rk]] <- NA_character_
      joined_old_new_taxa <- tax_df %>%
        dplyr::select(dplyr::any_of(c("Sequence", ranks)))
      
      # Bootstraps
      if (!is.null(boot_mat)) {
        boot_df <- as.data.frame(boot_mat) %>%
          tibble::rownames_to_column("Sequence")
        names(boot_df)[-1] <- paste0("boot.", names(boot_df)[-1])
        
        # Unir con taxonomía completa (todas las columnas taxonómicas)
        tax_for_boot <- joined_old_new_taxa %>%
          dplyr::select(Sequence, Kingdom, Phylum, Class, Order, Family, Genus, Species)
        
        # Unir bootstraps y taxonomía por Sequence
        ready_new_taxa <- dplyr::left_join(boot_df, tax_for_boot, by = "Sequence")
      } else {
        ready_new_taxa <- tibble::tibble(Sequence = joined_old_new_taxa$Sequence)
      }
      
      # Chequeos rápidos
      if (any(duplicated(joined_old_new_taxa$Sequence))) {
        warning("Secuencias duplicadas en tabla taxonómica.")
      }
      pct_unassigned <- vapply(ranks, function(rk) {
        mean(is.na(joined_old_new_taxa[[rk]]) | joined_old_new_taxa[[rk]] == "")
      }, numeric(1))
      readr::write_csv(
        tibble::tibble(level = ranks, pct_unassigned = pct_unassigned),
        file.path(output_location, "logs",
                  paste0(run_name, "_", primer.data$locus_shorthand[i], "_pct_unassigned.csv"))
      )
      
      # Guardar CSVs solicitados
      readr::write_csv(joined_old_new_taxa, taxonomy_file)
      readr::write_csv(ready_new_taxa, bootstrap_file)
      message("Taxonomía escrita en:\n- ", taxonomy_file, "\n- ", bootstrap_file)
    }
    
    ### Cargar sample_meta desde metadata/metadata.csv (usa sample_id) --------------
    meta_path <- if (exists("metadata_location")) {
      file.path(metadata_location, "metadata.csv")
    } else {
      here::here("metadata", "metadata.csv")
    }
    if (!file.exists(meta_path)) stop("No se encontró metadata.csv en: ", meta_path)
    
    sample_meta <- readr::read_csv(meta_path, show_col_types = FALSE)
    
    # Detectar columna ID
    id_col <- dplyr::case_when(
      "sample_id"   %in% names(sample_meta) ~ "sample_id",
      "Sample_name" %in% names(sample_meta) ~ "Sample_name",
      TRUE ~ NA_character_
    )
    if (is.na(id_col)) {
      stop("metadata.csv debe tener la columna 'sample_id' (o 'Sample_name').")
    }
    
    # Normalizar IDs y poner como rownames
    sample_meta <- sample_meta %>%
      dplyr::mutate(!!id_col := trimws(as.character(.data[[id_col]]))) %>%
      dplyr::rename(SampleID = !!id_col) %>%
      tibble::column_to_rownames("SampleID")
    
    # Alinear con la matriz de abundancias (otu_mat) y validar
    otu_mat <- if (exists("cleaned.seqtab.nochim")) cleaned.seqtab.nochim else seqtab.nochim
    storage.mode(otu_mat) <- "integer"
    
    smp_in_otu <- rownames(otu_mat)
    faltan <- setdiff(smp_in_otu, rownames(sample_meta))
    sobrantes <- setdiff(rownames(sample_meta), smp_in_otu)
    if (length(faltan))   warning("Faltan en metadata (sample_id): ", paste(faltan, collapse = ", "))
    if (length(sobrantes)) message("Sobrantes en metadata (no usados): ", paste(sobrantes, collapse = ", "))
    
    sample_meta <- sample_meta[smp_in_otu, , drop = FALSE]
    
    # Construir tax_mat desde joined_old_new_taxa
    tax_mat <- joined_old_new_taxa %>%
      dplyr::mutate(Sequence = as.character(Sequence)) %>%
      tibble::column_to_rownames("Sequence") %>%
      as.matrix()
    
    # Guardar insumos para phyloseq (+ track y tablas) en un único .RData
    rdata_path <- file.path(output_location, "rdata_output",
                            paste0(run_name,"_", primer.data$locus_shorthand[i], "_phyloseq_inputs.RData"))
    dir.create(dirname(rdata_path), recursive = TRUE, showWarnings = FALSE)
    
    # Asegurar tipos y coherencia
    if (!all(colnames(otu_mat) %in% rownames(tax_mat))) {
      warning("ASVs en otu_mat sin fila en tax_mat.")
    }
    storage.mode(otu_mat) <- "integer"
    if (!all(rownames(otu_mat) %in% rownames(sample_meta))) {
      warning("Muestras en otu_mat sin metadata.")
    }
    sample_meta <- sample_meta[rownames(otu_mat), , drop = FALSE]
    
    track        <- track      # ya creado
    seqtab_raw   <- seqtab.nochim
    seqtab_clean <- cleaned.seqtab.nochim
    taxonomy_df  <- joined_old_new_taxa
    bootstrap_df <- ready_new_taxa
    params <- list(
      run_name   = run_name,
      locus      = primer.data$locus_shorthand[i],
      truncLen_F = where_trim_all_Fs,
      truncLen_R = where_trim_all_Rs,
      min_len    = min_len,
      max_len    = max_len,
      date       = as.character(Sys.time())
    )
    
    save(otu_mat, tax_mat, sample_meta, track,
         seqtab_raw, seqtab_clean,
         taxonomy_df, bootstrap_df, params,
         file = rdata_path, compress = "xz")
    
    message("Guardado .RData para phyloseq (con track): ", rdata_path)
    message("Fin del procesamiento para locus ", primer.data$locus_shorthand[i], ": ", Sys.time())    

