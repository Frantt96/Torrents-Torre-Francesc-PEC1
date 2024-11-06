library(SummarizedExperiment)
library(dplyr)
# Instalación BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# Instalacion del paquete POMA
BiocManager::install("POMA")
library(pheatmap)
library(POMA)
library(magrittr)

# Cargamos los datos de estudio

data <- read.csv("human_cachexia.csv", header = TRUE, stringsAsFactors = FALSE)  # Ajusta según tus datos

# Separamos los datos en matriz de expresión (metabolitos) y metadatos (ID de paciente y pérdida de músculo)
# Extraemos las columnas de metabolitos como matriz de expresión
expression_data <- as.matrix(data[, 3:ncol(data)])
rownames(expression_data) <- data$Patient.ID

# Creamos los metadatos para las filas y columnas
# - Filas: (ID y condición de pérdida muscular)
# - Columnas: nombres de los metabolitos
row_metadata <- data.frame(Patient_ID = data$Patient.ID, Muscle_loss = data$Muscle.loss)
col_metadata <- data.frame(Metabolite = colnames(expression_data))

# Creamos el objeto SummarizedExperiment
se <- SummarizedExperiment(assays = list(counts = expression_data),
                                  rowData = row_metadata,
                                  colData = col_metadata)

# Resumen del objeto para verificar
se


# Normalizamos los datos de se
normalized <- se %>% 
  PomaNorm(method = "log_pareto")

# Diferentes gráficos
PomaBoxplots(normalized)
PomaDensity(normalized)


# Grafico heatmap
library(pheatmap)
pheatmap(assay(normalized), scale = "row", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", main = "Heatmap de Expresión")

# Gráfico boxplot

boxplot(assay(normalized), main = "Boxplot de los Metabolitos", las = 2, cex.axis = 0.7)


# Gráfico correlación
cor_matrix <- cor(assay(normalized))
heatmap(cor_matrix, main = "Matriz de Correlación de Metabolitos")

