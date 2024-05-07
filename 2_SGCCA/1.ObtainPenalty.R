################### GANADOR = 0.44 RNAseq ########
RNAseqnormalized <- read.delim("~/Documents/INMEGEN/multiomics/RNAseqnormalized.tsv")
X <- RNAseqnormalized[, -ncol(RNAseqnormalized)]  # Selecciona todas las columnas excepto la última como características
y <- RNAseqnormalized[, ncol(RNAseqnormalized)]   # La última columna es la variable de respuesta
set.seed(123)  # Para reproducibilidad
train_indices <- sample(1:nrow(RNAseqnormalized), 0.8 * nrow(RNAseqnormalized))  # Índices para entrenamiento (80% de los datos)
X_train <- X[train_indices, ]
y_train <- y[train_indices]
X_test <- X[-train_indices, ]
y_test <- y[-train_indices]
#install.packages("glmnet")  # Si aún no está instalado
library(glmnet)

# Ajustar modelo Elastic Net con validación cruzada para seleccionar la mejor penalización
X_train <- as.matrix(X_train)
cv_model <- cv.glmnet(X_train, y_train, alpha = 0.5)  # alpha = 0.5 para Elastic Net

best_penalty <- cv_model$lambda.min  # Lambda correspondiente al mínimo error de validación cruzada
X_test<-as.matrix(X_test)
predictions <- predict(cv_model, newx = X_test, s = best_penalty)
lambda_selected <- 187.4545
lambda_min <- min(cv_model$lambda)  # Valor mínimo de lambda seleccionado por cv.glmnet
lambda_max <- max(cv_model$lambda)  # Valor máximo de lambda seleccionado por cv.glmnet

# Normalizar el valor de lambda en el rango de 0 a 1
penalty_normalized <- (log(lambda_selected) - log(lambda_min)) / (log(lambda_max) - log(lambda_min))
################ visualizando el mejor parametro

# Define una secuencia de tamaños de conjunto de datos de entrenamiento
train_sizes <- seq(100, nrow(X_train), by = 200)

# Almacena los errores cuadráticos medios para cada tamaño de conjunto de datos de entrenamiento
mse <- numeric(length(train_sizes))

# Ajusta el modelo y evalúa el rendimiento para cada tamaño de conjunto de datos de entrenamiento
for (i in seq_along(train_sizes)) {
  subset_X_train <- X_train[1:train_sizes[i], ]
  subset_y_train <- y_train[1:train_sizes[i]]
  
  # Ajusta el modelo Elastic Net
  model <- glmnet(subset_X_train, subset_y_train, alpha = 0.5, lambda = lambda_selected)
  
  # Realiza predicciones en el conjunto de prueba
  predictions <- predict(model, newx = X_test)
  
  # Calcula el error cuadrático medio
  mse[i] <- mean((predictions - y_test)^2)
}

# Grafica la curva de aprendizaje
plot(train_sizes, mse, type = "b", xlab = "Tamaño del conjunto de datos de entrenamiento", ylab = "Error cuadrático medio", main = "Curva de aprendizaje")
######################### optimo de Methy GANADOR 0.24 ############
library(data.table)

# Abre el archivo usando fread
methyM <- fread("~/Documents/INMEGEN/multiomics/methyM.tsv")
methyM=as.matrix(methyM[,2:ncol(methyM)],rownames=methyM$V1)
X <- methyM[, -ncol(methyM)]  # Selecciona todas las columnas excepto la última como características
y <- methyM[, ncol(methyM)]   # La última columna es la variable de respuesta
set.seed(123)  # Para reproducibilidad
train_indices <- sample(1:nrow(methyM), 0.8 * nrow(methyM))  # Índices para entrenamiento (80% de los datos)
X_train <- X[train_indices, ]
y_train <- y[train_indices]
X_test <- X[-train_indices, ]
y_test <- y[-train_indices]
#install.packages("glmnet")  # Si aún no está instalado
library(glmnet)

# Ajustar modelo Elastic Net con validación cruzada para seleccionar la mejor penalización
X_train <- as.matrix(X_train)
cv_model <- cv.glmnet(X_train, y_train, alpha = 0.5)  # alpha = 0.5 para Elastic Net

best_penalty <- cv_model$lambda.min  # Lambda correspondiente al mínimo error de validación cruzada
X_test<-as.matrix(X_test)
predictions <- predict(cv_model, newx = X_test, s = best_penalty)
lambda_selected <- 0.00578
lambda_min <- min(cv_model$lambda)  # Valor mínimo de lambda seleccionado por cv.glmnet
lambda_max <- max(cv_model$lambda)  # Valor máximo de lambda seleccionado por cv.glmnet

# Normalizar el valor de lambda en el rango de 0 a 1
penalty_normalized <- (log(lambda_selected) - log(lambda_min)) / (log(lambda_max) - log(lambda_min))
################ visualizando el mejor parametro

# Define una secuencia de tamaños de conjunto de datos de entrenamiento
train_sizes <- seq(100, nrow(X_train), by = 200)

# Almacena los errores cuadráticos medios para cada tamaño de conjunto de datos de entrenamiento
mse <- numeric(length(train_sizes))

# Ajusta el modelo y evalúa el rendimiento para cada tamaño de conjunto de datos de entrenamiento
for (i in seq_along(train_sizes)) {
  subset_X_train <- X_train[1:train_sizes[i], ]
  subset_y_train <- y_train[1:train_sizes[i]]
  
  # Ajusta el modelo Elastic Net
  model <- glmnet(subset_X_train, subset_y_train, alpha = 0.5, lambda = lambda_selected)
  
  # Realiza predicciones en el conjunto de prueba
  predictions <- predict(model, newx = X_test)
  
  # Calcula el error cuadrático medio
  mse[i] <- mean((predictions - y_test)^2)
}

# Grafica la curva de aprendizaje
plot(train_sizes, mse, type = "b", xlab = "Tamaño del conjunto de datos de entrenamiento", ylab = "Error cuadrático medio", main = "Curva de aprendizaje")

#########################PENALIZACION PARA miRNAS GANADOR 0.25 #####
## ENET no ajusta correctamente por lo que ridge es una mejor opcion o en menor medida, lasso.
miRNAseqNormi <- read.delim("~/Documents/INMEGEN/multiomics/miRNAseqNormi.tsv")
X <- miRNAseqNormi[, -ncol(miRNAseqNormi)]  # Selecciona todas las columnas excepto la última como características
y <- miRNAseqNormi[, ncol(miRNAseqNormi)]   # La última columna es la variable de respuesta
set.seed(123)  # Para reproducibilidad
train_indices <- sample(1:nrow(miRNAseqNormi), 0.8 * nrow(miRNAseqNormi))  # Índices para entrenamiento (80% de los datos)
X_train <- X[train_indices, ]
y_train <- y[train_indices]
X_test <- X[-train_indices, ]
y_test <- y[-train_indices]
#install.packages("glmnet")  # Si aún no está instalado
library(glmnet)
X_test<- as.matrix(X_test)

## con ridge
X_train <- as.matrix(X_train)
alpha0 <- cv.glmnet(X_train, y_train, alpha = 0, type.measure = "mse")  # alpha = 0 para ridge
alpha0.predict<- predict(alpha0,s=alpha0$lambda.1se,newx = X_test)
mean((y_test-alpha0.predict)^2)
#[1] 491197087

## con lasso
alpha1 <- cv.glmnet(X_train, y_train, alpha = 1, type.measure = "mse")  # alpha = 1 para lasso
alpha1.predict<- predict(alpha1,s=alpha1$lambda.1se,newx = X_test)
mean((y_test-alpha1.predict)^2)
#[1] 195562034 <- en teoria este es menor asi que mejor

## con ENET
alpha0.5 <- cv.glmnet(X_train, y_train, alpha = 0.5, type.measure = "mse")  # alpha = 0.5 para enet
alpha0.5.predict<- predict(alpha1,s=alpha1$lambda.1se,newx = X_test)
mean((y_test-alpha0.5.predict)^2)

#[1] 195562034 <- lo mismo que lasso

list.of.fits <- list()
for (i in 0:10) {
  fit.name <- paste0("alpha", i/10)
  list.of.fits[[fit.name]] <-
    cv.glmnet(X_train, y_train, type.measure="mse", alpha=i/10, 
              family="gaussian")
}
results <- data.frame()
for (i in 0:10) {
  fit.name <- paste0("alpha", i/10)
  predicted <- 
    predict(list.of.fits[[fit.name]], 
            s=list.of.fits[[fit.name]]$lambda.1se, newx=X_test)
  mse <- mean((y_test - predicted)^2)
  temp <- data.frame(alpha=i/10, mse=mse, fit.name=fit.name)
  results <- rbind(results, temp)
}

results

list.of.fits <- list()
# Iterar sobre los valores de alpha de 0.25 a 0.35 con un paso de 0.01
for (i in seq(0.25, 0.35, 0.01)) {
  fit.name <- paste0("alpha", i)
  list.of.fits[[fit.name]] <- cv.glmnet(X_train, y_train, type.measure="mse", alpha=i, family="gaussian") }
 
results1 <- data.frame()
 
 # Iterar sobre los valores de alpha de 0.25 a 0.35 con un paso de 0.01
 for (i in seq(0.25, 0.35, 0.01)) {
         fit.name <- paste0("alpha", i)
         predicted <- predict(list.of.fits[[fit.name]], s=list.of.fits[[fit.name]]$lambda.1se, newx=X_test)
         mse <- mean((y_test - predicted)^2)
         temp <- data.frame(alpha=i, mse=mse, fit.name=fit.name)
         results1 <- rbind(results1, temp)}
 
 # Renombrar el dataframe resultante como results1
results1
#alpha      mse  fit.name
# 0.25   49682216 alpha0.25 <seleccionado para miRNAs


