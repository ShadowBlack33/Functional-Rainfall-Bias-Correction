#Script: Proyecto final - Apartado estadístico

#Limpiar área de trabajo
rm(list=ls())

#Cargar librerías
#Nota: en el ítem 1 las librerías solo se usan para comprobar los resultados paso a paso
library(fda)
library(fda.usc)     #install.packages("fda.usc") la primera vez
library(data.table)
library(readxl)
library(terra)

#Semilla
set.seed(2026)

#Directorio de trabajo
#Nota: el script usa rutas relativas. Abrir R en la carpeta "Proyecto final"

#Verificación
getwd()
dir()

#Los datos deben estar en esa carpeta
stopifnot(file.exists(c("P_1980_1990.csv", "P_1991_2003.csv", "CVC_Wheather_Stations_Long_Lat.csv",
                        "CatalogoEstacionesElevacion.xlsx", "ImagenesSatelitePentadal.zip")))

#Carpetas para guardar las figuras y las tablas
dir.create("Figuras_imputacion", showWarnings = FALSE)
dir.create("Salidas_imputacion", showWarnings = FALSE)

#Fondo blanco para que los png no queden transparentes
windows.options(bg = "white")






#Bloque 1. Carga de las estaciones (precipitación diaria)

#Precipitación diaria 1980-2003 en dos archivos
#STATION_ID se lee como texto: tiene 10 dígitos y como número entero se desborda
P_1980_1990 <- fread("P_1980_1990.csv", colClasses = c(STATION_ID = "character"))
P_1991_2003 <- fread("P_1991_2003.csv", colClasses = c(STATION_ID = "character"))

Estaciones_diaria <- rbind(P_1980_1990, P_1991_2003)
rm(P_1980_1990, P_1991_2003)

#Periodo del proyecto
Estaciones_diaria <- Estaciones_diaria[YEAR >= 1981 & YEAR <= 2003]

#Fecha y día del año
Estaciones_diaria[, Fecha := as.Date(paste(YEAR, MONTH, DAY, sep = "-"))]
Estaciones_diaria[, DiaAño := as.numeric(format(Fecha, "%j"))]

head(Estaciones_diaria)
dim(Estaciones_diaria)

#Verificación
sum(is.na(Estaciones_diaria$Fecha))                                #fechas inválidas: 0
sum(duplicated(Estaciones_diaria, by = c("STATION_ID", "Fecha")))  #duplicados: 0
uniqueN(Estaciones_diaria$STATION_ID)                              #estaciones
uniqueN(Estaciones_diaria, by = c("STATION_ID", "YEAR"))           #pares estación-año
summary(Estaciones_diaria$VALUE)


#Catálogo de estaciones con coordenadas (CVC). El archivo está en Latin-1
Catalogo_CVC <- read.csv("CVC_Wheather_Stations_Long_Lat.csv", fileEncoding = "latin1",
                         colClasses = c(CODIGO = "character"))
Catalogo_CVC <- as.data.table(Catalogo_CVC)

#Catálogo de estaciones con altitud
Catalogo_Alt <- read_excel("CatalogoEstacionesElevacion.xlsx", sheet = 1,
                           col_types = c("text", "text", "text", "text", "numeric"))
Catalogo_Alt <- as.data.table(Catalogo_Alt)
Catalogo_Alt[, cuenca := trimws(cuenca)]       #una celda trae un espacio al final

dim(Catalogo_CVC)
dim(Catalogo_Alt)

#Cruce de identificadores entre los tres archivos
ID_lluvia <- unique(Estaciones_diaria$STATION_ID)
c(Con_lluvia            = length(ID_lluvia),
  Con_coordenadas       = sum(ID_lluvia %in% Catalogo_CVC$CODIGO),
  Con_altitud           = sum(ID_lluvia %in% Catalogo_Alt$codigo),
  Con_coordenadas_y_alt = sum(ID_lluvia %in% Catalogo_CVC$CODIGO & ID_lluvia %in% Catalogo_Alt$codigo))






#Bloque 2. Carga del satélite (precipitación por péntada)

#El zip trae archivos ocultos de macOS ("._ppt198101.tif"). El patrón exige que el
#nombre empiece por "ppt" para dejarlos por fuera
carpeta_sat  <- "ImagenesSatelitePentadal"

#Descomprimir las imágenes si la carpeta todavía no existe
if(!dir.exists(carpeta_sat)) unzip("ImagenesSatelitePentadal.zip")

archivos_sat <- list.files(carpeta_sat, pattern = "^ppt[0-9]{6}\\.tif$", full.names = TRUE)

#Año y péntada salen del nombre: pptAAAAPP.tif
Info_sat <- data.table(Archivo = archivos_sat,
                       YEAR    = as.numeric(substr(basename(archivos_sat), 4, 7)),
                       Pentada = as.numeric(substr(basename(archivos_sat), 8, 9)))

#Periodo del proyecto, ordenado por año y péntada
Info_sat <- Info_sat[YEAR >= 1981 & YEAR <= 2003][order(YEAR, Pentada)]

#Verificación: 23 años x 72 péntadas = 1656 imágenes
nrow(Info_sat)
table(table(Info_sat$YEAR))

#Pila de imágenes: una capa por péntada
Satelite <- rast(Info_sat$Archivo)
names(Satelite) <- paste(Info_sat$YEAR, sprintf("%02d", Info_sat$Pentada), sep = "_")

Satelite
res(Satelite)
ext(Satelite)
crs(Satelite, describe = TRUE)$code

#Verificación: sin negativos ni vacíos
sum(is.na(values(Satelite)))
range(values(Satelite))






#Bloque 3. Dominio anual común para las dos fuentes

#Dominio: t en [0, 365], medido en días desde el inicio del año
#Cada registro se ubica en el punto medio del periodo que cubre: el día d en t = d - 0.5
#y la péntada en la mitad de sus días
rangval_Anual = c(0, 365)

#Unidad de las curvas: lluvia por periodo, en mm por día
#Estaciones: el registro diario. Satélite: total de la péntada dividido por sus días (mm_dia)

#Estaciones: todos los años quedan con 365 días
#El 29 de febrero se suma al 28 (más adelante, después de depurar) para no cambiar el
#total anual: la péntada 12 del satélite también incluye ese día
Estaciones_diaria[, Bisiesto := (YEAR %% 4 == 0)]
Estaciones_diaria[, Feb29    := (MONTH == 2 & DAY == 29)]
Estaciones_diaria[, Dia := DiaAño]
Estaciones_diaria[Bisiesto == TRUE & DiaAño >= 60, Dia := DiaAño - 1]

#Verificación
sum(Estaciones_diaria$Feb29)                #registros de 29 de febrero
range(Estaciones_diaria$Dia)                #1 365


#Satélite: calendario de las 72 péntadas sobre el año de 365 días
#Supuesto: 6 péntadas por mes; las cinco primeras tienen 5 días y la sexta va del día 26
#al fin de mes (convención del producto CHIRPS). Confirmar con el docente
DiasMes   = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
InicioMes = cumsum(c(0, DiasMes[-12]))

Pentadas <- data.table(Pentada = 1:72)
Pentadas[, Mes    := (Pentada - 1) %/% 6 + 1]
Pentadas[, Orden  := (Pentada - 1) %% 6 + 1]
Pentadas[, DiaIni := InicioMes[Mes] + 5*(Orden - 1) + 1]
Pentadas[, DiaFin := ifelse(Orden < 6, InicioMes[Mes] + 5*Orden, InicioMes[Mes] + DiasMes[Mes])]
Pentadas[, nDias  := DiaFin - DiaIni + 1]
Pentadas[, t      := (DiaIni - 1 + DiaFin)/2]

Pentadas

#Verificación
sum(Pentadas$nDias)                         #365
table(Pentadas$nDias)                       #3 a 6 días por péntada






#Bloque 4. Emparejamiento estación - píxel (paso 3.1)

#Regla de asignación: a cada estación le corresponde la celda de la malla que la contiene
#Una estación fuera de la malla no tiene píxel y no entra a la comparación
Emparejamiento <- Catalogo_CVC[CODIGO %in% ID_lluvia,
                               .(STATION_ID = CODIGO, Nombre = STATION_NAME, Lon, Lat = lat)]
Emparejamiento <- merge(Emparejamiento, Catalogo_Alt[, .(STATION_ID = codigo, Altitud = alt)],
                        by = "STATION_ID", all.x = TRUE)

Emparejamiento[, pixel_id := cellFromXY(Satelite, cbind(Lon, Lat))]
Emparejamiento[, Fila     := rowFromCell(Satelite, pixel_id)]
Emparejamiento[, Columna  := colFromCell(Satelite, pixel_id)]

#Centro del píxel y distancia de la estación a ese centro (km)
Centro <- xyFromCell(Satelite, Emparejamiento$pixel_id)
Emparejamiento[, Lon_pixel := Centro[, 1]]
Emparejamiento[, Lat_pixel := Centro[, 2]]

#Distancia entre dos puntos en km (aproximación plana, suficiente a esta escala)
Distancia_km <- function(lon1, lat1, lon2, lat2){
  a  = pi/180
  dx = (lon1 - lon2)*a*cos((lat1 + lat2)*a/2)
  dy = (lat1 - lat2)*a
  6371*sqrt(dx^2 + dy^2)
}
Emparejamiento[, Dist_centro_km := Distancia_km(Lon, Lat, Lon_pixel, Lat_pixel)]

#Estaciones que comparten píxel
Emparejamiento[!is.na(pixel_id), Estaciones_en_pixel := .N, by = pixel_id]
Emparejamiento[, Regla := ifelse(is.na(pixel_id), "Fuera de la malla: sin pixel",
                                 "Celda de 0.05 grados que contiene la estacion")]


#Altitud: en esta versión solo se usa la del catálogo, sin DEM
#Así esta versión cambia una sola cosa frente a la principal (el relleno de días del paso 5b) y se
#puede ver qué efecto tiene. Las estaciones sin altitud en el catálogo quedan sin altitud
Emparejamiento[, Altitud_catalogo := Altitud]

#Fuente de la altitud de cada estación
Emparejamiento[, Fuente_altitud := ifelse(!is.na(Altitud_catalogo), "Catalogo", "Sin altitud")]
table(Emparejamiento[!is.na(pixel_id)]$Fuente_altitud)

Emparejamiento <- Emparejamiento[order(STATION_ID)]
head(Emparejamiento)

#Resumen del emparejamiento
c(Con_coordenadas    = nrow(Emparejamiento),
  Dentro_de_la_malla = sum(!is.na(Emparejamiento$pixel_id)),
  Fuera_de_la_malla  = sum(is.na(Emparejamiento$pixel_id)),
  Pixeles_distintos  = uniqueN(Emparejamiento[!is.na(pixel_id)]$pixel_id),
  Con_altitud_catalogo = sum(!is.na(Emparejamiento$pixel_id) & !is.na(Emparejamiento$Altitud_catalogo)),
  Con_altitud        = sum(!is.na(Emparejamiento$pixel_id) & !is.na(Emparejamiento$Altitud)))

#Estaciones fuera de la malla
Emparejamiento[is.na(pixel_id), .(STATION_ID, Nombre, Lon, Lat)]

#Píxeles con más de una estación
table(Emparejamiento[!is.na(pixel_id), .N, by = pixel_id]$N)
Emparejamiento[Estaciones_en_pixel > 1, .(pixel_id, STATION_ID, Nombre, Altitud)][order(pixel_id)]

#Mapa de las estaciones sobre la malla del satélite
windows()
plot(mean(Satelite[[1:72]])*72, main = "Estaciones sobre la malla del satélite (lluvia de 1981, mm)")
points(Emparejamiento$Lon, Emparejamiento$Lat, pch = 19, cex = 0.6,
       col = ifelse(is.na(Emparejamiento$pixel_id), "red", "black"))


#Series del satélite en los píxeles emparejados
Pixeles <- sort(unique(Emparejamiento[!is.na(pixel_id)]$pixel_id))
Sat_px  <- as.matrix(Satelite[Pixeles])     #píxeles en las filas, péntadas en las columnas
dim(Sat_px)

Satelite_pentada <- data.table(pixel_id = rep(Pixeles, times = nlyr(Satelite)),
                               YEAR     = rep(Info_sat$YEAR, each = length(Pixeles)),
                               Pentada  = rep(Info_sat$Pentada, each = length(Pixeles)),
                               VALUE    = as.vector(Sat_px))
Satelite_pentada <- merge(Satelite_pentada, Pentadas[, .(Pentada, nDias, t)], by = "Pentada")
Satelite_pentada[, mm_dia := VALUE/nDias]
setkey(Satelite_pentada, pixel_id, YEAR, Pentada)

#Verificación: una fila por píxel, año y péntada
nrow(Satelite_pentada) == length(Pixeles)*23*72
summary(Satelite_pentada[, .(Total_anual = sum(VALUE)), by = .(pixel_id, YEAR)]$Total_anual)






#Bloque 5. Depuración de las estaciones y pares estación - año que entran (subcriterio 4.1)

#Parámetros de la depuración (propuesta para acordar con el grupo)
#Ninguno se eligió mirando los años de prueba (2001-2003)
umbral_imposible = 500      #mm en un día: por encima se toma como error de digitación
umbral_total     = 100      #mm: tamaño mínimo de un total acumulado registrado en un día
max_registros    = 5        #un mes con tan pocos registros y valores altos trae acumulados
n_vecinos        = 5        #estaciones más cercanas para contrastar un mes en cero
radio_vecinos    = 30       #km
umbral_vecinos   = 30       #mm: mediana mensual de los vecinos para declarar un cero falso
cobertura_min    = 0.95     #fracción de días con dato que necesita un par estación - año
min_años         = 4        #años válidos que necesita una estación (criterio de los scripts de clase)

#Tabla de trazabilidad: estaciones, pares estación - año y registros con dato en cada paso
Conteo <- function(Paso, Datos){
  Datos = Datos[!is.na(VALUE)]
  data.table(Paso = Paso, Estaciones = uniqueN(Datos$STATION_ID),
             Pares = uniqueN(Datos, by = c("STATION_ID", "YEAR")), Registros = nrow(Datos))
}

Depurada <- copy(Estaciones_diaria)
Depurada[, Bandera := "ok"]
Trazabilidad <- Conteo("0. Registros originales 1981-2003", Depurada)


#Resumen por estación y mes
Depurada[, DiasDelMes := DiasMes[MONTH] + (MONTH == 2 & Bisiesto)]
Resumen_mes <- function(Datos){
  Datos[, .(n        = sum(!is.na(VALUE)),
            nLluvia  = sum(VALUE > 0, na.rm = TRUE),
            Maximo   = max(c(VALUE, 0), na.rm = TRUE),
            Total    = sum(VALUE, na.rm = TRUE),
            Completo = sum(!is.na(VALUE)) == DiasDelMes[1]),
        by = .(STATION_ID, YEAR, MONTH)]
}
Mensual <- Resumen_mes(Depurada)


#Paso 1: totales acumulados registrados en un solo día
#Se marca el mes cuando toda la lluvia está en un único día de umbral_total o más, o
#cuando el mes tiene max_registros o menos y alguno llega a umbral_total
Acumulados <- Mensual[Maximo >= umbral_total & (nLluvia == 1 | n <= max_registros)]
nrow(Acumulados)
Acumulados[order(-Maximo)][1:10]
table(Acumulados$STATION_ID)

Depurada <- merge(Depurada, Acumulados[, .(STATION_ID, YEAR, MONTH, Acum = TRUE)],
                  by = c("STATION_ID", "YEAR", "MONTH"), all.x = TRUE)
Depurada[Acum == TRUE & !is.na(VALUE), `:=`(VALUE = NA, Bandera = "total_acumulado")]
Depurada[, Acum := NULL]
Trazabilidad <- rbind(Trazabilidad, Conteo("1. Sin totales acumulados en un día", Depurada))


#Paso 2: valores imposibles o errores de digitación
#Va después del paso 1 para no contar como imposible un total mensual registrado en un día
Depurada[VALUE < 0 | VALUE > umbral_imposible, .(STATION_ID, Fecha, VALUE)]
Depurada[VALUE < 0 | VALUE > umbral_imposible, `:=`(VALUE = NA, Bandera = "imposible")]
Trazabilidad <- rbind(Trazabilidad, Conteo("2. Sin valores imposibles", Depurada))


#Paso 3: estaciones con coordenadas y dentro de la malla del satélite
#Sin píxel no hay curva satelital para comparar ni para la red
Depurada <- Depurada[STATION_ID %in% Emparejamiento[!is.na(pixel_id)]$STATION_ID]
Trazabilidad <- rbind(Trazabilidad, Conteo("3. Con coordenadas y dentro de la malla", Depurada))


#Paso 4: ceros que en realidad son datos faltantes
#Un mes completo en cero se contrasta con las n_vecinos estaciones más cercanas (hasta
#radio_vecinos km) que tienen ese mes completo. Si la mediana de sus totales llega a
#umbral_vecinos, el cero no es sequía: es un faltante
Mensual <- Resumen_mes(Depurada)
Ceros   <- Mensual[Completo == TRUE & Total == 0]
nrow(Ceros)
table(Ceros$MONTH)

Emp_ok   <- Emparejamiento[!is.na(pixel_id)]
Dist_est <- outer(1:nrow(Emp_ok), 1:nrow(Emp_ok),
                  function(i, j) Distancia_km(Emp_ok$Lon[i], Emp_ok$Lat[i], Emp_ok$Lon[j], Emp_ok$Lat[j]))
dimnames(Dist_est) <- list(Emp_ok$STATION_ID, Emp_ok$STATION_ID)

Total_mes <- dcast(Mensual[Completo == TRUE], YEAR + MONTH ~ STATION_ID, value.var = "Total")

Mediana_vecinos <- function(id, año, mes){
  d        = sort(Dist_est[id, ])
  vecinos  = names(d)[d > 0 & d <= radio_vecinos]
  vecinos  = vecinos[vecinos %in% names(Total_mes)]
  totales  = unlist(Total_mes[YEAR == año & MONTH == mes, ..vecinos])
  totales  = totales[!is.na(totales)]
  if(length(totales) == 0) return(NA_real_)
  median(totales[1:min(n_vecinos, length(totales))])
}
Ceros[, Med_vecinos := mapply(Mediana_vecinos, STATION_ID, YEAR, MONTH)]

summary(Ceros$Med_vecinos)
Ceros_falsos <- Ceros[!is.na(Med_vecinos) & Med_vecinos >= umbral_vecinos]
nrow(Ceros_falsos)
table(Ceros_falsos$MONTH)

Depurada <- merge(Depurada, Ceros_falsos[, .(STATION_ID, YEAR, MONTH, CeroFalso = TRUE)],
                  by = c("STATION_ID", "YEAR", "MONTH"), all.x = TRUE)
Depurada[CeroFalso == TRUE, `:=`(VALUE = NA, Bandera = "cero_falso")]
Depurada[, CeroFalso := NULL]
Trazabilidad <- rbind(Trazabilidad, Conteo("4. Sin ceros falsos", Depurada))


#Paso 5: calendario completo de 365 días por par estación - año
#El día que no viene como fila queda como faltante. El 29 de febrero se suma al 28
Feb_29   <- Depurada[Feb29 == TRUE, .(STATION_ID, YEAR, Valor29 = VALUE)]
Depurada <- merge(Depurada[Feb29 == FALSE], Feb_29, by = c("STATION_ID", "YEAR"), all.x = TRUE)
Depurada[Dia == 59 & !is.na(Valor29), VALUE := ifelse(is.na(VALUE), Valor29, VALUE + Valor29)]

Calendario <- unique(Depurada[, .(STATION_ID, YEAR)])[, .(Dia = 1:365), by = .(STATION_ID, YEAR)]
Depurada   <- merge(Calendario, Depurada[, .(STATION_ID, YEAR, Dia, VALUE, Bandera)],
                    by = c("STATION_ID", "YEAR", "Dia"), all.x = TRUE)
Depurada[is.na(Bandera), Bandera := "faltante"]
Depurada[, t := Dia - 0.5]

#Verificación: 365 filas por par
table(Depurada[, .N, by = .(STATION_ID, YEAR)]$N)
table(Depurada$Bandera)
Trazabilidad <- rbind(Trazabilidad, Conteo("5. Calendario de 365 días", Depurada))


#Paso 5b: relleno de días faltantes con Kriging universal
#El profesor dijo que con las coordenadas se puede imputar con otras estaciones, pero con cuidado,
#porque la altitud pesa (a mayor altitud, más ruido). En la clase de geoestadística también dejó
#claro que los faltantes espaciales no se rellenan con promedios simples. Por eso se usa el Kriging
#universal de la clase del 31 de agosto (Ejemplo4_Geoestadística.R), un día a la vez:
#  Z(s) = b0 + b1*x + b2*y + b3*altitud + e(s)
#la tendencia se ajusta por MCO con las estaciones que tienen dato ese día y el residual e(s) se
#interpola con Kriging. Límites para no inventar datos:
#- solo pares que ya tienen por su cuenta al menos cobertura_imputar de los días
#- no se rellena por encima de alt_max_imputar m ni en los años de prueba (2001-2003)
#- solo entran datos originales de ese día, nunca días ya rellenados
#- la lluvia no puede ser negativa: una predicción negativa se deja en 0
cobertura_imputar = 0.80
alt_max_imputar   = 2000     #m
min_est_dia       = 20       #estaciones con dato que necesita un día para hacer el Kriging

Alt_est = setNames(Emparejamiento$Altitud, Emparejamiento$STATION_ID)

#Coordenadas planas en km (la misma proyección local de la clase)
#Las estaciones sin altitud en el catálogo no pueden entrar a la tendencia
Emp_krig <- Emparejamiento[!is.na(Altitud), .(STATION_ID, Lon, Lat, Altitud)]
lat0 = mean(Emp_krig$Lat)
Emp_krig[, x_km := (Lon - mean(Lon)) * 111.320 * cos(lat0 * pi / 180)]
Emp_krig[, y_km := (Lat - lat0) * 110.574]
nrow(Emp_krig)

#Valores originales: un par estación - año por fila y un día por columna
Ancha   <- dcast(Depurada, STATION_ID + YEAR ~ Dia, value.var = "VALUE")
Clave   = paste(Ancha$STATION_ID, Ancha$YEAR)
Valores = as.matrix(Ancha[, -c(1, 2)])
rownames(Valores) = Clave

#Estaciones con altitud y con dato en un día
Datos_dia <- function(año, dia){
  f = which(Ancha$YEAR == año & Ancha$STATION_ID %in% Emp_krig$STATION_ID)
  d = Emp_krig[match(Ancha$STATION_ID[f], STATION_ID)]
  d[, z := Valores[f, dia]]
  d[!is.na(z)]
}

#Semivariograma de los residuales de la tendencia
#Cada día tiene su escala (un día seco casi no varía y uno de aguacero sí), por eso los residuales
#se dividen por su desviación estándar antes de juntar los días. Los pesos del Kriging no cambian
#con la escala, solo con la forma del semivariograma
Dias_semi = data.table(YEAR = sample(1981:2000, 600, replace = TRUE), Dia = sample(1:365, 600, replace = TRUE))
Pares_semi <- rbindlist(lapply(1:nrow(Dias_semi), function(k){
  d = Datos_dia(Dias_semi$YEAR[k], Dias_semi$Dia[k])
  if(nrow(d) < min_est_dia || sd(d$z) == 0) return(NULL)
  r = residuals(lm(z ~ x_km + y_km + Altitud, data = d))
  r = r / sd(r)
  D = as.matrix(dist(d[, .(x_km, y_km)]))
  data.table(h = D[upper.tri(D)], gamma = (outer(r, r, "-")^2 / 2)[upper.tri(D)])
}))
#Las dos correcciones de la clase: solo hasta la mitad de la distancia máxima y bins con 30 pares o más
Pares_semi <- Pares_semi[h <= max(h) / 2]
Pares_semi[, Bin := cut(h, breaks = 20)]
Semivariograma <- Pares_semi[, .(h = mean(h), gamma = mean(gamma), N = .N), by = Bin][N >= 30][order(h)]
Semivariograma

#Modelo exponencial ajustado con varios arranques de phi (L-BFGS-B, como en la clase)
sce_exponencial <- function(theta, h, gamma_emp){
  sum((gamma_emp - (theta[1] + theta[2] * (1 - exp(-h / theta[3]))))^2)
}
Arranques = quantile(Semivariograma$h, probs = c(0.1, 0.25, 0.5, 0.75, 1))
Ajustes <- lapply(Arranques, function(phi0){
  optim(par = c(0, 1, phi0), fn = sce_exponencial, h = Semivariograma$h, gamma_emp = Semivariograma$gamma,
        method = "L-BFGS-B", lower = c(0, 0, 1), upper = c(2, 5, 3 * max(Semivariograma$h)))
})
sapply(Ajustes, function(a) a$value)              #si son parecidos, el ajuste es estable
Ajuste_semi = Ajustes[[which.min(sapply(Ajustes, function(a) a$value))]]
c0_hat   = unname(Ajuste_semi$par[1])             #nugget
c1_hat   = unname(Ajuste_semi$par[2])             #sill - nugget
phi_hat  = unname(Ajuste_semi$par[3])
sill_hat = c0_hat + c1_hat
c(nugget = c0_hat, sill = sill_hat, rango_practico_km = 3 * phi_hat)

#Covarianza a partir del semivariograma, con el nugget solo en la diagonal (C(h) = c1*exp(-h/phi) si h > 0)
C_hat <- function(h) ifelse(h == 0, sill_hat, c1_hat * exp(-h / phi_hat))

#Kriging del residual en uno o varios puntos nuevos, con la restricción sum(lambda) = 1
Kriging_residual <- function(s0, coords_xy, resid){
  n_loc = nrow(coords_xy)
  A = rbind(cbind(C_hat(as.matrix(dist(coords_xy))), 1), c(rep(1, n_loc), 0))
  b = sapply(1:nrow(s0), function(j) C_hat(sqrt((coords_xy[, 1] - s0[j, 1])^2 + (coords_xy[, 2] - s0[j, 2])^2)))
  b = rbind(matrix(b, nrow = n_loc), 1)
  pesos = solve(A, b)[1:n_loc, , drop = FALSE]
  as.vector(t(pesos) %*% resid)
}

#Predicción de un día: tendencia (X*beta) + residual interpolado
Predecir_dia <- function(d, nuevos){
  m = lm(z ~ x_km + y_km + Altitud, data = d)
  tendencia = predict(m, newdata = nuevos)
  if(sd(residuals(m)) == 0) return(pmax(0, unname(tendencia)))
  pmax(0, unname(tendencia) + Kriging_residual(as.matrix(nuevos[, .(x_km, y_km)]),
                                               as.matrix(d[, .(x_km, y_km)]), residuals(m)))
}

#Pares que se pueden rellenar y sus días faltantes
Propia <- Depurada[, .(Cob_propia = mean(!is.na(VALUE))), by = .(STATION_ID, YEAR)]
Candidatos <- Propia[Cob_propia >= cobertura_imputar & Cob_propia < 1 & YEAR <= 2000 &
                     !is.na(Alt_est[STATION_ID]) & Alt_est[STATION_ID] <= alt_max_imputar]
nrow(Candidatos)
Faltantes <- rbindlist(lapply(1:nrow(Candidatos), function(k){
  y = Valores[paste(Candidatos$STATION_ID[k], Candidatos$YEAR[k]), ]
  data.table(STATION_ID = Candidatos$STATION_ID[k], YEAR = Candidatos$YEAR[k], Dia = which(is.na(y)))
}))
Faltantes <- merge(Faltantes, Emp_krig[, .(STATION_ID, x_km, y_km, Altitud)], by = "STATION_ID")

#Un Kriging por día: las estaciones que faltan ese día se predicen con el mismo sistema
Dias_falta = unique(Faltantes[, .(YEAR, Dia)])
Imputados <- rbindlist(lapply(1:nrow(Dias_falta), function(k){
  año = Dias_falta$YEAR[k]
  dia = Dias_falta$Dia[k]
  d = Datos_dia(año, dia)
  if(nrow(d) < min_est_dia) return(NULL)
  nuevos = Faltantes[YEAR == año & Dia == dia]
  data.table(STATION_ID = nuevos$STATION_ID, YEAR = año, Dia = dia,
             Imputado = Predecir_dia(d, nuevos), Estaciones = nrow(d))
}))
c(Dias_rellenados = nrow(Imputados), Pares = uniqueN(Imputados, by = c("STATION_ID", "YEAR")),
  Estaciones = uniqueN(Imputados$STATION_ID))
summary(Imputados$Estaciones)                    #estaciones con dato que entraron a cada Kriging

Depurada <- merge(Depurada, Imputados[, .(STATION_ID, YEAR, Dia, Imputado)],
                  by = c("STATION_ID", "YEAR", "Dia"), all.x = TRUE)
Depurada[!is.na(Imputado), `:=`(VALUE = Imputado, Bandera = "imputado")]
Depurada[, Imputado := NULL]
table(Depurada$Bandera)
Trazabilidad <- rbind(Trazabilidad, Conteo("5b. Días rellenados con Kriging", Depurada))


#Paso 6: cobertura mínima por par estación - año
Cobertura <- Depurada[, .(Dias_con_dato = sum(!is.na(VALUE)), Dias_imputados = sum(Bandera == "imputado")),
                      by = .(STATION_ID, YEAR)]
Cobertura[, Cobertura := Dias_con_dato/365]

#Pares que quedarían con otros criterios de cobertura
sapply(c(0.80, 0.90, 0.95, 0.98, 1), function(p) sum(Cobertura$Cobertura >= p))

windows()
hist(Cobertura$Cobertura, breaks = 50, main = "Cobertura por par estación - año",
     xlab = "Fracción de días con dato", ylab = "Pares")
abline(v = cobertura_min, col = "red", lty = 2)

Pares <- Cobertura[Cobertura >= cobertura_min]
Depurada <- merge(Depurada, Pares[, .(STATION_ID, YEAR)], by = c("STATION_ID", "YEAR"))
Trazabilidad <- rbind(Trazabilidad, Conteo("6. Cobertura mínima por par", Depurada))


#Paso 7: mínimo de años válidos por estación
#Solo se cuentan los años de entrenamiento y validación (hasta 2000): ninguna decisión de la
#depuración usa los años de prueba (2001-2003)
Años_estacion <- Pares[YEAR <= 2000, .(Años = .N), by = STATION_ID]
table(Años_estacion$Años)

Pares    <- Pares[STATION_ID %in% Años_estacion[Años >= min_años]$STATION_ID]
Depurada <- Depurada[STATION_ID %in% Pares$STATION_ID]
Trazabilidad <- rbind(Trazabilidad, Conteo("7. Mínimo de años por estación", Depurada))

Trazabilidad

#Informativo: pares finales con altitud del catálogo
Pares <- merge(Pares, Emparejamiento[, .(STATION_ID, pixel_id, Altitud, Fuente_altitud)], by = "STATION_ID")
table(Pares$Fuente_altitud)
c(Pares = nrow(Pares), Con_altitud = sum(!is.na(Pares$Altitud)),
  Estaciones = uniqueN(Pares$STATION_ID), Estaciones_con_altitud = uniqueN(Pares[!is.na(Altitud)]$STATION_ID))


#Partición temporal (las estaciones reservadas por el docente se separan cuando se conozcan)
Pares[, Conjunto := ifelse(YEAR <= 1995, "entrenamiento", ifelse(YEAR <= 2000, "validacion", "prueba"))]
setkey(Pares, STATION_ID, YEAR)
table(Pares$Conjunto)
Pares[, .(Pares = .N, Estaciones = uniqueN(STATION_ID)), by = Conjunto]

#Años por estación
table(Pares[, .N, by = STATION_ID]$N)






#Bloque 6. Archivos de intercambio

dir.create("Salidas_imputacion", showWarnings = FALSE)

#Tabla estación - píxel
fwrite(Emparejamiento, "Salidas_imputacion/emparejamiento.csv")

#Trazabilidad de la depuración
fwrite(Trazabilidad, "Salidas_imputacion/trazabilidad.csv")

#Pares estación - año que entran, con su conjunto
fwrite(Pares[, .(STATION_ID, YEAR, pixel_id, Altitud, Fuente_altitud, Dias_con_dato, Dias_imputados, Conjunto)],
       "Salidas_imputacion/particion.csv")

#Pares que entran gracias al relleno (con menos del 95 % de días propios)
Pares[, .(Pares = .N, Con_relleno = sum(Dias_con_dato - Dias_imputados < cobertura_min*365),
          Dias_imputados = sum(Dias_imputados))]

#Registros diarios depurados de los pares que entran
Estaciones_depuradas <- merge(Depurada, Emparejamiento[, .(STATION_ID, Lon, Lat, Altitud, pixel_id)],
                              by = "STATION_ID")
setkey(Estaciones_depuradas, STATION_ID, YEAR, Dia)
fwrite(Estaciones_depuradas[, .(station_id = STATION_ID, lon = Lon, lat = Lat, altitud = Altitud,
                                pixel_id, year = YEAR, day = Dia, t, value_mm = VALUE, bandera = Bandera)],
       "Salidas_imputacion/estaciones_depuradas.csv")

#Series del satélite en los píxeles emparejados
fwrite(Satelite_pentada[, .(pixel_id, year = YEAR, pentada = Pentada, t, n_dias = nDias,
                            value_mm = VALUE, mm_dia)], "Salidas_imputacion/satelite_pentada.csv")






#Bloque 7. ONI y fase ENSO de cada año

#Índice Oceánico de El Niño (ONI): media móvil de 3 meses de la anomalía de la temperatura del mar
#en la región Niño 3.4. Valores: NOAA Climate Prediction Center, tabla "Cold & Warm Episodes by
#Season" (ERSSTv5), https://cpc.ncep.noaa.gov/products/analysis_monitoring/enso/oni/v5/
#(consultada el 14 de noviembre de 2026). Una fila por año y una columna por trimestre móvil
#Criterio: el que adopta el IDEAM para Colombia en el "Protocolo para la vigilancia de condiciones
#de fenómeno ENOS en Colombia" (código SMA-PC001, versión 05, 2026): El Niño si el ONI es >= +0.5 °C
#y La Niña si es <= -0.5 °C durante cinco trimestres móviles consecutivos
Trimestres = c("DEF", "EFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDE")
ONI <- matrix(c(-0.3, -0.5, -0.5, -0.4, -0.3, -0.3, -0.3, -0.2, -0.2, -0.1, -0.2, -0.1,     #1981
                 0.0,  0.1,  0.2,  0.5,  0.7,  0.7,  0.8,  1.1,  1.6,  2.0,  2.2,  2.2,     #1982
                 2.2,  1.9,  1.5,  1.3,  1.1,  0.7,  0.3, -0.1, -0.5, -0.8, -1.0, -0.9,     #1983
                -0.6, -0.4, -0.3, -0.4, -0.5, -0.4, -0.3, -0.2, -0.2, -0.6, -0.9, -1.1,     #1984
                -1.0, -0.8, -0.8, -0.8, -0.8, -0.6, -0.5, -0.5, -0.4, -0.3, -0.3, -0.4,     #1985
                -0.5, -0.5, -0.3, -0.2, -0.1,  0.0,  0.2,  0.4,  0.7,  0.9,  1.1,  1.2,     #1986
                 1.2,  1.2,  1.1,  0.9,  1.0,  1.2,  1.5,  1.7,  1.6,  1.5,  1.3,  1.1,     #1987
                 0.8,  0.5,  0.1, -0.3, -0.9, -1.3, -1.3, -1.1, -1.2, -1.5, -1.8, -1.8,     #1988
                -1.7, -1.4, -1.1, -0.8, -0.6, -0.4, -0.3, -0.3, -0.2, -0.2, -0.2, -0.1,     #1989
                 0.1,  0.2,  0.3,  0.3,  0.3,  0.3,  0.3,  0.4,  0.4,  0.3,  0.4,  0.4,     #1990
                 0.4,  0.3,  0.2,  0.3,  0.5,  0.6,  0.7,  0.6,  0.6,  0.8,  1.2,  1.5,     #1991
                 1.7,  1.6,  1.5,  1.3,  1.1,  0.7,  0.4,  0.1, -0.1, -0.2, -0.3, -0.1,     #1992
                 0.1,  0.3,  0.5,  0.7,  0.7,  0.6,  0.3,  0.3,  0.2,  0.1,  0.0,  0.1,     #1993
                 0.1,  0.1,  0.2,  0.3,  0.4,  0.4,  0.4,  0.4,  0.6,  0.7,  1.0,  1.1,     #1994
                 1.0,  0.7,  0.5,  0.3,  0.1,  0.0, -0.2, -0.5, -0.8, -1.0, -1.0, -1.0,     #1995
                -0.9, -0.8, -0.6, -0.4, -0.3, -0.3, -0.3, -0.3, -0.4, -0.4, -0.4, -0.5,     #1996
                -0.5, -0.4, -0.1,  0.3,  0.8,  1.2,  1.6,  1.9,  2.1,  2.3,  2.4,  2.4,     #1997
                 2.2,  1.9,  1.4,  1.0,  0.5, -0.1, -0.8, -1.1, -1.3, -1.4, -1.5, -1.6,     #1998
                -1.5, -1.3, -1.1, -1.0, -1.0, -1.0, -1.1, -1.1, -1.2, -1.3, -1.5, -1.7,     #1999
                -1.7, -1.4, -1.1, -0.8, -0.7, -0.6, -0.6, -0.5, -0.5, -0.6, -0.7, -0.7,     #2000
                -0.7, -0.5, -0.4, -0.3, -0.3, -0.1, -0.1, -0.1, -0.2, -0.3, -0.3, -0.3,     #2001
                -0.1,  0.0,  0.1,  0.2,  0.4,  0.7,  0.8,  0.9,  1.0,  1.2,  1.3,  1.1,     #2002
                 0.9,  0.6,  0.4,  0.0, -0.3, -0.2,  0.1,  0.2,  0.3,  0.3,  0.4,  0.4),    #2003
              nrow = 23, ncol = 12, byrow = TRUE, dimnames = list(1981:2003, Trimestres))

umbral_ONI     = 0.5
min_trimestres = 5

#Episodios: rachas de al menos 5 trimestres seguidos por encima (o por debajo) del umbral
#Se recorre la serie entera para no cortar los episodios que pasan de un año al otro
#(la serie empieza y termina sin episodio en curso: DEF de 1981 = -0.3 y NDE de 2003 = 0.4)
ONI_serie = as.vector(t(ONI))                       #276 trimestres en orden
Episodio <- function(Condicion){
  r = rle(Condicion)
  rep(r$values & r$lengths >= min_trimestres, r$lengths)
}
Trim_Niño = matrix(Episodio(ONI_serie >= umbral_ONI), nrow = 23, byrow = TRUE)
Trim_Niña = matrix(Episodio(ONI_serie <= -umbral_ONI), nrow = 23, byrow = TRUE)

#Fase del año: la del episodio que ocupa al menos la mitad de sus trimestres (6 de 12)
ONI_anual <- data.table(YEAR = 1981:2003, ONI = round(rowMeans(ONI), 2),
                        Trim_Niño = rowSums(Trim_Niño), Trim_Niña = rowSums(Trim_Niña))
ONI_anual[, Fase := ifelse(Trim_Niño >= 6, "El Niño", ifelse(Trim_Niña >= 6, "La Niña", "Neutro"))]
ONI_anual
table(ONI_anual$Fase)

#Años con trimestres de los dos episodios (1998 empieza en El Niño y termina en La Niña)
ONI_anual[Trim_Niño > 0 & Trim_Niña > 0]

fwrite(cbind(ONI_anual, Fuente = "ONI: NOAA CPC (ERSSTv5). Criterio: IDEAM, SMA-PC001"), "Salidas_imputacion/oni_anual.csv")






#Malla de prueba para verificar la base del paso 1.1 (72 péntadas, como en la clase del 19 de octubre)
TiempoEval_Prec = seq(1, 72)
rangval_Prec    = c(1, 72)






#Ítem 1. Construcción de datos funcionales paso a paso (sin librerías)

#Paso 1.1: Matriz de la base Phi con B-splines (recursión de Cox-de Boor)
#TiempoEval: tiempos donde se evalúa la base
#rangval:    extremos del dominio c(a, b)
#K:          número de funciones base
#norden:     orden del B-spline (cúbico = 4)
Base_Bspline <- function(TiempoEval, rangval, K, norden){
  n = length(TiempoEval)

  #1. Nodos interiores igualmente espaciados (K - norden nodos)
  nInteriores = K - norden
  Interiores  = seq(rangval[1], rangval[2], length.out = nInteriores + 2)[-c(1, nInteriores + 2)]

  #2. Vector de nodos completo: cada extremo repetido norden veces (longitud K + norden)
  Nodos  = c(rep(rangval[1], norden), Interiores, rep(rangval[2], norden))
  nNodos = length(Nodos)

  #3. Orden 1: indicadoras de cada intervalo (K + norden - 1 columnas)
  #B_k,1(t) = 1 si Nodos[k] <= t < Nodos[k + 1], y 0 en otro caso
  B = matrix(0, n, nNodos - 1)
  for(k in 1:(nNodos - 1)){
    B[, k] = as.numeric(TiempoEval >= Nodos[k] & TiempoEval < Nodos[k + 1])
  }
  #El extremo derecho del dominio se asigna al último intervalo no vacío (columna K)
  B[TiempoEval == rangval[2], K] = 1

  #4. Recursión de orden 2 hasta norden (en cada nivel se pierde una columna)
  #B_k,m(t) = (t - Nodos[k])/(Nodos[k+m-1] - Nodos[k]) * B_k,m-1(t)
  #         + (Nodos[k+m] - t)/(Nodos[k+m] - Nodos[k+1]) * B_k+1,m-1(t)
  #Con nodos repetidos el denominador es cero y ese término vale cero
  if(norden > 1){
    for(m in 2:norden){
      Bm = matrix(0, n, nNodos - m)
      for(k in 1:(nNodos - m)){
        d1 = Nodos[k + m - 1] - Nodos[k]
        d2 = Nodos[k + m] - Nodos[k + 1]
        Termino1 = 0
        Termino2 = 0
        if(d1 > 0) Termino1 = (TiempoEval - Nodos[k])/d1*B[, k]
        if(d2 > 0) Termino2 = (Nodos[k + m] - TiempoEval)/d2*B[, k + 1]
        Bm[, k] = Termino1 + Termino2
      }
      B = Bm
    }
  }

  #5. Devolver la matriz de n x K
  B
}


#Caso de prueba 1: precipitación, B-splines cúbicos (orden 4) con K = 20
K_prueba      = 20
norden_prueba = 4

Phi <- Base_Bspline(TiempoEval_Prec, rangval_Prec, K_prueba, norden_prueba)

#Comprobaciones
dim(Phi)                                    #72 x 20
range(rowSums(Phi))                         #1 1 (partición de la unidad)
min(Phi)                                    #mayor o igual que 0
max(rowSums(Phi > 1e-12))                   #a lo sumo 4 valores no nulos por fila
c(Phi[1, 1], Phi[nrow(Phi), ncol(Phi)])     #1 1

#Comparación con la librería "fda": la diferencia debe ser cero
BasesPrec = create.bspline.basis(rangval_Prec, nbasis = K_prueba, norder = norden_prueba)
Phi_fda   = eval.basis(TiempoEval_Prec, BasesPrec)
max(abs(Phi - Phi_fda))

#Gráfico de las funciones base
windows()
matplot(TiempoEval_Prec, Phi, type = "l", lty = 1,
        xlab = "Tiempo (pentadas)", ylab = "Funciones base", main = "B-splines de orden 4 con K = 20")


#Caso de prueba 2: 24 registros, B-splines de orden 2 con K = 14 (caso del SO2)
TiempoEval_SO2 = seq(1, 24)
rangval_SO2    = c(1, 24)

Phi_SO2     <- Base_Bspline(TiempoEval_SO2, rangval_SO2, 14, 2)
Phi_SO2_fda = eval.basis(TiempoEval_SO2, create.bspline.basis(rangval_SO2, nbasis = 14, norder = 2))
max(abs(Phi_SO2 - Phi_SO2_fda))






#Datos para construir las curvas (pares que entran después de la depuración)

#Estaciones: un par estación - año por fila y un día por columna (mm por día). NA = sin dato
Y_est_dt   <- dcast(Depurada, STATION_ID + YEAR ~ Dia, value.var = "VALUE")
Curvas_est <- Y_est_dt[, .(STATION_ID, YEAR)]
Curvas_est <- merge(Curvas_est, Pares[, .(STATION_ID, YEAR, pixel_id, Conjunto)],
                    by = c("STATION_ID", "YEAR"), sort = FALSE)
Y_est      <- as.matrix(Y_est_dt[, -c(1, 2)])
Tiempo_est = (1:365) - 0.5

#Satélite: un píxel - año por fila y una péntada por columna (mm por día)
Pixel_año  <- unique(Pares[, .(pixel_id, YEAR, Conjunto)])
X_sat_dt   <- dcast(merge(Satelite_pentada, Pixel_año[, .(pixel_id, YEAR)], by = c("pixel_id", "YEAR")),
                    pixel_id + YEAR ~ Pentada, value.var = "mm_dia")
Curvas_sat <- merge(X_sat_dt[, .(pixel_id, YEAR)], Pixel_año, by = c("pixel_id", "YEAR"), sort = FALSE)
X_sat      <- as.matrix(X_sat_dt[, -c(1, 2)])
Tiempo_sat = Pentadas$t

#N: número de curvas; n: registros por curva
dim(Y_est)                                  #2870 x 365
dim(X_sat)                                  #píxel - año x 72
sum(rowSums(is.na(Y_est)) == 0)             #curvas de estación sin días faltantes
all(Curvas_est$STATION_ID == Y_est_dt$STATION_ID & Curvas_est$YEAR == Y_est_dt$YEAR)
all(Curvas_sat$pixel_id == X_sat_dt$pixel_id & Curvas_sat$YEAR == X_sat_dt$YEAR)

#Orden de la base: 4 (cúbica), para que la primera y la segunda derivada sean continuas
norden = 4

#Curvas de entrenamiento (1981-1995): con ellas se eligen K y los pesos, sin usar los años de prueba
Entrena_est = Curvas_est$Conjunto == "entrenamiento"
Entrena_sat = Curvas_sat$Conjunto == "entrenamiento"
c(Estaciones = sum(Entrena_est), Satelite = sum(Entrena_sat))

#Par de ejemplo: una estación y su píxel en un año con los 365 días
Ejemplo_est = which(Curvas_est$Conjunto == "entrenamiento" & rowSums(is.na(Y_est)) == 0)[1]
Ejemplo_sat = which(Curvas_sat$pixel_id == Curvas_est$pixel_id[Ejemplo_est] &
                    Curvas_sat$YEAR == Curvas_est$YEAR[Ejemplo_est])
Curvas_est[Ejemplo_est]
Curvas_sat[Ejemplo_sat]






#Paso 1.2: Coeficientes por mínimos cuadrados

#Una curva, paso a paso: y = Phi c + e
K_ejemplo = 8
y   = Y_est[Ejemplo_est, ]
Phi = Base_Bspline(Tiempo_est, rangval_Anual, K_ejemplo, norden)

PtP   = t(Phi) %*% Phi                      #Phi'Phi (K x K)
Pty   = t(Phi) %*% y                        #Phi'y   (K x 1)
c_est = solve(PtP, Pty)                     #ecuaciones normales: Phi'Phi c = Phi'y
c_est

#Matriz que lleva los datos a los coeficientes: C = (Phi'Phi)^(-1) Phi'  ("y2cMap" en fda)
y2cMap = solve(PtP, t(Phi))
max(abs(y2cMap %*% y - c_est))              #cero

#Comparación con la librería "fda": la diferencia debe ser cero
Base_fda = create.bspline.basis(rangval_Anual, nbasis = K_ejemplo, norder = norden)
c_fda    = smooth.basis(Tiempo_est, y, Base_fda)$fd$coefs
max(abs(c_est - c_fda))

#Curva estimada en una malla fina
Tiempo_fino = seq(0, 365, length.out = 731)
Phi_fino    = Base_Bspline(Tiempo_fino, rangval_Anual, K_ejemplo, norden)

windows()
plot(Tiempo_est, y, pch = 20, col = "gray60", xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)",
     main = paste("Estación", Curvas_est$STATION_ID[Ejemplo_est], "en", Curvas_est$YEAR[Ejemplo_est], "- K =", K_ejemplo))
lines(Tiempo_fino, Phi_fino %*% c_est, lwd = 2, col = "blue")


#Todas las curvas de una fuente
#Y:      matriz con una curva por fila (NA = sin dato)
#Tiempo: tiempos de las columnas de Y
#Pesos:  un peso por tiempo. Con pesos iguales es mínimos cuadrados ordinarios;
#        con pesos distintos es mínimos cuadrados ponderados (paso 1.5)
#Cada curva se ajusta con sus tiempos observados. Las curvas completas comparten Phi
Ajustar_curvas <- function(Y, Tiempo, K, norden, Pesos = rep(1, ncol(Y))){
  Phi = Base_Bspline(Tiempo, rangval_Anual, K, norden)
  N   = nrow(Y)
  Coef = matrix(NA, N, K)
  Res  = matrix(NA, N, ncol(Y))
  n = SCE = SCE_pond = VC = rep(NA, N)

  for(i in 1:N){
    obs   = !is.na(Y[i, ])
    Phi_i = Phi[obs, , drop = FALSE]
    w     = Pesos[obs]
    C_i   = solve(t(Phi_i) %*% (w*Phi_i), t(Phi_i*w))     #(Phi'W Phi)^(-1) Phi'W
    c_i   = C_i %*% Y[i, obs]
    e_i   = Y[i, obs] - Phi_i %*% c_i
    s_ii  = rowSums(Phi_i*t(C_i))                         #diagonal de S

    Coef[i, ]   = c_i
    Res[i, obs] = e_i
    n[i]        = sum(obs)
    SCE[i]      = sum(e_i^2)
    SCE_pond[i] = sum(w*e_i^2)
    VC[i]       = mean((e_i/(1 - s_ii))^2)
  }
  list(Coef = Coef, Res = Res, n = n, SCE = SCE, SCE_pond = SCE_pond, VC = VC, Phi = Phi)
}

#Verificación con la curva de ejemplo
Ajuste_prueba <- Ajustar_curvas(Y_est[Ejemplo_est, , drop = FALSE], Tiempo_est, K_ejemplo, norden)
max(abs(Ajuste_prueba$Coef[1, ] - c_est))   #cero






#Paso 1.3: Matriz de suavizado S y su verificación

#S = Phi (Phi'Phi)^(-1) Phi' lleva los datos a los valores ajustados: y_hat = S y
S = Phi %*% y2cMap
dim(S)

#Verificación numérica de una matriz de suavizado
Verificar_S <- function(Fuente, Metodo, S, K){
  data.table(Fuente = Fuente, Metodo = Metodo, n = nrow(S), K = K,
             Simetria     = max(abs(S - t(S))),           #S = S'
             Idempotencia = max(abs(S %*% S - S)),        #S S = S
             Traza        = sum(diag(S)))                 #tr(S) = K
}
Matriz_S <- function(Tiempo, K, norden, Pesos = rep(1, length(Tiempo))){
  Phi = Base_Bspline(Tiempo, rangval_Anual, K, norden)
  Phi %*% solve(t(Phi) %*% (Pesos*Phi), t(Phi*Pesos))
}

#Tabla para varios K en las dos fuentes
Tabla_S <- rbind(Verificar_S("Estaciones", "MCO", Matriz_S(Tiempo_est, 8, norden), 8),
                 Verificar_S("Estaciones", "MCO", Matriz_S(Tiempo_est, 20, norden), 20),
                 Verificar_S("Estaciones", "MCO", Matriz_S(Tiempo_est, 40, norden), 40),
                 Verificar_S("Satelite", "MCO", Matriz_S(Tiempo_sat, 8, norden), 8),
                 Verificar_S("Satelite", "MCO", Matriz_S(Tiempo_sat, 20, norden), 20),
                 Verificar_S("Satelite", "MCO", Matriz_S(Tiempo_sat, 40, norden), 40))
Tabla_S






#Paso 1.4: Selección de K con ECM, VC y VCG para cada fuente

#Para cada K y cada curva:
#ECM de ajuste = SCE/n
#VC  = promedio de [(y - y_hat)/(1 - s_ii)]^2   (dejar uno por fuera, con el atajo de s_ii)
#VCG = n SCE/(n - K)^2
#Los tres se promedian sobre las curvas de entrenamiento de la fuente
Criterios_K <- function(Y, Tiempo, Malla_K, norden){
  Tabla = NULL
  for(K in Malla_K){
    A = Ajustar_curvas(Y, Tiempo, K, norden)
    Tabla = rbind(Tabla, data.table(K = K, ECM = mean(A$SCE/A$n), VC = mean(A$VC),
                                    VCG = mean(A$n*A$SCE/(A$n - K)^2)))
  }
  Tabla
}

#En las estaciones K llega hasta 20: así los nodos quedan a más de 18 días, el hueco más largo
#que permite la cobertura del 95 %, y ninguna función base se queda sin datos
Malla_K_est = 4:20
Malla_K_sat = 4:30

Criterios_est <- Criterios_K(Y_est[Entrena_est, ], Tiempo_est, Malla_K_est, norden)
Criterios_sat <- Criterios_K(X_sat[Entrena_sat, ], Tiempo_sat, Malla_K_sat, norden)

Criterios_est
Criterios_sat

#K que minimiza cada criterio
Minimos_K <- rbind(data.table(Fuente = "Estaciones", ECM = Criterios_est[which.min(ECM)]$K,
                              VC = Criterios_est[which.min(VC)]$K, VCG = Criterios_est[which.min(VCG)]$K),
                   data.table(Fuente = "Satelite", ECM = Criterios_sat[which.min(ECM)]$K,
                              VC = Criterios_sat[which.min(VC)]$K, VCG = Criterios_sat[which.min(VCG)]$K))
Minimos_K

#K elegido: el mínimo de la VCG de cada fuente
K_est = Criterios_est[which.min(VCG)]$K
K_sat = Criterios_sat[which.min(VCG)]$K
c(K_est = K_est, K_sat = K_sat)

#Gráfica de los tres criterios por fuente
Graficar_criterio <- function(Tabla, Criterio, Titulo){
  plot(Tabla$K, Tabla[[Criterio]], type = "b", pch = 19, xlab = "K", ylab = Criterio, main = Titulo)
  abline(v = Tabla$K[which.min(Tabla[[Criterio]])], col = "red", lty = 2)
}
windows(width = 12, height = 7)
par(mfrow = c(2, 3))
Graficar_criterio(Criterios_est, "ECM", "Estaciones: ECM de ajuste")
Graficar_criterio(Criterios_est, "VC",  "Estaciones: VC")
Graficar_criterio(Criterios_est, "VCG", "Estaciones: VCG")
Graficar_criterio(Criterios_sat, "ECM", "Satélite: ECM de ajuste")
Graficar_criterio(Criterios_sat[K <= 20], "VC",  "Satélite: VC (K hasta 20)")
Graficar_criterio(Criterios_sat, "VCG", "Satélite: VCG")
par(mfrow = c(1, 1))
dev.print(pdf, "Figuras_imputacion/Fig_Seleccion_K.pdf")

#Contraste: mínimos con todas las curvas (1981-2003)
Criterios_est_todas <- Criterios_K(Y_est, Tiempo_est, 4:15, norden)
Criterios_sat_todas <- Criterios_K(X_sat, Tiempo_sat, 4:15, norden)
c(K_est_todas = Criterios_est_todas[which.min(VCG)]$K, K_sat_todas = Criterios_sat_todas[which.min(VCG)]$K)

#Verificación de S con el K elegido de cada fuente
Tabla_S <- rbind(Tabla_S,
                 Verificar_S("Estaciones", "MCO, K elegido", Matriz_S(Tiempo_est, K_est, norden), K_est),
                 Verificar_S("Satelite", "MCO, K elegido", Matriz_S(Tiempo_sat, K_sat, norden), K_sat))
Tabla_S






#Paso 1.5: Covarianza residual

#Ajuste por mínimos cuadrados ordinarios con el K elegido
MCO_est <- Ajustar_curvas(Y_est, Tiempo_est, K_est, norden)
MCO_sat <- Ajustar_curvas(X_sat, Tiempo_sat, K_sat, norden)

#Modelo estándar: Var[y] = sigma^2 I, con sigma^2 = SCE/(n - K) en cada curva
Sigma2_est = MCO_est$SCE/(MCO_est$n - K_est)
Sigma2_sat = MCO_sat$SCE/(MCO_sat$n - K_sat)
summary(Sigma2_est)
summary(Sigma2_sat)

#¿Varianza constante? Varianza de los residuos en cada tiempo, entre las curvas de entrenamiento
#Es la diagonal de Sigma = (1/N) suma de e_j e_j'
Var_t_est = colMeans(MCO_est$Res[Entrena_est, ]^2, na.rm = TRUE)
Var_t_sat = colMeans(MCO_sat$Res[Entrena_sat, ]^2, na.rm = TRUE)
c(Estaciones = max(Var_t_est)/min(Var_t_est), Satelite = max(Var_t_sat)/min(Var_t_sat))   #razón máx/mín

#¿Autocorrelación? Correlación de rezago 1 y estadístico de Durbin-Watson de cada curva
Autocorrelacion <- function(Res){
  Res = Res[rowSums(is.na(Res)) == 0, ]
  n   = ncol(Res)
  r1  = rowSums(Res[, -1]*Res[, -n])/rowSums(Res^2)
  DW  = rowSums((Res[, -1] - Res[, -n])^2)/rowSums(Res^2)
  c(Rezago_1 = mean(r1), Durbin_Watson = mean(DW))
}
rbind(Estaciones = Autocorrelacion(MCO_est$Res[Entrena_est, ]),
      Satelite   = Autocorrelacion(MCO_sat$Res[Entrena_sat, ]))

windows(width = 12, height = 5)
par(mfrow = c(1, 2))
plot(Tiempo_est, Var_t_est, type = "h", col = "gray40", xlab = "Tiempo (días)", ylab = "Varianza residual",
     main = "Estaciones: varianza de los residuos por día")
plot(Tiempo_sat, Var_t_sat, type = "b", pch = 19, xlab = "Tiempo (días)", ylab = "Varianza residual",
     main = "Satélite: varianza de los residuos por péntada")
par(mfrow = c(1, 1))
dev.print(pdf, "Figuras_imputacion/Fig_Varianza_Residual.pdf")

#La varianza cambia a lo largo del año y la autocorrelación es baja (Durbin-Watson cerca de 2),
#así que Sigma se toma diagonal, Sigma = sigma^2 V con V = diag(Var_t)
#Prueba con mínimos cuadrados ponderados, W = V^(-1): c = (Phi'W Phi)^(-1) Phi'W y
Pesos_est = 1/Var_t_est
Pesos_sat = 1/Var_t_sat

MCP_est <- Ajustar_curvas(Y_est, Tiempo_est, K_est, norden, Pesos_est)
MCP_sat <- Ajustar_curvas(X_sat, Tiempo_sat, K_sat, norden, Pesos_sat)

#Varianza residual ponderada de cada curva: sigma^2 = e'W e/(n - K)
Sigma2p_est = MCP_est$SCE_pond/(MCP_est$n - K_est)
Sigma2p_sat = MCP_sat$SCE_pond/(MCP_sat$n - K_sat)
summary(Sigma2p_est)
summary(Sigma2p_sat)

#Total anual de las curvas con cada método frente al observado (mm)
#Con la malla diaria, el total de una curva en mm/día es la suma de sus 365 valores
Phi_dia_est = Base_Bspline((1:365) - 0.5, rangval_Anual, K_est, norden)
Phi_dia_sat = Base_Bspline((1:365) - 0.5, rangval_Anual, K_sat, norden)
Obs_est = rowSums(Y_est, na.rm = TRUE)*365/rowSums(!is.na(Y_est))
Obs_sat = as.vector(X_sat %*% Pentadas$nDias)
Totales_metodo <- data.table(Fuente    = c("Estaciones", "Satelite"),
                             Observado = c(mean(Obs_est), mean(Obs_sat)),
                             MCO       = c(mean(rowSums(MCO_est$Coef %*% t(Phi_dia_est))), mean(rowSums(MCO_sat$Coef %*% t(Phi_dia_sat)))),
                             MCP       = c(mean(rowSums(MCP_est$Coef %*% t(Phi_dia_est))), mean(rowSums(MCP_sat$Coef %*% t(Phi_dia_sat)))))
Totales_metodo[, `:=`(Sesgo_MCO_pct = round(100*(MCO/Observado - 1), 1), Sesgo_MCP_pct = round(100*(MCP/Observado - 1), 1))]
Totales_metodo
fwrite(Totales_metodo, "Salidas_imputacion/totales_por_metodo.csv")

#Decisión: coeficientes de mínimos cuadrados ordinarios
#Los ponderados bajan la lluvia anual de las estaciones: dan poco peso a los días de las temporadas
#de lluvia, que es donde la varianza es alta, y la curva se pega a los días secos. Con MCO los
#B-splines suman 1, los residuos suman cero y el total se conserva. Con varianza no constante los
#coeficientes de MCO siguen siendo insesgados; lo que cambia es su varianza, y eso se corrige en
#las bandas con C Sigma C' (paso 1.6)
Coef_est = MCO_est$Coef
Coef_sat = MCO_sat$Coef
dim(Coef_est)
dim(Coef_sat)

#Cambio en las curvas entre los dos métodos (mm/día) y valores ajustados negativos con MCO
c(Estaciones = mean(abs(MCP_est$Coef %*% t(MCP_est$Phi) - MCO_est$Coef %*% t(MCO_est$Phi))),
  Satelite   = mean(abs(MCP_sat$Coef %*% t(MCP_sat$Phi) - MCO_sat$Coef %*% t(MCO_sat$Phi))))
c(Estaciones = mean(Coef_est %*% t(MCO_est$Phi) < 0), Satelite = mean(Coef_sat %*% t(MCO_sat$Phi) < 0))

#Verificación de S ponderada: S = Phi (Phi'W Phi)^(-1) Phi'W
#Sigue siendo idempotente y con traza K, pero ya no es simétrica
Tabla_S <- rbind(Tabla_S,
                 Verificar_S("Estaciones", "MCP, K elegido", Matriz_S(Tiempo_est, K_est, norden, Pesos_est), K_est),
                 Verificar_S("Satelite", "MCP, K elegido", Matriz_S(Tiempo_sat, K_sat, norden, Pesos_sat), K_sat))
Tabla_S
fwrite(Tabla_S, "Salidas_imputacion/verificacion_S.csv")






#Paso 1.6: Bandas de confianza punto a punto

#Coeficientes de MCO: c = C y con C = (Phi'Phi)^(-1) Phi'
#Var[c] = C Sigma C', con Sigma la covarianza residual estimada con los residuos de las curvas de
#entrenamiento, Sigma = (1/N) suma e_j e_j'. Como la autocorrelación es baja se usa su diagonal, diag(Var_t)
#Var[x(t)] = phi(t)' Var[c] phi(t)   y   banda: x(t) +- 2 raíz(Var[x(t)])
#Para comparar se calcula también la banda con varianza constante, Var[c] = sigma^2 (Phi'Phi)^(-1)
Banda_curva <- function(y, Tiempo, K, norden, Var_t, Tiempo_fino){
  obs    = !is.na(y)
  Phi    = Base_Bspline(Tiempo, rangval_Anual, K, norden)[obs, , drop = FALSE]
  C      = solve(t(Phi) %*% Phi, t(Phi))
  c_hat  = C %*% y[obs]
  Var_c  = C %*% diag(Var_t[obs]) %*% t(C)                 #C Sigma C'
  e      = y[obs] - Phi %*% c_hat
  Sigma2 = sum(e^2)/(sum(obs) - K)
  Var_c_cte = Sigma2*solve(t(Phi) %*% Phi)
  Phi_f  = Base_Bspline(Tiempo_fino, rangval_Anual, K, norden)
  x_hat  = as.vector(Phi_f %*% c_hat)
  EE     = sqrt(rowSums((Phi_f %*% Var_c)*Phi_f))
  EE_cte = sqrt(rowSums((Phi_f %*% Var_c_cte)*Phi_f))
  list(x = x_hat, EE = EE, Inferior = x_hat - 2*EE, Superior = x_hat + 2*EE,
       Inferior_cte = x_hat - 2*EE_cte, Superior_cte = x_hat + 2*EE_cte,
       Comprobacion = max(abs(C %*% (Sigma2*diag(sum(obs))) %*% t(C) - Var_c_cte)))
}

Banda_est <- Banda_curva(Y_est[Ejemplo_est, ], Tiempo_est, K_est, norden, Var_t_est, Tiempo_fino)
Banda_sat <- Banda_curva(X_sat[Ejemplo_sat, ], Tiempo_sat, K_sat, norden, Var_t_sat, Tiempo_fino)

#Comprobación: con Sigma = sigma^2 I, C Sigma C' es sigma^2 (Phi'Phi)^(-1)
c(Estaciones = Banda_est$Comprobacion, Satelite = Banda_sat$Comprobacion)

#Ancho de la banda (4 EE) a lo largo del año
summary(4*Banda_est$EE)
summary(4*Banda_sat$EE)

windows(width = 12, height = 5)
par(mfrow = c(1, 2))
plot(Tiempo_est, Y_est[Ejemplo_est, ], pch = 20, col = "gray70", xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)",
     main = paste("Estación", Curvas_est$STATION_ID[Ejemplo_est], "en", Curvas_est$YEAR[Ejemplo_est], "- K =", K_est))
lines(Tiempo_fino, Banda_est$x, lwd = 2, col = "blue")
lines(Tiempo_fino, Banda_est$Inferior, lty = 2, col = "blue")
lines(Tiempo_fino, Banda_est$Superior, lty = 2, col = "blue")
lines(Tiempo_fino, Banda_est$Inferior_cte, lty = 3, col = "gray40")
lines(Tiempo_fino, Banda_est$Superior_cte, lty = 3, col = "gray40")
legend("topright", c("Curva (MCO)", "Banda con C Sigma C'", "Banda con varianza constante"),
       col = c("blue", "blue", "gray40"), lty = c(1, 2, 3), bty = "n", cex = 0.8)
plot(Tiempo_sat, X_sat[Ejemplo_sat, ], pch = 20, col = "gray70", xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)",
     main = paste("Satélite, píxel", Curvas_sat$pixel_id[Ejemplo_sat], "en", Curvas_sat$YEAR[Ejemplo_sat], "- K =", K_sat))
lines(Tiempo_fino, Banda_sat$x, lwd = 2, col = "red")
lines(Tiempo_fino, Banda_sat$Inferior, lty = 2, col = "red")
lines(Tiempo_fino, Banda_sat$Superior, lty = 2, col = "red")
lines(Tiempo_fino, Banda_sat$Inferior_cte, lty = 3, col = "gray40")
lines(Tiempo_fino, Banda_sat$Superior_cte, lty = 3, col = "gray40")
par(mfrow = c(1, 1))
dev.print(pdf, "Figuras_imputacion/Fig_Bandas.pdf")

#Lectura de las bandas
#Son punto a punto: en cada t fijo, el intervalo cubre el valor de la curva x(t) con cerca del
#95 % de confianza. No cubren la curva completa a la vez (no son bandas simultáneas)
#Miden la precisión de la curva estimada, no la dispersión de los datos: por eso la mayoría de
#los registros diarios queda por fuera. Solo incluyen varianza: si K es muy pequeño, el sesgo
#puede dejar la curva verdadera por fuera






#Paso 1.7: Curvas de ambas fuentes sobre el mismo dominio anual

#Se usan B-splines cúbicos y no Fourier: la lluvia de un año no es periódica (el 31 de diciembre no
#tiene que empalmar con el 1 de enero) y cada base solo afecta un tramo del año, así los dos picos se
#ajustan por separado. Además el enunciado usa B-splines para la red

#Malla común: los 365 puntos medios de los días, t = d - 0.5, en el dominio [0, 365]
Tiempo_comun  = Tiempo_est
Phi_comun_est = Base_Bspline(Tiempo_comun, rangval_Anual, K_est, norden)
Phi_comun_sat = Base_Bspline(Tiempo_comun, rangval_Anual, K_sat, norden)

#Curvas evaluadas en la malla común, en mm/día (una curva por fila)
#Se usan los coeficientes finales (mínimos cuadrados ordinarios) de todos los pares que entran
Curvas_comun_est = Coef_est %*% t(Phi_comun_est)
Curvas_comun_sat = Coef_sat %*% t(Phi_comun_sat)
dim(Curvas_comun_est)                       #2870 x 365
dim(Curvas_comun_sat)                       #2508 x 365

#Cada par estación - año con la curva de su píxel en el mismo año (se usa en el ítem 3)
Indice_sat     = match(paste(Curvas_est$pixel_id, Curvas_est$YEAR), paste(Curvas_sat$pixel_id, Curvas_sat$YEAR))
sum(is.na(Indice_sat))                      #0: todos los pares tienen curva satelital
Curvas_par_sat = Curvas_comun_sat[Indice_sat, ]
dim(Curvas_par_sat)                         #2870 x 365

#Verificación de unidades: total anual de la curva (suma de los 365 valores) Vs total observado
#En las estaciones el observado se escala por los días con dato
Total_curva_est = rowSums(Curvas_comun_est)
Total_obs_est   = rowSums(Y_est, na.rm = TRUE)*365/rowSums(!is.na(Y_est))
Total_curva_sat = rowSums(Curvas_comun_sat)
Total_obs_sat   = as.vector(X_sat %*% Pentadas$nDias)
summary(100*(Total_curva_est/Total_obs_est - 1))     #diferencia relativa (%)
summary(100*(Total_curva_sat/Total_obs_sat - 1))

#Nota: el satélite no tiene datos antes de t = 2.5 ni después de t = 362.5, ahí la curva se extrapola


#Figuras finales del ítem 1

#Funciones base con el K elegido de cada fuente
Phi_fino_est = Base_Bspline(Tiempo_fino, rangval_Anual, K_est, norden)
Phi_fino_sat = Base_Bspline(Tiempo_fino, rangval_Anual, K_sat, norden)

windows(width = 12, height = 4.5)
par(mfrow = c(1, 2))
matplot(Tiempo_fino, Phi_fino_est, type = "l", lty = 1, xlab = "Tiempo (días)", ylab = "Funciones base",
        main = paste("Estaciones: B-splines cúbicos con K =", K_est))
matplot(Tiempo_fino, Phi_fino_sat, type = "l", lty = 1, xlab = "Tiempo (días)", ylab = "Funciones base",
        main = paste("Satélite: B-splines cúbicos con K =", K_sat))
par(mfrow = c(1, 1))
dev.print(pdf, "Figuras_imputacion/Fig_Base_Phi.pdf")

#Todas las curvas de cada fuente en el mismo dominio y un par de ejemplo
#Se guarda en png: con miles de curvas un pdf queda muy pesado
Rango_y = c(min(0, quantile(Curvas_comun_est, 0.001)), quantile(Curvas_comun_est, 0.999))

windows(width = 15, height = 5)
par(mfrow = c(1, 3))
matplot(Tiempo_comun, t(Curvas_comun_est), type = "l", lty = 1, col = rgb(0, 0, 1, 0.03), ylim = Rango_y,
        xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)",
        main = paste("Estaciones:", nrow(Curvas_comun_est), "curvas estación - año, K =", K_est))
abline(h = 0, lty = 3)
matplot(Tiempo_comun, t(Curvas_comun_sat), type = "l", lty = 1, col = rgb(1, 0, 0, 0.03), ylim = Rango_y,
        xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)",
        main = paste("Satélite:", nrow(Curvas_comun_sat), "curvas píxel - año, K =", K_sat))
abline(h = 0, lty = 3)
plot(Tiempo_est, Y_est[Ejemplo_est, ], pch = 20, col = "gray70", ylim = c(0, 30),
     xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)",
     main = paste("Estación", Curvas_est$STATION_ID[Ejemplo_est], "y su píxel en", Curvas_est$YEAR[Ejemplo_est]))
points(Tiempo_sat, X_sat[Ejemplo_sat, ], pch = 17, col = "red")
lines(Tiempo_comun, Curvas_comun_est[Ejemplo_est, ], lwd = 2, col = "blue")
lines(Tiempo_comun, Curvas_comun_sat[Ejemplo_sat, ], lwd = 2, col = "red")
legend("top", c("Estación (diaria)", "Satélite (péntada)", "Curva de la estación", "Curva del satélite"),
       pch = c(20, 17, NA, NA), lty = c(NA, NA, 1, 1), lwd = c(NA, NA, 2, 2),
       col = c("gray70", "red", "blue", "red"), bty = "n", cex = 0.85, ncol = 2)
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Curvas_Dominio_Comun.png", width = 1500, height = 500, res = 110)


#Verificación del relleno (paso 5b): validación cruzada dejando una estación por fuera, como en la clase
#En días con dato se saca una estación, se predice con las demás y se compara con lo que midió.
#Se compara contra dos métodos más simples, para ver si el Kriging sí aporta:
#- el promedio de las vecinas parecidas (a 30 km o menos y con menos de 300 m de diferencia de altitud)
#- k vecinos más cercanos con peso por el inverso de la distancia al cuadrado (IDW, k = 5)
Vecinas <- lapply(setNames(rownames(Dist_est), rownames(Dist_est)), function(id){
  if(is.na(Alt_est[id])) return(character(0))
  d = Dist_est[id, ]
  v = names(d)[d > 0 & d <= radio_vecinos]
  v[!is.na(Alt_est[v]) & abs(Alt_est[v] - Alt_est[id]) <= 300]
})
KNN_idw <- function(d, nuevo, k = 5){
  dist_k = sqrt((d$x_km - nuevo$x_km)^2 + (d$y_km - nuevo$y_km)^2)
  o = order(dist_k)[1:k]
  sum(d$z[o] / dist_k[o]^2) / sum(1 / dist_k[o]^2)
}
Dias_val  = data.table(YEAR = sample(1981:2000, 150, replace = TRUE), Dia = sample(1:365, 150, replace = TRUE))
Elegibles = Emp_krig[Altitud <= alt_max_imputar, STATION_ID]
Validacion_relleno <- rbindlist(lapply(1:nrow(Dias_val), function(k){
  d = Datos_dia(Dias_val$YEAR[k], Dias_val$Dia[k])
  if(nrow(d) < min_est_dia) return(NULL)
  rbindlist(lapply(which(d$STATION_ID %in% Elegibles), function(i){
    v = intersect(Vecinas[[d$STATION_ID[i]]], d$STATION_ID)
    data.table(STATION_ID = d$STATION_ID[i], Altitud = d$Altitud[i], Real = d$z[i],
               Kriging = Predecir_dia(d[-i], d[i]),
               KNN = KNN_idw(d[-i], d[i]),
               Vecinas = if(length(v) > 0) mean(d$z[match(v, d$STATION_ID)]) else NA_real_)
  }))
}))
Validacion_relleno <- Validacion_relleno[!is.na(Vecinas)]          #los mismos casos para los dos métodos
Validacion_relleno[, Banda := ifelse(Altitud < 1200, "Menos de 1200 m", "1200 a 2000 m")]
Resumen_relleno <- Validacion_relleno[, .(Casos = .N,
                     RMSE_kriging = round(sqrt(mean((Real - Kriging)^2)), 2),
                     RMSE_knn     = round(sqrt(mean((Real - KNN)^2)), 2),
                     RMSE_vecinas = round(sqrt(mean((Real - Vecinas)^2)), 2),
                     MAE_kriging  = round(mean(abs(Real - Kriging)), 2),
                     MAE_knn      = round(mean(abs(Real - KNN)), 2),
                     MAE_vecinas  = round(mean(abs(Real - Vecinas)), 2),
                     Cor_kriging  = round(cor(Real, Kriging), 2),
                     Cor_knn      = round(cor(Real, KNN), 2),
                     Cor_vecinas  = round(cor(Real, Vecinas), 2),
                     Sesgo_kriging = round(mean(Kriging - Real), 2),
                     Sesgo_knn     = round(mean(KNN - Real), 2),
                     Sesgo_vecinas = round(mean(Vecinas - Real), 2)), by = Banda]
Resumen_relleno
fwrite(Resumen_relleno, "Salidas_imputacion/verificacion_relleno.csv")
fwrite(Semivariograma[, .(h, gamma, N)], "Salidas_imputacion/semivariograma_relleno.csv")

windows(width = 10, height = 5)
par(mfrow = c(1, 2))
curva_h = seq(0, max(Semivariograma$h), length.out = 100)
plot(Semivariograma$h, Semivariograma$gamma, pch = 19, col = "navy", ylim = c(0, max(Semivariograma$gamma)),
     xlab = "Distancia (km)", ylab = "Semivarianza (residual estandarizado)",
     main = "Semivariograma de la lluvia diaria")
lines(curva_h, c0_hat + c1_hat * (1 - exp(-curva_h / phi_hat)), col = "darkorange", lwd = 2)
legend("bottomright", c("Empírico", "Exponencial ajustado"), pch = c(19, NA), lty = c(NA, 1),
       lwd = c(NA, 2), col = c("navy", "darkorange"), bty = "n")
plot(Validacion_relleno$Real, Validacion_relleno$Real - Validacion_relleno$Kriging, pch = 20,
     col = rgb(0, 0, 0, 0.15), xlab = "Lluvia observada (mm/día)", ylab = "Error (observado - Kriging)",
     main = "Validación cruzada (LOO) del Kriging")
abline(h = 0, col = "steelblue", lty = 2, lwd = 2)
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Verificacion_Relleno.png", width = 1100, height = 550, res = 110)


#Verificación del ítem 1 contra la rúbrica
Resumen_Item1 <- data.table(
  Criterio  = c("Phi y coeficientes por mínimos cuadrados, sin librerías",
                "Selección de K con ECM, VC y VCG por fuente",
                "Matriz S idempotente y con traza K",
                "Covarianza residual y bandas punto a punto",
                "Ambas fuentes sobre el mismo dominio anual"),
  Evidencia = c(sprintf("Phi contra eval.basis: %.1e. Coeficientes contra smooth.basis: %.1e",
                        max(abs(Base_Bspline(TiempoEval_Prec, rangval_Prec, K_prueba, norden_prueba) - Phi_fda)),
                        max(abs(c_est - c_fda))),
                sprintf("Estaciones: K = %d por VCG (VC = %d). Satélite: K = %d por VCG (VC = %d). Fig_Seleccion_K.pdf",
                        K_est, Minimos_K[Fuente == "Estaciones"]$VC, K_sat, Minimos_K[Fuente == "Satelite"]$VC),
                sprintf("Idempotencia máxima %.1e y diferencia máxima entre la traza y K de %.1e. Salidas_imputacion/verificacion_S.csv",
                        max(Tabla_S$Idempotencia), max(abs(Tabla_S$Traza - Tabla_S$K))),
                sprintf("Varianza residual máx/mín de %.1f (estaciones) y %.1f (satélite), autocorrelación baja: MCO con bandas C Sigma C'. Fig_Varianza_Residual.pdf y Fig_Bandas.pdf",
                        max(Var_t_est)/min(Var_t_est), max(Var_t_sat)/min(Var_t_sat)),
                sprintf("%d curvas de estación y %d de satélite en t de 0 a 365 días, en mm/día. Fig_Curvas_Dominio_Comun.png",
                        nrow(Curvas_comun_est), nrow(Curvas_comun_sat))))
Resumen_Item1
fwrite(Resumen_Item1, "Salidas_imputacion/resumen_item1.csv")






#Ítem 2. Análisis exploratorio funcional

#Se usan todos los pares que pasan la depuración en 1981-2003 (conteos en Salidas_imputacion/trazabilidad.csv)
#y las curvas del satélite de esos mismos píxeles y años
Curvas_por_año <- rbind(Estaciones = table(Curvas_est$YEAR), Satelite = table(Curvas_sat$YEAR))
Curvas_por_año
range(Curvas_por_año["Estaciones", ])       #curvas por año: mínimo y máximo

#Estaciones: 2870 curvas. Satélite: 2508 (las estaciones del mismo píxel comparten curva)






#Paso 2.1: Media funcional

#Media funcional paso a paso: x_barra(t) = c_barra' phi(t)
CoefMedia_est = colMeans(Coef_est)
CoefMedia_sat = colMeans(Coef_sat)
Media_est = as.vector(Phi_comun_est %*% CoefMedia_est)
Media_sat = as.vector(Phi_comun_sat %*% CoefMedia_sat)

#Comprobación con el promedio de las curvas evaluadas
max(abs(Media_est - colMeans(Curvas_comun_est)))
max(abs(Media_sat - colMeans(Curvas_comun_sat)))

#Comprobación con la función mean.fd() de "fda": la diferencia debe ser cero
Base_est_fda = create.bspline.basis(rangval_Anual, nbasis = K_est, norder = norden)
Base_sat_fda = create.bspline.basis(rangval_Anual, nbasis = K_sat, norder = norden)
Est_fd = fd(t(Coef_est), Base_est_fda)
Sat_fd = fd(t(Coef_sat), Base_sat_fda)

MediaEst_fda = as.vector(eval.fd(Tiempo_comun, mean.fd(Est_fd)))
MediaSat_fda = as.vector(eval.fd(Tiempo_comun, mean.fd(Sat_fd)))
max(abs(Media_est - MediaEst_fda))
max(abs(Media_sat - MediaSat_fda))






#Paso 2.2: Varianza, desviación estándar y coeficiente de variación funcional

#Var(t) = (N-1)^{-1} sum_i [x_i(t) - x_barra(t)]^2
N_est = nrow(Curvas_comun_est)
N_sat = nrow(Curvas_comun_sat)
Var_est = colSums(sweep(Curvas_comun_est, 2, Media_est)^2)/(N_est - 1)
Var_sat = colSums(sweep(Curvas_comun_sat, 2, Media_sat)^2)/(N_sat - 1)

#Con los coeficientes: Var(t) = phi(t)' Sigma_c phi(t)
SigmaC_est = cov(Coef_est)
SigmaC_sat = cov(Coef_sat)
max(abs(Var_est - rowSums((Phi_comun_est %*% SigmaC_est)*Phi_comun_est)))
max(abs(Var_sat - rowSums((Phi_comun_sat %*% SigmaC_sat)*Phi_comun_sat)))

#Comprobación con var.fd(): la varianza es la diagonal de la covarianza
max(abs(Var_est - diag(eval.bifd(Tiempo_comun, Tiempo_comun, var.fd(Est_fd)))))
max(abs(Var_sat - diag(eval.bifd(Tiempo_comun, Tiempo_comun, var.fd(Sat_fd)))))

#Desviación estándar y coeficiente de variación funcional
DE_est = sqrt(Var_est)
DE_sat = sqrt(Var_sat)
CV_est = DE_est/Media_est
CV_sat = DE_sat/Media_sat


#Día del año a fecha (año de 365 días)
Dia_a_fecha <- function(t) format(as.Date(floor(t), origin = "2019-01-01"), "%d-%b")

#Picos, mínimo de mitad de año y total anual de la media
Resumen_fuente <- function(Fuente, Media, DE){
  Primera = Tiempo_comun < 182.5
  p1 = which.max(ifelse(Primera, Media, -Inf))
  p2 = which.max(ifelse(Primera, -Inf, Media))
  m  = (p1:p2)[which.min(Media[p1:p2])]
  data.table(Fuente = Fuente, Total_anual_mm = round(sum(Media)),
             Pico_1 = Dia_a_fecha(Tiempo_comun[p1]), Pico_1_mm_dia = round(Media[p1], 2),
             Minimo = Dia_a_fecha(Tiempo_comun[m]),  Minimo_mm_dia = round(Media[m], 2),
             Pico_2 = Dia_a_fecha(Tiempo_comun[p2]), Pico_2_mm_dia = round(Media[p2], 2),
             DE_media = round(mean(DE), 2), CV_medio = round(mean(DE/Media), 3))
}
Resumen_Media <- rbind(Resumen_fuente("Estaciones", Media_est, DE_est),
                       Resumen_fuente("Satelite", Media_sat, DE_sat))
Resumen_Media
fwrite(Resumen_Media, "Salidas_imputacion/resumen_media_DE.csv")

#Media, DE y CV de cada fuente en la malla común
fwrite(data.table(t = Tiempo_comun, Fecha = Dia_a_fecha(Tiempo_comun),
                  Media_est, DE_est, CV_est, Media_sat, DE_sat, CV_sat), "Salidas_imputacion/media_DE_funcional.csv")


#Gráfico de la media +- 2 desviaciones estándar (es dispersión, no intervalo de confianza)
Rango_2DE = range(Media_est - 2*DE_est, Media_est + 2*DE_est, Media_sat - 2*DE_sat, Media_sat + 2*DE_sat)

windows(width = 12, height = 9)
par(mfrow = c(2, 2))
matplot(Tiempo_comun, t(Curvas_comun_est), type = "l", lty = 1, col = rgb(0.5, 0.5, 0.5, 0.03), ylim = Rango_2DE,
        xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)", main = "Estaciones: media y media ± 2 DE")
lines(Tiempo_comun, Media_est, lwd = 3, col = "blue")
lines(Tiempo_comun, Media_est + 2*DE_est, lwd = 2, lty = 2, col = "blue")
lines(Tiempo_comun, Media_est - 2*DE_est, lwd = 2, lty = 2, col = "blue")
abline(h = 0, lty = 3)
matplot(Tiempo_comun, t(Curvas_comun_sat), type = "l", lty = 1, col = rgb(0.5, 0.5, 0.5, 0.03), ylim = Rango_2DE,
        xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)", main = "Satélite: media y media ± 2 DE")
lines(Tiempo_comun, Media_sat, lwd = 3, col = "red")
lines(Tiempo_comun, Media_sat + 2*DE_sat, lwd = 2, lty = 2, col = "red")
lines(Tiempo_comun, Media_sat - 2*DE_sat, lwd = 2, lty = 2, col = "red")
abline(h = 0, lty = 3)
plot(Tiempo_comun, Media_est, type = "l", lwd = 3, col = "blue", ylim = range(0, Media_est, Media_sat),
     xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)", main = "Media funcional")
lines(Tiempo_comun, Media_sat, lwd = 3, lty = 2, col = "red")
legend("bottomright", c("Estaciones", "Satélite"), col = c("blue", "red"), lty = c(1, 2), lwd = 3, bty = "n")
plot(Tiempo_comun, DE_est, type = "l", lwd = 3, col = "blue", ylim = range(0, DE_est, DE_sat),
     xlab = "Tiempo (días)", ylab = "mm/día", main = "Desviación estándar funcional")
lines(Tiempo_comun, DE_sat, lwd = 3, lty = 2, col = "red")
legend("bottomright", c("Estaciones", "Satélite"), col = c("blue", "red"), lty = c(1, 2), lwd = 3, bty = "n")
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Media_DE.png", width = 1400, height = 1050, res = 110)






#Paso 2.3: Covarianza y correlación funcional

#Covarianza funcional paso a paso: v(s,t) = phi(s)' Sigma_c phi(t)
Cov_est = Phi_comun_est %*% SigmaC_est %*% t(Phi_comun_est)        #365 x 365
Cov_sat = Phi_comun_sat %*% SigmaC_sat %*% t(Phi_comun_sat)

#Correlación funcional: r(s,t) = v(s,t)/raíz(v(s,s) v(t,t))
Cor_est = Cov_est/sqrt(outer(diag(Cov_est), diag(Cov_est)))
Cor_sat = Cov_sat/sqrt(outer(diag(Cov_sat), diag(Cov_sat)))

#La diagonal de la covarianza es la varianza funcional
max(abs(diag(Cov_est) - Var_est))
max(abs(diag(Cov_sat) - Var_sat))

#Comprobación con var.fd() y cor.fd(): la diferencia debe ser cero
max(abs(Cov_est - eval.bifd(Tiempo_comun, Tiempo_comun, var.fd(Est_fd))))
max(abs(Cov_sat - eval.bifd(Tiempo_comun, Tiempo_comun, var.fd(Sat_fd))))
max(abs(Cor_est - cor.fd(Tiempo_comun, Est_fd)))
max(abs(Cor_sat - cor.fd(Tiempo_comun, Sat_fd)))


#Correlación promedio según los días de separación
Cor_rezago <- function(Cor, Rezagos){
  sapply(Rezagos, function(h) mean(Cor[cbind(1:(365 - h), (1 + h):365)]))
}
Rezagos = c(15, 30, 60, 90, 180)
Memoria <- rbind(Estaciones = Cor_rezago(Cor_est, Rezagos), Satelite = Cor_rezago(Cor_sat, Rezagos))
colnames(Memoria) <- paste(Rezagos, "dias")
round(Memoria, 3)

#Posición de los dos picos y del mínimo de la media
Indices_media <- function(Media){
  Primera = Tiempo_comun < 182.5
  p1 = which.max(ifelse(Primera, Media, -Inf))
  p2 = which.max(ifelse(Primera, -Inf, Media))
  m  = (p1:p2)[which.min(Media[p1:p2])]
  c(Pico_1 = p1, Minimo = m, Pico_2 = p2)
}
i_est = Indices_media(Media_est)
i_sat = Indices_media(Media_sat)

#Correlación entre los picos y el mínimo
Cor_momentos <- rbind(Estaciones = c(Pico1_Pico2  = Cor_est[i_est["Pico_1"], i_est["Pico_2"]],
                                     Pico1_Minimo = Cor_est[i_est["Pico_1"], i_est["Minimo"]],
                                     Minimo_Pico2 = Cor_est[i_est["Minimo"], i_est["Pico_2"]]),
                      Satelite   = c(Pico1_Pico2  = Cor_sat[i_sat["Pico_1"], i_sat["Pico_2"]],
                                     Pico1_Minimo = Cor_sat[i_sat["Pico_1"], i_sat["Minimo"]],
                                     Minimo_Pico2 = Cor_sat[i_sat["Minimo"], i_sat["Pico_2"]]))
round(Cor_momentos, 3)

#Mapa de calor con curvas de nivel
Rango_cov = range(Cov_est, Cov_sat)
Paleta_cov = hcl.colors(40, "YlOrRd", rev = TRUE)
#Colores suaves: el gris, el azul y el verde - café se leen mejor que el magenta o el rosado
Paleta_cor = hcl.colors(40, "Green-Brown")

windows(width = 12, height = 10)
par(mfrow = c(2, 2))
image(Tiempo_comun, Tiempo_comun, Cov_est, zlim = Rango_cov, col = Paleta_cov,
      xlab = "s (días)", ylab = "t (días)", main = "Estaciones: covarianza (mm/día)²")
contour(Tiempo_comun, Tiempo_comun, Cov_est, nlevels = 8, labcex = 0.6, add = TRUE)
image(Tiempo_comun, Tiempo_comun, Cov_sat, zlim = Rango_cov, col = Paleta_cov,
      xlab = "s (días)", ylab = "t (días)", main = "Satélite: covarianza (mm/día)²")
contour(Tiempo_comun, Tiempo_comun, Cov_sat, nlevels = 8, labcex = 0.6, add = TRUE)
image(Tiempo_comun, Tiempo_comun, Cor_est, zlim = c(-1, 1), col = Paleta_cor,
      xlab = "s (días)", ylab = "t (días)", main = "Estaciones: correlación")
contour(Tiempo_comun, Tiempo_comun, Cor_est, levels = seq(-1, 1, 0.2), labcex = 0.6, add = TRUE)
image(Tiempo_comun, Tiempo_comun, Cor_sat, zlim = c(-1, 1), col = Paleta_cor,
      xlab = "s (días)", ylab = "t (días)", main = "Satélite: correlación")
contour(Tiempo_comun, Tiempo_comun, Cor_sat, levels = seq(-1, 1, 0.2), labcex = 0.6, add = TRUE)
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Covarianza_Correlacion.png", width = 1300, height = 1100, res = 110)






#Paso 2.4: Derivadas e inicio y fin de las temporadas de lluvia

#D x(t) = c' D phi(t) (mm/día por día). Con eval.fd(..., 1) se deriva la base directamente
DMedia_est = as.vector(eval.fd(Tiempo_comun, mean.fd(Est_fd), 1))
DMedia_sat = as.vector(eval.fd(Tiempo_comun, mean.fd(Sat_fd), 1))

#La derivada de la media es igual a la media de las derivadas
max(abs(rowMeans(eval.fd(Tiempo_comun, Est_fd, 1)) - DMedia_est))
max(abs(rowMeans(eval.fd(Tiempo_comun, Sat_fd, 1)) - DMedia_sat))

#deriv.fd() vuelve a proyectar la derivada en la base, por eso no da exacto
max(abs(eval.fd(Tiempo_comun, deriv.fd(mean.fd(Est_fd), 1)) - DMedia_est))

#Comprobación con diferencias centrales
h = 1e-3
D_numerica <- function(K, CoefMedia){
  as.vector((Base_Bspline(Tiempo_comun + h, rangval_Anual, K, norden) -
             Base_Bspline(Tiempo_comun - h, rangval_Anual, K, norden)) %*% CoefMedia)/(2*h)
}
max(abs(D_numerica(K_est, CoefMedia_est) - DMedia_est))
max(abs(D_numerica(K_sat, CoefMedia_sat) - DMedia_sat))


#Inicio de cada temporada: máximo de la derivada antes del pico. Fin: mínimo después del pico
#No se buscan en los 10 primeros ni en los 10 últimos días (borde de la base)
#También se da el criterio de nivel: días en que la media está por encima de su promedio anual
Temporadas <- function(Fuente, Media, DMedia){
  i  = Indices_media(Media)
  t1 = which(Tiempo_comun >= 10)[1]
  t2 = max(which(Tiempo_comun <= 355))
  v0 = (t1:i["Pico_1"])[which.min(Media[t1:i["Pico_1"]])]              #valle antes del primer pico
  Tramos = list(Subida_1 = v0:i["Pico_1"],     Bajada_1 = i["Pico_1"]:i["Minimo"],
                Subida_2 = i["Minimo"]:i["Pico_2"], Bajada_2 = i["Pico_2"]:t2)
  Ini = c(Tramos$Subida_1[which.max(DMedia[Tramos$Subida_1])], Tramos$Subida_2[which.max(DMedia[Tramos$Subida_2])])
  Fin = c(Tramos$Bajada_1[which.min(DMedia[Tramos$Bajada_1])], Tramos$Bajada_2[which.min(DMedia[Tramos$Bajada_2])])
  Nivel = mean(Media)
  Ini_n = c(Tramos$Subida_1[Media[Tramos$Subida_1] >= Nivel][1], Tramos$Subida_2[Media[Tramos$Subida_2] >= Nivel][1])
  Fin_n = c(Tramos$Bajada_1[Media[Tramos$Bajada_1] < Nivel][1], Tramos$Bajada_2[Media[Tramos$Bajada_2] < Nivel][1])
  data.table(Fuente = Fuente, Temporada = c("Primera", "Segunda"),
             Inicio = Dia_a_fecha(Tiempo_comun[Ini]),
             Pico   = Dia_a_fecha(Tiempo_comun[i[c("Pico_1", "Pico_2")]]),
             Fin    = Dia_a_fecha(Tiempo_comun[Fin]),
             Duracion_dias  = Tiempo_comun[Fin] - Tiempo_comun[Ini],
             Subida_max     = round(DMedia[Ini], 3),
             Bajada_max     = round(DMedia[Fin], 3),
             Derivada_pico  = signif(DMedia[i[c("Pico_1", "Pico_2")]], 2),
             Fin_en_borde   = Fin == t2,
             Inicio_nivel   = Dia_a_fecha(Tiempo_comun[Ini_n]),
             Fin_nivel      = Dia_a_fecha(Tiempo_comun[Fin_n]),
             Nivel_mm_dia   = round(Nivel, 2),
             t_inicio = Tiempo_comun[Ini], t_fin = Tiempo_comun[Fin])
}
Tabla_Temporadas <- rbind(Temporadas("Estaciones", Media_est, DMedia_est),
                          Temporadas("Satelite", Media_sat, DMedia_sat))
Tabla_Temporadas
fwrite(Tabla_Temporadas, "Salidas_imputacion/temporadas_lluvia.csv")

windows(width = 12, height = 5)
par(mfrow = c(1, 2))
plot(Tiempo_comun, Media_est, type = "l", lwd = 3, col = "blue", ylim = range(0, Media_est, Media_sat),
     xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)", main = "Media funcional e inicio y fin de cada temporada")
lines(Tiempo_comun, Media_sat, lwd = 3, lty = 2, col = "red")
abline(v = Tabla_Temporadas[Fuente == "Estaciones", c(t_inicio, t_fin)], col = "blue", lty = 3)
abline(v = Tabla_Temporadas[Fuente == "Satelite", c(t_inicio, t_fin)], col = "red", lty = 3)
legend("bottomright", c("Estaciones", "Satélite"), col = c("blue", "red"), lty = c(1, 2), lwd = 3, bty = "n")
plot(Tiempo_comun, DMedia_est, type = "l", lwd = 3, col = "blue", ylim = range(DMedia_est, DMedia_sat),
     xlab = "Tiempo (días)", ylab = "mm/día por día", main = "Derivada de la media")
lines(Tiempo_comun, DMedia_sat, lwd = 3, lty = 2, col = "red")
abline(h = 0, lty = 3)
abline(v = c(10, 355), col = "gray60", lty = 2)
legend("bottomleft", c("Estaciones", "Satélite"), col = c("blue", "red"), lty = c(1, 2), lwd = 3, bty = "n")
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Derivadas.png", width = 1400, height = 580, res = 110)






#Paso 2.5: Profundidad funcional y boxplot funcional

#Profundidad de Fraiman-Muniz con depth.FM()
Est_fdata  = fdata(Curvas_comun_est, Tiempo_comun, rangval_Anual)
Sat_fdata  = fdata(Curvas_comun_sat, Tiempo_comun, rangval_Anual)
ProfFM_est = depth.FM(Est_fdata, trim = 0.25)
ProfFM_sat = depth.FM(Sat_fdata, trim = 0.25)

#Fraiman-Muniz paso a paso: promedio en t de 1 - |1/2 - F_N,t(x_i(t))|
FM_paso <- function(X){
  N = nrow(X)
  F = apply(X, 2, function(v) rank(v, ties.method = "max")/N)
  rowMeans(1 - abs(0.5 - F))
}
FM_est = FM_paso(Curvas_comun_est)
FM_sat = FM_paso(Curvas_comun_sat)

#depth.FM() da otra escala, se compara el orden de las curvas
cor(FM_est, ProfFM_est$dep, method = "spearman")
cor(FM_sat, ProfFM_sat$dep, method = "spearman")


#Mediana funcional (curva más profunda) y media recortada al 25 %
Mediana_est  = as.vector(ProfFM_est$median$data)
Mediana_sat  = as.vector(ProfFM_sat$median$data)
MediaRec_est = as.vector(ProfFM_est$mtrim$data)
MediaRec_sat = as.vector(ProfFM_sat$mtrim$data)

#Par de la mediana y curvas en la media recortada
Curvas_est[ProfFM_est$lmed]
Curvas_sat[ProfFM_sat$lmed]
c(Estaciones = length(ProfFM_est$ltrim), Satelite = length(ProfFM_sat$ltrim))

windows(width = 12, height = 5)
par(mfrow = c(1, 2))
matplot(Tiempo_comun, t(Curvas_comun_est), type = "l", lty = 1, col = rgb(0.5, 0.5, 0.5, 0.03),
        ylim = range(0, Media_est, Mediana_est, Media_sat, Mediana_sat)*c(1, 1.6),
        xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)", main = "Estaciones: medidas de centralidad")
lines(Tiempo_comun, Media_est, lwd = 3, col = "blue")
lines(Tiempo_comun, Mediana_est, lwd = 3, col = "red")
lines(Tiempo_comun, MediaRec_est, lwd = 3, lty = 2, col = "darkgreen")
legend("topleft", c("Media", "Mediana (FM)", "Media recortada 25 %"),
       col = c("blue", "red", "darkgreen"), lty = c(1, 1, 2), lwd = 3, bty = "n")
matplot(Tiempo_comun, t(Curvas_comun_sat), type = "l", lty = 1, col = rgb(0.5, 0.5, 0.5, 0.03),
        ylim = range(0, Media_est, Mediana_est, Media_sat, Mediana_sat)*c(1, 1.6),
        xlab = "Tiempo (días)", ylab = "Precipitación (mm/día)", main = "Satélite: medidas de centralidad")
lines(Tiempo_comun, Media_sat, lwd = 3, col = "blue")
lines(Tiempo_comun, Mediana_sat, lwd = 3, col = "red")
lines(Tiempo_comun, MediaRec_sat, lwd = 3, lty = 2, col = "darkgreen")
legend("topleft", c("Media", "Mediana (FM)", "Media recortada 25 %"),
       col = c("blue", "red", "darkgreen"), lty = c(1, 1, 2), lwd = 3, bty = "n")
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Centralidad.png", width = 1400, height = 580, res = 110)


#MBD paso a paso: [(r - 1)(N - r) + (N - 1)]/C(N, 2) en cada t, con r el rango de x_i(t)
MBD_paso <- function(X){
  N = nrow(X)
  R = apply(X, 2, rank)
  (rowMeans((R - 1)*(N - R)) + N - 1)/choose(N, 2)
}
MBD_est = MBD_paso(Curvas_comun_est)
MBD_sat = MBD_paso(Curvas_comun_sat)


#Boxplot funcional con la profundidad de bandas modificada (MBD)
#El eje y llega al percentil 99.9: con todo el rango las atípicas aplastan la región central
Rango_box <- function(X) c(min(0, quantile(X, 0.001)), quantile(X, 0.999))
windows(width = 12, height = 5)
par(mfrow = c(1, 2))
Box_est = fbplot(t(Curvas_comun_est), Tiempo_comun, method = "MBD", color = "gray80", outliercol = "firebrick", barcol = "steelblue", xlab = "Tiempo (días)",
                 ylab = "Precipitación (mm/día)", xlim = rangval_Anual, ylim = Rango_box(Curvas_comun_est))
title("Estaciones: boxplot funcional (MBD)")
Box_sat = fbplot(t(Curvas_comun_sat), Tiempo_comun, method = "MBD", color = "gray80", outliercol = "firebrick", barcol = "steelblue", xlab = "Tiempo (días)",
                 ylab = "Precipitación (mm/día)", xlim = rangval_Anual, ylim = Rango_box(Curvas_comun_sat))
title("Satélite: boxplot funcional (MBD)")
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Boxplot_Funcional.png", width = 1400, height = 580, res = 110)

#Comprobación: la diferencia debe ser cero
max(abs(Box_est$depth - MBD_est))
max(abs(Box_sat$depth - MBD_sat))


#Atípicas paso a paso: región central del 50 % más profundo y bigotes a 1.5 veces su ancho
Atipicas_paso <- function(X, Prof){
  Centro     = X[order(Prof, decreasing = TRUE)[1:ceiling(nrow(X)/2)], ]
  Inf_c      = apply(Centro, 2, min)
  Sup_c      = apply(Centro, 2, max)
  Bigote_inf = Inf_c - 1.5*(Sup_c - Inf_c)
  Bigote_sup = Sup_c + 1.5*(Sup_c - Inf_c)
  Arriba = rowSums(sweep(X, 2, Bigote_sup) >= 0)
  Abajo  = rowSums(sweep(X, 2, Bigote_inf) <= 0)
  list(Indices = which(Arriba + Abajo > 0), Arriba = Arriba, Abajo = Abajo)
}
Atip_est = Atipicas_paso(Curvas_comun_est, MBD_est)
Atip_sat = Atipicas_paso(Curvas_comun_sat, MBD_sat)

#Comprobación con fbplot()$outpoint
c(Estaciones = setequal(Atip_est$Indices, Box_est$outpoint), Satelite = setequal(Atip_sat$Indices, Box_sat$outpoint))
c(Estaciones = length(Atip_est$Indices), Satelite = length(Atip_sat$Indices))


#Curva más profunda con MBD y con FM
cbind(Profundidad = c("MBD", "FM"), Curvas_est[c(which.max(MBD_est), which.max(FM_est))])
cbind(Profundidad = c("MBD", "FM"), Curvas_sat[c(which.max(MBD_sat), which.max(FM_sat))])
#Concordancia entre FM y MBD
c(Estaciones = cor(FM_est, MBD_est, method = "spearman"), Satelite = cor(FM_sat, MBD_sat, method = "spearman"))

#Profundidad de cada curva
fwrite(cbind(Curvas_est[, .(STATION_ID, YEAR, pixel_id, Conjunto)], FM = FM_est, MBD = MBD_est),
       "Salidas_imputacion/profundidad_estaciones.csv")
fwrite(cbind(Curvas_sat[, .(pixel_id, YEAR, Conjunto)], FM = FM_sat, MBD = MBD_sat),
       "Salidas_imputacion/profundidad_satelite.csv")


#Curvas atípicas por estación y año
Atipicas_est <- Curvas_est[Atip_est$Indices, .(STATION_ID, YEAR, pixel_id, Conjunto)]
Atipicas_est[, `:=`(Dias_arriba    = Atip_est$Arriba[Atip_est$Indices],
                    Dias_abajo     = Atip_est$Abajo[Atip_est$Indices],
                    Total_anual_mm = round(Total_curva_est[Atip_est$Indices]),
                    MBD            = round(MBD_est[Atip_est$Indices], 4))]
Atipicas_est <- merge(Atipicas_est, Emparejamiento[, .(STATION_ID, Nombre, Altitud, Lon, Lat)], by = "STATION_ID")
Atipicas_est <- merge(Atipicas_est, Catalogo_CVC[, .(STATION_ID = CODIGO, Cuenca = BASIN)], by = "STATION_ID")
setorder(Atipicas_est, YEAR, STATION_ID)
fwrite(Atipicas_est, "Salidas_imputacion/atipicas_estaciones.csv")

#Atípicas por año
Por_año <- merge(Curvas_est[, .(Curvas = .N), by = YEAR], Atipicas_est[, .(Atipicas = .N), by = YEAR],
                 by = "YEAR", all.x = TRUE)
Por_año[is.na(Atipicas), Atipicas := 0]
Por_año[, Porcentaje := round(100*Atipicas/Curvas, 1)]
Por_año

#Atípicas por estación
Por_estacion <- merge(Atipicas_est[, .(Años_atipicos = .N), by = .(STATION_ID, Nombre, Cuenca, Altitud)],
                      Curvas_est[, .(Años = .N), by = STATION_ID], by = "STATION_ID")
Por_estacion[order(-Años_atipicos)]

#Atípicas del satélite por píxel y año
Atipicas_sat <- Curvas_sat[Atip_sat$Indices, .(pixel_id, YEAR, Conjunto)]
Atipicas_sat[, `:=`(Dias_arriba    = Atip_sat$Arriba[Atip_sat$Indices],
                    Dias_abajo     = Atip_sat$Abajo[Atip_sat$Indices],
                    Total_anual_mm = round(Total_curva_sat[Atip_sat$Indices]))]
Atipicas_sat <- merge(Atipicas_sat, Emparejamiento[STATION_ID %in% Curvas_est$STATION_ID, .(Estaciones_del_pixel = paste(Nombre, collapse = ", ")),
                                                   by = pixel_id], by = "pixel_id")
setorder(Atipicas_sat, YEAR, pixel_id)
fwrite(Atipicas_sat, "Salidas_imputacion/atipicas_satelite.csv")
Atipicas_sat[, .N, by = YEAR][order(YEAR)]
Atipicas_sat[, .N, by = .(pixel_id, Estaciones_del_pixel)][order(-N)]


#Anomalías: cada curva menos la media de su estación (o píxel)
#Las atípicas son casi siempre las mismas estaciones, entonces se repite el boxplot con las anomalías
Anom_est = Curvas_comun_est - rowsum(Curvas_comun_est, Curvas_est$STATION_ID)[Curvas_est$STATION_ID, ]/
           as.vector(table(Curvas_est$STATION_ID)[Curvas_est$STATION_ID])
Anom_sat = Curvas_comun_sat - rowsum(Curvas_comun_sat, Curvas_sat$pixel_id)[as.character(Curvas_sat$pixel_id), ]/
           as.vector(table(Curvas_sat$pixel_id)[as.character(Curvas_sat$pixel_id)])

#Comprobación: la anomalía media de cada estación es cero
max(abs(rowsum(Anom_est, Curvas_est$STATION_ID)))
max(abs(rowsum(Anom_sat, Curvas_sat$pixel_id)))

windows(width = 12, height = 5)
par(mfrow = c(1, 2))
BoxA_est = fbplot(t(Anom_est), Tiempo_comun, method = "MBD", color = "gray80", outliercol = "firebrick", barcol = "steelblue", xlab = "Tiempo (días)",
                  ylab = "Anomalía (mm/día)", xlim = rangval_Anual, ylim = quantile(Anom_est, c(0.001, 0.999)))
title("Estaciones: anomalías respecto a su estación")
BoxA_sat = fbplot(t(Anom_sat), Tiempo_comun, method = "MBD", color = "gray80", outliercol = "firebrick", barcol = "steelblue", xlab = "Tiempo (días)",
                  ylab = "Anomalía (mm/día)", xlim = rangval_Anual, ylim = quantile(Anom_sat, c(0.001, 0.999)))
title("Satélite: anomalías respecto a su píxel")
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Boxplot_Anomalias.png", width = 1400, height = 580, res = 110)

#Anomalía del total anual (%) respecto al promedio de la estación o píxel
Anual_est <- data.table(Curvas_est, Total = Total_curva_est, Atipica = seq_len(N_est) %in% BoxA_est$outpoint)
Anual_est[, Anomalia := 100*(Total/mean(Total) - 1), by = STATION_ID]
Anual_sat <- data.table(Curvas_sat, Total = Total_curva_sat, Atipica = seq_len(N_sat) %in% BoxA_sat$outpoint)
Anual_sat[, Anomalia := 100*(Total/mean(Total) - 1), by = pixel_id]

#Resumen por año
Años_ENSO <- merge(Anual_est[, .(Anomalia_est = round(median(Anomalia), 1), Atipicas_est = sum(Atipica),
                                 Curvas_est = .N), by = YEAR],
                   Anual_sat[, .(Anomalia_sat = round(median(Anomalia), 1), Atipicas_sat = sum(Atipica),
                                 Curvas_sat = .N), by = YEAR], by = "YEAR")
#Fase ENSO de cada año según el ONI (Bloque 7)
Años_ENSO <- merge(Años_ENSO, ONI_anual[, .(YEAR, ONI, ENSO = Fase)], by = "YEAR")
Años_ENSO
Años_ENSO[, .(Años = .N, Anomalia_est = round(mean(Anomalia_est), 1), Anomalia_sat = round(mean(Anomalia_sat), 1),
              Atipicas_est = sum(Atipicas_est), Atipicas_sat = sum(Atipicas_sat)), by = ENSO]
fwrite(Años_ENSO, "Salidas_imputacion/anomalias_por_año.csv")

windows(width = 10, height = 5)
Barras = t(as.matrix(Años_ENSO[, .(Anomalia_est, Anomalia_sat)]))
Pos = barplot(Barras, beside = TRUE, plot = FALSE)
Rango_anom = range(Barras) + c(-3, 3)
plot(NA, xlim = range(Pos) + c(-1, 1), ylim = Rango_anom, xaxt = "n", xlab = "", ylab = "Anomalía del total anual (%)",
     main = "Anomalía anual mediana por fuente (fondo: azul La Niña, gris El Niño)")
Fondo = Años_ENSO$ENSO != "Neutro"
rect(colMeans(Pos)[Fondo] - 1.5, Rango_anom[1] - 10, colMeans(Pos)[Fondo] + 1.5, Rango_anom[2] + 10,
     col = ifelse(Años_ENSO$ENSO[Fondo] == "La Niña", rgb(0, 0, 1, 0.12), rgb(0, 0, 0, 0.12)), border = NA)
barplot(Barras, beside = TRUE, names.arg = Años_ENSO$YEAR, col = c("steelblue", "darkorange"), las = 2, add = TRUE, axes = FALSE)
abline(h = 0)
legend("topleft", c("Estaciones", "Satélite"), fill = c("steelblue", "darkorange"), bty = "n")
dev.print(png, "Figuras_imputacion/Fig_Anomalias_ENSO.png", width = 1200, height = 600, res = 110)






#Ítem 3. Comparación funcional por fuente

#Paso 3.2: Diferencia d(t) = estación - satélite (mm/día)
#Cada par contra la curva de su píxel en el mismo año
D_dif  = Curvas_comun_est - Curvas_par_sat           #2870 x 365
N_par  = nrow(D_dif)
MediaD = colMeans(D_dif)
EE_D   = apply(D_dif, 2, sd)/sqrt(N_par)             #error estándar de la media
InfD   = MediaD - 2*EE_D
SupD   = MediaD + 2*EE_D

#Diferencia acumulada en el año (mm)
sum(MediaD)

#Sensibilidad al método: la misma diferencia con los coeficientes ponderados (prueba del paso 1.5)
#y con los totales observados (mm por año)
MediaD_MCP = colMeans(MCP_est$Coef %*% t(Phi_comun_est) - (MCP_sat$Coef %*% t(Phi_comun_sat))[Indice_sat, ])
c(Ordinarios = sum(MediaD), Ponderados = sum(MediaD_MCP), Observado = mean(Total_obs_est - Total_obs_sat[Indice_sat]))

#Días donde la banda no contiene el cero
Intervalos <- function(Indicador){
  r   = rle(Indicador)
  Fin = cumsum(r$lengths)
  Ini = Fin - r$lengths + 1
  data.table(Desde = Dia_a_fecha(Tiempo_comun[Ini[r$values]]), Hasta = Dia_a_fecha(Tiempo_comun[Fin[r$values]]),
             Dias = r$lengths[r$values], t_desde = Tiempo_comun[Ini[r$values]], t_hasta = Tiempo_comun[Fin[r$values]])
}
Intervalos_D <- rbind(cbind(Signo = "Estación por debajo del satélite", Intervalos(SupD < 0)),
                      cbind(Signo = "Estación por encima del satélite", Intervalos(InfD > 0)))
setorder(Intervalos_D, t_desde)
Intervalos_D
fwrite(Intervalos_D, "Salidas_imputacion/intervalos_diferencia.csv")
fwrite(data.table(t = Tiempo_comun, Fecha = Dia_a_fecha(Tiempo_comun), Media_d = MediaD, EE = EE_D,
                  Inferior = InfD, Superior = SupD), "Salidas_imputacion/diferencia_funcional.csv")

#Nota: la banda es punto a punto y supone curvas independientes (no lo son del todo: hay estaciones
#repetidas en varios años y píxeles compartidos), así que el EE real es algo mayor

windows(width = 12, height = 5)
par(mfrow = c(1, 2))
matplot(Tiempo_comun, t(D_dif), type = "l", lty = 1, col = rgb(0.5, 0.5, 0.5, 0.03),
        ylim = quantile(D_dif, c(0.001, 0.999)), xlab = "Tiempo (días)", ylab = "Estación - satélite (mm/día)",
        main = "Diferencia por par estación - año")
lines(Tiempo_comun, MediaD, lwd = 3, col = "darkred")
abline(h = 0, lty = 3)
plot(Tiempo_comun, MediaD, type = "n", ylim = range(InfD, SupD, 0), xlab = "Tiempo (días)", ylab = "mm/día",
     main = "Diferencia media ± 2 EE")
rect(Intervalos_D$t_desde - 0.5, par("usr")[3], Intervalos_D$t_hasta + 0.5, par("usr")[4],
     col = ifelse(Intervalos_D$Signo == "Estación por debajo del satélite", rgb(0, 0, 0, 0.12), rgb(0, 0, 1, 0.12)),
     border = NA)
lines(Tiempo_comun, MediaD, lwd = 3, col = "darkred")
lines(Tiempo_comun, InfD, lty = 2, col = "darkred")
lines(Tiempo_comun, SupD, lty = 2, col = "darkred")
lines(Tiempo_comun, MediaD_MCP, lwd = 2, lty = 3, col = "gray30")
abline(h = 0, lty = 3)
legend("topright", c("Ordinarios (curvas finales)", "Ponderados (prueba)"), col = c("darkred", "gray30"), lty = c(1, 3),
       lwd = c(3, 2), bty = "n", cex = 0.85)
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Diferencia.png", width = 1400, height = 580, res = 110)






#Paso 3.3: Correlación cruzada funcional y comparación de las covarianzas

#Curvas emparejadas: X = satélite y Y = estación, una por cada par estación - año
Sat_par_fd = fd(t(Coef_sat[Indice_sat, ]), Base_sat_fda)

#Covarianza cruzada paso a paso
#Cov_XY(s,t) = (N-1)^{-1} sum_i [x_i(s) - x_barra(s)][y_i(t) - y_barra(t)]
X_centrados = sweep(Curvas_par_sat, 2, colMeans(Curvas_par_sat))
Y_centrados = sweep(Curvas_comun_est, 2, Media_est)
CovCruzada  = t(X_centrados) %*% Y_centrados/(N_par - 1)       #filas: s del satélite; columnas: t de la estación
DE_par_sat  = apply(Curvas_par_sat, 2, sd)
CorCruzada  = CovCruzada/outer(DE_par_sat, DE_est)

#Comprobación con var.fd() y cor.fd(): la diferencia debe ser cero
max(abs(CovCruzada - eval.bifd(Tiempo_comun, Tiempo_comun, var.fd(Sat_par_fd, Est_fd))))
max(abs(CorCruzada - cor.fd(Tiempo_comun, Sat_par_fd, Tiempo_comun, Est_fd)))

#La correlación cruzada no es simétrica
max(abs(CorCruzada - t(CorCruzada)))

#Correlación en la diagonal: r_XY(t,t)
CorDiag = diag(CorCruzada)
summary(CorDiag)
c(Minimo = Dia_a_fecha(Tiempo_comun[which.min(CorDiag)]), Maximo = Dia_a_fecha(Tiempo_comun[which.max(CorDiag)]))

#Diagonal con las anomalías del paso 2.5 (cada curva menos la media de su estación o píxel)
#Así se quita el efecto del lugar y queda si el satélite sigue los años lluviosos y secos de cada sitio
Anom_par_sat = Anom_sat[Indice_sat, ]
CorDiag_anom = sapply(1:365, function(j) cor(Anom_par_sat[, j], Anom_est[, j]))
summary(CorDiag_anom)

windows(width = 12, height = 5.5)
par(mfrow = c(1, 2))
image(Tiempo_comun, Tiempo_comun, CorCruzada, zlim = c(-1, 1), col = Paleta_cor,
      xlab = "s: satélite (días)", ylab = "t: estaciones (días)", main = "Correlación cruzada")
contour(Tiempo_comun, Tiempo_comun, CorCruzada, levels = (-10:10)/10, labcex = 0.6, add = TRUE)
abline(0, 1, lty = 3)
plot(Tiempo_comun, CorDiag, type = "l", lwd = 3, col = "darkblue", ylim = c(-0.2, 1),
     xlab = "Tiempo (días)", ylab = expression(r[XY](t,t)), main = "Diagonal de la correlación cruzada")
lines(Tiempo_comun, CorDiag_anom, lwd = 3, lty = 2, col = "darkorange")
abline(h = 0, lty = 3)
legend("topleft", c("Curvas", "Anomalías respecto al sitio"), col = c("darkblue", "darkorange"),
       lty = c(1, 2), lwd = 3, bty = "n")
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Correlacion_Cruzada.png", width = 1400, height = 640, res = 110)

fwrite(data.table(t = Tiempo_comun, Fecha = Dia_a_fecha(Tiempo_comun), Cor_diagonal = CorDiag,
                  Cor_diagonal_anomalias = CorDiag_anom), "Salidas_imputacion/correlacion_cruzada_diagonal.csv")


#Comparación de la covarianza de cada fuente
#Se usan las curvas emparejadas para que las dos fuentes tengan los mismos 2870 pares
Cov_par_sat = t(X_centrados) %*% X_centrados/(N_par - 1)
Cor_par_sat = Cov_par_sat/outer(DE_par_sat, DE_par_sat)

#Varianza total (integral de la varianza), desviación estándar y correlación a 30, 90 y 180 días
Comparacion_Cov <- data.table(Fuente         = c("Estaciones", "Satelite"),
                              Varianza_total = round(c(sum(diag(Cov_est)), sum(diag(Cov_par_sat)))),
                              DE_media       = round(c(mean(DE_est), mean(DE_par_sat)), 2),
                              DE_minima      = round(c(min(DE_est), min(DE_par_sat)), 2),
                              DE_maxima      = round(c(max(DE_est), max(DE_par_sat)), 2),
                              rbind(Cor_rezago(Cor_est, c(30, 90, 180)), Cor_rezago(Cor_par_sat, c(30, 90, 180))))
setnames(Comparacion_Cov, 6:8, c("Cor_30_dias", "Cor_90_dias", "Cor_180_dias"))
Comparacion_Cov
fwrite(Comparacion_Cov, "Salidas_imputacion/comparacion_covarianza.csv")

#Razón de varianzas y de desviaciones estándar (satélite/estaciones)
sum(diag(Cov_par_sat))/sum(diag(Cov_est))
summary(DE_par_sat/DE_est)

#Diferencia relativa entre las dos superficies
sqrt(sum((Cov_est - Cov_par_sat)^2)/sum(Cov_est^2))
mean(abs(Cor_est - Cor_par_sat))

windows(width = 15, height = 5)
par(mfrow = c(1, 3))
plot(Tiempo_comun, DE_est, type = "l", lwd = 3, col = "blue", ylim = c(0, max(DE_est)),
     xlab = "Tiempo (días)", ylab = "mm/día", main = "Desviación estándar de cada fuente")
lines(Tiempo_comun, DE_par_sat, lwd = 3, lty = 2, col = "red")
legend("top", c("Estaciones", "Satélite"), col = c("blue", "red"), lty = c(1, 2), lwd = 3, bty = "n")
image(Tiempo_comun, Tiempo_comun, Cov_est - Cov_par_sat, zlim = c(-1, 1)*max(abs(Cov_est - Cov_par_sat)), col = Paleta_cor,
      xlab = "s (días)", ylab = "t (días)", main = "Covarianza: estaciones - satélite")
contour(Tiempo_comun, Tiempo_comun, Cov_est - Cov_par_sat, nlevels = 8, labcex = 0.6, add = TRUE)
image(Tiempo_comun, Tiempo_comun, Cor_est - Cor_par_sat, zlim = c(-0.6, 0.6), col = Paleta_cor,
      xlab = "s (días)", ylab = "t (días)", main = "Correlación: estaciones - satélite")
contour(Tiempo_comun, Tiempo_comun, Cor_est - Cor_par_sat, levels = (-6:6)/10, labcex = 0.6, add = TRUE)
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Comparacion_Covarianza.png", width = 1600, height = 540, res = 110)






#Paso 3.4: Comparación de la diferencia por subregión

#Grupos: se usa la subregión y no la altitud porque todas las estaciones tienen cuenca y
#coordenadas, mientras que la altitud falta en varias
#Pacífico: cuencas de la vertiente del Pacífico (Anchicayá, Dagua, Naya, Yurumanguí y Calima)
#Norte y Centro - sur: valle del río Cauca y sus laderas, separados en la latitud 4° N
#Garrapatas queda en el Norte: su río drena al San Juan, pero sus estaciones están en el cañón
#interior del norte del departamento y llueve en ellas como en el valle, no como en la costa
#(abajo se revisa con los datos y se prueba qué pasa si se cambia al Pacífico)
Cuencas_Pacifico = c("Anchicaya", "Dagua", "Naya", "Yurumangui", "Calima")
Estaciones_D <- merge(Emparejamiento[STATION_ID %in% Curvas_est$STATION_ID, .(STATION_ID, Nombre, Lon, Lat, Altitud)],
                      Catalogo_CVC[, .(STATION_ID = CODIGO, Cuenca = BASIN)], by = "STATION_ID")
Estaciones_D[, Subregion := ifelse(Cuenca %in% Cuencas_Pacifico, "Pacifico", ifelse(Lat >= 4, "Norte", "Centro - sur"))]
table(Estaciones_D$Subregion)

#Subregión de cada par estación - año
Subregion_par = Estaciones_D$Subregion[match(Curvas_est$STATION_ID, Estaciones_D$STATION_ID)]
table(Subregion_par)

#Diferencia media, banda y días sin el cero en cada grupo de curvas
Diferencia_grupos <- function(Grupo){
  Niveles = sort(unique(Grupo))
  Media = sapply(Niveles, function(g) colMeans(D_dif[Grupo == g, ]))                               #365 x grupos
  EE    = sapply(Niveles, function(g) apply(D_dif[Grupo == g, ], 2, sd)/sqrt(sum(Grupo == g)))
  Tabla = data.table(Grupo = Niveles, Pares = as.vector(table(Grupo)[Niveles]),
                     Estaciones    = sapply(Niveles, function(g) uniqueN(Curvas_est$STATION_ID[Grupo == g])),
                     Estacion_mm   = round(sapply(Niveles, function(g) mean(Total_curva_est[Grupo == g]))),
                     Satelite_mm   = round(sapply(Niveles, function(g) mean(Total_curva_sat[Indice_sat][Grupo == g]))),
                     Diferencia_mm = round(colSums(Media)),
                     Dias_debajo   = colSums(Media + 2*EE < 0),
                     Dias_encima   = colSums(Media - 2*EE > 0))
  Tabla[, Diferencia_pct := round(100*Diferencia_mm/Satelite_mm, 1)]
  Tramos = rbindlist(lapply(Niveles, function(g){
    Debajo = Intervalos(Media[, g] + 2*EE[, g] < 0)[, Signo := "Estación por debajo del satélite"]
    Encima = Intervalos(Media[, g] - 2*EE[, g] > 0)[, Signo := "Estación por encima del satélite"]
    rbind(Debajo, Encima)[, Grupo := g]
  }))
  list(Niveles = Niveles, Media = Media, EE = EE, Tabla = Tabla, Tramos = Tramos[order(Grupo, t_desde)])
}

#Gráfico de la diferencia media +- 2 EE de cada grupo
Graficar_grupos <- function(Dif, Colores, Titulo){
  plot(Tiempo_comun, Dif$Media[, 1], type = "n", ylim = range(Dif$Media - 2*Dif$EE, Dif$Media + 2*Dif$EE, 0),
       xlab = "Tiempo (días)", ylab = "Estación - satélite (mm/día)", main = Titulo)
  for(j in 1:length(Dif$Niveles)){
    lines(Tiempo_comun, Dif$Media[, j], lwd = 3, col = Colores[j])
    lines(Tiempo_comun, Dif$Media[, j] - 2*Dif$EE[, j], lty = 2, col = Colores[j])
    lines(Tiempo_comun, Dif$Media[, j] + 2*Dif$EE[, j], lty = 2, col = Colores[j])
  }
  abline(h = 0, lty = 3)
  legend("top", Dif$Niveles, col = Colores, lwd = 3, bty = "n", horiz = TRUE)
}

Dif_Subregion <- Diferencia_grupos(Subregion_par)
Dif_Subregion$Tabla
Dif_Subregion$Tramos[, .(Grupo, Signo, Desde, Hasta, Dias)]
fwrite(Dif_Subregion$Tabla, "Salidas_imputacion/diferencia_subregion.csv")
fwrite(Dif_Subregion$Tramos, "Salidas_imputacion/intervalos_diferencia_subregion.csv")


#ANOVA funcional por permutación, como en la clase de inferencia funcional (5 de noviembre)
#X trae una curva por unidad (filas). Bajo H0 las etiquetas de grupo son intercambiables, SST(t) no
#cambia al permutar y basta recalcular SSR(t). En la malla diaria cada día pesa lo mismo
F_permutacion <- function(X, Grupo, B = 5000){
  n = nrow(X)
  Z_grupo = sapply(sort(unique(Grupo)), function(g) as.numeric(Grupo == g))      #indicadoras n x k
  k = ncol(Z_grupo)
  Media_global = colMeans(X)
  SST_t = colSums(sweep(X, 2, Media_global)^2)
  SSR_t <- function(Z){
    Medias = t(Z) %*% X/colSums(Z)                                               #k x 365
    colSums(colSums(Z)*sweep(Medias, 2, Media_global)^2)
  }
  SSR_obs = SSR_t(Z_grupo)
  F_t = (SSR_obs/(k - 1))/((SST_t - SSR_obs)/(n - k))
  F_n = (sum(SSR_obs)/(k - 1))/(sum(SST_t - SSR_obs)/(n - k))
  SSR_nulos = sapply(1:B, function(b) SSR_t(Z_grupo[sample(n), ]))               #365 x B
  F_t_nulos = (SSR_nulos/(k - 1))/((SST_t - SSR_nulos)/(n - k))
  F_n_nulos = (colSums(SSR_nulos)/(k - 1))/(colSums(SST_t - SSR_nulos)/(n - k))
  Fmax_nulos = apply(F_t_nulos, 2, max)
  #Valores p por permutación (se suma 1 por la asignación observada)
  p_perm <- function(nulos, obs) (1 + sum(nulos >= obs - 1e-12))/(B + 1)
  c_maximo  = quantile(Fmax_nulos, 0.95)
  p_puntual = (1 + rowSums(F_t_nulos >= F_t - 1e-12))/(B + 1)
  list(F_t = F_t, c_maximo = c_maximo, Unidades = colSums(Z_grupo),
       p_F = p_perm(F_n_nulos, F_n), p_GPF = p_perm(colSums(F_t_nulos), sum(F_t)), p_Fmax = p_perm(Fmax_nulos, max(F_t)),
       Tramos_puntual = Intervalos(p_puntual < 0.05), Tramos = Intervalos(F_t > c_maximo))
}

#Las unidades son las estaciones: una curva de diferencia por estación (promedio de sus años)
#Así no se cuentan como independientes los años de una misma estación
D_estacion = rowsum(D_dif, Curvas_est$STATION_ID)/as.vector(table(Curvas_est$STATION_ID))
dim(D_estacion)                             #162 x 365
Subregion_estacion = Estaciones_D$Subregion[match(rownames(D_estacion), Estaciones_D$STATION_ID)]

Prueba_Subregion <- F_permutacion(D_estacion, Subregion_estacion)
Prueba_Subregion$Unidades
c(F_max = max(Prueba_Subregion$F_t), c_maximo = unname(Prueba_Subregion$c_maximo),
  p_F = Prueba_Subregion$p_F, p_GPF = Prueba_Subregion$p_GPF, p_Fmax = Prueba_Subregion$p_Fmax)
Prueba_Subregion$Tramos[, .(Desde, Hasta, Dias)]


#Diferencia de cada estación frente a su lluvia y a su altitud (mm por año)
Estaciones_D <- merge(Estaciones_D,
                      data.table(STATION_ID = Curvas_est$STATION_ID, E = Total_curva_est, S = Total_curva_sat[Indice_sat])[
                        , .(Años = .N, Estacion_mm = round(mean(E)), Satelite_mm = round(mean(S))), by = STATION_ID],
                      by = "STATION_ID")
Estaciones_D[, Diferencia_mm := Estacion_mm - Satelite_mm]
Estaciones_D[order(Diferencia_mm)][c(1:5, (.N - 4):.N), .(Nombre, Cuenca, Subregion, Altitud, Estacion_mm, Satelite_mm, Diferencia_mm)]
fwrite(Estaciones_D, "Salidas_imputacion/diferencia_por_estacion.csv")

#Garrapatas frente a las cuencas del Pacífico: lluvia anual media de sus estaciones
Estaciones_D[Cuenca == "Garrapatas" | Subregion == "Pacifico",
             .(Estaciones = .N, Estacion_mm = round(mean(Estacion_mm)), Satelite_mm = round(mean(Satelite_mm))), by = Cuenca]

#Sensibilidad: la misma prueba con Garrapatas en el Pacífico
Cuenca_estacion = Estaciones_D$Cuenca[match(rownames(D_estacion), Estaciones_D$STATION_ID)]
Prueba_Subregion_G <- F_permutacion(D_estacion, ifelse(Cuenca_estacion == "Garrapatas", "Pacifico", Subregion_estacion))
c(p_F = Prueba_Subregion_G$p_F, p_GPF = Prueba_Subregion_G$p_GPF, p_Fmax = Prueba_Subregion_G$p_Fmax)

#El satélite comprime el rango de la lluvia: recta satélite ~ estación con pendiente menor que 1
Recta_fuentes = lm(Satelite_mm ~ Estacion_mm, data = Estaciones_D)
coef(Recta_fuentes)
cor(Estaciones_D$Estacion_mm, Estaciones_D$Satelite_mm)

#Con la altitud (solo las estaciones que la tienen)
sum(!is.na(Estaciones_D$Altitud))
cor(Estaciones_D$Diferencia_mm, Estaciones_D$Altitud, use = "complete.obs")
cor(Estaciones_D[Subregion != "Pacifico"]$Diferencia_mm, Estaciones_D[Subregion != "Pacifico"]$Altitud, use = "complete.obs")

Colores_sub = c("darkgreen", "darkorange", "purple")       #Centro - sur, Norte, Pacifico
windows(width = 15, height = 5)
par(mfrow = c(1, 3))
Graficar_grupos(Dif_Subregion, Colores_sub, "Diferencia media ± 2 EE por subregión")
plot(Tiempo_comun, Prueba_Subregion$F_t, type = "l", lwd = 2, col = "red", xlab = "Tiempo (días)", ylab = "F(t)",
     main = "Prueba F funcional entre subregiones", ylim = c(0, max(Prueba_Subregion$F_t)*1.2))
abline(h = Prueba_Subregion$c_maximo, lty = 2)
legend("topleft", c("Estadístico observado", "Valor crítico máximo (0.05)"), col = c("red", "black"),
       lty = c(1, 2), lwd = c(2, 1), bty = "n")
plot(Estaciones_D$Estacion_mm, Estaciones_D$Satelite_mm, pch = 19,
     col = Colores_sub[match(Estaciones_D$Subregion, Dif_Subregion$Niveles)],
     xlab = "Estación (mm por año)", ylab = "Satélite (mm por año)", main = "Lluvia anual media por estación")
abline(0, 1, lty = 2)
abline(Recta_fuentes, col = "red")
legend("bottomright", c(Dif_Subregion$Niveles, "Igualdad", "Recta ajustada"), col = c(Colores_sub, "black", "red"),
       pch = c(19, 19, 19, NA, NA), lty = c(NA, NA, NA, 2, 1), bty = "n")
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Diferencia_Subregion.png", width = 1600, height = 540, res = 110)






#Paso 3.5: Comparación de la diferencia por fase ENSO

#Fase de cada par según el ONI del año (Bloque 7)
Fase_par = ONI_anual$Fase[match(Curvas_est$YEAR, ONI_anual$YEAR)]
table(Fase_par)

Dif_ENSO <- Diferencia_grupos(Fase_par)
Dif_ENSO$Tabla
Dif_ENSO$Tramos[, .(Grupo, Signo, Desde, Hasta, Dias)]
fwrite(Dif_ENSO$Tabla, "Salidas_imputacion/diferencia_ENSO.csv")
fwrite(Dif_ENSO$Tramos, "Salidas_imputacion/intervalos_diferencia_ENSO.csv")

#Las unidades son los años: una curva de diferencia por año (promedio de sus estaciones)
D_año = rowsum(D_dif, Curvas_est$YEAR)/as.vector(table(Curvas_est$YEAR))
dim(D_año)                                  #23 x 365
Fase_año = ONI_anual$Fase[match(rownames(D_año), ONI_anual$YEAR)]

Prueba_ENSO <- F_permutacion(D_año, Fase_año)
Prueba_ENSO$Unidades
c(F_max = max(Prueba_ENSO$F_t), c_maximo = unname(Prueba_ENSO$c_maximo),
  p_F = Prueba_ENSO$p_F, p_GPF = Prueba_ENSO$p_GPF, p_Fmax = Prueba_ENSO$p_Fmax)
Prueba_ENSO$Tramos_puntual[, .(Desde, Hasta, Dias)]     #sin corrección por multiplicidad
Prueba_ENSO$Tramos[, .(Desde, Hasta, Dias)]             #criterio del máximo

#Sensibilidad: 1998 como neutro (empieza en El Niño y termina en La Niña)
Prueba_ENSO_98 <- F_permutacion(D_año, replace(Fase_año, rownames(D_año) == "1998", "Neutro"))
c(p_F = Prueba_ENSO_98$p_F, p_GPF = Prueba_ENSO_98$p_GPF, p_Fmax = Prueba_ENSO_98$p_Fmax)

#Solo los años El Niño contra los años La Niña
Extremos = Fase_año != "Neutro"
Prueba_Niño_Niña <- F_permutacion(D_año[Extremos, ], Fase_año[Extremos])
c(F_max = max(Prueba_Niño_Niña$F_t), c_maximo = unname(Prueba_Niño_Niña$c_maximo),
  p_F = Prueba_Niño_Niña$p_F, p_GPF = Prueba_Niño_Niña$p_GPF, p_Fmax = Prueba_Niño_Niña$p_Fmax)
Prueba_Niño_Niña$Tramos[, .(Desde, Hasta, Dias)]

#Sin depender del umbral: ONI del año frente a la diferencia anual media (mm) y a la lluvia de cada fuente
Diferencia_año <- merge(ONI_anual, data.table(YEAR = as.integer(rownames(D_año)), Diferencia_mm = round(rowSums(D_año))),
                        by = "YEAR")
Diferencia_año <- merge(Diferencia_año, Años_ENSO[, .(YEAR, Anomalia_est, Anomalia_sat)], by = "YEAR")
Diferencia_año
c(Diferencia = cor(Diferencia_año$ONI, Diferencia_año$Diferencia_mm),
  Estaciones = cor(Diferencia_año$ONI, Diferencia_año$Anomalia_est),
  Satelite   = cor(Diferencia_año$ONI, Diferencia_año$Anomalia_sat))
fwrite(Diferencia_año, "Salidas_imputacion/diferencia_por_año.csv")

Colores_ENSO = c("red", "blue", "gray40")                  #El Niño, La Niña, Neutro
windows(width = 15, height = 5)
par(mfrow = c(1, 3))
Graficar_grupos(Dif_ENSO, Colores_ENSO, "Diferencia media ± 2 EE por fase ENSO")
plot(Tiempo_comun, Prueba_ENSO$F_t, type = "l", lwd = 2, col = "red", xlab = "Tiempo (días)", ylab = "F(t)",
     main = "Prueba F funcional entre fases", ylim = c(0, max(Prueba_ENSO$F_t, Prueba_ENSO$c_maximo)*1.2))
abline(h = Prueba_ENSO$c_maximo, lty = 2)
legend("topleft", c("Estadístico observado", "Valor crítico máximo (0.05)"), col = c("red", "black"),
       lty = c(1, 2), lwd = c(2, 1), bty = "n")
plot(Diferencia_año$ONI, Diferencia_año$Diferencia_mm, type = "n", xlab = "ONI promedio del año",
     ylab = "Estación - satélite (mm por año)", main = "Diferencia anual frente al ONI")
text(Diferencia_año$ONI, Diferencia_año$Diferencia_mm, substr(Diferencia_año$YEAR, 3, 4),
     col = Colores_ENSO[match(Diferencia_año$Fase, Dif_ENSO$Niveles)])
abline(h = 0, lty = 3)
abline(lm(Diferencia_mm ~ ONI, data = Diferencia_año), col = "black")
legend("topright", Dif_ENSO$Niveles, text.col = Colores_ENSO, bty = "n")
par(mfrow = c(1, 1))
dev.print(png, "Figuras_imputacion/Fig_Diferencia_ENSO.png", width = 1600, height = 540, res = 110)






#Verificación de los ítems 2 y 3 contra la rúbrica
Resumen_Item2 <- data.table(
  Criterio  = c("Media funcional paso a paso y verificada con la librería",
                "Varianza y desviación estándar funcional",
                "Covarianza y correlación funcional de cada fuente",
                "Derivadas e inicio y fin de las temporadas de lluvia",
                "Profundidad (FM y MBD), boxplot funcional y atípicas por estación y año",
                "Todo el periodo 1981-2003"),
  Evidencia = c(sprintf("Diferencia con mean.fd(): %.1e (estaciones) y %.1e (satélite). Fig_Media_DE.png",
                        max(abs(Media_est - MediaEst_fda)), max(abs(Media_sat - MediaSat_fda))),
                sprintf("DE media de %.2f (estaciones) y %.2f (satélite) mm/día. Salidas_imputacion/media_DE_funcional.csv",
                        mean(DE_est), mean(DE_sat)),
                sprintf("Diferencia con var.fd(): %.1e y con cor.fd(): %.1e. Fig_Covarianza_Correlacion.png",
                        max(abs(Cov_est - eval.bifd(Tiempo_comun, Tiempo_comun, var.fd(Est_fd)))),
                        max(abs(Cor_est - cor.fd(Tiempo_comun, Est_fd)))),
                sprintf("Estaciones: %s a %s y %s a %s. Salidas_imputacion/temporadas_lluvia.csv y Fig_Derivadas.png",
                        Tabla_Temporadas$Inicio[1], Tabla_Temporadas$Fin[1], Tabla_Temporadas$Inicio[2], Tabla_Temporadas$Fin_nivel[2]),
                sprintf("MBD contra fbplot(): %.1e. Atípicas: %d de estación y %d de satélite. Salidas_imputacion/atipicas_estaciones.csv",
                        max(abs(Box_est$depth - MBD_est)), length(Atip_est$Indices), length(Atip_sat$Indices)),
                sprintf("%d curvas de estación y %d de satélite; entre %d y %d curvas de estación por año",
                        N_est, N_sat, min(Curvas_por_año["Estaciones", ]), max(Curvas_por_año["Estaciones", ]))))
Resumen_Item2
fwrite(Resumen_Item2, "Salidas_imputacion/resumen_item2.csv")

Resumen_Item3 <- data.table(
  Criterio  = c("Emparejamiento estación - píxel documentado",
                "Diferencia d(t) = estación - satélite con banda e intervalos sin el cero",
                "Correlación cruzada funcional y su diagonal",
                "Comparación de la covarianza de cada fuente",
                "Comparación por subregión",
                "Comparación por fase ENSO (ONI de la NOAA, criterio del IDEAM)"),
  Evidencia = c(sprintf("%d estaciones en %d píxeles, %d comparten píxel. Salidas_imputacion/emparejamiento.csv",
                        sum(!is.na(Emparejamiento$pixel_id)), uniqueN(Emparejamiento[!is.na(pixel_id)]$pixel_id),
                        sum(Emparejamiento$Estaciones_en_pixel > 1, na.rm = TRUE)),
                sprintf("Diferencia anual de %.0f mm; %d días con la estación por debajo y %d por encima. Fig_Diferencia.png",
                        sum(MediaD), sum(SupD < 0), sum(InfD > 0)),
                sprintf("Diferencia con cor.fd(): %.1e. Diagonal entre %.2f y %.2f. Fig_Correlacion_Cruzada.png",
                        max(abs(CorCruzada - cor.fd(Tiempo_comun, Sat_par_fd, Tiempo_comun, Est_fd))), min(CorDiag), max(CorDiag)),
                sprintf("Varianza total del satélite: %.0f %% de la de estaciones. Fig_Comparacion_Covarianza.png",
                        100*sum(diag(Cov_par_sat))/sum(diag(Cov_est))),
                sprintf("%s. Valor p: %.4f (F), %.4f (GPF) y %.4f (F máximo). Fig_Diferencia_Subregion.png",
                        paste(Dif_Subregion$Tabla$Grupo, Dif_Subregion$Tabla$Diferencia_mm, "mm", collapse = "; "),
                        Prueba_Subregion$p_F, Prueba_Subregion$p_GPF, Prueba_Subregion$p_Fmax),
                sprintf("%s. Valor p: %.3f (F), %.3f (GPF) y %.3f (F máximo). Correlación con el ONI: %.2f. Fig_Diferencia_ENSO.png",
                        paste(Dif_ENSO$Tabla$Grupo, Dif_ENSO$Tabla$Diferencia_mm, "mm", collapse = "; "),
                        Prueba_ENSO$p_F, Prueba_ENSO$p_GPF, Prueba_ENSO$p_Fmax, cor(Diferencia_año$ONI, Diferencia_año$Diferencia_mm))))
Resumen_Item3
fwrite(Resumen_Item3, "Salidas_imputacion/resumen_item3.csv")






#Archivos para Python: base B-spline y coeficientes de cada curva

#Nodos interiores de cada base (los mismos de Base_Bspline)
Nodos_interiores <- function(K) seq(rangval_Anual[1], rangval_Anual[2], length.out = K - norden + 2)[-c(1, K - norden + 2)]

writeLines(c("{",
             sprintf('  "orden": %d,', norden),
             sprintf('  "dominio": [%s],', paste(rangval_Anual, collapse = ", ")),
             sprintf('  "estaciones": {"K": %d, "nodos_interiores": [%s]},', K_est, paste(Nodos_interiores(K_est), collapse = ", ")),
             sprintf('  "satelite": {"K": %d, "nodos_interiores": [%s]}', K_sat, paste(Nodos_interiores(K_sat), collapse = ", ")),
             "}"), "Salidas_imputacion/base_bspline.json")

Coef_est_dt <- data.table(Curvas_est[, .(STATION_ID, YEAR, pixel_id, Conjunto)], Coef_est)
Coef_sat_dt <- data.table(Curvas_sat[, .(pixel_id, YEAR, Conjunto)], Coef_sat)
setnames(Coef_est_dt, 5:ncol(Coef_est_dt), paste0("c", 1:K_est))
setnames(Coef_sat_dt, 4:ncol(Coef_sat_dt), paste0("c", 1:K_sat))
fwrite(Coef_est_dt, "Salidas_imputacion/coeficientes_estaciones.csv")
fwrite(Coef_sat_dt, "Salidas_imputacion/coeficientes_satelite.csv")

#Curva de control para verificar la base en Python: primera curva del satélite evaluada en la malla diaria
fwrite(data.table(t = Tiempo_comun, valor = Curvas_comun_sat[1, ]), "Salidas_imputacion/curva_control_satelite.csv")






#Ítem 4.4. Prueba de Wilcoxon pareada por estación

#Orden de trabajo: este script, luego Proyecto_Final_IA.py y otra vez este bloque, que lee las
#métricas por curva que calcula Python (promedio de las 3 semillas de cada cabeza)
#Si Python todavía no se ha corrido, el bloque se salta
Archivo_metricas = "Salidas_IA/metricas_por_curva.csv"

if(file.exists(Archivo_metricas)){
  Metricas <- fread(Archivo_metricas, colClasses = c(STATION_ID = "character"))

  #Cada métrica se promedia por estación: las curvas de una misma estación no son independientes
  #Sesgo: valor absoluto del sesgo del acumulado anual (%). Picos: error medio de los dos picos (días)
  Metricas_est <- Metricas[, .(L2_rel    = mean(L2_rel),
                               KGE       = mean(KGE),
                               Sesgo_abs = mean(abs(sesgo_pct)),
                               Err_picos = mean((err_pico1 + err_pico2)/2)),
                           by = .(conjunto, modelo, STATION_ID)]

  #Las mismas comparaciones del bootstrap de Python
  Comparaciones = list(c("kan", "linear"), c("mlp", "linear"), c("mlp", "kan"),
                       c("linear", "satelite_sin_corregir"), c("kan", "satelite_sin_corregir"),
                       c("mlp", "satelite_sin_corregir"), c("linear", "climatologia_estacion"),
                       c("kan", "climatologia_estacion"), c("mlp", "climatologia_estacion"))
  Metricas_prueba = c("L2_rel", "KGE", "Sesgo_abs", "Err_picos")
  Menor_es_mejor  = c(L2_rel = TRUE, KGE = FALSE, Sesgo_abs = TRUE, Err_picos = TRUE)

  #Prueba de rangos con signo de Wilcoxon sobre las diferencias por estación (modelo A - modelo B)
  Wilcoxon <- rbindlist(lapply(c("validacion", "prueba"), function(Conj){
    rbindlist(lapply(Metricas_prueba, function(Met){
      Tabla = dcast(Metricas_est[conjunto == Conj], STATION_ID ~ modelo, value.var = Met)
      rbindlist(lapply(Comparaciones, function(Par){
        d = Tabla[[Par[1]]] - Tabla[[Par[2]]]
        d = d[!is.na(d)]
        Prueba = wilcox.test(d, mu = 0, exact = FALSE)
        data.table(Conjunto = Conj, Metrica = Met, Comparacion = paste(Par[1], "-", Par[2]),
                   Estaciones = length(d), Mediana_dif = median(d), Media_dif = mean(d),
                   Estadistico_V = unname(Prueba$statistic), Valor_p = Prueba$p.value)
      }))
    }))
  }))

  #Ajuste de Holm por las nueve comparaciones de cada conjunto y métrica, como en las
  #comparaciones por pares del ANOVA funcional
  Wilcoxon[, Valor_p_Holm := p.adjust(Valor_p, method = "holm"), by = .(Conjunto, Metrica)]
  Wilcoxon[, Significativa := Valor_p_Holm < 0.05]

  #Qué modelo es mejor cuando la diferencia es significativa (según el signo de la mediana)
  Wilcoxon[, Mejor := ifelse(!Significativa, "Ninguno",
                             ifelse((Mediana_dif < 0) == Menor_es_mejor[Metrica],
                                    sub(" - .*", "", Comparacion), sub(".* - ", "", Comparacion)))]

  print(Wilcoxon[Metrica == "L2_rel"])
  print(Wilcoxon[Conjunto == "prueba", .(Metrica, Comparacion, Mediana_dif = round(Mediana_dif, 4),
                                   Valor_p_Holm = signif(Valor_p_Holm, 3), Mejor)])
  fwrite(Wilcoxon, "Salidas_imputacion/wilcoxon_por_estacion.csv")
} else {
  message("Falta ", Archivo_metricas, ": correr Proyecto_Final_IA.py y después este bloque")
}






#Índice: dónde sale cada figura y cada tabla
Indice_salidas <- rbind(
  data.table(Paso = "Bloque 4", Archivo = "Salidas_imputacion/emparejamiento.csv"),
  data.table(Paso = "Bloque 5", Archivo = c("Salidas_imputacion/trazabilidad.csv", "Salidas_imputacion/verificacion_relleno.csv", "Salidas_imputacion/semivariograma_relleno.csv",
                                         "Figuras_imputacion/Fig_Verificacion_Relleno.png")),
  data.table(Paso = "Bloque 6", Archivo = c("Salidas_imputacion/particion.csv", "Salidas_imputacion/estaciones_depuradas.csv", "Salidas_imputacion/satelite_pentada.csv")),
  data.table(Paso = "Bloque 7", Archivo = "Salidas_imputacion/oni_anual.csv"),
  data.table(Paso = "1.3",      Archivo = "Salidas_imputacion/verificacion_S.csv"),
  data.table(Paso = "1.4",      Archivo = "Figuras_imputacion/Fig_Seleccion_K.pdf"),
  data.table(Paso = "1.5",      Archivo = c("Figuras_imputacion/Fig_Varianza_Residual.pdf", "Salidas_imputacion/totales_por_metodo.csv")),
  data.table(Paso = "1.6",      Archivo = "Figuras_imputacion/Fig_Bandas.pdf"),
  data.table(Paso = "1.7",      Archivo = c("Figuras_imputacion/Fig_Base_Phi.pdf", "Figuras_imputacion/Fig_Curvas_Dominio_Comun.png", "Salidas_imputacion/resumen_item1.csv")),
  data.table(Paso = "2.2",      Archivo = c("Figuras_imputacion/Fig_Media_DE.png", "Salidas_imputacion/resumen_media_DE.csv", "Salidas_imputacion/media_DE_funcional.csv")),
  data.table(Paso = "2.3",      Archivo = "Figuras_imputacion/Fig_Covarianza_Correlacion.png"),
  data.table(Paso = "2.4",      Archivo = c("Figuras_imputacion/Fig_Derivadas.png", "Salidas_imputacion/temporadas_lluvia.csv")),
  data.table(Paso = "2.5",      Archivo = c("Figuras_imputacion/Fig_Centralidad.png", "Figuras_imputacion/Fig_Boxplot_Funcional.png", "Figuras_imputacion/Fig_Boxplot_Anomalias.png",
                                            "Figuras_imputacion/Fig_Anomalias_ENSO.png", "Salidas_imputacion/profundidad_estaciones.csv", "Salidas_imputacion/profundidad_satelite.csv",
                                            "Salidas_imputacion/atipicas_estaciones.csv", "Salidas_imputacion/atipicas_satelite.csv", "Salidas_imputacion/anomalias_por_año.csv")),
  data.table(Paso = "3.2",      Archivo = c("Figuras_imputacion/Fig_Diferencia.png", "Salidas_imputacion/diferencia_funcional.csv", "Salidas_imputacion/intervalos_diferencia.csv")),
  data.table(Paso = "3.3",      Archivo = c("Figuras_imputacion/Fig_Correlacion_Cruzada.png", "Figuras_imputacion/Fig_Comparacion_Covarianza.png",
                                            "Salidas_imputacion/correlacion_cruzada_diagonal.csv", "Salidas_imputacion/comparacion_covarianza.csv")),
  data.table(Paso = "3.4",      Archivo = c("Figuras_imputacion/Fig_Diferencia_Subregion.png", "Salidas_imputacion/diferencia_subregion.csv",
                                            "Salidas_imputacion/intervalos_diferencia_subregion.csv", "Salidas_imputacion/diferencia_por_estacion.csv")),
  data.table(Paso = "3.5",      Archivo = c("Figuras_imputacion/Fig_Diferencia_ENSO.png", "Salidas_imputacion/diferencia_ENSO.csv",
                                            "Salidas_imputacion/intervalos_diferencia_ENSO.csv", "Salidas_imputacion/diferencia_por_año.csv")),
  data.table(Paso = "Cierre",   Archivo = c("Salidas_imputacion/resumen_item2.csv", "Salidas_imputacion/resumen_item3.csv", "Salidas_imputacion/base_bspline.json",
                                            "Salidas_imputacion/coeficientes_estaciones.csv", "Salidas_imputacion/coeficientes_satelite.csv",
                                            "Salidas_imputacion/curva_control_satelite.csv")),
  data.table(Paso = "4.4",      Archivo = "Salidas_imputacion/wilcoxon_por_estacion.csv"))
Indice_salidas[, Existe := file.exists(Archivo)]
Indice_salidas
fwrite(Indice_salidas, "Salidas_imputacion/indice_figuras_tablas.csv")






#Versiones de R y de los paquetes
sessionInfo()
writeLines(capture.output(sessionInfo()), "Salidas_imputacion/versiones_paquetes.txt")
