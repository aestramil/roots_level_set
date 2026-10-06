# ============================================================
# EXPERIMENTO 2D - VERSION LIMPIA PARA TESIS
#
# Metodos en la tabla:
#   1. Uniforme puro (iid + banda)
#   2. RW isotropico + MH
#   3. Copula adaptativa + MH
#   4. Copula + MH + KNNmin aumentado (K = 5, 10, 100)
#   5. GP Random Straddle
#
# Figuras (4 paneles por funcion):
#   1. RW isotropico + MH
#   2. Copula adaptativa + MH
#   3. Copula + MH + KNNmin aumentado(K_grafico)
#   4. GP Random Straddle
#
## KNN local con preservacion espacial:
#   - busca, para cada estado MH, el mejor vecino entre sus K vecinos
#     segun g(x)=|f(x)-c|;
#   - solo incorpora el vecino si mejora g;
#   - NO reemplaza ni elimina el punto MH original;
#   - la estimacion final es la UNION de los puntos MH originales y
#     los vecinos mejorados.
#
# Consecuencia: con la metrica de cobertura usada aqui, agregar puntos
# no puede disminuir la cobertura respecto de Copula + MH.
# ============================================================

rm(list = ls())

# ============================================================
# 0. PAQUETES
# ============================================================

library(VineCopula)
library(FNN)
library(DiceKriging)

# ============================================================
# 1. FUNCIONES DE PRUEBA
# ============================================================

ff1 <- list(
  f = function(x) {
    1 - x[1]^2 - sin(x[2]^2)
  },
  L = -2,
  U = 2,
  c = 0,
  nombre = "ff1"
)

ff2 <- list(
  f = function(x) {
    (0.25 - (x[1] - 0.5)^2 - x[2]^2) *
      (0.25 - (x[1] + 0.5)^2 - x[2]^2)
  },
  L = -1,
  U = 1.5,
  c = 0,
  nombre = "ff2"
)

ff3 <- list(
  f = function(x) {
    1 - x[1]^2 - x[1] * sin(x[2]^2)
  },
  L = -2,
  U = 2,
  c = 0,
  nombre = "ff3"
)

ff4 <- list(
  f = function(x) {
    sin(10 * x[1]) +
      cos(4 * x[2]) -
      cos(3 * x[1] * x[2])
  },
  L = 0,
  U = 2,
  c = 0,
  nombre = "ff4"
)

funciones <- list(ff1, ff2, ff3, ff4)

# ============================================================
# 2. PARAMETROS
# ============================================================

# Presupuesto de los metodos MCMC:
# 500 evaluaciones iniciales + 4500 propuestas = 5000 evaluaciones.
N_f <- 500
B_total <- 5000
B_mh <- B_total - N_f

# Copula adaptativa
sigma_g <- 0.10
eta_rho <- 0.02
rho_max <- 1 - eta_rho

# Mezcla de propuesta copular
w_cop <- 0.90
w_unif <- 0.10

# Target MH
# pi_tau(x) proporcional a exp(-|f(x)-c|/tau_mh)
tau_mh <- 0.05
burn_frac <- 0.10

# RW isotropico
# Se interpreta como fraccion del ancho U-L.
sigma_rw_rel <- 0.05

# KNNmin aumentado
K_valores <- c(5, 10, 100)
K_grafico <- 10

# KNN local: solo se reemplazan puntos espacialmente redundantes.
# radio_local se define con la distancia al k_local-esimo vecino de X_MH.
k_local_proteccion <- 5
min_vecinos_proteccion <- 4
lambda_radio <- 1.0

# Metricas
# Las distancias se calculan tras normalizar el dominio a [0,1]^2.
delta_geom <- 0.01
ngrid <- 500
n_arco <- 6000

# GP Random Straddle
N_gp <- 200
N_init_gp <- 30
N_cand_gp <- 3000
N_pred_gp <- 20000

# ============================================================
# 3. CURVA VERDADERA / CURVA DE REFERENCIA
#    SOLO PARA EVALUACION, NO FORMA PARTE DE LOS ALGORITMOS
# ============================================================

extraer_curva <- function(f, L, U, c0, ngrid = 500) {
  
  xs <- seq(L, U, length.out = ngrid)
  ys <- seq(L, U, length.out = ngrid)
  
  Z <- outer(
    xs,
    ys,
    Vectorize(function(x, y) f(c(x, y)))
  )
  
  cls <- contourLines(xs, ys, Z, levels = c0)
  
  if (length(cls) == 0) {
    stop("No se encontro la curva de nivel.")
  }
  
  df <- do.call(
    rbind,
    lapply(seq_along(cls), function(j) {
      data.frame(
        x = cls[[j]]$x,
        y = cls[[j]]$y,
        componente = j
      )
    })
  )
  
  P <- as.matrix(df[, c("x", "y")])
  colnames(P) <- c("x1", "x2")
  
  list(df = df, P = P)
}

# ============================================================
# 4. REMUESTREO DE LA CURVA POR LONGITUD DE ARCO
# ============================================================

remuestrear_arco <- function(df, n = 6000) {
  
  componentes <- split(df, df$componente)
  
  longitudes <- sapply(componentes, function(dd) {
    if (nrow(dd) < 2) return(0)
    sum(sqrt(diff(dd$x)^2 + diff(dd$y)^2))
  })
  
  total <- sum(longitudes)
  if (total <= 0) stop("Longitud total de la curva igual a cero.")
  
  n_comp <- pmax(2, round(n * longitudes / total))
  salida <- list()
  
  for (j in seq_along(componentes)) {
    
    dd <- componentes[[j]]
    if (nrow(dd) < 2) next
    
    ds <- sqrt(diff(dd$x)^2 + diff(dd$y)^2)
    s <- c(0, cumsum(ds))
    if (max(s) == 0) next
    
    ss <- seq(0, max(s), length.out = n_comp[j])
    
    xx <- approx(s, dd$x, xout = ss, ties = "ordered")$y
    yy <- approx(s, dd$y, xout = ss, ties = "ordered")$y
    
    salida[[length(salida) + 1]] <- cbind(x1 = xx, x2 = yy)
  }
  
  do.call(rbind, salida)
}

# ============================================================
# 5. METRICAS GEOMETRICAS
#
# precision:
#   proporcion de puntos estimados a distancia <= delta de la curva.
# cobertura:
#   proporcion de la curva a distancia <= delta de la estimacion.
# F1:
#   media armonica de precision y cobertura.
# ============================================================

metricas_geom <- function(Xhat, Lref, L, U, delta = 0.01) {
  
  Xhat <- as.matrix(Xhat)
  Lref <- as.matrix(Lref)
  
  if (nrow(Xhat) == 0 || nrow(Lref) == 0) {
    return(c(
      precision = NA_real_,
      cobertura = NA_real_,
      F1 = NA_real_
    ))
  }
  
  # Normalizacion a [0,1]^2
  Xn <- (Xhat - L) / (U - L)
  Ln <- (Lref - L) / (U - L)
  
  # estimacion -> curva verdadera
  d_est_true <- FNN::get.knnx(
    data = Ln,
    query = Xn,
    k = 1
  )$nn.dist[, 1]
  
  # curva verdadera -> estimacion
  d_true_est <- FNN::get.knnx(
    data = Xn,
    query = Ln,
    k = 1
  )$nn.dist[, 1]
  
  precision <- mean(d_est_true <= delta)
  cobertura <- mean(d_true_est <= delta)
  
  if (precision + cobertura == 0) {
    F1 <- 0
  } else {
    F1 <- 2 * precision * cobertura / (precision + cobertura)
  }
  
  c(
    precision = precision,
    cobertura = cobertura,
    F1 = F1
  )
}

# ============================================================
# 6. HERRAMIENTAS DE LA PROPUESTA COPULAR ADAPTATIVA
# ============================================================

a_unif <- function(x, L, U) {
  u <- (x - L) / (U - L)
  pmin(pmax(u, 1e-10), 1 - 1e-10)
}

rho_adaptativo <- function(g,
                           sigma = sigma_g,
                           rho_lim = rho_max) {
  rho <- rho_lim * exp(-(g / sigma)^2)
  pmin(pmax(rho, 0), rho_lim)
}

dens_cop_gauss <- function(u, v, rho) {
  
  if (!is.finite(rho) || abs(rho) < 1e-12) return(1)
  
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  v <- pmin(pmax(v, 1e-10), 1 - 1e-10)
  
  val <- VineCopula::BiCopPDF(
    u1 = u,
    u2 = v,
    family = 1,
    par = rho
  )
  
  if (!is.finite(val) || val < 0) return(0)
  val
}

rcond_gauss <- function(u, rho) {
  
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  
  if (!is.finite(rho) || abs(rho) < 1e-12) {
    return(runif(1))
  }
  
  z1 <- qnorm(u)
  z2 <- rnorm(
    1,
    mean = rho * z1,
    sd = sqrt(max(1 - rho^2, 1e-12))
  )
  
  v <- pnorm(z2)
  pmin(pmax(v, 1e-10), 1 - 1e-10)
}

# q_C(y|x) = prod_j c_rho(x)(u_xj,u_yj)/(U-L)^2
q_cop_cond <- function(y, x, gx, L, U) {
  
  if (
    any(!is.finite(y)) ||
    any(y < L) ||
    any(y > U)
  ) {
    return(0)
  }
  
  rho_x <- rho_adaptativo(gx)
  ux <- a_unif(x, L, U)
  uy <- a_unif(y, L, U)
  
  c1 <- dens_cop_gauss(ux[1], uy[1], rho_x)
  c2 <- dens_cop_gauss(ux[2], uy[2], rho_x)
  
  val <- c1 * c2 / (U - L)^2
  
  if (!is.finite(val) || val < 0) return(0)
  val
}

r_cop_cond <- function(x, gx, L, U) {
  
  rho_x <- rho_adaptativo(gx)
  ux <- a_unif(x, L, U)
  
  v1 <- rcond_gauss(ux[1], rho_x)
  v2 <- rcond_gauss(ux[2], rho_x)
  
  c(
    L + (U - L) * v1,
    L + (U - L) * v2
  )
}

q_unif_2d <- function(y, L, U) {
  
  if (
    any(!is.finite(y)) ||
    any(y < L) ||
    any(y > U)
  ) {
    return(0)
  }
  
  1 / (U - L)^2
}

r_unif_2d <- function(L, U) {
  c(
    runif(1, L, U),
    runif(1, L, U)
  )
}

# q(y|x) = w_cop q_C(y|x) + w_unif Unif([L,U]^2)
q_total <- function(y, x, gx, L, U) {
  
  w_cop * q_cop_cond(
    y = y,
    x = x,
    gx = gx,
    L = L,
    U = U
  ) +
    w_unif * q_unif_2d(
      y = y,
      L = L,
      U = U
    )
}

r_total <- function(x, gx, L, U) {
  
  if (runif(1) < w_cop) {
    r_cop_cond(x = x, gx = gx, L = L, U = U)
  } else {
    r_unif_2d(L = L, U = U)
  }
}

# ============================================================
# 7. KNNmin LOCAL CON PRESERVACION ESPACIAL
#
# Idea:
#   - para cada punto MH x_i, calculamos un radio local a partir de la
#     distancia a sus vecinos dentro de la propia nube MH;
#   - solo refinamos x_i si su zona esta suficientemente representada
#     (hay redundancia local);
#   - entre los K vecinos de TODAS las evaluaciones, buscamos el de menor
#     g=|f-c|, pero exigimos que permanezca dentro del radio local;
#   - si no existe una mejora local admisible, x_i se conserva.
#
# A diferencia de la version aumentada, aqui se reemplaza un punto solo
# cuando es espacialmente redundante. Esto permite mejorar precision sin
# desplazar puntos aislados que son importantes para la cobertura.
# ============================================================

refinar_KNN_min_local <- function(
    X,
    gX,
    Xref,
    gref,
    K,
    k_local = k_local_proteccion,
    min_vecinos = min_vecinos_proteccion,
    lambda = lambda_radio,
    tol_self = 1e-12
) {
  
  X <- as.matrix(X)
  Xref <- as.matrix(Xref)
  gX <- as.numeric(gX)
  Xout <- X
  
  # Geometria local de la nube MH. Pedimos k_local + 1 porque el primer
  # vecino es el propio punto (distancia cero).
  kk <- min(k_local + 1, nrow(X))
  kn_local <- FNN::get.knnx(data = X, query = X, k = kk)
  
  # Radio local: distancia al vecino k_local (sin contar el propio punto).
  if (kk <= 1) return(Xout)
  radio_local <- lambda * kn_local$nn.dist[, kk]
  
  # Vecinos candidatos entre TODAS las evaluaciones reales.
  k_busqueda <- min(K + 10, nrow(Xref))
  kn_ref <- FNN::get.knnx(data = Xref, query = X, k = k_busqueda)
  
  reemplazado <- logical(nrow(X))
  
  for (i in seq_len(nrow(X))) {
    
    # Numero de otros estados MH dentro del radio local.
    dloc <- kn_local$nn.dist[i, ]
    n_vecinos <- sum(dloc > tol_self & dloc <= radio_local[i])
    
    # Si el punto es poco representado, se protege para preservar cobertura.
    if (n_vecinos < min_vecinos) next
    
    ids <- kn_ref$nn.index[i, ]
    dist <- kn_ref$nn.dist[i, ]
    
    # Quitamos el propio punto y exigimos cercania geometrica.
    keep <- dist > tol_self & dist <= radio_local[i]
    ids <- ids[keep]
    dist <- dist[keep]
    
    if (length(ids) == 0) next
    
    # Nos quedamos como maximo con K candidatos admisibles.
    ord <- order(dist)
    ids <- ids[ord][seq_len(min(K, length(ids)))]
    
    id_best <- ids[which.min(gref[ids])]
    
    # Reemplazo SOLO si mejora |f-c|.
    if (is.finite(gref[id_best]) && is.finite(gX[i]) &&
        gref[id_best] < gX[i]) {
      Xout[i, ] <- Xref[id_best, ]
      reemplazado[i] <- TRUE
    }
  }
  
  colnames(Xout) <- c("x1", "x2")
  attr(Xout, "reemplazado") <- reemplazado
  attr(Xout, "prop_reemplazada") <- mean(reemplazado)
  Xout
}

# ============================================================
# 8. METODO PRINCIPAL: COPULA ADAPTATIVA + MH + KNNmin LOCAL
# ============================================================

run_metodo <- function(problema, seed = 123) {
  
  set.seed(seed)
  
  f <- problema$f
  L <- problema$L
  U <- problema$U
  c0 <- problema$c
  
  tiempo_inicio <- proc.time()[["elapsed"]]
  
  # ----------------------------------------------------------
  # 8.1. MUESTRA INICIAL
  # ----------------------------------------------------------
  
  X0 <- cbind(
    runif(N_f, L, U),
    runif(N_f, L, U)
  )
  colnames(X0) <- c("x1", "x2")
  
  f0 <- apply(X0, 1, f)
  g0 <- abs(f0 - c0)
  
  # ----------------------------------------------------------
  # 8.2. INICIALIZACION MH
  # ----------------------------------------------------------
  
  ibest <- which.min(g0)
  x <- as.numeric(X0[ibest, ])
  fx <- f0[ibest]
  gx <- g0[ibest]
  
  # ----------------------------------------------------------
  # 8.3. OBJETOS DE LA CADENA
  # ----------------------------------------------------------
  
  chain <- matrix(
    NA_real_,
    nrow = B_mh + 1,
    ncol = 2
  )
  colnames(chain) <- c("x1", "x2")
  chain[1, ] <- x
  
  accepted <- logical(B_mh)
  
  propuestas <- matrix(
    NA_real_,
    nrow = B_mh,
    ncol = 2
  )
  colnames(propuestas) <- c("x1", "x2")
  
  f_propuestas <- numeric(B_mh)
  g_propuestas <- numeric(B_mh)
  rho_chain <- numeric(B_mh)
  
  # ----------------------------------------------------------
  # 8.4. METROPOLIS-HASTINGS
  # ----------------------------------------------------------
  
  for (t in seq_len(B_mh)) {
    
    rho_chain[t] <- rho_adaptativo(gx)
    
    # 90% copula condicional adaptativa + 10% uniforme
    y <- r_total(
      x = x,
      gx = gx,
      L = L,
      U = U
    )
    
    # UNA evaluacion real de f
    fy <- f(y)
    gy <- abs(fy - c0)
    
    propuestas[t, ] <- y
    f_propuestas[t] <- fy
    g_propuestas[t] <- gy
    
    # Target
    log_ratio_target <- -(gy - gx) / tau_mh
    
    # Hastings: q(y|x) usa rho(x), q(x|y) usa rho(y)
    q_y_given_x <- max(
      q_total(
        y = y,
        x = x,
        gx = gx,
        L = L,
        U = U
      ),
      1e-300
    )
    
    q_x_given_y <- max(
      q_total(
        y = x,
        x = y,
        gx = gy,
        L = L,
        U = U
      ),
      1e-300
    )
    
    log_alpha <- log_ratio_target +
      log(q_x_given_y) -
      log(q_y_given_x)
    
    if (log(runif(1)) < min(0, log_alpha)) {
      x <- y
      fx <- fy
      gx <- gy
      accepted[t] <- TRUE
    }
    
    chain[t + 1, ] <- x
  }
  
  # ----------------------------------------------------------
  # 8.5. BURN-IN
  # ----------------------------------------------------------
  
  burn <- floor(burn_frac * nrow(chain))
  
  X_mh <- chain[
    (burn + 1):nrow(chain),
    ,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  # 8.6. TODAS LAS EVALUACIONES REALES
  # ----------------------------------------------------------
  
  Xref <- rbind(X0, propuestas)
  fref <- c(f0, f_propuestas)
  gref <- c(g0, g_propuestas)
  
  # Recuperamos g para los estados MH sin nuevas evaluaciones.
  idx_mh_en_ref <- FNN::get.knnx(
    data = Xref,
    query = X_mh,
    k = 1
  )$nn.index[, 1]
  
  g_mh <- gref[idx_mh_en_ref]
  
  # ----------------------------------------------------------
  # 8.7. KNNmin LOCAL
  # ----------------------------------------------------------
  
  X5_min <- refinar_KNN_min_local(
    X = X_mh,
    gX = g_mh,
    Xref = Xref,
    gref = gref,
    K = 5
  )
  
  X10_min <- refinar_KNN_min_local(
    X = X_mh,
    gX = g_mh,
    Xref = Xref,
    gref = gref,
    K = 10
  )
  
  X100_min <- refinar_KNN_min_local(
    X = X_mh,
    gX = g_mh,
    Xref = Xref,
    gref = gref,
    K = 100
  )
  
  tiempo_total <- proc.time()[["elapsed"]] - tiempo_inicio
  
  # ----------------------------------------------------------
  # 8.8. CURVA VERDADERA SOLO PARA EVALUACION
  # ----------------------------------------------------------
  
  curva <- extraer_curva(
    f = f,
    L = L,
    U = U,
    c0 = c0,
    ngrid = ngrid
  )
  
  Larc <- remuestrear_arco(
    df = curva$df,
    n = n_arco
  )
  
  # ----------------------------------------------------------
  # 8.9. METRICAS
  # ----------------------------------------------------------
  
  met_mh <- metricas_geom(
    Xhat = X_mh,
    Lref = Larc,
    L = L,
    U = U,
    delta = delta_geom
  )
  
  met_5_min <- metricas_geom(
    Xhat = X5_min,
    Lref = Larc,
    L = L,
    U = U,
    delta = delta_geom
  )
  
  met_10_min <- metricas_geom(
    Xhat = X10_min,
    Lref = Larc,
    L = L,
    U = U,
    delta = delta_geom
  )
  
  met_100_min <- metricas_geom(
    Xhat = X100_min,
    Lref = Larc,
    L = L,
    U = U,
    delta = delta_geom
  )
  
  metricas <- rbind(
    "Copula adaptativa + MH" = met_mh,
    "Copula + MH + KNNmin(5)" = met_5_min,
    "Copula + MH + KNNmin(10)" = met_10_min,
    "Copula + MH + KNNmin(100)" = met_100_min
  )
  
  list(
    X_mh = X_mh,
    X5_min = X5_min,
    X10_min = X10_min,
    X100_min = X100_min,
    
    prop_reemplazada_K5 = attr(X5_min, "prop_reemplazada"),
    prop_reemplazada_K10 = attr(X10_min, "prop_reemplazada"),
    prop_reemplazada_K100 = attr(X100_min, "prop_reemplazada"),
    
    Xref = Xref,
    curva = curva,
    Larc = Larc,
    metricas = metricas,
    aceptacion = mean(accepted),
    N_f_real = length(fref),
    tiempo = tiempo_total,
    rho_medio = mean(rho_chain),
    rho_mediana = median(rho_chain),
    L = L,
    U = U
  )
}

# ============================================================
# 9. BASELINE: MUESTREO UNIFORME PURO
#
# Generamos B_total puntos iid uniformes y construimos la estimacion con
# aquellos que caen en la banda |f(x)-c| < delta_uniforme.
# Para que el baseline tenga deliberadamente alta cobertura y baja precision,
# usamos una banda mas ancha que delta_geom. IMPORTANTE: este parametro debe
# fijarse antes de mirar las metricas finales; no debe ajustarse por funcion.
# ============================================================

run_uniforme_puro <- function(
    problema,
    seed = 123,
    delta_uniforme = 3 * delta_geom
) {
  
  set.seed(seed)
  f <- problema$f
  L <- problema$L
  U <- problema$U
  c0 <- problema$c
  tiempo_inicio <- proc.time()[["elapsed"]]
  
  Xall <- cbind(runif(B_total, L, U), runif(B_total, L, U))
  colnames(Xall) <- c("x1", "x2")
  fall <- apply(Xall, 1, f)
  gall <- abs(fall - c0)
  
  # La banda se expresa en escala de f. Si quedara vacia, conservamos los
  # puntos con menor |f-c| para que las metricas sigan siendo definibles.
  keep <- gall < delta_uniforme
  if (!any(keep)) {
    nkeep <- min(200, nrow(Xall))
    keep_ids <- order(gall)[seq_len(nkeep)]
    Xhat <- Xall[keep_ids, , drop = FALSE]
  } else {
    Xhat <- Xall[keep, , drop = FALSE]
  }
  
  list(
    Xhat = Xhat,
    Xall = Xall,
    delta_uniforme = delta_uniforme,
    N_f_real = B_total,
    tiempo = proc.time()[["elapsed"]] - tiempo_inicio
  )
}

# ============================================================
# 10. BASELINE MCMC: RANDOM WALK ISOTROPICO + MH
#
# Propuesta:
#   y = x + sigma_RW Z,  Z ~ N_2(0,I_2).
# Si y cae fuera del hipercubo, se rechaza automaticamente.
# La propuesta gaussiana original es simetrica, por lo que para puntos
# interiores el cociente de propuestas se cancela.
# ============================================================

run_rw_mh <- function(
    problema,
    seed = 123,
    sigma_rw_rel = sigma_rw_rel
) {
  
  set.seed(seed)
  
  f <- problema$f
  L <- problema$L
  U <- problema$U
  c0 <- problema$c
  
  sigma_rw <- sigma_rw_rel * (U - L)
  
  tiempo_inicio <- proc.time()[["elapsed"]]
  
  X0 <- cbind(
    runif(N_f, L, U),
    runif(N_f, L, U)
  )
  colnames(X0) <- c("x1", "x2")
  
  f0 <- apply(X0, 1, f)
  g0 <- abs(f0 - c0)
  
  ibest <- which.min(g0)
  x <- as.numeric(X0[ibest, ])
  gx <- g0[ibest]
  
  chain <- matrix(NA_real_, nrow = B_mh + 1, ncol = 2)
  colnames(chain) <- c("x1", "x2")
  chain[1, ] <- x
  
  accepted <- logical(B_mh)
  n_eval_propuestas <- 0L
  
  for (t in seq_len(B_mh)) {
    
    y <- x + rnorm(2, mean = 0, sd = sigma_rw)
    
    # Si cae fuera, rechazo automatico SIN evaluar f.
    if (any(y < L) || any(y > U)) {
      chain[t + 1, ] <- x
      next
    }
    
    fy <- f(y)
    gy <- abs(fy - c0)
    n_eval_propuestas <- n_eval_propuestas + 1L
    
    log_alpha <- -(gy - gx) / tau_mh
    
    if (log(runif(1)) < min(0, log_alpha)) {
      x <- y
      gx <- gy
      accepted[t] <- TRUE
    }
    
    chain[t + 1, ] <- x
  }
  
  burn <- floor(burn_frac * nrow(chain))
  
  Xhat <- chain[
    (burn + 1):nrow(chain),
    ,
    drop = FALSE
  ]
  
  list(
    Xhat = Xhat,
    aceptacion = mean(accepted),
    N_f_real = N_f + n_eval_propuestas,
    N_propuestas = B_mh,
    tiempo = proc.time()[["elapsed"]] - tiempo_inicio,
    sigma_rw = sigma_rw
  )
}

# ============================================================
# 11. GP RANDOM STRADDLE
# ============================================================

ajustar_gp <- function(X, y) {
  
  DiceKriging::km(
    design = data.frame(
      x1 = X[, 1],
      x2 = X[, 2]
    ),
    response = y,
    covtype = "gauss",
    nugget.estim = TRUE,
    control = list(trace = FALSE)
  )
}

predecir_gp <- function(gp, X) {
  
  pr <- predict(
    gp,
    newdata = data.frame(
      x1 = X[, 1],
      x2 = X[, 2]
    ),
    type = "UK",
    checkNames = FALSE
  )
  
  list(
    mean = as.numeric(pr$mean),
    sd = pmax(as.numeric(pr$sd), 1e-8)
  )
}

run_gp_straddle <- function(
    problema,
    seed = 123,
    N = N_gp
) {
  
  set.seed(seed)
  
  f <- problema$f
  L <- problema$L
  U <- problema$U
  c0 <- problema$c
  
  tiempo_inicio <- proc.time()[["elapsed"]]
  
  # Diseno inicial
  n0 <- min(N_init_gp, N)
  
  X <- cbind(
    runif(n0, L, U),
    runif(n0, L, U)
  )
  colnames(X) <- c("x1", "x2")
  
  y <- apply(X, 1, f)
  
  # Straddle secuencial
  if (N > n0) {
    
    for (t in seq_len(N - n0)) {
      
      gp <- ajustar_gp(X = X, y = y)
      
      Xcand <- cbind(
        runif(N_cand_gp, L, U),
        runif(N_cand_gp, L, U)
      )
      colnames(Xcand) <- c("x1", "x2")
      
      pred <- predecir_gp(gp = gp, X = Xcand)
      
      score <- 1.96 * pred$sd - abs(pred$mean - c0)
      ibest <- which.max(score)
      
      xnew <- Xcand[ibest, , drop = FALSE]
      ynew <- f(as.numeric(xnew[1, ]))
      
      X <- rbind(X, xnew)
      y <- c(y, ynew)
    }
  }
  
  gp <- ajustar_gp(X = X, y = y)
  
  # Puntos baratos para aproximar el level set con la media predictiva.
  Xpred <- cbind(
    runif(N_pred_gp, L, U),
    runif(N_pred_gp, L, U)
  )
  colnames(Xpred) <- c("x1", "x2")
  
  pred <- predecir_gp(gp = gp, X = Xpred)
  error_gp <- abs(pred$mean - c0)
  
  n_hat <- min(N, nrow(Xpred))
  ids_gp <- order(error_gp)[seq_len(n_hat)]
  
  X_gp <- Xpred[ids_gp, , drop = FALSE]
  
  list(
    X_gp = X_gp,
    X_real = X,
    N_f_real = nrow(X),
    tiempo = proc.time()[["elapsed"]] - tiempo_inicio
  )
}

# ============================================================
# 12. GRAFICOS: SOLO CUATRO METODOS
#
# Rojo: metodos MCMC.
# Azul: GP.
# Negro: curva de referencia.
# ============================================================

graficar_metodos <- function(
    rr,
    rr_rw,
    rgp,
    nombre,
    K_grafico = 10
) {
  
  par(
    mfrow = c(1, 4),
    mar = c(4, 4, 3, 1)
  )
  
  X_knn <- switch(
    as.character(K_grafico),
    "5" = rr$X5_min,
    "10" = rr$X10_min,
    "100" = rr$X100_min,
    stop("K_grafico debe ser 5, 10 o 100.")
  )
  
  conjuntos <- list(
    rr_rw$Xhat,
    rr$X_mh,
    X_knn,
    rgp$X_gp
  )
  
  titulos <- c(
    "RW isotropico + MH",
    "Copula adaptativa + MH",
    paste0("Copula + MH + KNNmin(", K_grafico, ")"),
    "GP Random Straddle"
  )
  
  colores <- c("red", "red", "red", "blue")
  tamanios <- c(0.45, 0.45, 0.45, 0.80)
  transparencias <- c(0.60, 0.60, 0.60, 1.00)
  
  for (k in seq_along(conjuntos)) {
    
    plot(
      rr$curva$P,
      type = "p",
      pch = 16,
      cex = 0.10,
      col = "black",
      asp = 1,
      xlim = c(rr$L, rr$U),
      ylim = c(rr$L, rr$U),
      xlab = expression(x[1]),
      ylab = expression(x[2]),
      main = paste0(nombre, ": ", titulos[k])
    )
    
    points(
      conjuntos[[k]],
      pch = 16,
      cex = tamanios[k],
      col = adjustcolor(
        colores[k],
        alpha.f = transparencias[k]
      )
    )
  }
  
  par(mfrow = c(1, 1))
}

# ============================================================
# 13. EXPERIMENTO COMPLETO
# ============================================================

tabla_final <- data.frame()
resultados <- list()

for (j in seq_along(funciones)) {
  
  problema <- funciones[[j]]
  seed_j <- 123 + 100 * j
  
  cat(
    "\n\n========================================\n",
    "FUNCION: ", problema$nombre, "\n",
    "========================================\n",
    sep = ""
  )
  
  # ----------------------------------------------------------
  # 13.1. Copula adaptativa + MH + KNNmin local
  # ----------------------------------------------------------
  
  rr <- run_metodo(
    problema = problema,
    seed = seed_j
  )
  
  # ----------------------------------------------------------
  # 13.2. Uniforme puro
  # ----------------------------------------------------------
  
  runif_base <- run_uniforme_puro(
    problema = problema,
    seed = seed_j + 1,
    delta_uniforme = 3 * delta_geom
  )
  
  # ----------------------------------------------------------
  # 13.3. RW isotropico + MH
  # ----------------------------------------------------------
  
  rr_rw <- run_rw_mh(
    problema = problema,
    seed = seed_j + 2,
    sigma_rw_rel = sigma_rw_rel
  )
  
  # ----------------------------------------------------------
  # 13.4. GP Random Straddle
  # ----------------------------------------------------------
  
  rgp <- run_gp_straddle(
    problema = problema,
    seed = seed_j + 3,
    N = N_gp
  )
  
  # ----------------------------------------------------------
  # 13.5. Metricas de los baselines
  # ----------------------------------------------------------
  
  met_unif <- metricas_geom(
    Xhat = runif_base$Xhat,
    Lref = rr$Larc,
    L = problema$L,
    U = problema$U,
    delta = delta_geom
  )
  
  met_rw <- metricas_geom(
    Xhat = rr_rw$Xhat,
    Lref = rr$Larc,
    L = problema$L,
    U = problema$U,
    delta = delta_geom
  )
  
  met_gp <- metricas_geom(
    Xhat = rgp$X_gp,
    Lref = rr$Larc,
    L = problema$L,
    U = problema$U,
    delta = delta_geom
  )
  
  # ----------------------------------------------------------
  # 13.6. Tabla: Uniforme puro
  # ----------------------------------------------------------
  
  tabla_final <- rbind(
    tabla_final,
    data.frame(
      funcion = problema$nombre,
      metodo = "Uniforme",
      precision = unname(met_unif["precision"]),
      cobertura = unname(met_unif["cobertura"]),
      F1 = unname(met_unif["F1"]),
      aceptacion_MH = NA_real_,
      N_f = runif_base$N_f_real,
      tiempo_seg = runif_base$tiempo,
      stringsAsFactors = FALSE
    )
  )
  
  # ----------------------------------------------------------
  # 13.7. Tabla: RW isotropico + MH
  # ----------------------------------------------------------
  
  tabla_final <- rbind(
    tabla_final,
    data.frame(
      funcion = problema$nombre,
      metodo = "RW isotropico + MH",
      precision = unname(met_rw["precision"]),
      cobertura = unname(met_rw["cobertura"]),
      F1 = unname(met_rw["F1"]),
      aceptacion_MH = rr_rw$aceptacion,
      N_f = rr_rw$N_f_real,
      tiempo_seg = rr_rw$tiempo,
      stringsAsFactors = FALSE
    )
  )
  
  # ----------------------------------------------------------
  # 13.8. Tabla: Copula + MH y KNNmin(K)
  # ----------------------------------------------------------
  
  for (k in seq_len(nrow(rr$metricas))) {
    
    nombre_metodo <- rownames(rr$metricas)[k]
    
    tabla_final <- rbind(
      tabla_final,
      data.frame(
        funcion = problema$nombre,
        metodo = nombre_metodo,
        precision = rr$metricas[k, "precision"],
        cobertura = rr$metricas[k, "cobertura"],
        F1 = rr$metricas[k, "F1"],
        aceptacion_MH = rr$aceptacion,
        N_f = rr$N_f_real,
        tiempo_seg = rr$tiempo,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # ----------------------------------------------------------
  # 13.9. Tabla: GP
  # ----------------------------------------------------------
  
  tabla_final <- rbind(
    tabla_final,
    data.frame(
      funcion = problema$nombre,
      metodo = "GP Random Straddle",
      precision = unname(met_gp["precision"]),
      cobertura = unname(met_gp["cobertura"]),
      F1 = unname(met_gp["F1"]),
      aceptacion_MH = NA_real_,
      N_f = rgp$N_f_real,
      tiempo_seg = rgp$tiempo,
      stringsAsFactors = FALSE
    )
  )
  
  # ----------------------------------------------------------
  # 13.10. Consola
  # ----------------------------------------------------------
  
  cat("\n--- METRICAS ---\n")
  
  tabla_funcion <- tabla_final[
    tabla_final$funcion == problema$nombre,
    ,
    drop = FALSE
  ]
  
  print(
    transform(
      tabla_funcion,
      precision = round(precision, 4),
      cobertura = round(cobertura, 4),
      F1 = round(F1, 4),
      aceptacion_MH = round(aceptacion_MH, 4),
      tiempo_seg = round(tiempo_seg, 2)
    ),
    row.names = FALSE
  )
  
  cat(
    "\nProporcion de estados reemplazados por KNN local:",
    "\n  K=5   :", sprintf("%.4f", rr$prop_reemplazada_K5),
    "\n  K=10  :", sprintf("%.4f", rr$prop_reemplazada_K10),
    "\n  K=100 :", sprintf("%.4f", rr$prop_reemplazada_K100),
    "\n"
  )
  
  cat(
    "\nAceptacion Copula + MH =", round(rr$aceptacion, 4),
    "\nAceptacion RW + MH =", round(rr_rw$aceptacion, 4),
    "\n"
  )
  
  # ----------------------------------------------------------
  # 13.11. Graficos: SOLO 4 PANELES
  # ----------------------------------------------------------
  
  graficar_metodos(
    rr = rr,
    rr_rw = rr_rw,
    rgp = rgp,
    nombre = problema$nombre,
    K_grafico = K_grafico
  )
  
  # Guardamos resultados por si se quieren inspeccionar despues.
  resultados[[problema$nombre]] <- list(
    copula = rr,
    uniforme = runif_base,
    rw = rr_rw,
    gp = rgp,
    met_uniforme = met_unif,
    met_rw = met_rw,
    met_gp = met_gp
  )
}

# ============================================================
# 14. TABLA FINAL
# ============================================================

tabla_slides <- tabla_final

tabla_slides$precision <- round(tabla_slides$precision, 4)
tabla_slides$cobertura <- round(tabla_slides$cobertura, 4)
tabla_slides$F1 <- round(tabla_slides$F1, 4)
tabla_slides$aceptacion_MH <- round(tabla_slides$aceptacion_MH, 4)
tabla_slides$tiempo_seg <- round(tabla_slides$tiempo_seg, 2)

cat(
  "\n\n========================================\n",
  "TABLA FINAL\n",
  "========================================\n\n",
  sep = ""
)

print(tabla_slides, row.names = FALSE)

# ============================================================
# 15. TABLAS POR FUNCION
# ============================================================

for (nombre_funcion in unique(tabla_slides$funcion)) {
  
  cat(
    "\n\n==================== ",
    nombre_funcion,
    " ====================\n",
    sep = ""
  )
  
  print(
    tabla_slides[
      tabla_slides$funcion == nombre_funcion,
      ,
      drop = FALSE
    ],
    row.names = FALSE
  )
}

# ============================================================
# 16. RESUMEN DEL EXPERIMENTO
# ============================================================

cat(
  "\n\n========================================\n",
  "PARAMETROS PRINCIPALES\n",
  "========================================\n",
  "Evaluaciones nominales MCMC: ", B_total, "\n",
  "Evaluaciones reales GP: ", N_gp, "\n",
  "Puntos iniciales MCMC: ", N_f, "\n",
  "Propuestas MH: ", B_mh, "\n",
  "sigma_g: ", sigma_g, "\n",
  "rho maximo: ", rho_max, "\n",
  "Peso copula: ", w_cop, "\n",
  "Peso uniforme en propuesta copular: ", w_unif, "\n",
  "tau MH: ", tau_mh, "\n",
  "Burn-in: ", 100 * burn_frac, "%\n",
  "sigma RW relativo al ancho: ", sigma_rw_rel, "\n",
  "Banda uniforme |f-c| < delta: ", 3 * delta_geom, "\n",
  "K usados: ", paste(K_valores, collapse = ", "), "\n",
  "K mostrado en figuras: ", K_grafico, "\n",
  "delta geometrico: ", delta_geom, "\n",
  "========================================\n",
  "FIN DEL EXPERIMENTO\n",
  "========================================\n",
  sep = ""
)
