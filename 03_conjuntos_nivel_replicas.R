# ============================================================
# EXPERIMENTO 2D - VERSION LIMPIA PARA TESIS
#
# Metodos en la tabla:
#   1. Uniforme puro (iid + banda)
#   2. RW isotropico + MH
#   3. Copula adaptativa + MH
#   4. Copula + MH + KNN sin reutilizacion
#   5. GP Random Straddle
#
# Figuras (4 paneles por funcion):
#   1. RW isotropico + MH
#   2. Copula adaptativa + MH
#   3. Copula + MH + KNN sin reutilizacion
#   4. GP Random Straddle
#
## KNN-min sin reutilizacion:
#   - procesa los estados MH desde mayor a menor g(x)=|f(x)-c|;
#   - para cada estado busca sus K vecinos entre TODAS las evaluaciones
#     ya realizadas (incluidas propuestas MH aceptadas y rechazadas);
#   - elige, entre los candidatos aun no utilizados, el de menor g;
#   - reemplaza el estado MH solo si el candidato mejora g;
#   - cada punto evaluado puede utilizarse como reemplazo a lo sumo una vez.
#
# Es un postprocesamiento heuristico: no modifica la cadena MH ni agrega
# evaluaciones de f; aprovecha el historial de evaluaciones ya pagadas.
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

# KNN-min sin reutilizacion de candidatos
# Cada punto evaluado puede usarse como reemplazo a lo sumo una vez.
# Los estados MH se procesan desde mayor a menor g=|f-c|.
K_knn <- 100

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

# Repeticiones Monte Carlo del experimento completo
# Para una prueba rapida usar nn <- 2; para la corrida final usar nn <- 30.
nn <- 30
seed_base <- 20260929

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
# 7. KNN-MIN SIN REUTILIZACION DE CANDIDATOS
#
# Para cada estado X_i de la cadena MH:
#   1. los estados se procesan de mayor a menor g(X_i)=|f(X_i)-c|;
#   2. se buscan sus K vecinos mas proximos entre TODAS las
#      evaluaciones reales ya realizadas;
#   3. entre esos vecinos se descartan los candidatos que ya fueron
#      utilizados previamente como reemplazo;
#   4. se elige el candidato disponible con menor g;
#   5. X_i se reemplaza solamente si ese candidato mejora g(X_i).
#
# Cada fila de Xref puede utilizarse como reemplazo a lo sumo una vez.
# No usa la curva verdadera y no agrega evaluaciones de f.
# ============================================================

refinar_KNN_sin_reutilizacion <- function(
    X,
    gX,
    Xref,
    gref,
    K = 10,
    tol_self = 1e-12
) {
  X <- as.matrix(X)
  Xref <- as.matrix(Xref)
  gX <- as.numeric(gX)
  gref <- as.numeric(gref)

  n <- nrow(X)
  nr <- nrow(Xref)
  Xout <- X

  reemplazado <- logical(n)
  usado_ref <- logical(nr)
  id_reemplazo <- rep(NA_integer_, n)

  if (n == 0L || nr == 0L) {
    attr(Xout, "prop_reemplazada") <- 0
    attr(Xout, "n_candidatos_usados") <- 0L
    return(Xout)
  }

  # Pedimos K+1 porque normalmente X_i aparece exactamente en Xref
  # y luego se elimina mediante tol_self.
  kk <- min(K + 1L, nr)
  kn <- FNN::get.knnx(data = Xref, query = X, k = kk)

  # Primero intentamos corregir los estados con mayor discrepancia.
  orden <- order(gX, decreasing = TRUE, na.last = NA)

  for (i in orden) {

    ids <- kn$nn.index[i, ]
    ds  <- kn$nn.dist[i, ]

    # Quitar coincidencias exactas con el estado actual.
    keep <- (ds > tol_self)
    ids <- ids[keep]

    if (length(ids) == 0L) next

    # Mantener como maximo K vecinos distintos.
    ids <- head(ids, K)

    # No reutilizar candidatos ya asignados como reemplazo.
    ids_disp <- ids[!usado_ref[ids]]

    if (length(ids_disp) == 0L) next

    # Mejor candidato funcional entre los disponibles.
    id_best <- ids_disp[which.min(gref[ids_disp])]

    if (is.finite(gref[id_best]) &&
        is.finite(gX[i]) &&
        gref[id_best] < gX[i]) {

      Xout[i, ] <- Xref[id_best, ]
      reemplazado[i] <- TRUE
      usado_ref[id_best] <- TRUE
      id_reemplazo[i] <- id_best
    }
  }

  colnames(Xout) <- c("x1", "x2")

  attr(Xout, "prop_reemplazada") <- mean(reemplazado)
  attr(Xout, "reemplazado") <- reemplazado
  attr(Xout, "id_reemplazo") <- id_reemplazo
  attr(Xout, "n_candidatos_usados") <- sum(usado_ref)

  Xout
}

# ============================================================
# 8. METODO PRINCIPAL: COPULA ADAPTATIVA + MH + KNN SIN REUTILIZACION
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
  # 8.7. KNN-MIN SIN REUTILIZACION DE CANDIDATOS
  # ----------------------------------------------------------

  X_knn_unico <- refinar_KNN_sin_reutilizacion(
    X = X_mh,
    gX = g_mh,
    Xref = Xref,
    gref = gref,
    K = K_knn
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
    Xhat = X_mh, Lref = Larc, L = L, U = U, delta = delta_geom
  )

  met_knn_unico <- metricas_geom(
    Xhat = X_knn_unico, Lref = Larc, L = L, U = U, delta = delta_geom
  )

  metricas <- rbind(
    "Copula adaptativa + MH" = met_mh,
    "Copula + MH + KNN sin reutilizacion" = met_knn_unico
  )

  list(
    X_mh = X_mh,
    X_knn_unico = X_knn_unico,
    prop_reemplazada_knn = attr(X_knn_unico, "prop_reemplazada"),
    n_candidatos_usados_knn = attr(X_knn_unico, "n_candidatos_usados"),

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
# Se usa una banda fija, comun a todas las funciones. Este parametro debe
# fijarse antes de comparar los resultados y no ajustarse por funcion.
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
    nombre
) {

  par(
    mfrow = c(2,2),
    mar = c(4, 4, 3, 1)
  )

  conjuntos <- list(
    rr_rw$Xhat,
    rr$X_mh,
    rr$X_knn_unico,
    rgp$X_gp
  )

  titulos <- c(
    "RW isotropico + MH",
    "Copula adaptativa + MH",
    paste0("Copula + MH + KNN unico (K=", K_knn, ")"),
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
# 13. EXPERIMENTO REPETIDO nn VECES
#
# Se guardan los resultados de cada replica y luego se reportan
# media y desvio estandar de precision, cobertura, F1, aceptacion,
# numero de evaluaciones y tiempo.
#
# IMPORTANTE:
#   - para verificar que todo corre: nn <- 2
#   - para la corrida final nocturna: nn <- 30
# ============================================================

tabla_replicas <- data.frame()
resultados_primera_rep <- list()

for (b in seq_len(nn)) {

  cat(
    "\n\n############################################################\n",
    "REPLICA ", b, " / ", nn, "\n",
    "############################################################\n",
    sep = ""
  )

  for (j in seq_along(funciones)) {

    problema <- funciones[[j]]

    # Semillas reproducibles y distintas por replica, funcion y metodo.
    seed_j <- seed_base + 10000L * b + 100L * j

    cat(
      "\n--- ", problema$nombre,
      " | replica ", b, " / ", nn, " ---\n",
      sep = ""
    )

    # ----------------------------------------------------------
    # 13.1. Copula adaptativa + MH + KNN sin reutilizacion
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
      seed = seed_j + 1L,
      delta_uniforme = 3 * delta_geom
    )

    # ----------------------------------------------------------
    # 13.3. RW isotropico + MH
    # ----------------------------------------------------------
    rr_rw <- run_rw_mh(
      problema = problema,
      seed = seed_j + 2L,
      sigma_rw_rel = sigma_rw_rel
    )

    # ----------------------------------------------------------
    # 13.4. GP Random Straddle
    # ----------------------------------------------------------
    rgp <- run_gp_straddle(
      problema = problema,
      seed = seed_j + 3L,
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
    # 13.6. Guardar una fila por metodo y replica
    # ----------------------------------------------------------
    agregar_fila <- function(metodo, met, aceptacion, N_f_real, tiempo) {
      data.frame(
        replica = b,
        funcion = problema$nombre,
        metodo = metodo,
        precision = unname(met["precision"]),
        cobertura = unname(met["cobertura"]),
        F1 = unname(met["F1"]),
        aceptacion_MH = aceptacion,
        N_f = N_f_real,
        tiempo_seg = tiempo,
        stringsAsFactors = FALSE
      )
    }

    tabla_replicas <- rbind(
      tabla_replicas,
      agregar_fila(
        "Uniforme", met_unif, NA_real_,
        runif_base$N_f_real, runif_base$tiempo
      ),
      agregar_fila(
        "RW isotropico + MH", met_rw, rr_rw$aceptacion,
        rr_rw$N_f_real, rr_rw$tiempo
      ),
      agregar_fila(
        "Copula adaptativa + MH",
        rr$metricas["Copula adaptativa + MH", ],
        rr$aceptacion, rr$N_f_real, rr$tiempo
      ),
      agregar_fila(
        "Copula + MH + KNN sin reutilizacion",
        rr$metricas["Copula + MH + KNN sin reutilizacion", ],
        rr$aceptacion, rr$N_f_real, rr$tiempo
      ),
      agregar_fila(
        "GP Random Straddle", met_gp, NA_real_,
        rgp$N_f_real, rgp$tiempo
      )
    )

    # ----------------------------------------------------------
    # 13.7. Diagnostico breve en consola
    # ----------------------------------------------------------
    cat(
      "   F1 Copula+MH = ",
      sprintf("%.4f", rr$metricas["Copula adaptativa + MH", "F1"]),
      " | F1 KNN = ",
      sprintf("%.4f", rr$metricas["Copula + MH + KNN sin reutilizacion", "F1"]),
      " | prop. reemplazada = ",
      sprintf("%.4f", rr$prop_reemplazada_knn),
      "\n",
      sep = ""
    )

    # Guardamos y graficamos solo la primera replica para no abrir
    # 30 veces las mismas ventanas graficas durante la corrida nocturna.
    if (b == 1L) {
      resultados_primera_rep[[problema$nombre]] <- list(
        copula = rr,
        uniforme = runif_base,
        rw = rr_rw,
        gp = rgp,
        met_uniforme = met_unif,
        met_rw = met_rw,
        met_gp = met_gp
      )

      graficar_metodos(
        rr = rr,
        rr_rw = rr_rw,
        rgp = rgp,
        nombre = problema$nombre
      )
    }
  }

  # Guardado incremental: si la corrida se interrumpe, no se pierde
  # lo ya calculado.
  write.csv(
    tabla_replicas,
    file = "resultados_replicas_parcial.csv",
    row.names = FALSE
  )
}

# ============================================================
# 14. RESUMEN: MEDIA Y DESVIO ESTANDAR
# ============================================================

media_na <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

desvio_na <- function(x) {
  z <- x[is.finite(x)]
  if (length(z) <= 1L) return(NA_real_)
  sd(z)
}

claves <- unique(tabla_replicas[, c("funcion", "metodo")])

tabla_resumen <- do.call(
  rbind,
  lapply(seq_len(nrow(claves)), function(i) {

    fun_i <- claves$funcion[i]
    met_i <- claves$metodo[i]

    dd <- tabla_replicas[
      tabla_replicas$funcion == fun_i &
        tabla_replicas$metodo == met_i,
      ,
      drop = FALSE
    ]

    data.frame(
      funcion = fun_i,
      metodo = met_i,
      nn = nrow(dd),

      precision_media = media_na(dd$precision),
      precision_sd = desvio_na(dd$precision),

      cobertura_media = media_na(dd$cobertura),
      cobertura_sd = desvio_na(dd$cobertura),

      F1_media = media_na(dd$F1),
      F1_sd = desvio_na(dd$F1),

      aceptacion_media = media_na(dd$aceptacion_MH),
      aceptacion_sd = desvio_na(dd$aceptacion_MH),

      N_f_media = media_na(dd$N_f),
      N_f_sd = desvio_na(dd$N_f),

      tiempo_media = media_na(dd$tiempo_seg),
      tiempo_sd = desvio_na(dd$tiempo_seg),

      stringsAsFactors = FALSE
    )
  })
)

# Orden estable de metodos en las tablas.
orden_metodos <- c(
  "Uniforme",
  "RW isotropico + MH",
  "Copula adaptativa + MH",
  "Copula + MH + KNN sin reutilizacion",
  "GP Random Straddle"
)

tabla_resumen$metodo <- factor(
  tabla_resumen$metodo,
  levels = orden_metodos
)

tabla_resumen <- tabla_resumen[
  order(tabla_resumen$funcion, tabla_resumen$metodo),
  ,
  drop = FALSE
]

tabla_resumen$metodo <- as.character(tabla_resumen$metodo)

# Version redondeada para consola / tesis.
tabla_resumen_print <- tabla_resumen
cols_num <- setdiff(
  names(tabla_resumen_print),
  c("funcion", "metodo", "nn")
)
tabla_resumen_print[cols_num] <- lapply(
  tabla_resumen_print[cols_num],
  function(x) round(x, 4)
)

cat(
  "\n\n========================================\n",
  "RESUMEN FINAL: MEDIA Y DESVIO (nn = ", nn, ")\n",
  "========================================\n\n",
  sep = ""
)

print(tabla_resumen_print, row.names = FALSE)

# Tabla compacta con formato media (sd) para precision, cobertura y F1.
fmt_media_sd <- function(mu, sig, dig = 4) {
  if (!is.finite(mu)) return("--")
  if (!is.finite(sig)) return(sprintf(paste0("%.", dig, "f"), mu))
  sprintf(
    paste0("%.", dig, "f (%. ", dig, "f)"),
    mu, sig
  )
}

# La funcion anterior se redefine sin espacios en el formato para evitar
# diferencias entre versiones de R.
fmt_media_sd <- function(mu, sig, dig = 4) {
  if (!is.finite(mu)) return("--")
  if (!is.finite(sig)) return(formatC(mu, format = "f", digits = dig))
  paste0(
    formatC(mu, format = "f", digits = dig),
    " (",
    formatC(sig, format = "f", digits = dig),
    ")"
  )
}

tabla_compacta <- data.frame(
  funcion = tabla_resumen$funcion,
  metodo = tabla_resumen$metodo,
  precision = mapply(
    fmt_media_sd,
    tabla_resumen$precision_media,
    tabla_resumen$precision_sd
  ),
  cobertura = mapply(
    fmt_media_sd,
    tabla_resumen$cobertura_media,
    tabla_resumen$cobertura_sd
  ),
  F1 = mapply(
    fmt_media_sd,
    tabla_resumen$F1_media,
    tabla_resumen$F1_sd
  ),
  aceptacion_MH = mapply(
    fmt_media_sd,
    tabla_resumen$aceptacion_media,
    tabla_resumen$aceptacion_sd
  ),
  stringsAsFactors = FALSE
)

cat(
  "\n\n========================================\n",
  "TABLA COMPACTA: media (desvio)\n",
  "========================================\n\n",
  sep = ""
)

print(tabla_compacta, row.names = FALSE)

# ============================================================
# 15. GUARDAR RESULTADOS
# ============================================================

write.csv(
  tabla_replicas,
  file = "resultados_30_replicas_crudos.csv",
  row.names = FALSE
)

write.csv(
  tabla_resumen,
  file = "resultados_30_replicas_resumen.csv",
  row.names = FALSE
)

write.csv(
  tabla_compacta,
  file = "resultados_30_replicas_media_sd.csv",
  row.names = FALSE
)

saveRDS(
  list(
    parametros = list(
      nn = nn,
      seed_base = seed_base,
      N_f = N_f,
      B_total = B_total,
      B_mh = B_mh,
      sigma_g = sigma_g,
      eta_rho = eta_rho,
      w_cop = w_cop,
      w_unif = w_unif,
      tau_mh = tau_mh,
      burn_frac = burn_frac,
      sigma_rw_rel = sigma_rw_rel,
      K_knn = K_knn,
      delta_geom = delta_geom,
      N_gp = N_gp
    ),
    replicas = tabla_replicas,
    resumen = tabla_resumen,
    primera_rep = resultados_primera_rep
  ),
  file = "experimento_2D_30_replicas.rds"
)

cat(
  "\n========================================\n",
  "ARCHIVOS GUARDADOS\n",
  "========================================\n",
  "resultados_30_replicas_crudos.csv\n",
  "resultados_30_replicas_resumen.csv\n",
  "resultados_30_replicas_media_sd.csv\n",
  "experimento_2D_30_replicas.rds\n",
  "========================================\n",
  sep = ""
)
