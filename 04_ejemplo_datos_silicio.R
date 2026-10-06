rm(list = ls())

# ============================================================
# LINGOTES DE SILICIO — MÉTODO DEFINITIVO
#
# 1. Interpolación bilineal de lifetime
# 2. Cópula gaussiana condicional adaptativa + MH
# 3. MH + KNN(5)
# 4. MH + KNN(10)
# 5. MH + KNN(100)
# 6. GP Random Straddle
# 7. Precisión, cobertura y F1
# 8. Gráficos comparativos
#
# Nivel objetivo: c = 230
# ============================================================

library(FNN)
library(ggplot2)

# GP
if (!requireNamespace("DiceKriging", quietly = TRUE)) {
  stop("Instalá DiceKriging con install.packages('DiceKriging')")
}

set.seed(123)

# ============================================================
# 1. CARGAR DATOS
# ============================================================

datos <- read.table(
  "data4.txt",
  header = FALSE
)

colnames(datos) <- c("x1", "x2", "lifetime")

c_level <- 230

cat("\n============================================\n")
cat("LINGOTES DE SILICIO\n")
cat("============================================\n")
cat("Observaciones:", nrow(datos), "\n")
cat("Rango x1:", range(datos$x1), "\n")
cat("Rango x2:", range(datos$x2), "\n")
cat("Rango lifetime:", range(datos$lifetime), "\n")
cat("Nivel objetivo:", c_level, "\n")


# ============================================================
# 2. CONSTRUIR GRILLA OBSERVADA
# ============================================================

x1_vals <- sort(unique(datos$x1))
x2_vals <- sort(unique(datos$x2))

n1 <- length(x1_vals)
n2 <- length(x2_vals)

Z <- matrix(
  NA_real_,
  nrow = n1,
  ncol = n2
)

for (i in seq_len(nrow(datos))) {
  
  ii <- match(datos$x1[i], x1_vals)
  jj <- match(datos$x2[i], x2_vals)
  
  Z[ii, jj] <- datos$lifetime[i]
}

if (anyNA(Z)) {
  stop("La grilla tiene combinaciones (x1,x2) sin observar.")
}


# ============================================================
# 3. DOMINIO
# ============================================================

L1 <- min(x1_vals)
U1 <- max(x1_vals)

L2 <- min(x2_vals)
U2 <- max(x2_vals)

range1 <- U1 - L1
range2 <- U2 - L2

cat("\nDominio:\n")
cat("x1 en [", L1, ",", U1, "]\n")
cat("x2 en [", L2, ",", U2, "]\n")


# ============================================================
# 4. INTERPOLACIÓN BILINEAL
# ============================================================

f_interp <- function(x) {
  
  x1 <- x[1]
  x2 <- x[2]
  
  if (
    x1 < L1 || x1 > U1 ||
    x2 < L2 || x2 > U2
  ) {
    return(NA_real_)
  }
  
  i <- findInterval(x1, x1_vals)
  j <- findInterval(x2, x2_vals)
  
  if (i <= 0) i <- 1
  if (i >= n1) i <- n1 - 1
  
  if (j <= 0) j <- 1
  if (j >= n2) j <- n2 - 1
  
  xa <- x1_vals[i]
  xb <- x1_vals[i + 1]
  
  ya <- x2_vals[j]
  yb <- x2_vals[j + 1]
  
  z11 <- Z[i,     j]
  z21 <- Z[i + 1, j]
  z12 <- Z[i,     j + 1]
  z22 <- Z[i + 1, j + 1]
  
  tx <- (x1 - xa)/(xb - xa)
  ty <- (x2 - ya)/(yb - ya)
  
  z <-
    (1 - tx)*(1 - ty)*z11 +
    tx*(1 - ty)*z21 +
    (1 - tx)*ty*z12 +
    tx*ty*z22
  
  as.numeric(z)
}


# ============================================================
# 5. EVALUACIÓN VECTORIAL
# ============================================================

eval_f <- function(X) {
  
  X <- as.matrix(X)
  
  apply(
    X,
    1,
    f_interp
  )
}


# ============================================================
# 6. PARÁMETROS GENERALES
# ============================================================

# ---------- Cópula + MH ----------

N_f       <- 500L
B_total   <- 5000L
B_mh      <- B_total - N_f

sigma_g   <- 30
eta       <- 0.02

w_cop     <- 0.90
w_unif    <- 0.10

tau_mh    <- 5

burn_frac <- 0.10

# ---------- KNN ----------

K_values <- c(5, 10, 100)

lambda_knn <- 5
eps_dist   <- 1e-6
tol_self   <- 1e-12

# ---------- GP ----------

N_gp       <- 200L
N_init_gp  <- 30L
N_cand_gp  <- 3000L
N_pred_gp  <- 20000L

# ---------- Métricas ----------

# Las distancias se calculan después de normalizar a [0,1]^2.
delta_metric <- 0.03

# Número aproximado de puntos de referencia
M_ref <- 5000L


# ============================================================
# 7. NORMALIZACIÓN DEL DOMINIO
# ============================================================

clamp01 <- function(u) {
  pmin(pmax(u, 1e-10), 1 - 1e-10)
}

to_unit <- function(X) {
  
  X <- as.matrix(X)
  
  cbind(
    (X[, 1] - L1)/range1,
    (X[, 2] - L2)/range2
  )
}

from_unit <- function(U) {
  
  U <- as.matrix(U)
  
  cbind(
    L1 + range1*U[, 1],
    L2 + range2*U[, 2]
  )
}


# ============================================================
# 8. DEPENDENCIA ADAPTATIVA
#
# rho(x) = (1-eta) exp[-(g(x)/sigma_g)^2]
# ============================================================

rho_fun <- function(fx) {
  
  g <- abs(fx - c_level)
  
  rho <-
    (1 - eta) *
    exp(-(g/sigma_g)^2)
  
  pmin(
    pmax(rho, 0),
    1 - eta
  )
}


# ============================================================
# 9. CÓPULA GAUSSIANA
#
# Densidad c_rho(u,v)
# ============================================================

log_cop_gaussian <- function(u, v, rho) {
  
  u <- clamp01(u)
  v <- clamp01(v)
  
  rho <- pmin(
    pmax(rho, -0.999999),
    0.999999
  )
  
  z1 <- qnorm(u)
  z2 <- qnorm(v)
  
  -0.5*log(1 - rho^2) +
    (
      2*rho*z1*z2 -
        rho^2*(z1^2 + z2^2)
    ) /
    (
      2*(1 - rho^2)
    )
}


# ============================================================
# 10. SIMULACIÓN CONDICIONAL GAUSSIANA
#
# V | U=u
# ============================================================

sim_gaussian_conditional <- function(u, rho) {
  
  u <- clamp01(u)
  
  rho <- pmin(
    pmax(rho, -0.999999),
    0.999999
  )
  
  z1 <- qnorm(u)
  
  z2 <- rnorm(
    1,
    mean = rho*z1,
    sd = sqrt(1 - rho^2)
  )
  
  clamp01(
    pnorm(z2)
  )
}


# ============================================================
# 11. DENSIDAD q_C(y|x)
#
# Cada coordenada propuesta se vincula con
# la correspondiente coordenada actual.
# ============================================================

log_q_copula <- function(y, x, fx) {
  
  ux <- c(
    clamp01((x[1] - L1)/range1),
    clamp01((x[2] - L2)/range2)
  )
  
  uy <- c(
    clamp01((y[1] - L1)/range1),
    clamp01((y[2] - L2)/range2)
  )
  
  rho <- rho_fun(fx)
  
  logc1 <- log_cop_gaussian(
    ux[1],
    uy[1],
    rho
  )
  
  logc2 <- log_cop_gaussian(
    ux[2],
    uy[2],
    rho
  )
  
  # Jacobiano de [0,1]^2 -> dominio original
  log_jac <- -log(range1) - log(range2)
  
  logc1 + logc2 + log_jac
}


# ============================================================
# 12. DENSIDAD DE LA PROPUESTA COMPLETA
#
# q(y|x) =
#   0.90 q_C(y|x) +
#   0.10 Unif(D)
# ============================================================

log_sum_exp2 <- function(a, b) {
  
  m <- max(a, b)
  
  m + log(
    exp(a - m) +
      exp(b - m)
  )
}


log_q_mix <- function(y, x, fx) {
  
  log_qc <- log_q_copula(
    y = y,
    x = x,
    fx = fx
  )
  
  log_qunif <-
    -log(range1) -
    log(range2)
  
  log_sum_exp2(
    log(w_cop) + log_qc,
    log(w_unif) + log_qunif
  )
}


# ============================================================
# 13. GENERAR PROPUESTA
# ============================================================

propose_point <- function(x, fx) {
  
  # Componente uniforme
  if (runif(1) < w_unif) {
    
    return(
      c(
        runif(1, L1, U1),
        runif(1, L2, U2)
      )
    )
  }
  
  # Componente copular
  rho <- rho_fun(fx)
  
  ux1 <- clamp01(
    (x[1] - L1)/range1
  )
  
  ux2 <- clamp01(
    (x[2] - L2)/range2
  )
  
  v1 <- sim_gaussian_conditional(
    ux1,
    rho
  )
  
  v2 <- sim_gaussian_conditional(
    ux2,
    rho
  )
  
  c(
    L1 + range1*v1,
    L2 + range2*v2
  )
}


# ============================================================
# 14. CÓPULA ADAPTATIVA + MH
# ============================================================

run_copula_MH <- function(seed = 123) {
  
  set.seed(seed)
  
  tiempo_ini <- proc.time()[3]
  
  # ----------------------------------------------------------
  # 14.1 Evaluaciones iniciales
  # ----------------------------------------------------------
  
  X0 <- cbind(
    runif(N_f, L1, U1),
    runif(N_f, L2, U2)
  )
  
  f0 <- eval_f(X0)
  g0 <- abs(f0 - c_level)
  
  ibest <- which.min(g0)
  
  x_curr <- X0[ibest, ]
  f_curr <- f0[ibest]
  
  # ----------------------------------------------------------
  # 14.2 Almacenamiento
  # ----------------------------------------------------------
  
  cadena <- matrix(
    NA_real_,
    nrow = B_mh + 1L,
    ncol = 2
  )
  
  colnames(cadena) <- c("x1", "x2")
  
  f_chain <- numeric(B_mh + 1L)
  
  cadena[1, ] <- x_curr
  f_chain[1]  <- f_curr
  
  accepted <- logical(B_mh)
  
  # Guardamos TODAS las propuestas evaluadas
  propuestas <- matrix(
    NA_real_,
    nrow = B_mh,
    ncol = 2
  )
  
  colnames(propuestas) <- c("x1", "x2")
  
  f_propuestas <- numeric(B_mh)
  
  # ----------------------------------------------------------
  # 14.3 MH
  # ----------------------------------------------------------
  
  for (b in seq_len(B_mh)) {
    
    y <- propose_point(
      x = x_curr,
      fx = f_curr
    )
    
    fy <- f_interp(y)
    
    propuestas[b, ] <- y
    f_propuestas[b] <- fy
    
    # Target
    log_pi_ratio <-
      -(
        abs(fy - c_level) -
          abs(f_curr - c_level)
      ) / tau_mh
    
    # Hastings
    log_q_fwd <- log_q_mix(
      y = y,
      x = x_curr,
      fx = f_curr
    )
    
    log_q_rev <- log_q_mix(
      y = x_curr,
      x = y,
      fx = fy
    )
    
    log_alpha <- min(
      0,
      log_pi_ratio +
        log_q_rev -
        log_q_fwd
    )
    
    if (log(runif(1)) < log_alpha) {
      
      x_curr <- y
      f_curr <- fy
      
      accepted[b] <- TRUE
    }
    
    cadena[b + 1L, ] <- x_curr
    f_chain[b + 1L]  <- f_curr
  }
  
  # ----------------------------------------------------------
  # 14.4 Burn-in
  # ----------------------------------------------------------
  
  burn <- floor(
    burn_frac*nrow(cadena)
  )
  
  keep <- seq.int(
    burn + 1L,
    nrow(cadena)
  )
  
  X_mh <- cadena[keep, , drop = FALSE]
  
  # ----------------------------------------------------------
  # 14.5 Banco de evaluaciones reales
  # ----------------------------------------------------------
  
  Xref <- rbind(
    X0,
    propuestas
  )
  
  fref <- c(
    f0,
    f_propuestas
  )
  
  gref <- abs(
    fref - c_level
  )
  
  tiempo <- proc.time()[3] - tiempo_ini
  
  list(
    X0 = X0,
    f0 = f0,
    cadena = cadena,
    f_chain = f_chain,
    X_mh = X_mh,
    Xref = Xref,
    fref = fref,
    gref = gref,
    propuestas = propuestas,
    f_propuestas = f_propuestas,
    acceptance = mean(accepted),
    N_f = length(fref),
    tiempo = tiempo
  )
}


# ============================================================
# 15. CORRER CÓPULA + MH
# ============================================================

cat("\nCorriendo Cópula adaptativa + MH...\n")

res_mh <- run_copula_MH(
  seed = 123
)

cat(
  "Aceptación MH:",
  round(res_mh$acceptance, 4),
  "\n"
)

cat(
  "Evaluaciones reales:",
  res_mh$N_f,
  "\n"
)

cat(
  "Tiempo:",
  round(res_mh$tiempo, 2),
  "s\n"
)


# ============================================================
# 16. REFINAMIENTO KNN PONDERADO — CORREGIDO
#
# w_i = exp(-g_i/lambda)/(d_i + eps)
# ============================================================

refinar_KNN <- function(
    X,
    Xref,
    gref,
    K,
    lambda = 5,
    eps = 1e-6,
    tol = 1e-12
) {
  
  X    <- as.matrix(X)
  Xref <- as.matrix(Xref)
  gref <- as.numeric(gref)
  
  # ----------------------------------------------------------
  # Controles
  # ----------------------------------------------------------
  
  if (ncol(X) != 2L || ncol(Xref) != 2L) {
    stop("X y Xref deben tener dos columnas.")
  }
  
  if (nrow(Xref) != length(gref)) {
    stop("nrow(Xref) debe coincidir con length(gref).")
  }
  
  if (K >= nrow(Xref)) {
    stop("K debe ser menor que el número de puntos de referencia.")
  }
  
  n <- nrow(X)
  
  Xout <- matrix(
    NA_real_,
    nrow = n,
    ncol = 2
  )
  
  colnames(Xout) <- c("x1", "x2")
  
  # Pedimos vecinos extra porque el propio X_t puede
  # aparecer dentro del banco Xref.
  k_busqueda <- min(
    K + 10L,
    nrow(Xref)
  )
  
  knn <- FNN::get.knnx(
    data  = Xref,
    query = X,
    k     = k_busqueda
  )
  
  # ----------------------------------------------------------
  # Refinar cada estado X_t
  # ----------------------------------------------------------
  
  for (i in seq_len(n)) {
    
    ids  <- knn$nn.index[i, ]
    dist <- knn$nn.dist[i, ]
    
    # Eliminar coincidencias exactas o numéricamente idénticas
    keep <- dist > tol
    
    ids  <- ids[keep]
    dist <- dist[keep]
    
    # --------------------------------------------------------
    # Si después de eliminar coincidencias quedan menos de K,
    # calculamos todas las distancias explícitamente.
    # --------------------------------------------------------
    
    if (length(ids) < K) {
      
      diferencias <- sweep(
        Xref,
        MARGIN = 2,
        STATS = X[i, ],
        FUN = "-"
      )
      
      d_all <- sqrt(
        rowSums(diferencias^2)
      )
      
      ord <- order(d_all)
      
      ord <- ord[
        d_all[ord] > tol
      ]
      
      if (length(ord) < K) {
        stop(
          paste0(
            "No hay suficientes vecinos distintos para K = ",
            K
          )
        )
      }
      
      ids <- ord[seq_len(K)]
      
      dist <- d_all[ids]
      
    } else {
      
      ids  <- ids[seq_len(K)]
      dist <- dist[seq_len(K)]
    }
    
    # --------------------------------------------------------
    # Pesos:
    #
    #   gref[ids] = |f(X_(i)) - c|
    #   dist      = ||X_(i) - X_t||
    # --------------------------------------------------------
    
    w <- exp(
      -gref[ids] / lambda
    ) / (
      dist + eps
    )
    
    # Seguridad numérica
    if (
      any(!is.finite(w)) ||
      sum(w) <= 0
    ) {
      
      w <- 1 / (
        dist + eps
      )
    }
    
    # --------------------------------------------------------
    # Promedio ponderado
    # --------------------------------------------------------
    
    Xout[i, ] <-
      colSums(
        Xref[ids, , drop = FALSE] * w
      ) /
      sum(w)
  }
  
  return(Xout)
}


# ============================================================
# 17. KNN(5), KNN(10), KNN(100)
# ============================================================

cat("\nAplicando KNN(5)...\n")

X_knn5 <- refinar_KNN(
  X      = res_mh$X_mh,
  Xref   = res_mh$Xref,
  gref   = res_mh$gref,
  K      = 5,
  lambda = lambda_knn,
  eps    = eps_dist,
  tol    = tol_self
)


cat("Aplicando KNN(10)...\n")

X_knn10 <- refinar_KNN(
  X      = res_mh$X_mh,
  Xref   = res_mh$Xref,
  gref   = res_mh$gref,
  K      = 10,
  lambda = lambda_knn,
  eps    = eps_dist,
  tol    = tol_self
)


cat("Aplicando KNN(100)...\n")

X_knn100 <- refinar_KNN(
  X      = res_mh$X_mh,
  Xref   = res_mh$Xref,
  gref   = res_mh$gref,
  K      = 100,
  lambda = lambda_knn,
  eps    = eps_dist,
  tol    = tol_self
)

cat("KNN finalizado correctamente.\n")

# ============================================================
# 17. KNN(5), KNN(10), KNN(100)
# ============================================================

cat("\nAplicando KNN...\n")

X_knn5 <- refinar_KNN(
  X = res_mh$X_mh,
  Xref = res_mh$Xref,
  gref = res_mh$gref,
  K = 5
)

X_knn10 <- refinar_KNN(
  X = res_mh$X_mh,
  Xref = res_mh$Xref,
  gref = res_mh$gref,
  K = 10
)

X_knn100 <- refinar_KNN(
  X = res_mh$X_mh,
  Xref = res_mh$Xref,
  gref = res_mh$gref,
  K = 100
)


# ============================================================
# 18. CONJUNTO DE NIVEL DE REFERENCIA
#
# Se extrae el contorno c=230 de una grilla fina.
# SOLO se utiliza para evaluar.
# ============================================================

cat("\nConstruyendo referencia L_c...\n")

n_ref1 <- 600L
n_ref2 <- 450L

gx <- seq(
  L1,
  U1,
  length.out = n_ref1
)

gy <- seq(
  L2,
  U2,
  length.out = n_ref2
)

# Evaluación rápida de la interpolación
grid_ref <- expand.grid(
  x1 = gx,
  x2 = gy
)

grid_ref$z <- eval_f(
  as.matrix(
    grid_ref[, c("x1", "x2")]
  )
)

Zref <- matrix(
  grid_ref$z,
  nrow = length(gx),
  ncol = length(gy)
)

cl <- contourLines(
  x = gx,
  y = gy,
  z = Zref,
  levels = c_level
)

if (length(cl) == 0L) {
  stop("No se encontró el nivel c=230.")
}

L_ref_raw <- do.call(
  rbind,
  lapply(
    cl,
    function(a) {
      cbind(
        x1 = a$x,
        x2 = a$y
      )
    }
  )
)

# Si hay demasiados puntos, remuestreamos aproximadamente
if (nrow(L_ref_raw) > M_ref) {
  
  ids <- round(
    seq(
      1,
      nrow(L_ref_raw),
      length.out = M_ref
    )
  )
  
  L_ref <- L_ref_raw[ids, , drop = FALSE]
  
} else {
  
  L_ref <- L_ref_raw
}

cat(
  "Puntos de referencia:",
  nrow(L_ref),
  "\n"
)


# ============================================================
# 19. MÉTRICAS
#
# Distancias en dominio normalizado [0,1]^2
# ============================================================

metricas_curva <- function(
    Xhat,
    Lref,
    delta = delta_metric
) {
  
  Xhat <- as.matrix(Xhat)
  Lref <- as.matrix(Lref)
  
  # Eliminar NA / infinitos
  Xhat <- Xhat[
    apply(Xhat, 1, function(z) all(is.finite(z))),
    ,
    drop = FALSE
  ]
  
  if (
    nrow(Xhat) == 0L ||
    nrow(Lref) == 0L
  ) {
    
    return(
      c(
        precision = NA_real_,
        cobertura = NA_real_,
        F1 = NA_real_
      )
    )
  }
  
  Xu <- to_unit(Xhat)
  Lu <- to_unit(Lref)
  
  # Distancia estimación -> referencia
  d_hat <- FNN::get.knnx(
    data = Lu,
    query = Xu,
    k = 1
  )$nn.dist[, 1]
  
  precision <- mean(
    d_hat <= delta
  )
  
  # Distancia referencia -> estimación
  d_ref <- FNN::get.knnx(
    data = Xu,
    query = Lu,
    k = 1
  )$nn.dist[, 1]
  
  cobertura <- mean(
    d_ref <= delta
  )
  
  F1 <-
    if (
      precision + cobertura > 0
    ) {
      2*precision*cobertura /
        (precision + cobertura)
    } else {
      0
    }
  
  c(
    precision = precision,
    cobertura = cobertura,
    F1 = F1
  )
}


# ============================================================
# 20. GP RANDOM STRADDLE
#
# 200 evaluaciones reales de f
# ============================================================

run_GP_straddle <- function(seed = 321) {
  
  set.seed(seed)
  
  tiempo_ini <- proc.time()[3]
  
  # ----------------------------------------------------------
  # Diseño inicial
  # ----------------------------------------------------------
  
  X_train <- cbind(
    runif(N_init_gp, L1, U1),
    runif(N_init_gp, L2, U2)
  )
  
  y_train <- eval_f(X_train)
  
  # Trabajamos en coordenadas normalizadas para el GP
  X_train_u <- to_unit(X_train)
  
  # ----------------------------------------------------------
  # Active learning
  # ----------------------------------------------------------
  
  while (nrow(X_train) < N_gp) {
    
    modelo <- DiceKriging::km(
      design = data.frame(
        x1 = X_train_u[, 1],
        x2 = X_train_u[, 2]
      ),
      response = y_train,
      covtype = "gauss",
      nugget.estim = TRUE,
      control = list(trace = FALSE)
    )
    
    Ucand <- cbind(
      runif(N_cand_gp),
      runif(N_cand_gp)
    )
    
    pred <- predict(
      modelo,
      newdata = data.frame(
        x1 = Ucand[, 1],
        x2 = Ucand[, 2]
      ),
      type = "UK",
      checkNames = FALSE
    )
    
    mu <- as.numeric(pred$mean)
    sd <- sqrt(
      pmax(
        as.numeric(pred$sd)^2,
        1e-12
      )
    )
    
    # Straddle:
    # mayor = más interesante
    score <-
      1.96*sd -
      abs(mu - c_level)
    
    ibest <- which.max(score)
    
    x_new <- from_unit(
      matrix(
        Ucand[ibest, ],
        nrow = 1
      )
    )[1, ]
    
    y_new <- f_interp(x_new)
    
    X_train <- rbind(
      X_train,
      x_new
    )
    
    y_train <- c(
      y_train,
      y_new
    )
    
    X_train_u <- rbind(
      X_train_u,
      Ucand[ibest, ]
    )
    
    if (
      nrow(X_train) %% 25 == 0
    ) {
      cat(
        "GP:",
        nrow(X_train),
        "/",
        N_gp,
        "\n"
      )
    }
  }
  
  # ----------------------------------------------------------
  # Modelo final
  # ----------------------------------------------------------
  
  modelo <- DiceKriging::km(
    design = data.frame(
      x1 = X_train_u[, 1],
      x2 = X_train_u[, 2]
    ),
    response = y_train,
    covtype = "gauss",
    nugget.estim = TRUE,
    control = list(trace = FALSE)
  )
  
  # ----------------------------------------------------------
  # Predicción densa
  # ----------------------------------------------------------
  
  n_side <- ceiling(
    sqrt(N_pred_gp)
  )
  
  u1_pred <- seq(
    0,
    1,
    length.out = n_side
  )
  
  u2_pred <- seq(
    0,
    1,
    length.out = n_side
  )
  
  Upred <- expand.grid(
    x1 = u1_pred,
    x2 = u2_pred
  )
  
  pred_final <- predict(
    modelo,
    newdata = Upred,
    type = "UK",
    checkNames = FALSE
  )
  
  mu_pred <- as.numeric(
    pred_final$mean
  )
  
  Zgp <- matrix(
    mu_pred,
    nrow = length(u1_pred),
    ncol = length(u2_pred)
  )
  
  cl_gp <- contourLines(
    x = u1_pred,
    y = u2_pred,
    z = Zgp,
    levels = c_level
  )
  
  if (length(cl_gp) == 0L) {
    
    X_gp <- matrix(
      numeric(0),
      ncol = 2
    )
    
  } else {
    
    X_gp_u <- do.call(
      rbind,
      lapply(
        cl_gp,
        function(a) {
          cbind(
            a$x,
            a$y
          )
        }
      )
    )
    
    X_gp <- from_unit(
      X_gp_u
    )
    
    colnames(X_gp) <- c(
      "x1",
      "x2"
    )
  }
  
  tiempo <- proc.time()[3] - tiempo_ini
  
  list(
    X = X_gp,
    X_train = X_train,
    y_train = y_train,
    modelo = modelo,
    tiempo = tiempo,
    N_f = nrow(X_train)
  )
}


# ============================================================
# 21. CORRER GP
# ============================================================

cat("\n============================================\n")
cat("GP RANDOM STRADDLE\n")
cat("============================================\n")

res_gp <- run_GP_straddle(
  seed = 321
)

cat(
  "Evaluaciones reales GP:",
  res_gp$N_f,
  "\n"
)

cat(
  "Tiempo GP:",
  round(res_gp$tiempo, 2),
  "s\n"
)


# ============================================================
# 22. MÉTRICAS DE TODOS LOS MÉTODOS
# ============================================================

m_mh <- metricas_curva(
  res_mh$X_mh,
  L_ref
)

m_k5 <- metricas_curva(
  X_knn5,
  L_ref
)

m_k10 <- metricas_curva(
  X_knn10,
  L_ref
)

m_k100 <- metricas_curva(
  X_knn100,
  L_ref
)

m_gp <- metricas_curva(
  res_gp$X,
  L_ref
)


# ============================================================
# 23. TABLA COMPARATIVA
# ============================================================

tabla_resultados <- data.frame(
  
  metodo = c(
    "Copula adaptativa + MH",
    "MH + KNN(5)",
    "MH + KNN(10)",
    "MH + KNN(100)",
    "GP Random Straddle"
  ),
  
  precision = c(
    m_mh["precision"],
    m_k5["precision"],
    m_k10["precision"],
    m_k100["precision"],
    m_gp["precision"]
  ),
  
  cobertura = c(
    m_mh["cobertura"],
    m_k5["cobertura"],
    m_k10["cobertura"],
    m_k100["cobertura"],
    m_gp["cobertura"]
  ),
  
  F1 = c(
    m_mh["F1"],
    m_k5["F1"],
    m_k10["F1"],
    m_k100["F1"],
    m_gp["F1"]
  ),
  
  aceptacion_MH = c(
    res_mh$acceptance,
    res_mh$acceptance,
    res_mh$acceptance,
    res_mh$acceptance,
    NA
  ),
  
  N_f = c(
    res_mh$N_f,
    res_mh$N_f,
    res_mh$N_f,
    res_mh$N_f,
    res_gp$N_f
  ),
  
  tiempo_seg = c(
    res_mh$tiempo,
    res_mh$tiempo,
    res_mh$tiempo,
    res_mh$tiempo,
    res_gp$tiempo
  )
)

tabla_resultados[
  c(
    "precision",
    "cobertura",
    "F1",
    "aceptacion_MH",
    "tiempo_seg"
  )
] <-
  lapply(
    tabla_resultados[
      c(
        "precision",
        "cobertura",
        "F1",
        "aceptacion_MH",
        "tiempo_seg"
      )
    ],
    round,
    4
  )

cat("\n============================================\n")
cat("RESULTADOS — LINGOTES\n")
cat("============================================\n")

print(
  tabla_resultados,
  row.names = FALSE
)


# ============================================================
# 24. DATA FRAMES PARA GRÁFICOS
# ============================================================

df_ref <- data.frame(
  x1 = L_ref[, 1],
  x2 = L_ref[, 2]
)

df_mh <- data.frame(
  x1 = res_mh$X_mh[, 1],
  x2 = res_mh$X_mh[, 2]
)

df_k5 <- data.frame(
  x1 = X_knn5[, 1],
  x2 = X_knn5[, 2]
)

df_k10 <- data.frame(
  x1 = X_knn10[, 1],
  x2 = X_knn10[, 2]
)

df_k100 <- data.frame(
  x1 = X_knn100[, 1],
  x2 = X_knn100[, 2]
)

df_gp <- data.frame(
  x1 = res_gp$X[, 1],
  x2 = res_gp$X[, 2]
)

df_gp_real <- data.frame(
  x1 = res_gp$X_train[, 1],
  x2 = res_gp$X_train[, 2]
)


# ============================================================
# 25. FUNCIÓN DE GRÁFICO
# ============================================================

plot_metodo <- function(
    df_est,
    titulo,
    color_est = "red",
    alpha_est = 0.65,
    size_est = 0.55
) {
  
  ggplot() +
    
    geom_point(
      data = df_ref,
      aes(
        x = x1,
        y = x2
      ),
      color = "black",
      size = 0.35,
      alpha = 0.75
    ) +
    
    geom_point(
      data = df_est,
      aes(
        x = x1,
        y = x2
      ),
      color = color_est,
      size = size_est,
      alpha = alpha_est
    ) +
    
    coord_equal() +
    
    labs(
      x = expression(x[1]),
      y = expression(x[2]),
      title = titulo
    ) +
    
    theme_bw(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      )
    )
}


# ============================================================
# 26. GRÁFICOS INDIVIDUALES
# ============================================================

p_mh <- plot_metodo(
  df_mh,
  "Lingotes: Cópula adaptativa + MH",
  color_est = "red"
)

p_k5 <- plot_metodo(
  df_k5,
  "Lingotes: MH + KNN(5)",
  color_est = "red"
)

p_k10 <- plot_metodo(
  df_k10,
  "Lingotes: MH + KNN(10)",
  color_est = "red"
)

p_k100 <- plot_metodo(
  df_k100,
  "Lingotes: MH + KNN(100)",
  color_est = "red"
)

p_gp <- plot_metodo(
  df_gp,
  "Lingotes: GP Random Straddle",
  color_est = "blue",
  alpha_est = 0.8,
  size_est = 0.7
)

p_gp_real <- plot_metodo(
  df_gp_real,
  "Lingotes: Evaluaciones reales de f (GP)",
  color_est = "darkgreen",
  alpha_est = 0.8,
  size_est = 1
)

print(p_mh)
print(p_k5)
print(p_k10)
print(p_k100)
print(p_gp)
print(p_gp_real)


# ============================================================
# 27. FIGURA COMPARATIVA 2 x 3
# ============================================================

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Instalá patchwork con install.packages('patchwork')")
}

library(patchwork)

p_comparacion <-
  (
    p_mh |
      p_k5 |
      p_k10
  ) /
  (
    p_k100 |
      p_gp |
      p_gp_real
  )

print(
  p_comparacion
)


# ============================================================
# 28. HEATMAP DE LIFETIME + NIVEL c=230
# ============================================================

p_heatmap <- ggplot(
  datos,
  aes(
    x = x1,
    y = x2,
    z = lifetime
  )
) +
  
  geom_raster(
    aes(
      fill = lifetime
    )
  ) +
  
  geom_contour(
    breaks = c_level,
    color = "white",
    linewidth = 1
  ) +
  
  coord_equal() +
  
  scale_fill_viridis_c(
    option = "turbo",
    name = "Lifetime"
  ) +
  
  labs(
    x = expression(x[1]),
    y = expression(x[2]),
    title = "Lingote de silicio",
    subtitle = "Nivel objetivo c = 230"
  ) +
  
  theme_bw(
    base_size = 14
  )

print(
  p_heatmap
)


# ============================================================
# 29. HEATMAP + KNN(100)
# ============================================================

p_heatmap_k100 <- ggplot(
  datos,
  aes(
    x = x1,
    y = x2,
    z = lifetime
  )
) +
  
  geom_raster(
    aes(
      fill = lifetime
    )
  ) +
  
  geom_contour(
    breaks = c_level,
    color = "white",
    linewidth = 1
  ) +
  
  geom_point(
    data = df_k100,
    aes(
      x = x1,
      y = x2
    ),
    inherit.aes = FALSE,
    color = "magenta",
    size = 0.6,
    alpha = 0.65
  ) +
  
  coord_equal() +
  
  scale_fill_viridis_c(
    option = "turbo",
    name = "Lifetime"
  ) +
  
  labs(
    x = expression(x[1]),
    y = expression(x[2]),
    title = "Lingote: MH + KNN(100)",
    subtitle = "Nivel objetivo c = 230"
  ) +
  
  theme_bw(
    base_size = 14
  )

print(
  p_heatmap_k100
)


# ============================================================
# 30. RESUMEN
# ============================================================

cat("\n============================================\n")
cat("RESUMEN FINAL\n")
cat("============================================\n")

cat(
  "Nivel objetivo:",
  c_level,
  "\n"
)

cat(
  "Evaluaciones Cópula + MH:",
  res_mh$N_f,
  "\n"
)

cat(
  "Aceptación MH:",
  round(res_mh$acceptance, 4),
  "\n"
)

cat(
  "Evaluaciones GP:",
  res_gp$N_f,
  "\n"
)

cat(
  "delta métricas:",
  delta_metric,
  "\n"
)

cat("\n")

print(
  tabla_resultados,
  row.names = FALSE
)

cat("============================================\n")