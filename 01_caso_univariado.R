rm(list = ls())

library(VineCopula)
library(rootSolve)
library(reshape2)

set.seed(34)

## ───────────────────────────────────────────────────────────────
## Funciones de prueba
## ───────────────────────────────────────────────────────────────
ff <- list(
  F = function(x) cos(1 + x^2),
  L = -2,
  U = 4,
  x_true = c(
    -sqrt(pi * (0 + 0.5) - 1),
    -sqrt(pi * (1 + 0.5) - 1),
    sqrt(pi * (2 + 0.5) - 1),
    sqrt(pi * (3 + 0.5) - 1),
    sqrt(pi * (4 + 0.5) - 1),
    sqrt(pi * (0 + 0.5) - 1),
    sqrt(pi * (1 + 0.5) - 1)
  )
)

#gg <- list(
#  F = function(x, aaa = 1, bbb = 1/2, ccc = -1) 2 * (x - aaa) * (x - bbb) * (x - ccc),
#  L = -2,
#  U = 3,
#  x_true = c(-1, 1, 0.5)
#)

hh <- list(
  F = function(x, r = 0.05, K = 10, bet = 0.1, alf = 1) r * x * (1 - x / K) - bet * x^2 / (x^2 + alf^2),
  L = -1,
  U = 9,
  x_true = c(0, 2, 4 - sqrt(11), 4 + sqrt(11))
)

ii <- list(
  F = function(x, bb0 = -3, bb1 = 5, c0 = 2, r0 = 0) {
    (exp(c0 - r0) * (bb0 + x) * (bb1 - x)) / ((bb0 + x + 1) * (bb1 - x - 1)) - 1 
  },
  L = 1,
  U = 6,
  x_true = c(
    (3 - 4 * exp(2) + sqrt(1 - exp(2) + exp(4))) / (1 - exp(2)),
    (-3 + 4 * exp(2) + sqrt(1 - exp(2) + exp(4))) / (exp(2) - 1)
  )
)


ww <- list(
  # W_n(x) = ∏_{k=1}^n (x - k)
  F = function(x, n = 20) {
    x <- as.numeric(x)
    sapply(x, function(xx) prod(xx - seq_len(n)))
  },
  L = 0.5,          # intervalo cómodo para n=20
  U = 20+0.5,         # así quedan dentro 1,2,...,20
  x_true = 1:20     # raíces verdaderas de W_20(x)
)

## 3. Ejemplo donde uniroot falla

ss <- list(
  F = function(x) (x+2)^2*(x-1)^2*(x-3)^2,
  L = -4,
  U = 5,
  x_true = c(-2,1,3)
)



## 4. Ejemplo adicional: raíces cercanas
##
## Este ejemplo permite estudiar si el procedimiento distingue raíces
## próximas sin introducir un parámetro de distancia entre clusters.
## En particular, las raíces 0.50 y 0.55 están separadas solamente
## por 0.05 unidades.
raices_cercanas <- list(
  F = function(x) {
    2 * (x + 1) * (x - 0.50) * (x - 0.55)
  },
  L = -2,
  U = 2,
  x_true = c(-1, 0.50, 0.55)
)


## ───────────────────────────────────────────────────────────────
## Copulas
## ───────────────────────────────────────────────────────────────

gaussiana <- list(
  phi = function(y, sig = sigg) { 
rho=min(exp(-abs(y/sig)^2),0.999999)
    return(rho)
  },
  nro = 1
)


gumbel <- list(
  phi = function(y,sig=sigg) 1 + 15.9*exp(- (abs(y)/sig)^2),
  nro = 4
)

frank <- list(
  phi = function(y, sig = sigg) {
    rho <- 34.9 * (1 - 2/pi * atan(abs(y)/sig))
    max(rho, 1e-13)
  },
  nro = 5
)


## ───────────────────────────────────────────────────────────────
## Seleccion de la Copula
## ───────────────────────────────────────────────────────────────

cop <- gaussiana# Cambiar por: gaussiana, gumbel, frank si se desea

## ───────────────────────────────────────────────────────────────
## Selección de función
## ───────────────────────────────────────────────────────────────
f1 <- ss         # Cambiar por ff, gg, hh, ii si se desea
f <- f1$F
L <- f1$L
U <- f1$U

## ───────────────────────────────────────────────────────────────
## Parámetros
## ───────────────────────────────────────────────────────────────
N         <- 5000
sigg      <- 3
tau       <- max(abs(f(seq(L,U,0.1)))[abs(f(seq(L,U,0.1))) < Inf])/(16*log(10))
epsilon   <- .1
zero_tol  <- 0.001

## ───────────────────────────────────────────────────────────────
## Inicialización
## ───────────────────────────────────────────────────────────────
puntos        <- numeric(N + 1)
correlaciones <- numeric(N)
puntos[1]     <- runif(1, L, U)
delta=1e-16

## ───────────────────────────────────────────────────────────────
## MCMC con Copula + MH + Exploración
## ───────────────────────────────────────────────────────────────

aplicar <- function(N, zero_tol) {
  
  puntos <- numeric(N + 1)
  puntos[1] <- runif(1, L, U)
  
  # Densidad del kernel de propuesta q(u_to | u_from)
  q_prop <- function(u_to, u_from, f_from) {
    
    # Lejos del conjunto de ceros:
    # propuesta uniforme en (0,1)
    if (abs(f_from) > zero_tol) {
      return(1)
    }
    
    # Cerca del conjunto de ceros:
    # mezcla epsilon * Uniforme +
    #        (1-epsilon) * propuesta condicional por cópula
    
    rho_from <- cop$phi(f_from, sigg)
    
    dens_cop <- BiCopPDF(
      u_to,
      u_from,
      family = cop$nro,
      par = rho_from
    )
    
    epsilon + (1 - epsilon) * dens_cop
  }
  
  
  for (j in seq_len(N)) {
    
    #--------------------------------------------------
    # Estado actual
    #--------------------------------------------------
    
    xj <- puntos[j]
    fx_curr <- f(xj)
    
    uj <- punif(xj, L, U)
    
    #--------------------------------------------------
    # Generación de la propuesta
    #--------------------------------------------------
    
    if (abs(fx_curr) > zero_tol) {
      
      # Lejos de la raíz: propuesta global uniforme
      u_next <- runif(1)
      
    } else {
      
      rho_curr <- cop$phi(fx_curr, sigg)
      
      if (runif(1) < epsilon) {
        
        # componente uniforme de la mezcla
        u_next <- runif(1)
        
      } else {
        
        # componente guiada por cópula
        u_next <- BiCopCondSim(
          1,
          cond.val = uj,
          cond.var = 2,
          BiCop(
            family = cop$nro,
            par = rho_curr
          )
        )
      }
    }
    
    x_next <- qunif(u_next, L, U)
    fx_next <- f(x_next)
    
    #--------------------------------------------------
    # Densidades forward y reverse
    #--------------------------------------------------
    
    q_fwd <- q_prop(
      u_to   = u_next,
      u_from = uj,
      f_from = fx_curr
    )
    
    q_rev <- q_prop(
      u_to   = uj,
      u_from = u_next,
      f_from = fx_next
    )
    
    #--------------------------------------------------
    # Razón de Metropolis-Hastings
    #--------------------------------------------------
    
    # π(x) proporcional a exp(-|f(x)| / tau)
    #
    # pi_next/pi_curr =
    # exp(-( |f(x_next)|-|f(xj)| ) / tau)
    #
    # Lo calculamos así para evitar underflow.
    
    log_pi_ratio <-
      -(abs(fx_next) - abs(fx_curr)) / tau
    
    if (
      is.finite(log_pi_ratio) &&
      is.finite(q_fwd) &&
      is.finite(q_rev) &&
      q_fwd > 0 &&
      q_rev > 0
    ) {
      
      log_hastings_ratio <-
        log_pi_ratio + log(q_rev) - log(q_fwd)
      
      alpha <- min(
        1,
        exp(min(0, log_hastings_ratio))
      )
      
    } else {
      
      alpha <- 0
    }
    
    #--------------------------------------------------
    # Aceptación / rechazo
    #--------------------------------------------------
    
    if (runif(1) < alpha) {
      puntos[j + 1] <- x_next
    } else {
      puntos[j + 1] <- xj
    }
  }
  
  puntos
}

#plot(density(lapply(seq(0.01,10,3000),aplicar, N=10000)[[1]] ))


# --------------------------------------------------
# Trayectoria MCMC
# --------------------------------------------------

puntos <- aplicar(N, zero_tol)

xu <- sort(unique(puntos))
zu <- abs(f(xu))

idx_min <- which(
  zu[2:(length(zu)-1)] < zu[1:(length(zu)-2)] &
    zu[2:(length(zu)-1)] < zu[3:length(zu)]
) + 1

root_estimates <- xu[idx_min]
print(root_estimates)


## ───────────────────────────────────────────────────────────────
## Comparación con raíces verdaderas
## ───────────────────────────────────────────────────────────────
(uni <- uniroot.all(f, c(L, U)))

cat("Suma de |f(x)| en raíces reales (uniroot):", sum(abs(f(uni))), "\n")
if (length(root_estimates) > 0) {
  cat("Suma de |f(x)| en estimadas:", sum(abs(f(root_estimates))), "\n")
}

## ───────────────────────────────────────────────────────────────
## Crear carpeta de salida
## ───────────────────────────────────────────────────────────────
#if (!dir.exists("plots")) dir.create("plots")

## ───────────────────────────────────────────────────────────────
## Guardar PDF estimada
## ───────────────────────────────────────────────────────────────
d <- density(puntos[1:N], from = L, to = U)


#png("plots/estimated_pdf.png", width = 800, height = 600)
plot(d, main = "PDF estimada de la cadena")
abline(v = f1$x_true, col = 4, lty = 1, lwd = 2)      # verdaderas
abline(v = uni, col = 2, lty = 3, lwd = 2)            # uniroot
if (length(root_estimates) > 0) {
  abline(v = root_estimates, col = 3, lty = 2, lwd = 2)  # clustering
}
legend("topleft",
       legend = c("Clustering", "Uniroot", "Raíces verdaderas"),
       col = c(3, 2, 4),
       lty = c(2, 3, 1),
       lwd = 2)
#dev.off()

## ───────────────────────────────────────────────────────────────
## Guardar traza de la cadena
## ───────────────────────────────────────────────────────────────
#png("plots/chain_trace.png", width = 800, height = 600)
plot(puntos, type = "l", main = "Evolución de la cadena", ylab = "x", xlab = "Iteración",ylim=c(L,U))
abline(h = f1$x_true, col = 4)
#dev.off()

xx <- seq(L, U, length.out = 2000)

# ── Data frame con la curva (sin Inf/NaN) ───────────
yy <- f(xx)
df_curve <- data.frame(x = xx, y = yy)
df_curve <- df_curve[is.finite(df_curve$y), ]  # quita asintotas

# ── Data frame con puntos ───────────────────────────
pts <- as.numeric(puntos)
df_pts <- data.frame(x = pts, y = f(pts))
df_pts <- df_pts[is.finite(df_pts$y), ]
df_pts$id <- seq_len(nrow(df_pts))

root_score <- function(true_roots,
                       est_roots,
                       match_tol = 0.05,
                       penalize_distance = TRUE) {
  
  true_roots <- sort(as.numeric(true_roots))
  est_roots  <- sort(as.numeric(est_roots))
  
  K <- length(true_roots)
  M <- length(est_roots)
  
  if (M == 0) {
    return(list(score = 0,
                precision = 0,
                recall = 0,
                f1 = 0,
                matched = 0))
  }
  
  pairs <- expand.grid(i = seq_len(K), j = seq_len(M))
  pairs$d <- abs(true_roots[pairs$i] - est_roots[pairs$j])
  pairs <- pairs[order(pairs$d), ]
  
  matched_true <- rep(FALSE, K)
  matched_est  <- rep(FALSE, M)
  match_d <- numeric(0)
  
  for (row in seq_len(nrow(pairs))) {
    i <- pairs$i[row]
    j <- pairs$j[row]
    d <- pairs$d[row]
    
    if (d <= match_tol && !matched_true[i] && !matched_est[j]) {
      matched_true[i] <- TRUE
      matched_est[j]  <- TRUE
      match_d <- c(match_d, d)
    }
  }
  
  n_matched <- sum(matched_true)
  
  precision <- n_matched / M
  recall    <- n_matched / K
  
  if (precision + recall == 0) {
    f1 <- 0
  } else {
    f1 <- 2 * precision * recall / (precision + recall)
  }
  
  score <- f1
  
  # penalización por distancia
  if (penalize_distance && length(match_d) > 0) {
    dist_penalty <- mean(pmin(match_d / match_tol, 1))
    score <- score * (1 - 0.5 * dist_penalty)
  }
  
  list(score = score,
       precision = precision,
       recall = recall,
       f1 = f1,
       matched = n_matched,
       spurious = M - n_matched,
       missed = K - n_matched)
}

res <- root_score(sort(f1$x_true),
                  sort(root_estimates),
                  match_tol = 0.1)

res

## ───────────────────────────────────────────────────────────────
## Figura: h(x) + densidad empírica
## Paneles horizontalmente alineados
## ───────────────────────────────────────────────────────────────

library(ggplot2)
library(patchwork)

make_tol_rects <- function(xgrid, fx, tol) {
  
  inside <- abs(fx) < tol
  
  if (!any(inside)) {
    return(
      data.frame(
        xmin = numeric(0),
        xmax = numeric(0)
      )
    )
  }
  
  r <- rle(inside)
  
  ends <- cumsum(r$lengths)
  
  starts <- c(
    1,
    head(ends, -1) + 1
  )
  
  idx <- which(r$values)
  
  data.frame(
    xmin = xgrid[starts[idx]],
    xmax = xgrid[ends[idx]]
  )
}

## ───────────────────────────────────────────────────────────────
## Parámetros gráficos
## ───────────────────────────────────────────────────────────────

xlim_plot <- c(L, U)
n_grid <- 4000

## ───────────────────────────────────────────────────────────────
## Curva h(x)
## ───────────────────────────────────────────────────────────────

xx_plot <- seq(
  xlim_plot[1],
  xlim_plot[2],
  length.out = n_grid
)

yy_plot <- f(xx_plot)

df_curve <- data.frame(
  x = xx_plot,
  y = yy_plot
)

df_curve <- df_curve[
  is.finite(df_curve$y),
  ,
  drop = FALSE
]

## ───────────────────────────────────────────────────────────────
## Regiones donde |h(x)| < zero_tol
## ───────────────────────────────────────────────────────────────

rects <- make_tol_rects(
  df_curve$x,
  df_curve$y,
  zero_tol
)

if (nrow(rects) > 0) {
  rects$ymin <- -Inf
  rects$ymax <- Inf
}

## ───────────────────────────────────────────────────────────────
## Raíces verdaderas
## ───────────────────────────────────────────────────────────────

roots_true <- sort(
  as.numeric(f1$x_true)
)

df_true <- data.frame(
  x = roots_true,
  y = 0
)

## ───────────────────────────────────────────────────────────────
## Raíces estimadas
## ───────────────────────────────────────────────────────────────

roots_est <- if (
  exists("root_estimates") &&
  length(root_estimates) > 0
) {
  
  sort(
    as.numeric(root_estimates)
  )
  
} else {
  
  numeric(0)
}

df_est <- data.frame(
  x = roots_est,
  y = 0
)

## ───────────────────────────────────────────────────────────────
## Trayectoria para densidad
## ───────────────────────────────────────────────────────────────

chain_x <- as.numeric(
  puntos[-1]
)

chain_x <- chain_x[
  is.finite(chain_x)
]

df_chain <- data.frame(
  x = chain_x
)

## ───────────────────────────────────────────────────────────────
## Límites eje y panel superior
## ───────────────────────────────────────────────────────────────

y_med <- median(
  df_curve$y,
  na.rm = TRUE
)

y_mad <- mad(
  df_curve$y,
  na.rm = TRUE
)

if (
  !is.finite(y_mad) ||
  y_mad == 0
) {
  
  y_mad <- sd(
    df_curve$y,
    na.rm = TRUE
  )
}

if (
  !is.finite(y_mad) ||
  y_mad == 0
) {
  
  y_mad <- 1
}

ymin_plot <- min(
  -4 * y_mad + y_med,
  -1.5 * abs(zero_tol)
)

ymax_plot <- max(
  4 * y_mad + y_med,
  1.5 * abs(zero_tol)
)

## ───────────────────────────────────────────────────────────────
## Escala x COMÚN
## ───────────────────────────────────────────────────────────────

x_scale_common <- scale_x_continuous(
  limits = xlim_plot,
  expand = c(0, 0)
)

## ───────────────────────────────────────────────────────────────
## Tema común para alinear paneles
## ───────────────────────────────────────────────────────────────

theme_common <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      linewidth = 0.25,
      colour = "grey88"
    ),
    axis.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      2, 8, 2, 8
    )
  )

## ───────────────────────────────────────────────────────────────
## Panel superior
## ───────────────────────────────────────────────────────────────

p_top <- ggplot(
  df_curve,
  aes(x = x, y = y)
) +
  
  ## Región de activación
  {
    if (nrow(rects) > 0)
      geom_rect(
        data = rects,
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax,
          fill = "Región de activación"
        ),
        inherit.aes = FALSE,
        alpha = 0.32
      )
  } +
  
  ## Curva h(x)
  geom_line(
    linewidth = 0.9,
    colour = "black"
  ) +
  
  ## y = 0
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey45"
  ) +
  
  ## + zero_tol
  geom_hline(
    aes(
      yintercept = zero_tol,
      linetype = "Umbral"
    ),
    linewidth = 0.35,
    colour = "firebrick2"
  ) +
  
  ## - zero_tol
  geom_hline(
    aes(
      yintercept = -zero_tol,
      linetype = "Umbral"
    ),
    linewidth = 0.35,
    colour = "firebrick2"
  ) +
  
  ## Primero: raíces estimadas
  {
    if (nrow(df_est) > 0)
      geom_point(
        data = df_est,
        aes(
          x = x,
          y = y,
          colour = "Raíces estimadas",
          shape = "Raíces estimadas"
        ),
        inherit.aes = FALSE,
        stroke = 1.3,
        size = 3.2
      )
  } +
  
  ## Después: raíces verdaderas
  ## así quedan dibujadas ENCIMA de las estimadas
  {
    if (nrow(df_true) > 0)
      geom_point(
        data = df_true,
        aes(
          x = x,
          y = y,
          colour = "Raíces verdaderas",
          shape = "Raíces verdaderas"
        ),
        inherit.aes = FALSE,
        size = 2.8
      )
  } +
  
  scale_fill_manual(
    name = NULL,
    values = c(
      "Región de activación" = "palegreen3"
    ),
    labels = c(
      "Región de activación" =
        expression("|h(x)| < zero_tol")
    )
  ) +
  
  scale_colour_manual(
    name = NULL,
    values = c(
      "Raíces verdaderas" = "darkorange2",
      "Raíces estimadas" = "deepskyblue3"
    )
  ) +
  
  scale_shape_manual(
    name = NULL,
    values = c(
      "Raíces verdaderas" = 16,
      "Raíces estimadas" = 3
    )
  ) +
  
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Umbral" = "dashed"
    ),
    labels = c(
      "Umbral" =
        expression("|h(x)| == zero_tol")
    )
  ) +
  
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(alpha = 0.32)
    ),
    linetype = guide_legend(
      order = 2,
      override.aes = list(
        colour = "firebrick2",
        linewidth = 0.6
      )
    ),
    colour = guide_legend(order = 3),
    shape = guide_legend(order = 3)
  ) +
  
  labs(
    x = NULL,
    y = expression(h(x))
  ) +
  
  x_scale_common +
  
  coord_cartesian(
    ylim = c(
      ymin_plot,
      ymax_plot
    ),
    expand = FALSE
  ) +
  
  theme_common +
  
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 9)
  )

## ───────────────────────────────────────────────────────────────
## Panel inferior
## ───────────────────────────────────────────────────────────────

p_bottom <- ggplot(
  df_chain,
  aes(x = x)
) +
  
  ## Región de activación
  {
    if (nrow(rects) > 0)
      geom_rect(
        data = rects,
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = -Inf,
          ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "palegreen3",
        alpha = 0.32
      )
  } +
  
  ## Densidad
  geom_density(
    fill = "grey35",
    alpha = 0.22,
    linewidth = 0.8,
    colour = "black"
  ) +
  
  ## Primero estimadas
  {
    if (nrow(df_est) > 0)
      geom_point(
        data = df_est,
        aes(
          x = x,
          y = 0
        ),
        inherit.aes = FALSE,
        shape = 3,
        stroke = 1.3,
        size = 3.2,
        colour = "deepskyblue3"
      )
  } +
  
  ## Después verdaderas
  ## para que queden encima
  {
    if (nrow(df_true) > 0)
      geom_point(
        data = df_true,
        aes(
          x = x,
          y = 0
        ),
        inherit.aes = FALSE,
        shape = 16,
        size = 2.8,
        colour = "darkorange2"
      )
  } +
  
  labs(
    x = "x",
    y = "Densidad"
  ) +
  
  x_scale_common +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.04, 0.05)
    )
  ) +
  
  coord_cartesian(
    expand = FALSE
  ) +
  
  theme_common +
  
  theme(
    legend.position = "none"
  )

## ───────────────────────────────────────────────────────────────
## Composición
## ───────────────────────────────────────────────────────────────

fig_bonita <- p_top / p_bottom +
  plot_layout(
    heights = c(1.15, 0.85)
  )

## ───────────────────────────────────────────────────────────────
## Mostrar
## ───────────────────────────────────────────────────────────────

print(fig_bonita)

## ───────────────────────────────────────────────────────────────
## Guardar
## ───────────────────────────────────────────────────────────────

if (!dir.exists("plots")) {
  dir.create("plots")
}

ggsave(
  filename = sprintf(
    "plots/fig_h_tol_%g.pdf",
    zero_tol
  ),
  plot = fig_bonita,
  width = 8,
  height = 5.4
)

library(ggplot2)
library(patchwork)

## ============================================================
## Datos
## ============================================================

df_chain <- data.frame(
  iter = seq_along(puntos),
  x    = as.numeric(puntos)
)

roots_true <- sort(
  as.numeric(f1$x_true)
)

df_roots <- data.frame(
  root = roots_true
)

## ============================================================
## Densidad empírica
## ============================================================

dens <- density(
  df_chain$x,
  from = L,
  to   = U,
  n    = 1000
)

df_dens <- data.frame(
  x       = dens$x,
  density = dens$y
)

## ============================================================
## Polígono para rellenar la densidad
## ============================================================

df_poly <- rbind(
  
  data.frame(
    x = L,
    density = 0
  ),
  
  df_dens,
  
  data.frame(
    x = U,
    density = 0
  )
)

## ============================================================
## PANEL IZQUIERDO: evolución de la cadena
## ============================================================

p_chain <- ggplot(
  df_chain,
  aes(x = iter, y = x)
) +
  
  geom_line(
    linewidth = 0.30,
    colour = "black"
  ) +
  
  ## Raíces verdaderas
  geom_hline(
    data = df_roots,
    aes(yintercept = root),
    colour = "deepskyblue3",
    linewidth = 0.65,
    linetype = "dashed"
  ) +
  
  ## Evitar solapamiento con el panel derecho
  scale_x_continuous(
    breaks = c(
      0,
      1000,
      2000,
      3000,
      4000
    ),
    expand = c(0, 0)
  ) +
  
  ## Misma escala vertical en ambos paneles
  scale_y_continuous(
    limits = c(L, U),
    breaks = pretty(
      c(L, U),
      n = 5
    ),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Iteración",
    y = "x",
    title = "Evolución de la cadena"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      linewidth = 0.25,
      colour = "grey88"
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    plot.margin = margin(
      5, 3, 5, 5
    )
  )

## ============================================================
## PANEL DERECHO: densidad empírica
## ============================================================

p_density <- ggplot() +
  
  ## ----------------------------------------------------------
## Área gris bajo la densidad
## ----------------------------------------------------------

geom_polygon(
  data = df_poly,
  aes(
    x = density,
    y = x,
    group = 1
  ),
  fill = "grey75",
  colour = NA,
  alpha = 0.55
) +
  
  ## ----------------------------------------------------------
## Curva de densidad
##
## IMPORTANTE:
## geom_path() respeta el orden de los puntos.
## geom_line() los reordena por density y genera el desastre.
## ----------------------------------------------------------

geom_path(
  data = df_dens,
  aes(
    x = density,
    y = x,
    group = 1
  ),
  colour = "black",
  linewidth = 0.85
) +
  
  ## ----------------------------------------------------------
## Raíces verdaderas
## ----------------------------------------------------------

geom_hline(
  data = df_roots,
  aes(yintercept = root),
  colour = "deepskyblue3",
  linewidth = 0.8,
  linetype = "dashed"
) +
  
  ## Misma escala vertical que la trayectoria
  scale_y_continuous(
    limits = c(L, U),
    breaks = pretty(
      c(L, U),
      n = 5
    ),
    expand = c(0, 0)
  ) +
  
  scale_x_continuous(
    breaks = pretty(
      c(
        0,
        max(df_dens$density)
      ),
      n = 4
    ),
    expand = expansion(
      mult = c(0, 0.03)
    )
  ) +
  
  labs(
    x = "Densidad",
    y = NULL,
    title = "Densidad empírica"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      linewidth = 0.25,
      colour = "grey88"
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    ## No repetir eje vertical
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    plot.margin = margin(
      5, 5, 5, 3
    )
  )

## ============================================================
## COMBINAR
## ============================================================

fig_chain_density <-
  p_chain + p_density +
  plot_layout(
    widths = c(3.3, 1)
  )

## ============================================================
## MOSTRAR
## ============================================================

print(fig_chain_density)

## ============================================================
## GUARDAR
## ============================================================

if (!dir.exists("plots")) {
  dir.create("plots")
}

ggsave(
  filename = "plots/trayectoria_densidad.pdf",
  plot = fig_chain_density,
  width = 9,
  height = 5
)