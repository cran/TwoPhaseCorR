# Register global variables for CRAN compliance
utils::globalVariables(c("CEF", "rho", "rho_ortho"))

#' Two-Phase Experimental Design Construction and Analysis
#'
#' Constructs and evaluates a two-phase experimental design using cyclic methods.
#' Calculates information matrices and Canonical Efficiency Factor (CEF) under
#' correlated error structures.
#'
#' @name TwoPhaseDesign
#' @param v Integer (>=3). Number of treatments in Phase II.
#' @param rho Numeric (-1 < rho < 1). Correlation coefficient.
#' @param plot Logical. If TRUE (default), generates a CEF plot using ggplot2.
#' @return A list containing the Phase I and Phase II layouts, combined layout,
#'         information matrices for treatment and interaction effects, and a
#'         table and plot of Canonical Efficiency Factors.
#'
#' @import Matrix MASS ggplot2
#' @references
#' McIntyre, G. A. (1955). \emph{Design and analysis of two-phase experiments}. Biometrics, 11(3), 324-334.
#' <doi:10.2307/3001770>
#' @examples
#' result <- TwoPhaseDesign(v = 3, rho = 0.1, plot = FALSE)
#' print(result$cef_table)
#'
#' @export
TwoPhaseDesign <- function(v, rho, plot = TRUE) {
  if (!is.numeric(v) || v < 3 || v != floor(v)) stop("`v` must be an integer >= 3.")
  if (!is.numeric(rho) || abs(rho) >= 1) stop("`rho` must be between -1 and 1.")

  # Phase I design (cyclic construction mod v+1)
  S <- t(1:v)
  d1 <- matrix(NA, nrow = v, ncol = v + 1)
  for (j in 0:v) d1[, j + 1] <- (S + j - 1) %% (v + 1) + 1
  d1 <- t(d1)

  # Phase II base design (cyclic construction mod v)
  Si <- t(1:(v - 1))
  base_d2 <- matrix(NA, nrow = v - 1, ncol = v)
  for (j in 0:(v - 1)) base_d2[, j + 1] <- (Si + j - 1) %% v + 1
  d2 <- matrix(NA, nrow = v - 1, ncol = v * (v + 1))
  for (i in 1:(v + 1)) d2[, ((i - 1) * v + 1):(i * v)] <- base_d2
  d2 <- t(d2)

  # Combined two-phase design
  d11 <- matrix(as.vector(t(d1)), ncol = 1)
  d3 <- cbind(d11, d2)

  v1 <- max(d1)
  v2 <- max(d2)
  tot_trt <- v1 * v2
  nObs <- nrow(d2) * ncol(d2)

  uniq1 <- unique(as.vector(d1))
  uniq2 <- unique(as.vector(d2))

  facto <- matrix(NA, nrow = tot_trt, ncol = 2)
  k <- 1
  for (i in 1:v1) {
    for (j in 1:v2) {
      facto[k, ] <- c(uniq1[i], uniq2[j])
      k <- k + 1
    }
  }

  m <- matrix(1, nrow = nObs, ncol = 1)
  block_new <- matrix(0, nrow = nObs, ncol = nrow(d2))
  obs <- 1
  for (i in 1:nrow(d2)) {
    for (j in 1:ncol(d2)) {
      block_new[obs, i] <- 1
      obs <- obs + 1
    }
  }

  trt_d12 <- matrix(0, nrow = nObs, ncol = tot_trt)
  obs <- 1
  for (i in 1:nrow(d3)) {
    p1 <- d3[i, 1]
    for (j in 2:ncol(d3)) {
      p2 <- d3[i, j]
      match <- which(facto[, 1] == p1 & facto[, 2] == p2)
      if (length(match) > 0) {
        trt_d12[obs, match] <- 1
        obs <- obs + 1
      }
    }
  }

  V <- matrix(rho, nrow = ncol(d2), ncol = ncol(d2))
  diag(V) <- 1
  Vd <- kronecker(diag(nrow(d2)), V)

  x1 <- trt_d12
  x2 <- cbind(block_new, m)
  C_mat <- round(t(x1) %*% ginv(Vd) %*% x1 -
                   t(x1) %*% ginv(Vd) %*% x2 %*%
                   ginv(t(x2) %*% ginv(Vd) %*% x2) %*%
                   t(x2) %*% ginv(Vd) %*% x1, 5)

  # Projection matrices
  p1 <- matrix(1, nrow = v1 - 1, ncol = v1)
  for (i in 1:(v1 - 1)) {
    p1[i, i + 1] <- -i
    if (i + 2 <= v1) p1[i, (i + 2):v1] <- p1[i, (i + 2):v1] - 1
  }

  p2 <- matrix(1, nrow = v2 - 1, ncol = v2)
  for (i in 1:(v2 - 1)) {
    p2[i, i + 1] <- -i
    if (i + 2 <= v2) p2[i, (i + 2):v2] <- p2[i, (i + 2):v2] - 1
  }

  p10 <- matrix(1, nrow = v1, ncol = 1)
  p20 <- matrix(1, nrow = v2, ncol = 1)
  p11 <- kronecker(p1, t(p20))
  p22 <- kronecker(t(p10), p2)
  p12 <- kronecker(p1, p2)

  Np11 <- p11 / sqrt(rowSums(p11^2))
  Np22 <- p22 / sqrt(rowSums(p22^2))
  Np12 <- p12 / sqrt(rowSums(p12^2))

  C_mat_trt1 <- round(Np11 %*% C_mat %*% t(Np11), 5)
  C_mat_trt2 <- round(Np22 %*% C_mat %*% t(Np22), 5)
  C_mat_trt1_trt2 <- round(Np12 %*% C_mat %*% t(Np12), 5)

  # Canonical Efficiency Factor Computation
  lambda2 <- v2 - 2
  r <- v2 - 1
  k2 <- v2 - 1
  rho_ortho <- (v2 - k2) / (k2 * (v2 - 1))
  upper_rho <- 1 - (lambda2 * v2) / (r * k2 * 2)
  lower_rho <- -0.50

  fine_rho_vals <- seq(lower_rho, upper_rho, length.out = 50)
  CEF_vals <- (lambda2 * v2) / (r * k2 * (1 - fine_rho_vals))
  cef_all <- data.frame(rho = fine_rho_vals, CEF = CEF_vals)
  cef_plot <- subset(cef_all, is.finite(CEF) & CEF > 0)

  cef_range <- subset(cef_plot, CEF >= 0.75 & CEF <= 1)
  if (nrow(cef_range) >= 10) {
    idx <- round(seq(1, nrow(cef_range), length.out = 10))
    cef_range <- cef_range[idx, ]
  }

  # ⬇️ Add v, k2, lambda2, rho_ortho to cef_table
  cef_table <- data.frame(
    v = v,
    k2 = k2,
    lambda2 = lambda2,
    rho_ortho = round(rho_ortho, 4),
    rho = cef_range$rho,
    CEF = cef_range$CEF
  )

  if (plot) {
    print(
      ggplot(cef_plot, aes(x = rho, y = CEF)) +
        geom_line(color = "blue") +
        geom_point(color = "blue", size = 1) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
        geom_vline(xintercept = rho_ortho, linetype = "dotted", color = "red") +
        annotate("text", x = max(cef_plot$rho), y = min(cef_plot$CEF, na.rm = TRUE),
                 label = paste("Design achieves orthogonal efficiency at rho =", round(rho_ortho, 2)),
                 hjust = 1, vjust = 0, color = "red", size = 4) +
        labs(title = "Canonical Efficiency Factor vs Intra-block Correlation",
             x = "Intra-block Correlation",
             y = "Canonical Efficiency Factor (CEF)") +
        theme_minimal()
    )
  }

  return(list(
    d1 = d1,
    d2 = d2,
    d3 = d3,
    C_mat_trt1 = C_mat_trt1,
    C_mat_trt2 = C_mat_trt2,
    C_mat_trt1_trt2 = C_mat_trt1_trt2,
    cef_plot = cef_plot,
    cef_table = cef_table
  ))
}
