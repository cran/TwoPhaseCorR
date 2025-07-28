# Register global variables for CRAN compliance
utils::globalVariables(c("Class", "Efficiency", "Rho"))

#' Construct Two-Phase Experimental Designs with Correlated Errors
#'
#' @description Constructs a two-phase experimental design, computes component information matrices, and evaluates the efficiency factor across intra-block correlations values.
#'
#' @name TwoPhaseDesign
#' @param v Integer (greater than or equal to 3). Number of treatments in Phase II.
#' @param rho Intra-block correlation coefficient. A numeric value in (-1, 1).
#' @param plot Logical. If TRUE, plots Efficiency Factor vs Intra-block correlation coefficient.
#' @param n_table Number of efficiency values to display in the output table.
#' @param tol Tolerance level for classifying Efficiency Factor as approximately equal to 1. Default is 1e-3.
#'
#' @return A list with design layouts, component information matrices, efficiency plot, summary efficiency table, and filtered efficiency table.
#'
#' @import Matrix dplyr ggplot2
#' @importFrom MASS ginv
#'
#' @references
#' McIntyre, G. A. (1955). Design and analysis of two-phase experiments. \emph{Biometrics}, 11(3), 324-334.
#' <doi:10.2307/3001770>
#'
#' @examples
#' result <- TwoPhaseDesign(v = 4, rho = 0.25, plot = FALSE)
#' print(result$eff_summary)
#'
#' @export
TwoPhaseDesign <- function(v, rho, plot = TRUE, n_table = 10, tol = 1e-3) {
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

  ##--- Efficiency Computation ---##
  efficiency_phi <- function(rho, v2, k2) {
    num <- v2 * (k2 - 1) * (1 + rho * k2 - rho)
    den <- k2 * ((1 - rho) * (v2 - 1) + rho * v2 * (k2 - 1))
    return(num / den)
  }

  k2 <- v2 - 1
  lambda <- v2 - 2
  rho_seq <- seq(-0.99, 0.99, by = 0.01)
  eff_vals <- efficiency_phi(rho_seq, v2, k2)

  # Classification with tolerance
  class_labels <- case_when(
    abs(eff_vals - 1) <= tol ~ "= 1",
    eff_vals < 1 - tol ~ "< 1",
    eff_vals > 1 + tol ~ "> 1"
  )

  df_tab <- data.frame(
    Rho = rho_seq,
    Efficiency = round(eff_vals, 4),
    Class = factor(class_labels, levels = c("< 1", "= 1", "> 1"))
  )

  eff_summary <- df_tab %>%
    group_by(Class) %>%
    summarise(
      Rho_Min = round(min(Rho), 2),
      Rho_Max = round(max(Rho), 2),
      .groups = "drop"
    ) %>%
    mutate(Interpretation = case_when(
      Class == "= 1" ~ "Exactly Efficient",
      Class == "< 1" ~ "At Par Efficient in the Class",
      Class == "> 1" ~ "Efficient in the Class"
    ))

  eff_table <- df_tab %>%
    filter(Efficiency >= 0.75, Efficiency < 1) %>%
    slice(round(seq(1, n(), length.out = n_table))) %>%
    transmute(
      v = v2, k = k2, lambda = lambda,
      rho = round(Rho, 3),
      E = round(Efficiency, 4)
    )

  ##--- Plot ---##
  eff_plot <- NULL
  if (plot) {
    df_plot <- data.frame(
      Rho = seq(-0.99, 0.99, length.out = 20),
      Efficiency = efficiency_phi(seq(-0.99, 0.99, length.out = 20), v2, k2)
    )

    eff_plot <- ggplot(df_plot, aes(x = Rho, y = Efficiency)) +
      geom_line(color = "blue", linewidth = 1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
      labs(
        title = bquote("Plot of Efficiency Factor (E) against "~rho~" for " ~ v == .(v2)),
        x = "Intra-Block Correlation ("~rho~")",
        y = "Efficiency Factor (E)"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
      )
  }
  # ----------- Compute Efficiency at Input rho -----------
  efficiency_at_rho <- efficiency_phi(rho, v2, k2)
  cat(sprintf("Efficiency Factor at rho = %.3f is E = %.4f\n", rho, efficiency_at_rho))

  ##--- Return Results ---##
  return(list(
    d1 = d1,
    d2 = d2,
    d3 = d3,
    C_mat_trt1 = C_mat_trt1,
    C_mat_trt2 = C_mat_trt2,
    C_mat_trt1_trt2 = C_mat_trt1_trt2,
    eff_plot = eff_plot,
    eff_table = eff_table,
    eff_summary = eff_summary
  ))
}
