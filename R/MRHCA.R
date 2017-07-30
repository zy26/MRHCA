MRHCA <- function(cor_M, MR_E, step_size0 = 50, p_sig0 = 0.05, Hit_score_cutoff = 100) {
  cor_c <- cor_M
  diag(cor_c) <- 0
  cor_c <- abs(cor_c)
  cor_r_rank1 <- apply(cor_c, 1, rank)
  cor_r_rank2 <- t(apply(cor_c, 2, rank))
  cor_r_rank1 <- nrow(cor_r_rank1) - cor_r_rank1 + 1
  cor_r_rank2 <- nrow(cor_r_rank2) - cor_r_rank2 + 1
  MR_M <- cor_r_rank1 * cor_r_rank2
  ddd <- MR_E
  null_growth_rate <- c()
  step_size <- step_size0 #parameter
  print("Generating Empirical Null Growth Rate")
  for (i in 1:nrow(ddd)) {
    x0 <- sqrt(sort(ddd[i,]))
    null_growth_rate <- rbind(null_growth_rate, calculate_growth_rate(x0, step = step_size))
  }
  print("Hub Identification")
  hub_info_all <- c()
  sig_N <- ncol(null_growth_rate) #parameter
  p_sig <- p_sig0 #parameter
  hub_sig_cut <- Hit_score_cutoff #parameter
  for (i in 1:nrow(MR_M)) {
    x0 <- sqrt(sort(MR_M[i,]))
    tg_growth_rate <- calculate_growth_rate(x0, step = 50)
    is_sig <- growth_rate_significance(tg_growth_rate, null_growth_rate, sig_N) < p_sig
    S_all <- growth_rate_scoring(is_sig)
    cluster_bound <- 0
    if_hub <- 0
    if ((max(S_all) > hub_sig_cut) & (sum(tg_growth_rate[1:50] > 0.7) == 0)) {
      cluster_bound <- max(which(S_all == max(S_all)))
      if_hub <- 1
    }
    hub_info <- c(if_hub, cluster_bound)
    hub_info_all <- rbind(hub_info_all, hub_info)
  }
  rownames(hub_info_all) <- rownames(MR_M)
  return(list(hub_info_all, MR_M))
}

growth_rate_significance <- function(tg_growth_rate, null_growth_rate, sig_N) {
  growth_rate_sig <- rep(0, sig_N)
  for (i in 1:sig_N) {
    growth_rate_sig[i] <- sum(tg_growth_rate[i] >= null_growth_rate[, i])
  }
  growth_rate_sig <- growth_rate_sig / nrow(null_growth_rate)
  return(growth_rate_sig)
}

calculate_growth_rate <- function(x0, step = 20) {
  gr_all <- rep(0, length(x0))
  for (i in 1:length(x0)) {
    sid <- max(1, i - step)
    eid <- min(length(x0), i + step)
    tg_id <- c(sid:eid)
    l <- eid - sid + 1
    gr_all[i] <- (x0[eid] - x0[sid]) / l
  }
  return(gr_all)
}

growth_rate_scoring <- function(is_sig, up_score = 1, down_score = -2, delete_top = 10) {
  S <- 0
  S_all <- c()
  for (i in 1:length(is_sig)) {
    if (i < delete_top) {
      S <- S + 0
      S_all <- c(S_all, S)
    }
    else {
      if (is_sig[i] == 1) {
        S <- S + up_score
        S_all <- c(S_all, S)
      }
      if (is_sig[i] == 0) {
        S <- S + down_score
        S_all <- c(S_all, S)
      }
    }
  }
  return(S_all)
}

Empirical_null_distribution_MR <- function(data, Rounds = 1000) {
  N <- nrow(data)
  aaa <- runif(N * Rounds, 0, 1)
  bbb1 <- rep(0, length(aaa))
  bbb2 <- rep(0, length(aaa))
  for (i in 1:length(aaa)) {
    bbb1[i] <- rbinom(1, N - 1, aaa[i])
    bbb2[i] <- rbinom(1, N - 1, aaa[i])
  }
  ccc <- (bbb1 + 1) * (bbb2 + 1)
  ddd <- matrix(ccc, Rounds, N)
  return(ddd)
}
