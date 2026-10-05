# Rev 20, step 9 (Major Comment 5, Minor 8): context models on the main spatial sample. Associations, not effects.
#  C1 within-state population (within-between decomposition), with and without the 2023 sampling-frequency attribute;
#     conditional ORs, population-averaged ORs (Zeger et al. 1988 approximation) and average marginal effects (AME)
#     of a tenfold density increase on the probability scale, integrating over the station random effect;
#  C2 population-averaged within-state ORs from logistic GLMs with state and year fixed effects and station-clustered
#     SEs, unweighted and with inverse-probability-of-location weights; leave-one-state-out; state-specific slopes;
#  C3 persistent class, conditional logit stratified by state;
#  C4 rainfall (SI): year-to-year exposure (per SD), climate (per SD), crossed station and rain-cell intercepts,
#     flagged cells excluded as sensitivity; lag x rainfall interaction;
#  C5 basin variance decomposition; C6 upstream population and discharge if the local-kit output is present.
suppressMessages({library(lme4); library(sandwich); library(lmtest); library(survival)})
B <- Sys.getenv("WQ_ROOT", unset = getwd()); O <- paste0(B, "/results/")   # run from the repository root
C <- read.csv(paste0(B, "/analysis/context_r20.csv")); S <- read.csv(paste0(B, "/analysis/stations_r20.csv"))
C$fYear <- factor(C$Year); C$freq <- factor(C$freq, levels = c("monthly", "quarterly", "yearly", "not reported"))
D <- C[C$spatial_main %in% c(TRUE, "True") & !is.na(C$lpop5) & !is.na(C$rain_clim), ]
sdc <- sd(tapply(D$rain_clim, D$station_uid, mean)); sde <- sd(D$rain_expo, na.rm = TRUE)
D$clim_s <- D$rain_clim / sdc; D$expo_s <- D$rain_expo / sde
st <- unique(D[, c("station_uid", "state", "lpop5", "clim_s")])
sm <- aggregate(cbind(lpop5, clim_s) ~ state, st, mean); names(sm)[2:3] <- c("pop_sm", "clim_sm")
D <- merge(D, sm, by = "state"); D$pop_w <- D$lpop5 - D$pop_sm; D$clim_w <- D$clim_s - D$clim_sm
D$basin_sub <- as.character(D$basin_sub); D$basin_sub[is.na(D$basin_sub)] <- paste0("nb_", D$station_uid[is.na(D$basin_sub)])
args <- commandArgs(trailingOnly = TRUE)
IND <- if (length(args)) args else c("ADV_BOD", "ADV_DO", "ADV_FC_cens", "ADV_PH", "ADV_MULTI_GE2_cens")
ci <- function(b, se) exp(b + c(0, -1.96, 1.96) * se)
out <- list(); k <- 0
add <- function(...) { k <<- k + 1; out[[k]] <<- data.frame(...) }
cat("spatial sample stations", length(unique(D$station_uid)), "station-years", nrow(D), "\n")
cat("SD rain exposure", sde, "-> x", exp(sde), "; SD climate", sdc, "-> x", exp(sdc), "\n")
for (ind in IND) {
  X <- D[!is.na(D[[ind]]), ]; X$y <- X[[ind]]
  f0 <- y ~ fYear + pop_w + pop_sm + clim_w + clim_sm + (1 | state) + (1 | basin_sub) + (1 | station_uid)
  m <- glmer(f0, X, family = binomial, control = glmerControl(optimizer = "bobyqa"))
  cf <- summary(m)$coefficients; vc <- as.data.frame(VarCorr(m)); s2 <- vc$vcov[vc$grp == "station_uid"]; s2all <- sum(vc$vcov)
  r <- ci(cf["pop_w", 1], cf["pop_w", 2])
  # AME of +1 log10 density, integrating over the station random effect (fixed effects + state and basin BLUPs kept)
  eta <- predict(m, re.form = ~ (1 | state) + (1 | basin_sub)); set.seed(1); u <- rnorm(200, 0, sqrt(s2))
  p0 <- sapply(u, function(z) plogis(eta + z)); p1 <- sapply(u, function(z) plogis(eta + cf["pop_w", 1] + z))
  ame <- mean(p1 - p0); pbar <- mean(p0)
  b_pa <- cf["pop_w", 1] / sqrt(1 + 0.346 * s2all)
  add(indicator = ind, model = "C1 within-state, GLMM (rev 19 specification)", term = "log10 density (within state)",
      OR = r[1], lo = r[2], hi = r[3], n = nrow(X), stations = length(unique(X$station_uid)),
      PA_OR_approx = exp(b_pa), AME_pp = 100 * ame, baseline_prob = pbar, station_var = s2)
  r <- ci(cf["clim_w", 1], cf["clim_w", 2])
  add(indicator = ind, model = "C1 within-state, GLMM (rev 19 specification)", term = "rain climate per SD (within state)",
      OR = r[1], lo = r[2], hi = r[3], n = nrow(X), stations = length(unique(X$station_uid)), PA_OR_approx = NA, AME_pp = NA, baseline_prob = NA, station_var = s2)
  mf <- glmer(update(f0, . ~ . + freq), X, family = binomial, control = glmerControl(optimizer = "bobyqa"))
  cf2 <- summary(mf)$coefficients; r <- ci(cf2["pop_w", 1], cf2["pop_w", 2])
  add(indicator = ind, model = "C1 + sampling frequency 2023", term = "log10 density (within state)",
      OR = r[1], lo = r[2], hi = r[3], n = nrow(X), stations = length(unique(X$station_uid)), PA_OR_approx = NA, AME_pp = NA, baseline_prob = NA, station_var = NA)
  # C2 population-averaged GLM, state + year FE, station-clustered SE
  g <- glm(y ~ fYear + factor(state) + lpop5 + clim_s, X, family = binomial)
  v <- vcovCL(g, cluster = ~station_uid); r <- ci(coef(g)["lpop5"], sqrt(v["lpop5", "lpop5"]))
  pr0 <- mean(predict(g, type = "response")); Xn <- X; Xn$lpop5 <- Xn$lpop5 + 1; pr1 <- mean(predict(g, Xn, type = "response"))
  add(indicator = ind, model = "C2 GLM state+year FE (population-averaged), unweighted", term = "log10 density",
      OR = r[1], lo = r[2], hi = r[3], n = nrow(X), stations = length(unique(X$station_uid)), PA_OR_approx = NA, AME_pp = 100 * (pr1 - pr0), baseline_prob = pr0, station_var = NA)
  if (!exists("W")) {
    # selection model for being in the main spatial sample, all Panel A_long stations
    Sx <- S[S$n_obs >= 2, ]; Sx$sel <- Sx$spatial_main %in% c(TRUE, "True")
    Sx$lbod <- log(pmax(Sx$bod_med, 0.1)); Sx$lbod[is.na(Sx$lbod)] <- median(Sx$lbod, na.rm = TRUE)
    Sx$lcond_na <- is.na(Sx$cond_med); Sx$lcond <- log(pmax(Sx$cond_med, 1)); Sx$lcond[Sx$lcond_na] <- median(Sx$lcond, na.rm = TRUE)
    Sx$cohort <- cut(Sx$first, c(2015, 2016, 2017, 2020, 2024))
    sel <- glm(sel ~ factor(state) + cohort + n_obs + lbod + lcond + lcond_na, Sx, family = binomial)
    Sx$ps <- fitted(sel); w <- mean(Sx$sel) / Sx$ps; w <- pmin(pmax(w, quantile(w[Sx$sel], .01)), quantile(w[Sx$sel], .99))
    W <- setNames(w, Sx$station_uid)
    write.csv(data.frame(station_uid = Sx$station_uid, selected = Sx$sel, ps = Sx$ps, w = w), paste0(O, "C0_selection_weights.csv"), row.names = FALSE)
  }
  X$w <- W[X$station_uid]
  gw <- suppressWarnings(glm(y ~ fYear + factor(state) + lpop5 + clim_s, X, family = binomial, weights = w))
  v <- vcovCL(gw, cluster = ~station_uid); r <- ci(coef(gw)["lpop5"], sqrt(v["lpop5", "lpop5"]))
  add(indicator = ind, model = "C2 GLM state+year FE, inverse-probability-of-location weights", term = "log10 density",
      OR = r[1], lo = r[2], hi = r[3], n = nrow(X), stations = length(unique(X$station_uid)), PA_OR_approx = NA, AME_pp = NA, baseline_prob = NA, station_var = NA)
  gf <- glm(y ~ fYear + factor(state) + lpop5 + clim_s + freq, X, family = binomial)
  v <- vcovCL(gf, cluster = ~station_uid); r <- ci(coef(gf)["lpop5"], sqrt(v["lpop5", "lpop5"]))
  add(indicator = ind, model = "C2 GLM state+year FE + sampling frequency 2023", term = "log10 density",
      OR = r[1], lo = r[2], hi = r[3], n = nrow(X), stations = length(unique(X$station_uid)), PA_OR_approx = NA, AME_pp = NA, baseline_prob = NA, station_var = NA)
  # leave one state out (states with >= 20 stations), GLMM within specification, nAGQ = 0
  big <- names(which(tapply(X$station_uid, X$state, function(s) length(unique(s))) >= 20))
  for (s in big) {
    Y <- X[X$state != s, ]
    ml <- tryCatch(glmer(f0, Y, family = binomial, nAGQ = 0), error = function(e) NULL)
    if (!is.null(ml)) { c3 <- summary(ml)$coefficients; r <- ci(c3["pop_w", 1], c3["pop_w", 2])
      add(indicator = ind, model = paste0("LOSO without ", s), term = "log10 density (within state)", OR = r[1], lo = r[2], hi = r[3],
          n = nrow(Y), stations = length(unique(Y$station_uid)), PA_OR_approx = NA, AME_pp = NA, baseline_prob = NA, station_var = NA) }
  }
  # state-specific slopes (population-averaged GLM within each state, >= 15 stations, outcome varies)
  for (s in unique(X$state)) {
    Y <- X[X$state == s, ]; ns <- length(unique(Y$station_uid))
    if (ns < 15 || length(unique(Y$y)) < 2 || mean(Y$y) < 0.03 || mean(Y$y) > 0.97) next
    gs <- tryCatch(suppressWarnings(glm(y ~ fYear + lpop5, Y, family = binomial)), error = function(e) NULL)
    if (is.null(gs) || is.na(coef(gs)["lpop5"])) next
    v <- vcovCL(gs, cluster = ~station_uid); se <- sqrt(v["lpop5", "lpop5"]); r <- ci(coef(gs)["lpop5"], se)
    add(indicator = ind, model = paste0("state slope: ", s), term = "log10 density", OR = r[1], lo = r[2], hi = r[3], n = nrow(Y), stations = ns,
        PA_OR_approx = NA, AME_pp = NA, baseline_prob = mean(Y$y), station_var = se)
  }
  # C3 persistent class, conditional logit by state (stations >= 4 years)
  Z <- aggregate(y ~ station_uid + state + lpop5 + clim_s, X, function(v) c(n = length(v), k = sum(v)))
  Z <- do.call(data.frame, Z); names(Z)[5:6] <- c("n", "k")
  # persistent per rev 19 rule needs the sequence; recompute from X
  pers <- tapply(seq_len(nrow(X)), X$station_uid, function(i) { o <- order(X$Year[i]); v <- X$y[i][o]; if (length(v) < 4) return(NA)
    rl <- rle(v); run <- max(c(0, rl$lengths[rl$values == 1])); as.numeric(mean(v) >= 0.6 & run >= 3) })
  Z$persistent <- pers[Z$station_uid]; Z <- Z[!is.na(Z$persistent), ]
  cl <- tryCatch(clogit(persistent ~ lpop5 + clim_s + strata(state), data = Z), error = function(e) NULL)
  if (!is.null(cl)) { c4 <- summary(cl)$coefficients; r <- ci(c4["lpop5", 1], c4["lpop5", 3])
    add(indicator = ind, model = "C3 persistent class, conditional logit by state", term = "log10 density", OR = r[1], lo = r[2], hi = r[3],
        n = nrow(Z), stations = nrow(Z), PA_OR_approx = NA, AME_pp = NA, baseline_prob = mean(Z$persistent), station_var = NA) }
  # C4 rainfall (SI)
  Xr <- X[!is.na(X$expo_s), ]
  for (spec in c("all cells", "flagged cells excluded")) {
    Y <- if (spec == "all cells") Xr else Xr[!(Xr$rain_flag %in% c(TRUE, "True")), ]
    mr <- tryCatch(glmer(y ~ fYear + expo_s + clim_s + (1 | station_uid) + (1 | cell) + (1 | state), Y, family = binomial,
                          control = glmerControl(optimizer = "bobyqa")), error = function(e) NULL)
    if (!is.null(mr)) { c5 <- summary(mr)$coefficients; r <- ci(c5["expo_s", 1], c5["expo_s", 2])
      add(indicator = ind, model = paste0("C4 rainfall, year-to-year exposure per SD, ", spec), term = "rain exposure per SD",
          OR = r[1], lo = r[2], hi = r[3], n = nrow(Y), stations = length(unique(Y$station_uid)), PA_OR_approx = NA, AME_pp = NA, baseline_prob = NA, station_var = NA) }
  }
  # C5 basin decomposition
  mb <- tryCatch(glmer(y ~ fYear + (1 | state) + (1 | basin_main) + (1 | basin_sub) + (1 | station_uid), X, family = binomial,
                       control = glmerControl(optimizer = "bobyqa")), error = function(e) NULL)
  if (!is.null(mb)) { vc <- as.data.frame(VarCorr(mb)); tot <- sum(vc$vcov)
    for (g in vc$grp) add(indicator = ind, model = "C5 variance share", term = g, OR = vc$vcov[vc$grp == g] / tot, lo = NA, hi = NA,
                          n = nrow(X), stations = length(unique(X$station_uid)), PA_OR_approx = NA, AME_pp = NA, baseline_prob = NA, station_var = NA) }
  # C6 upstream population and discharge (if available)
  if ("pop_upstream" %in% names(X)) {
    U <- X[!is.na(X$pop_upstream) & X$pop_upstream > 0 & !is.na(X$DIS_AV_CMS) & X$DIS_AV_CMS > 0, ]
    U$lup <- log10(U$pop_upstream); U$lq <- log10(U$DIS_AV_CMS); U$lppq <- log10(U$pop_upstream / U$DIS_AV_CMS)
    for (spec in list(c("lpop5"), c("lpop5", "lq"), c("lup", "lq"), c("lppq"), c("lpop5", "lppq"))) {
      fm <- as.formula(paste("y ~ fYear + factor(state) + clim_s +", paste(spec, collapse = " + ")))
      gu <- glm(fm, U, family = binomial); v <- vcovCL(gu, cluster = ~station_uid)
      for (tm in spec) { r <- ci(coef(gu)[tm], sqrt(v[tm, tm]))
        add(indicator = ind, model = paste0("C6 upstream GLM state+year FE: ", paste(spec, collapse = " + ")), term = tm,
            OR = r[1], lo = r[2], hi = r[3], n = nrow(U), stations = length(unique(U$station_uid)), PA_OR_approx = NA, AME_pp = NA, baseline_prob = NA, station_var = NA) }
    }
  }
  cat(ind, "done\n")
}
R <- do.call(rbind, out); write.csv(R, paste0(O, if (length(args) == 1) paste0("C1_part_", args[1], ".csv") else "C1_context_models.csv"), row.names = FALSE)
print(R[!grepl("LOSO|state slope", R$model), c("indicator", "model", "term", "OR", "lo", "hi", "PA_OR_approx", "AME_pp", "baseline_prob")])
