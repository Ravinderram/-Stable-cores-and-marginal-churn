# Rev 20, step 6 (Major Comment 3a/3b, Minor 7): mixed models with monitoring-regime covariates.
# M1  y ~ year + (1|station)                                        -> latent ICC
# M1f y ~ year + frequency_2023 + (1|station)                       -> ICC after the 2023 sampling-frequency attribute
# M2  y ~ lag + y0 + year + (1|station)                             -> short-run dependence (Wooldridge initial condition)
# M2c M2 + entry cohort (first observed year: 2016, 2017, 2018-2020, 2021-2024) + frequency_2023
# Gaps: the lag is defined only when the previous calendar year was observed; rows after a gap are dropped from M2
# (counts reported). y0 is the status in the station's first observed year, whatever the calendar year; entry cohort
# lets that initial condition differ by when the station entered.
# Fixed effort (3b): in the linked monthly subset, M1 refitted with status from all samples and from one sample per quarter
# (50 draws), same stations and years.
suppressMessages({library(lme4); library(blme)})
B <- Sys.getenv("WQ_ROOT", unset = getwd()); O <- paste0(B, "/results/")   # run from the repository root
D <- read.csv(paste0(B, "/analysis/analysis_r20.csv"))
D <- D[D$in_panel_A_long == "True" | D$in_panel_A_long == TRUE, ]
fr <- tapply(D$reported_frequency, D$station_uid, function(v) { v <- v[!is.na(v) & v != ""]; if (length(v)) names(sort(table(v), decreasing = TRUE))[1] else "not reported" })
D$freq <- factor(fr[D$station_uid], levels = c("monthly", "quarterly", "yearly", "not reported"))
first <- tapply(D$Year, D$station_uid, min)
D$cohort <- cut(first[D$station_uid], c(2015, 2016, 2017, 2020, 2024), labels = c("2016", "2017", "2018-2020", "2021-2024"))
D$fYear <- factor(D$Year)
icc <- function(m) { v <- as.numeric(VarCorr(m)[[1]]); v / (v + pi^2 / 3) }
fit <- function(f, dat, fc) {
  if (fc) bglmer(f, data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
  else glmer(f, data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
}
grab <- function(m, term) { cf <- summary(m)$coefficients; i <- grep(term, rownames(cf), fixed = TRUE)
  if (!length(i)) return(c(NA, NA, NA)); i <- i[1]; exp(cf[i, 1] + c(0, -1.96, 1.96) * cf[i, 2]) }
rows <- list()
for (ind in c("ADV_BOD", "ADV_DO", "ADV_PH", "ADV_FC_cens", "ADV_FC")) {
  fc <- grepl("FC", ind)
  X <- D[!is.na(D[[ind]]), c("station_uid", "Year", "fYear", "freq", "cohort", ind)]; names(X)[6] <- "y"
  X <- X[order(X$station_uid, X$Year), ]
  X$lag <- ave(X$y, X$station_uid, FUN = function(v) c(NA, head(v, -1)))
  X$lagyear <- ave(X$Year, X$station_uid, FUN = function(v) c(NA, head(v, -1)))
  X$y0 <- ave(X$y, X$station_uid, FUN = function(v) v[1])
  X$isfirst <- ave(X$Year, X$station_uid, FUN = function(v) v == min(v))
  m1 <- fit(y ~ fYear + (1 | station_uid), X, fc)
  m1f <- fit(y ~ fYear + freq + (1 | station_uid), X, fc)
  Y <- X[X$isfirst == 0 & !is.na(X$lagyear) & X$Year - X$lagyear == 1, ]
  n_gap <- sum(X$isfirst == 0 & X$Year - X$lagyear > 1, na.rm = TRUE)
  m2 <- fit(y ~ lag + y0 + fYear + (1 | station_uid), Y, fc)
  m2c <- fit(y ~ lag + y0 + fYear + cohort + freq + (1 | station_uid), Y, fc)
  r <- data.frame(indicator = ind, n_sy = nrow(X), stations = length(unique(X$station_uid)),
    ICC_M1 = icc(m1), ICC_M1_freq = icc(m1f),
    OR_quarterly_vs_monthly = grab(m1f, "freqquarterly")[1], OR_notreported_vs_monthly = grab(m1f, "freqnot reported")[1],
    M2_rows = nrow(Y), rows_after_gap_dropped = n_gap,
    lagOR_M2 = grab(m2, "lag")[1], lagOR_M2_lo = grab(m2, "lag")[2], lagOR_M2_hi = grab(m2, "lag")[3],
    lagOR_M2c = grab(m2c, "lag")[1], lagOR_M2c_lo = grab(m2c, "lag")[2], lagOR_M2c_hi = grab(m2c, "lag")[3],
    ICC_M2 = icc(m2), ICC_M2c = icc(m2c),
    cohort2021_OR = grab(m2c, "cohort2021-2024")[1], cohort2017_OR = grab(m2c, "cohort2017")[1])
  print(r); rows[[ind]] <- r
}
R <- do.call(rbind, rows); write.csv(R, paste0(O, "R3_models_regime.csv"), row.names = FALSE)

# ---- fixed-effort ICC (monthly subset)
Q <- read.csv(paste0(O, "M8_fixed_effort_status.csv"))
out <- list()
for (ind in c("BOD", "PH", "FC")) {
  Z <- Q[Q$indicator == ind, ]; Z$fYear <- factor(Z$year)
  ok <- function(v) length(unique(v)) > 1
  m_all <- tryCatch(glmer(status_all ~ fYear + (1 | station_uid), Z, family = binomial), error = function(e) NULL)
  ic <- c()
  for (k in 0:49) { Z$s <- Z[[paste0("q", k)]]
    m <- tryCatch(suppressWarnings(glmer(s ~ fYear + (1 | station_uid), Z, family = binomial)), error = function(e) NULL)
    if (!is.null(m)) ic <- c(ic, icc(m)) }
  out[[ind]] <- data.frame(indicator = ind, stations = length(unique(Z$station_uid)), station_years = nrow(Z),
    prevalence_all = mean(Z$status_all), prevalence_quarterly = mean(as.matrix(Z[, paste0("q", 0:49)])),
    ICC_all_samples = if (is.null(m_all)) NA else icc(m_all), ICC_quarterly_median = median(ic),
    ICC_quarterly_lo = quantile(ic, .1), ICC_quarterly_hi = quantile(ic, .9), draws = length(ic))
  print(out[[ind]])
}
write.csv(do.call(rbind, out), paste0(O, "R4_fixed_effort_icc.csv"), row.names = FALSE)
