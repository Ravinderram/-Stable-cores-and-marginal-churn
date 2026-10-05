# Rev 20, step 9b: leave-one-state-out for the population-averaged within-state GLM (state + year fixed effects,
# station-clustered SE), comparable with the full-sample GLM estimate; states with >= 20 stations in the sample.
suppressMessages({library(sandwich)})
B <- Sys.getenv("WQ_ROOT", unset = getwd()); O <- paste0(B, "/results/")   # run from the repository root
C <- read.csv(paste0(B, "/analysis/context_r20.csv")); C$fYear <- factor(C$Year)
D <- C[C$spatial_main %in% c(TRUE, "True") & !is.na(C$lpop5) & !is.na(C$rain_clim), ]
sdc <- sd(tapply(D$rain_clim, D$station_uid, mean)); D$clim_s <- D$rain_clim / sdc
out <- list()
for (ind in c("ADV_BOD", "ADV_DO", "ADV_FC_cens", "ADV_PH", "ADV_MULTI_GE2_cens")) {
  X <- D[!is.na(D[[ind]]), ]; X$y <- X[[ind]]
  full <- glm(y ~ fYear + factor(state) + lpop5 + clim_s, X, family = binomial)
  big <- names(which(tapply(X$station_uid, X$state, function(s) length(unique(s))) >= 20))
  for (s in c("none", big)) {
    Y <- if (s == "none") X else X[X$state != s, ]
    g <- suppressWarnings(glm(y ~ fYear + factor(state) + lpop5 + clim_s, Y, family = binomial))
    se <- sqrt(vcovCL(g, cluster = ~station_uid)["lpop5", "lpop5"])
    out[[length(out) + 1]] <- data.frame(indicator = ind, left_out = s, OR = exp(coef(g)["lpop5"]),
      lo = exp(coef(g)["lpop5"] - 1.96 * se), hi = exp(coef(g)["lpop5"] + 1.96 * se), stations = length(unique(Y$station_uid)))
  }
}
R <- do.call(rbind, out); write.csv(R, paste0(O, "C2_loso_glm.csv"), row.names = FALSE)
print(aggregate(OR ~ indicator, R[R$left_out != "none", ], function(v) c(min = min(v), max = max(v))))
print(R[R$left_out == "none", ])
