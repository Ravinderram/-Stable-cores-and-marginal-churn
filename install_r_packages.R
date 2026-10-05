# R packages for steps 06 and 09 (last run with R 4.3.3: lme4 1.1-35, blme 1.0-5, survival 3.5-8,
# sandwich 3.1-0, lmtest 0.9-40).
pkgs <- c("lme4", "blme", "survival", "sandwich", "lmtest")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
for (p in pkgs) cat(p, as.character(packageVersion(p)), "\n")
