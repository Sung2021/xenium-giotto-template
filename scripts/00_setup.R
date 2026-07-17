options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

if (!requireNamespace("Giotto", quietly = TRUE)) {
  pak::pkg_install("giotto-suite/Giotto")
}

library(Giotto)

if (!checkGiottoEnvironment()) {
  message("Installing the default Giotto Python environment. This is needed only once.")
  installGiottoEnvironment()
}

message("Giotto setup is complete.")
print(sessionInfo())
