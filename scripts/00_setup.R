options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

required_packages <- c("data.table", "ggplot2", "reticulate")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) pak::pkg_install(missing_packages)

if (!requireNamespace("Giotto", quietly = TRUE)) {
  pak::pkg_install("giotto-suite/Giotto")
}

# Xenium parquet files require Arrow with ZSTD support.
has_arrow <- requireNamespace("arrow", quietly = TRUE)
has_zstd <- has_arrow && isTRUE(arrow::arrow_info()$capabilities[["zstd"]])
if (!has_zstd) {
  Sys.setenv(ARROW_WITH_ZSTD = "ON")
  install.packages(
    "arrow",
    repos = c("https://apache.r-universe.dev", getOption("repos")),
    type = "source"
  )
}

library(Giotto)

if (!checkGiottoEnvironment()) {
  message("Installing the default Giotto Python environment. This is needed only once.")
  installGiottoEnvironment()
}

set_giotto_python_path(NULL)

python_packages_ready <- try(
  GiottoUtils::package_check(
    pkg_name = c("tifffile", "imagecodecs", "scrublet"),
    repository = c("pip:tifffile", "pip:imagecodecs", "pip:scrublet")
  ),
  silent = TRUE
)

if (!isTRUE(python_packages_ready)) {
  active_environment <- GiottoUtils::py_active_env()
  reticulate::conda_install(
    envname = active_environment,
    packages = c("tifffile", "imagecodecs", "scrublet"),
    pip = TRUE
  )
}

message("Giotto, Arrow/ZSTD, and optional image/reference dependencies are ready.")
print(sessionInfo())
