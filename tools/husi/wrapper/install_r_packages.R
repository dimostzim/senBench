# Installs the R packages of the hUSI wrapper into the library of the running R and writes their versions.
# Usage: Rscript --vanilla install_r_packages.R <lock file to write>

.libPaths(.Library)  # no library outside the environment

SNAPSHOT <- "https://packagemanager.posit.co/cran/2023-12-15"
# dplyr and mclust: loaded by sc/hUSI.R, versions as in the hUSI README. jsonlite: writes run_info.json.
PINNED <- c(dplyr = "1.1.4", mclust = "6.0.1", jsonlite = "1.8.8")
CPUS <- 8

lock_path <- commandArgs(trailingOnly = TRUE)
stopifnot(length(lock_path) == 1)

install.packages(names(PINNED), lib = .Library, repos = SNAPSHOT, type = "source", Ncpus = CPUS)

installed <- as.data.frame(installed.packages(lib.loc = .Library), stringsAsFactors = FALSE)
stopifnot(all(installed[names(PINNED), "Version"] == PINNED))

lock <- data.frame(
    package = installed$Package,
    version = installed$Version,
    source = ifelse(is.na(installed$Priority), SNAPSHOT, "r-base")
)
write.table(lock[order(lock$package), ], lock_path, sep = "\t", quote = FALSE, row.names = FALSE)
