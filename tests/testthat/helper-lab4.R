pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
r_files <- list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE)
for (f in r_files) sys.source(f, envir = parent.env(environment()))
