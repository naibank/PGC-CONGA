# Shared validation; sourced relative to the executable, not the job directory.
parse_args <- function(args, required) {
  if (!length(args) || length(args) %% 2L ||
      any(!grepl("^--[^-]", args[seq.int(1L, length(args), 2L)])) ||
      any(grepl("^--", args[seq.int(2L, length(args), 2L)])))
    stop("Arguments must be --name value pairs")
  keys <- sub("^--", "", args[seq.int(1L, length(args), 2L)])
  if (anyDuplicated(keys)) stop("Duplicate argument names")
  params <- as.list(setNames(args[seq.int(2L, length(args), 2L)], keys))
  if (length(setdiff(required, keys))) stop("Missing arguments: ", paste(setdiff(required, keys), collapse=", "))
  params
}
load_genesets <- function(path) {
  env <- new.env(parent=emptyenv())
  objects <- load(path, envir=env)
  if (length(objects) != 1L) stop("Expected exactly one object in geneset RData")
  gs <- env[[objects]]
  if (basename(path) == "human_gsData_goincludingIEA100-1000_pathways15-500_20211125.RData") {
    gs <- gs$gs2gene
    gs <- gs[lengths(gs) >= 50L & lengths(gs) <= 500L]
    names(gs) <- gsub("[:.-]", "_", names(gs))
  }
  names(gs) <- gsub("[ -]", ".", names(gs))
  if (!is.list(gs) || is.null(names(gs)) || anyNA(names(gs)) ||
      any(!nzchar(names(gs))) || anyDuplicated(names(gs)))
    stop("Genesets must be a uniquely named list after name normalization")
  gs
}
subset_genesets <- function(gs, from, end) {
  bounds <- suppressWarnings(as.numeric(c(from, end)))
  if (length(bounds) != 2L || anyNA(bounds) || any(!is.finite(bounds)) ||
      any(bounds != floor(bounds)) || bounds[1] < 1 || bounds[2] < bounds[1] || bounds[2] > length(gs))
    stop("Invalid geneset range")
  gs[seq.int(bounds[1], bounds[2])]
}
