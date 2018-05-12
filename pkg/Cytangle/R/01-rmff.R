### I tihnk this function shgould have been written for and
### included in Biobase.....
readMIAMEfromFile <- function(filename) {
  temp <- scan(filename, what="c", quiet=TRUE, sep="\n")
  listers <- c("samples", "hybridizations", "normControls",
               "preprocessing", "other")
  singles <- c("name", "lab", "contact", "title", "abstract",
               "url", "pubMedIds")
  args <- list()
  ## parse out the pieces
  for (N in 1:length(temp)) {
    ## split each line into tag and value
    foo <- strsplit(temp[N], ": ")[[1]]
    if (length(foo) < 2) stop("Lines must be in 'TAG: VALUE' format.\n")
    tag <- foo[1]
    value <- paste(foo, collapse = ": ") # in case this string is in the value
    ## make sure it's a legal tag
    if (tag %in% singles) {
      if (tag %in% names(args)) {
        warning("Repeated tag:", tag, "\n")
        value <- paste(args[[tag]], value)
      }
      args[[tag]] <- value
    } else if (tag %in% listers) {
      V <- args[[tag]] # must be a list
      L <- length(V)
      if (L == 0) {
        V <- list(value)
      } else {
        V[L+1] <- value
      }
      args[[tag]] <- V
    } else { # unknown
      stop("Unknown tag:", tag, "\n")
    }
  }
  do.call(MIAME, args)
}
