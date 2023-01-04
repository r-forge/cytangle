\name{LoopCircos-class}
\docType{class}
\alias{LoopCircos-class}
\alias{LoopCircos}
\alias{image,LoopCircos-method}
\alias{angleMeans}
\title{
The \code{LoopCircos} Class
}
\description{
  A \code{LoopCircos} object represents a one-dimensional cycle (a
  loop, or topogical eequivalent of a circle), along with a set of
  features that vary around the loop and can be plotted in a "circos"
  plot to help explsin why the loop exists. The \code{LoopCircos} class
  inherits from the \code{\link{Cycle}} class, so you can think of it as
  a cycle object with extra information (namely, the explanatory
  features).
}
\usage{
LoopCircos(cycle, angles, colors)
angleMeans(view, rips, cycle = NULL, dset, angleWidth = 20, incr = 15)
\S4method{image}{LoopCircos}(x, ...)
}
\arguments{
  \item{values}{A numeric vector.}
  \item{name}{A character vector of length one.}
  \item{colors}{A character vector of length at least two, containing
    the names or hexadecimal representations of colors used to create a
    color ramp to display the feature.}
  \item{meaning}{A character vector of length two containing the
    interpretations of the low and high extreme values of the
    feature. Note that this works perfectly well for binary factors
    represented as numeric 0-1 vectors.}
  \item{x}{A \code{LoopCircos} object.}
  \item{y}{A (layout) matrix with coordinates showing where to plot
    each point in the feature.}
  \item{view}{ A (layout) matrix with coordinates showing where to plot
    each point in the feature.}
  \item{pch}{A number; the graphical plotting character parameter \code{pch}.}
  \item{\dots}{The usual set of additional graphical parameters. }
}
\section{Slots}{
  \describe{
    \item{\code{angles}:}{}
    \item{\code{colors}:}{}
  }
}
\section{Methods}{
  \describe{
    \item{image(x, ...):}{
      Produce a circos plot of a \code{LoopCircos} object.
    }
  }
}
\value{
  The \code{LoopCircos} function constructs and returns an object of the
  \code{LoopCircos} class

  The \code{image} method returns (invisibly) the LoopCircos
  object that was its first argument.
}
\seealso{
  \code{\link{Cycle}}
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(CLL)
cyc1 <- Cycle(ripdiag, 1, 236, "forestgreen")
V <- cmdscale(daisydist)
poof <- angleMeans(V, ripdiag, cyc1@index, clinical)
poof[is.na(poof)] <- 0
angle.df <- poof[, c("mutation.status", "CatCD38", "CatB2M", "CatRAI",
                     "Massive.Splenomegaly", "Hypogammaglobulinemia")]
colorScheme <- list(c(M = "green", U = "magenta"),
                    c(Hi = "cyan", Lo ="red"),
                    c(Hi = "blue", Lo = "yellow"),
                    c(Hi = "#dddddd", Lo = "#111111"),
                    c(No = "#dddddd", Yes = "brown"),
                    c(No = "#dddddd", Yes = "purple"))
annote <- LoopCircos(cyc1, angle.df, colorScheme)
image(annote)
}
\keyword{ classes }
\keyword{ graphics }
