
# this function draws the annotation overlay for kinds of cells
overlay <- function(lwd=2) {
  opar <- par(lwd=lwd)
  on.exit( par(opar) )
  text(0.8, 0.035, "CD34+ CD38-")
  draw.ellipse(0.53, 0.17, 0.25, 0.09, -17, lwd=lwd)

  text(0.8, 0.48, "CD34+")
  polyvert <- ucoords[c(3, 291, 238, 443, 104, 329, 479, 449),] +
    0.02*matrix(c(-1, 0, 0, 1, 1,  0,  0,  -1,
                  1, 1, 1, 1, 0, -1, -1, -1), ncol=2)
  polygon(polyvert)

  text(-0.95, -0.9,  "* Promyelocytes", adj=0)
  text(-0.14, 0.45, "*")
  draw.ellipse(-0.095, 0.45, 0.03, 0.03, lwd=lwd)

  text(0.05, 0, "Platelets", adj=0.5)
  draw.ellipse(0.04, 0.14, 0.04, 0.13, 20, lwd=lwd)

  text(-0.53, 0.85, "Mature", adj=0.5)
  text(-0.53, 0.80, "granulocytes", adj=0.5)
  draw.ellipse(-0.18, 0.56, 0.045, 0.045, lwd=lwd)
  arrows(-0.53, 0.77, -0.23, 0.57, length=0.15, angle=20)

  text(-0.05, 0.95, "Early", adj=0.5)
  text(-0.05, 0.90, "granulocytes", adj=0.5)
  polyvert <- ucoords[c(197, 271, 221, 188, 48, 428, 170, 257, 314, 37),] +
    0.02*matrix(c(-1, -1, -1, -1, -1, 0, 1, 1, 1, 1,
                  -1, 0, 0, 0, 1, 1, 1, 0, 0, -1), ncol=2)
  polygon(polyvert)
  arrows(-0.05, 0.87, -0.05, 0.80, length=0.15, angle=20)

  text(0.55, -0.95, "Plasma cells", adj=0)
  draw.ellipse(0.4424, -0.95, 0.1, 0.1, lwd=lwd)

  text(0.55, -0.76, "Early B cells", adj=0)
  draw.ellipse(0.45, -0.76, 0.11, 0.085, 135, lwd=lwd)

  text(0.52, -0.56, "Mature B cells", adj=0)
  draw.ellipse(0.35, -0.58, 0.20, 0.11, 35, lwd=lwd)

  text(-0.72, -0.15, "Mature CD14+", adj=1)
  text(-0.72, -0.2, "Monocytes", adj=1)
  draw.ellipse(-0.4, -0.10, 0.32, 0.15, 7, lwd=lwd)

  text(-0.3, -0.45, "pDCs", adj=1)
  draw.ellipse(-0.21, -0.45, 0.12, 0.05, 32, lwd=lwd)

  text(-0.9, 0.77, "NK cells", adj=0.5)
  draw.ellipse(-0.84, 0.59, 0.22, 0.14, 160, lwd=lwd)

  text(-0.7, 0.2, "T cells", adj=0.5)
  polyvert <- ucoords[c(11, 418, 355, 409, 460, 340, 234, 78, 453),] +
    0.03*matrix(c(1, 1, 0, 0, -1, -1, -1, -1, 1,
                  1, -1, -1, -1, -1, 1, 0, 1, 1), ncol=2)
  polygon(polyvert)

  text(-0.35,-0.8, "Early monocytes", adj=0.5)
  polyvert <- ucoords[c(311, 349, 413, 184, 445, 266, 6, 416, 303, 96,
                        165, 478, 467, 356),] +
    0.02*matrix(c(-2, -1, -1, -1, -1, -1, -1, 1, 2, 2,  1,  1,  1, 0,
                  0,  1,  1,  0,  0,  0,  1, 1.5, 0, 0, -1, -1, -1, -1), ncol=2)
  polygon(polyvert)

  text(0.1, 0.85, "Erythroid cells", adj=0)
  polyvert <- ucoords[c(258, 415, 92, 403, 395, 53, 286, 354, 204),] +
    0.02*matrix(c( 0, -1, -1, -1, -1, -1, 1, 1, 1,
                  -1,  0,  0,  0,  0,  1, 1, 0, -1), ncol=2)
  polygon(polyvert)

  text(0.2, 0.06, "Basophils")
  polyvert <- ucoords[c(57, 334, 160, 160, 335, 339, 281, 307, 112),] +
    0.02*matrix(c(-1, -1, -1,  2, -1, 1, 1, 1, 1,
                  1,  0, -1, -1,  0, 0, 0, 0, 1), ncol=2)
  polygon(polyvert)
}

