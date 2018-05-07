### Copyright 2018, Kevin R. Coombes and R. B. McGee.
###
### Defines classes describing the design of experiment handled by the Cytangle package.

### Completely provisional; idea is to list the tihings we want to bundle
### so we can later figure out how to represent each one.
setClass("CytangleDesign",
         slots=c(
           experiments = "character", # e.g., "Tube A", "Tube B"
           markers = "character",     # e.g., "CD99". Also need kinds of markers
           groups = "character",      # e.g. "FLT3-ITD"
           samples = "character"      # e.g., "Nl6"
         ))
