## Copyright (C) 2022 Kevin R. Coombes, RB McGee, and Jake Reed


setClass("Feature",
         slots = c(name = "character",
                   values = "numeric",
                   colRamp = "function",
                   meaning = "character" # or data.frame?
                   ))

setClass("Cycle",
         slots = c(index = "matrix",
                   dimension = "numeric",
                   color = "character"))

setClass("LoopCircos",
         contains = "Cycle",
         slots = c(angles = "matrix",
                   colors = "list")
         )

