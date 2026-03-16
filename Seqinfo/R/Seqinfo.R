setClass("Seqinfo",
         slots = c(
           seqnames = "character",
           seqlengths = "numeric",
           isCircular = "logical",
           genome = "character"
         ),
         prototype = list(
           seqnames = character(),
           seqlengths = numeric(),
           isCircular = logical(),
           genome = character()
         )
)

Seqinfo <- function(seqnames = character(),
                    seqlengths = rep(NA_real_, length(seqnames)),
                    isCircular = rep(NA, length(seqnames)),
                    genome = rep(NA_character_, length(seqnames))) {
  new("Seqinfo",
      seqnames = seqnames,
      seqlengths = seqlengths,
      isCircular = isCircular,
      genome = genome)
}
