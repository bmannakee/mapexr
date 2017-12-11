# Utility functions. Not exported.


# A utility function to make pulling a list from a data_frame column easy
.pull <- function(x,y) {x[,if(is.name(substitute(y))) deparse(substitute(y)) else y, drop = FALSE][[1]]}
