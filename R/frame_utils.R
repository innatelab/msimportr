# Utilities for working with data frames
#
# Author: Alexey Stukalov
###############################################################################

# Expands a frame by the column that contains sep-separated list of entities.
# In the new frame there is a row for each entity.
#' @export
expand_collapsed <- function(df, collapsed_col, separated_col,
                             extra_cols=NULL, sep=fixed(";")) {
  exp_list <- str_split(df[[collapsed_col]], sep)
  exp_lengths <- sapply(exp_list, length)
  res <- df[rep.int(1:nrow(df), exp_lengths), c(collapsed_col, extra_cols)]
  res[[separated_col]] <- unlist(exp_list)
  res[,c(separated_col, collapsed_col, extra_cols)]
}
