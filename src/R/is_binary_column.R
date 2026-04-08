is_binary_column <- function(x) {
    vec <- as.vector(na.omit(x))
    check <- all(vec %in% c(0, 1))
    return(check)
}
