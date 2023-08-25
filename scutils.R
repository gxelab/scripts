#' Utils for analyzing single cell RNA-Seq dataset


#' read output_filtered.h5 produced by CellBender
#' return a sparseMatrix of counts (rows: features; columns, barcodes)
read_cellbender_h5 <- function(path_h5){
    cb_h5 <- H5File$new(path_h5, mode = 'r')
    # when creating sparseMatrix, parameter `i` is 1-based. That's why the `+1`.
    # For a sparse matrix, index slot (`m@i` or `h5[["matrix/indices"]]`) is 0-based.
    m <- Matrix::sparseMatrix(
        i = cb_h5[["matrix/indices"]][] + 1,
        p = cb_h5[["matrix/indptr"]][],
        x = cb_h5[["matrix/data"]][],
        dims = cb_h5[["matrix/shape"]][][],
        dimnames = list(
            cb_h5[["matrix/features/name"]][],
            cb_h5[["matrix/barcodes"]][]),
        repr = "C")
    cb_h5$close_all()
    return(m)
}
