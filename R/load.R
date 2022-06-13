#' Load example Seurat datasets for testing
#' 
#' @param ex Seurat example dataset to load, options include: "3k"
#' @param norm A boolean specifying whether to normalize RNA assay and cluster (default FALSE)
#' @return Seurat object of example dataset
#' @export
#' @examples  
#' obj = loadExample("3k")
#' 
loadExample = function(ex, norm = FALSE) {
  if (file.exists("data/pbmc3k/filtered_gene_bc_matrices/hg19/") == "FALSE") {
    print("data directory not found, make sure data/pbmc... is in your working directory")
    return(NULL)
  } else if (ex == "3k") {
    data = Seurat::Read10X(data.dir = "data/pbmc3k/filtered_gene_bc_matrices/hg19/")
  } else {
    print("dataset does not exist, available: 3k")
    return(NULL)
  } 
    obj = Seurat::CreateSeuratObject(counts = data, project = ex, min.cells = 3, min.features = 200)
    if (norm) {
      obj <- Seurat::NormalizeData(obj)
      obj <- Seurat::FindVariableFeatures(obj)
      obj <- Seurat::ScaleData(obj)
      obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj))
      obj <- Seurat::FindNeighbors(obj, dims = 1:10)
      obj <- Seurat::FindClusters(obj, resolution = 0.5)
      obj <- Seurat::RunUMAP(obj, dims = 1:10)
    } else {
      return(obj)
    }
}