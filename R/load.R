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
  if (ex == "3k") {
    data = readRDS("https://raw.githubusercontent.com/hmumme/scToolbox/main/data/pbmc3k_data.Rda")
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