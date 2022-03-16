## Cell Annotation Functions

#' Perform automatic SingleR annotation on Seurat scRNA-seq object
#' 
#' @param obj Seurat scRNA-seq object (default assay must be RNA or SCT, already normalized and clustered)
#' @param type A string specifying the type of annotation labels to produce (fine or broad labels)
#' @return A vector of Labels to add into Seurat Object
#' @export 
#' @examples
#' cellTypePredictions = scAnno(obj = seuratObj, type = "fine")
#' seuratObj[["cellAnnotation"]] = cellTypePredictions
scAnno <- function(obj, type) {
  hpca.se = celldex::HumanPrimaryCellAtlasData()
  counts = Seurat::GetAssayData(obj)
  if (type == "fine") {
    pred = SingleR::SingleR(test = counts, ref = hpca.se, assay.type.test = 1, labels = hpca.se$label.fine)
  } else {
    pred = SingleR::SingleR(test = counts, ref = hpca.se, assay.type.test = 1, labels = hpca.se$label.main)
  }
  return(pred)
}

#' Generate DotPlot with informative details to assist with cluster annotation
#'
#' @param obj Seurat scRNA-seq object (default assay must be RNA or SCT, already normalized)
#' @param clusters clusters of interest to provide dotplot (default is all)
#' @param type A string specifiying the type of annotation
#' @return A ggplot object of dotplot with annotations
#' @export
#' 
#' 
dotAnno = function(obj, clusters = "all", type = "all") {
  if (type == "all") {

  }
}
