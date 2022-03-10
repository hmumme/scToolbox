## Automatic Cell Annotation Functions

#' Perform SingleR annotation on Seurat scRNA-seq Object
#' 
#' @param obj Seurat scRNA-seq object (default assay must be RNA, already normalized)
#' @param type A string specifying the type of annotation labels to produce (fine or broad labels)
#' @return A vector of Labels to add into Seurat Object
#' @export 
#' @examples
#' cellTypePredictions = scAnno(obj = seuratObj, type = "fine")
#' seuratObj$cellAnnotation = cellTypePredictions
scAnno <- function(obj, type) {
  hpca.se = celldex::HumanPrimaryCellAtlasData()
  counts = GetAssayData(obj)
  if (type == "fine") {
    pred = SingleR::SingleR(test = counts, ref = hpca.se, assay.type.test = 1, labels = hpca.se$label.fine)
  } else {
    pred = SingleR::SingleR(test = counts, ref = hpca.se, assay.type.test = 1, labels = hpca.se$label)
  }
  return(pred)
}


##