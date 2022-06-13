## Cell Annotation Functions

#' Perform automatic SingleR annotation on Seurat scRNA-seq object
#' 
#' @param obj Seurat scRNA-seq object (default assay must be RNA or SCT, already normalized and clustered)
#' @param type A string specifying the type of annotation labels to produce (fine or main)
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
  } else if (type == "main") {
    pred = SingleR::SingleR(test = counts, ref = hpca.se, assay.type.test = 1, labels = hpca.se$label.main)
  } else {
    print("type not found ... available type options: fine, main")
  }
  return(pred)
}

#' Generate DotPlot with informative details to assist with cluster annotation
#'
#' @param obj Seurat scRNA-seq object (default assay must be RNA or SCT, already normalized)
#' @param clusters clusters of interest to provide dotplot (default is 'all')
#' @param type A string specifiying the type of annotation ('general' (default), 't')
#' @return A ggplot object of dotplot with annotations
#' @export
#' 
#' 
dotAnno = function(obj, clusters = "all", type = "general") {
  checkFile = TRUE
  if (type == "t") {
    checkFile = all(file.exists("data/markers/tMarkers.txt"), file.exists("data/markers/tColors.txt"))
    features = read.table("data/markers/tMarkers.txt", sep = "\t")$V1
    groups = factor(read.table("data/markers/tColors.txt", sep = "\t")$V1)
    c = helpColor(groups)
  } else if (type == "all") {
    checkFile = all(file.exists("data/markers/generalMarkers.txt"), file.exists("data/markers/generalColors.txt"))
    features = read.table("data/markers/generalMarkers.txt", sep = "\t")$V1
    groups = factor(read.table("data/markers/generalColors.txt", sep = "\t")$V1)
  }

  if (checkFile) {
    plot = Seurat::DotPlot(obj, features = features, group.by = "seurat_clusters")
    plot = plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.25, colour = c))
    return(plot)
  } else {
    print("data directory not found, make sure data/markers... is in your working directory")
  }
}

#' Helper Function to Generate Color Vector for DotPlot
#' 
helpColor = function(groups) {
  cols = RColorBrewer::brewer.pal(length(levels(groups)),"Set1")
  names(cols) = levels(groups)
  c = rep("a",length(groups))
  for (i in 1:length(groups)) {
    c[i] = cols[groups[i]]
  }
  return(c)
}
