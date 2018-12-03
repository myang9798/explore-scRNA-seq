#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

rm(list=ls())
library(shiny)
#library(viridis)
#library(scRNAseq)
library(ggplot2)
#library(gridExtra)
load("sce_count_matrix.RData")
load("expr_annotation.RData")
load("tsne_obj.RData")
names(expr_annotation) <- c("breed", "sample")
tsne.obj <- sce.logcounts.tsne.wopca
tsne.obj$sType <- rep(c("A", "B", "F", "I"), times=table(expr_annotation$sample))
gene_intersect <- rownames(exprMat)
#load("sce.RData")
#exprMat <- counts(sce)
#save(exprMat, file="sce_count_matrix.RData")
#expr_annotation <- data.frame(cbind(as.character(sce$batch_name), 
#                                    as.character(sce$breed_name)),
#                              stringsAsFactors = F)
#names(expr_annotation) <- c("sample", "breed")
#batch <- mapvalues(expr_annotation$sample, from=c("A", "B", "F", "I"),
#                   to=c("Mutant, Sample A", "Mutant, Sample B", 
#                     "Wild Type, Sample I", "Wild Type, Sample F"))
#breed <- mapvalues(expr_annotation$breed, from=c("MT", "WT"), 
#                   to=c("Mutant", "Wild Type"))
#expr_annotation <- data.frame(breed=breed, batch=batch, stringsAsFactors = F)
#save(expr_annotation, file="expr_annotation.RData")

#notes on function
#first arguments are tsne object, counts in expression data,
#annotations of cells (e.g. a matrix with sample and breed),
#facet_var=column from expression_annotation
#genes=1 or 2 genes for which to make the plot
#k=cutoff to use for each gene
#output_dir=which directory to output graphs to
sce_scatter <- function(tsne.obj, sam1="A",sam2="A",sam3="A",sam4="A",
                        expression_data=eData.a, 
                        expression_annotation=expression_annotation,
                        facet_var, genes, k=0,
                        output_dir){
  #tsne.obj
  #currently colors are chosen here.
  #one gene:
  #blue
  #two genes:
  #blue for gene1, orange for gene2, and purple for marking both genes
  if(length(genes) == 2){
    gene_cols <- c("#0092F1", "#FF7D27")
  }
  if(length(genes) == 1){
    gene_cols <- c("#0092F1")
  }
  
  pDat <- data.frame(
    tSne1 = tsne.obj@reducedDims$TSNE[, 1][which((tsne.obj$sType==sam1)|(tsne.obj$sType==sam2)|(tsne.obj$sType==sam3)|(tsne.obj$sType==sam4))], 
    tSne2 = tsne.obj@reducedDims$TSNE[, 2][which((tsne.obj$sType==sam1)|(tsne.obj$sType==sam2)|(tsne.obj$sType==sam3)|(tsne.obj$sType==sam4))])
  
  ncell <- dim(pDat)[1]
  
  if(length(genes) == 2){
    gene1 <- genes[1]
    gene2 <- genes[2]
    g1 <- expression_data[gene1, ]
    g2 <- expression_data[gene2, ]
    g_cat <- rep(NA, ncell)
    g_cat[g1 > k & g2 > k] = "Both"
    g_cat[g1 > k & g2 == k] = gene1
    g_cat[g1 == k & g2 > k] = gene2
    g_cat[g1 == k & g2 == k] = "Neither"
    expressed <- g_cat != "Neither"
  }
  
  if(length(genes) == 1){
    gene1 <- genes[1]
    g1 <- expression_data[gene1, ]
    g_cat <- rep(NA, ncell)
    g_cat[g1 > k] <- gene1
    #off_label is the label we use to say that the one gene is not above the
    #threshold
    off_label <- paste("No", gene1, sep=" ") 
    g_cat[g1 <= k] <- off_label
    expressed <- g_cat != off_label
  }
  pDat
  pDat_final <- data.frame(pDat, g_cat=g_cat, 
                           expressed=expressed, 
                           expression_annotation)
  if(length(genes) == 2)
    pDat_final$g_cat <- factor(pDat_final$g_cat, c("Both", gene1, gene2, "Neither"))
  if(length(genes) == 1)
    pDat_final$g_cat <- factor(pDat_final$g_cat, c(gene1, off_label))
  
  pDat_neither <- pDat_final[!pDat_final$expressed, ]  
  pDat_expressed <- pDat_final[pDat_final$expressed, ]  
  
  gglabs <- labs(x = "tSne Component 1", y = "tSne Component 2") 
  gg_settings <-  gglabs +
    theme_bw() +
    theme(
      legend.key = element_blank(),
      legend.background = element_rect(size = 0),
      legend.margin=margin(0,0,0,0, unit="cm"),
      legend.box.margin=margin(0, 0, 0, 0, unit="cm"),
      legend.box.background=element_blank(),
      legend.box.spacing=unit(0, "cm")) 
  if(length(genes) == 2){                
    scale_cols <- c("#6B30B6", gene_cols)
  }
  if(length(genes) == 1){                
    scale_cols <- gene_cols
  }
  
  #not sure how the nuts and bolts work, see ggplot website
  #https://ggplot2.tidyverse.org/reference/vars.html
  quo_facet_var <- enquo(facet_var)
  wrap_by <- function(...){
    facet_wrap(vars(...), nrow = 2)
  }
  
  basic_plot <- ggplot(pDat_neither, aes(x=tSne1, y=tSne2)) +
    geom_point(shape=1, color=c("#C0C0C0")) +
    geom_point(aes(x=tSne1, y=tSne2, color=g_cat), data=pDat_expressed, shape=16, alpha=.8,
               inherit.aes=F) +
    scale_color_manual(values=scale_cols) +
    facet_wrap(vars(!!quo_facet_var), nrow=2) +
    gg_settings +
    guides(color=
             guide_legend(
               title = NULL,
               override.aes=list(size=3.5)))
  gene_output_name <- paste(genes, collapse="_")
  gene_output_name <- paste(c(gene_output_name, ".pdf"), collapse="")
  # output_file <- paste(output_dir, gene_output_name, sep="/")
  
  #pdf(output_file)
  plot(basic_plot)
  dev.off()
  return(basic_plot)
}



# Define UI ----
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      h3("The plot of a single gene or a pair of genes"),
      h5("(If you want the plot for a single gene, please choose the same two genes below. Otherwise, choose two different ones)",style="color:purple"),
      selectizeInput("Gene1",
                     label="Choose a gene",
                     choice=gene_intersect
      ),
      selectizeInput("Gene2", 
                     label="Choose a gene",
                     choice=gene_intersect
      ),
      selectizeInput("Sample",
                     label  = "Choose the Sample(s)",
                     choice = c("All (Sample A, B, F, I)","Mutant (Sample A, B)","Wild Type (Sample F, I)",
                                "Mutant, Sample A","Mutant, Sample B","Wild Type, Sample F","Wild Type, Sample I")
      ),
      selectInput("Way",
                  label = "Choose the way of producing genes",
                  choice=c("tSne")
      )
      
    ),
    mainPanel(
      plotOutput("plot1")
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  # Combine the selected variables into a new data frame
  output$plot1<-renderPlot({
    
    sample.data<-switch(input$Sample,
                        "All (Sample A, B, F, I)"=c("Mutant, Sample A", "Mutant, Sample B", "Wild Type, Sample I", "Wild Type, Sample F"),
                        "Mutant (Sample A, B)" = c("Mutant, Sample A", "Mutant, Sample B"),
                        "Wild Type (Sample F, I)" = c("Wild Type Sample F","Wild Type Sample I"),
                        "Mutant, Sample A"= "Mutant, Sample A",
                        "Mutant, Sample B"= "Mutant, Sample B",
                        "Wild Type, Sample I"= "Wild Type, Sample I",
                        "Wild Type, Sample F"= "Wild Type, Sample F")
    way<-switch(input$Way,
                "tSne"=tsne.obj)
    
    #try subsetting larger data sets rather than building
    #up the data sets from eData.x
    if(input$Sample=="All (Sample A, B, F, I)"){
      expression.data=exprMat
      sam1="A";sam2="B";sam3="F";sam4="I"
      expr_anno=expr_annotation
    }else if(input$Sample=="Mutant (Sample A, B)"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Mutant, Sample A", "Mutant, Sample B")]
      sam1="A";sam2="B";sam3="A";sam4="B"
      expr_anno=subset(expr_annotation,sample=="Mutant, Sample A"|sample=="Mutant, Sample B")
    }else if(input$Sample=="Wild Type (Sample F, I)"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Wild Type, Sample F", 
                                                              "Wild Type, Sample I")]
      sam1="F";sam2="I";sam3="F";sam4="I"
      expr_anno=subset(expr_annotation,sample=="Wild Type, Sample F"|sample=="Wild Type, Sample I")
    }else if(input$Sample=="Mutant, Sample A"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Mutant, Sample A")]
      sam1="A";sam2="A";sam3="A";sam4="A"
      expr_anno=subset(expr_annotation,sample=="Mutant, Sample A")
    }else if(input$Sample=="Mutant, Sample B"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Mutant, Sample B")]
      sam1="B";sam2="B";sam3="B";sam4="B"
      expr_anno=subset(expr_annotation,sample=="Mutant, Sample B")
    }else if(input$Sample=="Wild Type, Sample F"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Wild Type, Sample F")]
      sam1="F";sam2="F";sam3="F";sam4="F"
      expr_anno=subset(expr_annotation,sample=="Wild Type, Sample F")
    }else if(input$Sample=="Wild Type, Sample I"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Wild Type, Sample I")]
      sam1="I";sam2="I";sam3="I";sam4="I"
      expr_anno=subset(expr_annotation,sample=="Wild Type, Sample I")
    }
    
    if(input$Gene1==input$Gene2){
      sce_scatter(tsne.obj = way,sam1=sam1,sam2=sam2,sam3=sam3,sam4=sam4,expression_data=expression.data,
                  expression_annotation=expr_anno, 
                  facet_var=sample, genes=input$Gene1)
    }else{
      sce_scatter(tsne.obj = way,sam1=sam1,sam2=sam2,sam3=sam3,sam4=sam4,expression_data=expression.data,
                  expression_annotation=expr_anno, 
                  facet_var=sample, genes=c(input$Gene1,input$Gene2))
    }
    
  })
  
}



# Run the application 
shinyApp(ui = ui, server = server)