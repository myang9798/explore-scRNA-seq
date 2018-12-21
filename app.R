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
library(rlang)
library(gridExtra)
#library(gridExtra)
load("sce_count_matrix.RData")#name:exprMat
load("expr_annotation.RData")#expre_annotation
load("tsne_obj.RData")
load("norm_count.RData")
load("raw_count.RData")
load("anno_dat.RData")
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
               override.aes=list(size=3.5))) +
    ggpubr::theme_pubclean()
  gene_output_name <- paste(genes, collapse="_")
  gene_output_name <- paste(c(gene_output_name, ".pdf"), collapse="")
  # output_file <- paste(output_dir, gene_output_name, sep="/")
  
  #pdf(output_file)
  #plot(basic_plot)
  #dev.off()
  return(basic_plot)
}
#The function of violin plot
#notes about the function
#-genes should be a two dimensional vector; to get a plot with only
# one gene, the two genes in genes should be the same.
#-as of now, facet_var should be either breed or batch if we apply some
#clustering algorithm later we could add it to anno_dat
sce_violin_plot <- function(norm=T, genes, subset=NULL, facet_var=NULL){
  if(norm){ 
    gene_cdat <- norm_count
  } else {
    gene_cdat <- raw_count
  }
  if(is.null(subset)){
    subset <- rep(T, dim(gene_cdat)[2])
  }
  
  ugenes <- unique(genes)
  if(length(ugenes) == 2){
    gene_cols <- c("#0092F1", "#FF7D27")
  }
  if(length(ugenes) == 1){
    gene_cols <- c("#0092F1")
  }
  
  gene_cdat <- gene_cdat[, subset] 
  gene_cdat <- gene_cdat[genes, ]
  if(genes[1] == genes[2]){
    gene_cdat <- gene_cdat[1, , drop=F]
  }
  genes <- unique(genes)
  genes <- factor(genes, genes)
  anno_subset <- anno_dat[subset, ]
  violin_g <- gene_cdat
  violin_g <- data.frame(GeneName=rownames(violin_g), violin_g)
  rownames(violin_g) <- NULL
  violin_g <- tidyr::gather(violin_g, cell_id, expression, -GeneName)
  violin_g$cell_id <- gsub("\\.", "-", violin_g$cell_id)
  violin_g <- merge(violin_g, anno_subset, by.x="cell_id", by.y="cell_id")
  violin_p <- ggplot(violin_g, aes(x=GeneName, y=expression,fill=GeneName)) + 
    geom_violin()
  violin_p <- violin_p + labs(x=" ",y="Gene Expression",fill=" ")
  if(!is.null(facet_var)){
    violin_p <- violin_p + facet_grid(eval(expr(~ !!ensym(facet_var)))) 
  }
  violin_p <- violin_p  + ggpubr::theme_pubclean()
  violin_p <- violin_p + scale_fill_manual(values=gene_cols)
  return(violin_p)
}


# Define UI ----
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      h3("Exploratory plots for scRNA-seq data"),
      h5("Choose either two different genes or the same gene for data display of a single gene.",style="color:purple"),
      selectizeInput("Gene1",
                     label="Choose a gene",
                     choice=gene_intersect
      ),
      selectizeInput("Gene2", 
                     label="Choose a gene",
                     choice=gene_intersect
      ),
      selectizeInput("Sample",
                     label  = "Choose the samples",
                     choice = c("All Samples: A, B, F, and I","Mutant: Sample A, B","Wild Type: Sample F, I",
                                "Mutant, Sample A","Mutant, Sample B","Wild Type, Sample F","Wild Type, Sample I")
      ),
      selectInput("Way",
                  label = "Choose the dimension reduction method",
                  choice=c("tSne")
      ),
      selectInput("facet",
                  label = "Choose a variable for stratifying the violin plot",
                  choice=c("Batch","Breed")
      ),
      radioButtons("datatype", "Data type for the violin plot",
                   c("Normalized Data", "Raw Data")
      )

    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Scatter Plot", plotOutput("plot1")), 
        tabPanel("Violin Plot", plotOutput("plot2")))
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  # Combine the selected variables into a new data frame
  output$plot1<-renderPlot({
    
    sample.data<-switch(input$Sample,
                        "All Samples: A, B, F, and I"=c("Mutant, Sample A", "Mutant, Sample B", "Wild Type, Sample I", "Wild Type, Sample F"),
                        "Mutant: Sample A, B" = c("Mutant, Sample A", "Mutant, Sample B"),
                        "Wild Type: Sample F, I" = c("Wild Type Sample F","Wild Type Sample I"),
                        "Mutant, Sample A"= "Mutant, Sample A",
                        "Mutant, Sample B"= "Mutant, Sample B",
                        "Wild Type, Sample I"= "Wild Type, Sample I",
                        "Wild Type, Sample F"= "Wild Type, Sample F")
    way<-switch(input$Way,
                "tSne"=tsne.obj)
    
    #try subsetting larger data sets rather than building
    #up the data sets from eData.x
    output$plot1=renderPlot({ 
    if(input$Sample=="All Samples: A, B, F, and I"){
      expression.data=exprMat
      sam1="A";sam2="B";sam3="F";sam4="I"
      expr_anno=expr_annotation
      subset=NULL
    }else if(input$Sample=="Mutant: Sample A, B"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Mutant, Sample A", "Mutant, Sample B")]
      sam1="A";sam2="B";sam3="A";sam4="B"
      expr_anno=subset(expr_annotation,sample=="Mutant, Sample A"|sample=="Mutant, Sample B")
      subset=which(anno_dat$Breed=="MT")
    }else if(input$Sample=="Wild Type: Sample F, I"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Wild Type, Sample F", 
                                                              "Wild Type, Sample I")]
      sam1="F";sam2="I";sam3="F";sam4="I"
      expr_anno=subset(expr_annotation,sample=="Wild Type, Sample F"|sample=="Wild Type, Sample I")
      subset=which(anno_dat$Breed=="WT")
    }else if(input$Sample=="Mutant, Sample A"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Mutant, Sample A")]
      sam1="A";sam2="A";sam3="A";sam4="A"
      expr_anno=subset(expr_annotation,sample=="Mutant, Sample A")
      subset=which(anno_dat$Batch=="A")
    }else if(input$Sample=="Mutant, Sample B"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Mutant, Sample B")]
      sam1="B";sam2="B";sam3="B";sam4="B"
      expr_anno=subset(expr_annotation,sample=="Mutant, Sample B")
      subset=which(anno_dat$Batch=="B")
    }else if(input$Sample=="Wild Type, Sample F"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Wild Type, Sample F")]
      sam1="F";sam2="F";sam3="F";sam4="F"
      expr_anno=subset(expr_annotation,sample=="Wild Type, Sample F")
      subset=which(anno_dat$Batch=="F")
    }else if(input$Sample=="Wild Type, Sample I"){
      expression.data=exprMat[, expr_annotation$sample %in% c("Wild Type, Sample I")]
      sam1="I";sam2="I";sam3="I";sam4="I"
      expr_anno=subset(expr_annotation,sample=="Wild Type, Sample I")
      subset=which(anno_dat$Batch=="I")
    }
    if(input$Gene1==input$Gene2){

     g1=sce_scatter(tsne.obj = way,sam1=sam1,sam2=sam2,sam3=sam3,sam4=sam4,expression_data=expression.data,
                  expression_annotation=expr_anno, 
                  facet_var=sample, genes=input$Gene1)
     g1<-g1+theme(
       text = element_text(size=15,face="bold",color="black"),
       axis.text=element_text(size=12,face="bold"),
       legend.text = element_text(size=12,face = "bold"),
       legend.title = element_text(size=15,face="bold"),
       axis.title.x =element_text( size=15, face="bold"),
       axis.title.y = element_text(size=15, face="bold"))

    }else{
      g1=sce_scatter(tsne.obj = way,sam1=sam1,sam2=sam2,sam3=sam3,sam4=sam4,expression_data=expression.data,
                  expression_annotation=expr_anno, 
                  facet_var=sample, genes=c(input$Gene1,input$Gene2))
      g1<-g1+theme(
        text = element_text(size=15,face="bold",color="black"),
        axis.text=element_text(size=12,face="bold"),
        legend.text = element_text(size=12,face = "bold"),
        legend.title = element_text(size=15,face="bold"),
        axis.title.x =element_text( size=15, face="bold"),
        axis.title.y = element_text(size=15, face="bold"))
    
    }
    grid.arrange(g1)
      })
    output$plot2=renderPlot({ 
      if(input$Sample=="All Samples: A, B, F, and I"){
        subset=NULL
      }else if(input$Sample=="Mutant: Sample A, B"){
        subset=which(anno_dat$Breed=="MT")
      }else if(input$Sample=="Wild Type: Sample F, I"){
        subset=which(anno_dat$Breed=="WT")
      }else if(input$Sample=="Mutant, Sample A"){
        subset=which(anno_dat$Batch=="A")
      }else if(input$Sample=="Mutant, Sample B"){
        subset=which(anno_dat$Batch=="B")
      }else if(input$Sample=="Wild Type, Sample F"){
        subset=which(anno_dat$Batch=="F")
      }else if(input$Sample=="Wild Type, Sample I"){
        subset=which(anno_dat$Batch=="I")
      }
      if(input$datatype=="Normalized Data"){
     g2=sce_violin_plot(genes=c(input$Gene1,input$Gene2), norm=T,subset=subset,
                       facet_var=input$facet)
     g2<-g2+theme(
       text = element_text(size=15,face="bold",color="black"),
       axis.text.x=element_text(size=11,face="bold",color="black"),
       axis.text.y = element_text(size=14,face="bold",color= "black"),
       legend.text = element_text(size=14,face = "bold"),
       legend.title = element_text(size=18,face="bold"),
       axis.title.x =element_text( size=18, face="bold",color="black"),
       axis.title.y = element_text( size=18, face="bold",color="black"))
      }else{
      g2=sce_violin_plot(genes=c(input$Gene1,input$Gene2), norm=F,subset=subset,
                           facet_var=input$facet)
      g2<-g2+theme(
        text = element_text(size=15,face="bold",color="black"),
        axis.text.x=element_text(size=11,face="bold",color="black"),
        axis.text.y = element_text(size=14,face="bold",color= "black"),
        legend.text = element_text(size=14,face = "bold"),
        legend.title = element_text(size=18,face="bold"),
        axis.title.x =element_text( size=18, face="bold",color="black"),
        axis.title.y = element_text( size=18, face="bold",color="black"))
      }
    grid.arrange(g2)
    })
    
  })
  
}



# Run the application 
shinyApp(ui = ui, server = server)