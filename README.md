---
title: "User manual for AtRTD2 expression analysis and visualization"
author: "Wenbin Guo"
date: "11 September 2017"
output:
  html_document:
    code_folding: show
    fig_caption: yes
    highlight: textmate
    theme: cerulean
    toc: yes
    toc_depth: 4
    toc_float:
      collapsed: no
      smooth_scroll: yes
  pdf_document:
    toc: yes
    toc_depth: '4'
  word_document:
    toc: yes
    toc_depth: '4'
---


<style>
body {
text-align: justify}
</style>
# **Shiny Dashboard to visualize AtRTD2 expression profiles**
## **Mode 1: Shiny Dashboard--Web App**
We have provide a graphical interface to visualize the gene and transcript expression profiles of RNA-seq data in <u>t</u>ranscripts <u>p</u>er <u>m</u>illion (TPM) for the paper “Rapid and dynamic alternative splicing impacts the Arabidopsis cold response”. The Dashboard can be directly accessed in the Shiny Apps webpage: https://wyguo.shinyapps.io/atrtd2_profile_app/. The research was conducted with BBSRC funding and is a collaboration between the University of Dundee, The James Hutton Institute and the University of Glasgow. The Shiny Dashboard was designed by Wenbin Guo.

Deep RNA-seq data was generated on a diel time-course of 5-week old Arabidopsis Col-0 plants transferred from 20°C to 4°C. Arabidopsis rosettes were harvested every three hours from the last day at 20°C, and the first and fourth day following transfer to 4°C, totaling 26 time-points. Three biological repeats generated over 9.5 billion raw paired-end reads giving around 366 M paired-end reads per time-point. Individual transcripts were quantified in TPM units using Salmon (Patro et al., 2017) and AtRTD2-QUASI as reference transcriptome (Zhang et al., 2017). 

This Dashboard includes three sidebar menus (left-hand side of the webpage):

1. **Introduction:** this page.
2. **Expression profiles:** graphical interface to visualize expression profiles. Once the Dashboard is open, it may take a few seconds to to load the expression data. The interface includes four box panels:
    + _Input a gene_: Users can type in a gene ID in the TAIR/Araport11 format, e.g. AT1G10600 (Lamesch et al., 2012; Cheng et al., 2017; Zhang et al., 2017) (use upper case letters). Then click the button "Click to plot", the plot of gene and transcript expression profiles will be shown in the “Profile plot” panel.
    + _Profile plot:_ The expression profile display box.
    + _Select a type:_ Select an image type (png or pdf) and click “Download” to save the profile plots to a local folder.
    + _Profile in TPM:_ A data table of expression values in TPM across the time-series is given. Users can type a gene ID to check individual gene/transcript information.
3. **Contact us:** contact details.

However, the free service to hold the Dashboard only allows 25 hours of active use per week. In case of not being able to access the service, the same data and Shiny Dashboard are also available for access for R users. The R source code of the Dashboard and a user manual are on Github: https://github.com/wyguo/AtRTD2_profiles.

## **Mode 2: Command lines**

The implementation also can be run using command lines.

### **Install dependency packages**
If any dependency R package is missing in you PC, please install accordingly.

```r
##--install missing packages
required.packages <- c("shiny","ggplot2","reshape2","gtable","shinydashboard",
                       "gtools","gridExtra","leaflet")
required.packages <- required.packages[-which(required.packages %in% 
                                                rownames(installed.packages()))]
if(length(required.packages)>0) install.packages(required.packages, dependencies=TRUE)
```

### **Download source codes and data**
The source codes and transcript TPM expression data are available at: 
https://github.com/wyguo/AtRTD2_profiles/tree/master/Dashboard. 
It includes three files: 

- data.exp.RData
- server.R
- ui.R 

Please save the files to a local folder. Subsequently, setting the R working directory to this folder.

### **Open the Dashboard in command line**
The Shiny Dashboard to visualize expression profiles in **Mode 1** can be opened in command line:

```r
library(shiny)
runApp()
```
We recommand to run the Shiny Dashboard in RStudio. It seems R Shiny App (designed for RStudio) are currently not working properly in console of R version < 3.4.1 (2017-04-21). Higher versions were not tested. 


### **Visualize expression with R code**

The expression profiles also can be visualized in R-code mode.

```r
# load library
library(ggplot2)
library(gridExtra)
# load source function
source('ui.R')
# load transcript level expression
load('data.exp.RData')
head(data.exp)

# make the expression start from the first column and the row names are the transcripts
rownames(data.exp) <- data.exp[,1]
data.exp <- data.exp[,-1]

# make the plot
## provide a gene name
gene='AT1G01020'
## call plot_ArRTD2 function
plots <- plot_AtRTD2(data.exp = data.exp,gene = gene,reps = rep(1:26,each=9),
            y.lab = 'TPM',marker.size = 3,legend.text.size = 11,
            error.type = 'stderr',plot.gene = T,error.bar = T,show.plot = T)
names(plots)
# [1] "profiles" "lightbar" "plots"  

# make the plot
## save plot to png
ggsave(paste0(gene,'.png'),grid.arrange(plots$plots$profiles,plots$plots$lightbar,                                        ncol=1,heights=c(6,1)),width = 10,height = 5.5,units = 'in')

## save plot to pdf
ggsave(paste0(gene,'.pdf'),grid.arrange(plots$plots$profiles,plots$plots$lightbar,
                                        ncol=1,heights=c(6,1)),width = 10,height = 5.5,
       units = 'in')
```

Parameters of function plot_AtRTD2: 

- **data.exp**: transcript level expresssion to make plot, such as TPM and read counts. The columns are the samples and rows are the transcripts names.
- **gene**: a gene name to make the plot. Transcripts expression data matching to the gene name will be used to make the plot.
- **y.lab**: y axis title, defualt is \code{y.lab='TPM'}.
- **legend.text.size**: legend text size.
- **marker.size**: the marker size passed to \code(geom_point(size=marker.size)) in \link{\code{ggplot2}}.
- **error.type**: the method to calculate the errors shown as error bars. Options are "stderr" for standard error (default) and "sd" for standard deviation.
- **plot.gene**: logical, to show the gene expression profile in the plot (TRUE) or not (FALSE).
- **error.bar**: logical, to show error bars (TRUE) or not (FALSE).
- **show.plot**: logical, to return a grahical plot (TRUE) or not (FALSE).
- **return** a list of plots. In the list, "profiles" is the line plot of profiles, "lightbar" is the lightbar and "plots" is a list after aligning the profiles and lightbar plots to consisent x coordinates.



# **Reference**
1. Calixto et al. (submitted)  Rapid and dynamic alternative splicing impacts the Arabidopsis cold response.
2. Cheng,C.Y. et al. (2017) Araport11: a complete reannotation of the Arabidopsis thaliana reference genome. Plant J., 89, 789–804.
3. Lamesch,P. et al. (2012) The Arabidopsis Information Resource (TAIR): Improved gene annotation and new tools. Nucleic Acids Res., 40.
4. Patro,R. et al. (2017) Salmon provides fast and bias-aware quantification of transcript expression. Nat. Methods, 14, 417–419.
5. Zhang,R. et al. (2017) A high quality Arabidopsis transcriptome for accurate transcript-level analysis of alternative splicing. Nucleic Acids Res.

# **Session Information**

```
## R version 3.4.0 (2017-04-21)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 7 x64 (build 7601) Service Pack 1
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United Kingdom.1252 
## [2] LC_CTYPE=English_United Kingdom.1252   
## [3] LC_MONETARY=English_United Kingdom.1252
## [4] LC_NUMERIC=C                           
## [5] LC_TIME=English_United Kingdom.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] knitr_1.16           leaflet_1.1.0        gridExtra_2.2.1     
## [4] ggplot2_2.2.1.9000   shinydashboard_0.6.0 shiny_1.0.5         
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.12      magrittr_1.5      munsell_0.4.3    
##  [4] colorspace_1.3-2  xtable_1.8-2      R6_2.2.1         
##  [7] rlang_0.1.1       stringr_1.2.0     plyr_1.8.4       
## [10] tools_3.4.0       grid_3.4.0        gtable_0.2.0     
## [13] gtools_3.5.0      htmltools_0.3.6   crosstalk_1.0.0  
## [16] yaml_2.1.14       lazyeval_0.2.0    digest_0.6.12    
## [19] tibble_1.3.3      htmlwidgets_0.9   evaluate_0.10    
## [22] mime_0.5          stringi_1.1.5     compiler_3.4.0   
## [25] scales_0.4.1.9002 jsonlite_1.5      httpuv_1.3.5
```
