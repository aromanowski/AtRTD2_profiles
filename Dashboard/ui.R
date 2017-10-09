## app.R ##
library(shinydashboard)
library(ggplot2)
library(gridExtra)
library(shiny)
library(leaflet)

#############################################################################
#   █▀▀ █▀▀█ █──█ █▀▀█ █▀▀ █▀▀ 　 █▀▀ █──█ █▀▀▄ █▀▀ ▀▀█▀▀ ─▀─ █▀▀█ █▀▀▄     #
#   ▀▀█ █──█ █──█ █▄▄▀ █── █▀▀ 　 █▀▀ █──█ █──█ █── ──█── ▀█▀ █──█ █──█     #
#   ▀▀▀ ▀▀▀▀ ─▀▀▀ ▀─▀▀ ▀▀▀ ▀▀▀ 　 ▀── ─▀▀▀ ▀──▀ ▀▀▀ ──▀── ▀▀▀ ▀▀▀▀ ▀──▀     #
#############################################################################

##------start
#' Plot expression profiles of RNA-seq data of AtRTD2
#' 
#' @param data.exp transcript level expresssion to make plot, such as TPM and read counts. The columns are the samples and rows are the transcripts names.
#' @param gene a gene name to make the plot. Subset of transcripts expression matching to the gene name will be used to make the plot.
#' @param reps a vector that indicates the replicates. In the AtRTD2 dataset, \code{reps=rep(1:26,each=9)}.
#' @param y.lab y axis title, defualt is \code{y.lab='TPM'}
#' @param legend.text.size legend text size
#' @param marker.size the marker size passed to \code(geom_point(size=marker.size)) in \link{\code{ggplot2}}.
#' @param error.type the method to calculate the errors shown as error bars. Options are "stderr" for standard error (default) and "sd" for standard deviation.
#' @param plot.gene logical, to show the gene expression profile in the plot (TRUE) or not (FALSE). 
#' @param error.bar logical, to show error bars (TRUE) or not (FALSE).
#' @param show.plot logical, to return a grahical plot (TRUE) or not (FALSE). 
#' @return a list of plots. In the list, "profiles" is the line plot of profiles, "lightbar" is the lightbar and "plots" is a list after aligning the profiles and lightbar plots to consisent x coordinates.
#' @examples 
#' data(AtRTD2.example)
#' plots <- plot_AtRTD2(data.exp=AtRTD2.example$TPM,gene='AT5G53180',y.lab='TPM',
#'                     reps=rep(1:26,each=9),error.type='stderr',plot.gene=T,error.bar = T)
#' names(plots)
#' grid.arrange(plots$profiles,plots$lightbar,ncol=1,heights=c(6,1))

plot_AtRTD2 <- function(data.exp,gene,
                        reps=rep(1:26,each=9),y.lab='TPM',marker.size=3,legend.text.size=13,text.size=16,
                        error.type='stderr',plot.gene=T,error.bar=T,show.plot=T){
  #####################################################################
  ##prepare plot data                                                
  #####################################################################
  
  trans.idx <- grep(pattern = gene,rownames(data.exp))
  expression.sum <- if(length(trans.idx)==1) sum(data.exp[trans.idx,]) else rowSums(data.exp[trans.idx,])
  trans.idx <- trans.idx[expression.sum>0]
  
  if(length(trans.idx)==0)
    stop(paste0(gene, ' is invalid'))
  
  data2plot<-data.exp[trans.idx,]
  trans <- rownames(data.exp)[trans.idx]
  if(length(trans.idx)==1)
    data2plot <- rbind(data2plot,data2plot) else data2plot <- rbind(colSums(data2plot),data2plot)
  rownames(data2plot) <- c(gene,trans)
  
  ##--mean and error
  mean2plot <- t(rowmean(t(data2plot),group = reps))
  sd2plot <- by(t(data2plot),INDICES = reps,FUN = function(x){
    apply(x,2,function(y){
      stderr(y,error.type = error.type)
    })
  })
  sd2plot <- do.call(cbind,sd2plot)
  sd2plot <- sd2plot[rownames(mean2plot),]
  
  ##--add null to time point 17.5
  mean2plot <- cbind(mean2plot,`17.5`=rep(NA,nrow(mean2plot)))
  sd2plot <- cbind(sd2plot,`17.5`=rep(NA,nrow(sd2plot)))
  
  data2plot <- data.frame(reshape2::melt(mean2plot),reshape2::melt(sd2plot)[,3])
  colnames(data2plot) <- c('Targets','time','mean','sd')
  
  
  ##--plot gene or not
  if(!plot.gene)
    data2plot <- data2plot[-which(data2plot$Targets==gene),]
  
  
  profiles <- ggplot(data2plot,aes(x=time,y=mean,color=Targets))+geom_line()+
    geom_point(size=marker.size,aes(fill=Targets,shape=Targets))+
    scale_shape_manual(name="Targets",values=c(25:0,25:0,rep(25,500)))+
    scale_x_continuous(limits=c(0.8,26.2),breaks=1:26)+
    geom_vline(xintercept =9,color='skyblue',linetype='dotdash',size=1.1)+
    # geom_segment(x=9,xend=10,y=max(na.omit(data2plot$mean+data2plot$sd)),
    #              yend=max(na.omit(data2plot$mean+data2plot$sd)),color='black',size = 0.4,
    #              arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "closed"))+
    annotate('text',x=9,y=1.02*max(na.omit(data2plot$mean+data2plot$sd)),
             label='Cold stress',color='black',size=5)+
    labs(x='Time',y=y.lab,title=gene)+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.text=element_text(color='black',size=legend.text.size),
          text = element_text(size=text.size),axis.text = element_text(size=text.size+1))
  if(error.bar)
    profiles <- profiles+
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.3,color='black',size=0.3)
  
  
  lightbar <- ggplot() +
    geom_rect(data = NULL,
              aes(xmin=c(1,9,18),xmax=c(5,13,22),ymin=c(0,0,0),ymax=c(1.5,1.5,1.5)),
              fill='black',color='black')+
    geom_rect(data = NULL,
              aes(xmin=c(5,13,22),xmax=c(9,17,26),ymin=c(0,0,0),ymax=c(1.5,1.5,1.5)),
              fill='white',color='black')+
    
    geom_text(data=NULL,aes(x=c(3,7,5,13,22),y=c(2,2,-0.5,-0.5,-0.5),
                            label=c('dark','light',paste(20,'^o~~C'),paste(4,'^o~C~~Day1'),paste(4,'^o~C~~Day4'))),
              parse=TRUE,size=5)+
    theme_bw()+
    scale_x_continuous(limits=c(0.8,26.2),breaks=1:26)+
    theme(axis.text  = element_blank(),panel.grid = element_blank(),axis.title = element_blank(),
          axis.ticks = element_blank(),panel.border = element_blank(),
          text = element_text(size=text.size))+
    coord_cartesian(ylim=c(-1,2.2))
  
  suppressWarnings(plots<- AlignPlots(profiles, lightbar))
  names(plots) <- c('profiles','lightbar')
  if(show.plot)
    grid.arrange(plots$profiles,plots$lightbar,ncol=1,heights=c(6,1))
  return(list(profiles=profiles,lightbar=lightbar,plots=plots))
}


#############################################################################
#               █▀▀ █──█ ─▀─ █▀▀▄ █──█ 　 █▀▀█ █▀▀█ █▀▀█                    #
#               ▀▀█ █▀▀█ ▀█▀ █──█ █▄▄█ 　 █▄▄█ █──█ █──█                    #
#               ▀▀▀ ▀──▀ ▀▀▀ ▀──▀ ▄▄▄█ 　 ▀──▀ █▀▀▀ █▀▀▀                    #
#############################################################################
##------start


AlignPlots <- function(...) {
  setTimeLimit(elapsed = Inf,cpu = Inf, transient = TRUE)
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(grid::unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  if(any(sapply(legends.widths,is.null)))
    max.legends.width <- NULL else max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (gtable::is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),"mm"))
    }
    x
  })
  plots.grobs.eq.widths.aligned
}



rowmean <- function (x, group, reorder = T, na.rm = T) {
  order.idx <- as.character(unique(group))
  if (reorder) 
    order.idx <- gtools::mixedsort(order.idx)
  counts <- table(group)[order.idx]
  sums <- rowsum(x, group = group)[order.idx, ]
  means <- (diag(1/counts)) %*% as.matrix(sums)
  rownames(means) <- order.idx
  if (na.rm) 
    means[is.na(means)] <- 0
  return(means)
}


stderr <- function(x,error.type='stderr') {
  if(error.type=='stderr')
    error<-sd(na.omit(x))/sqrt(length(na.omit(x)))
  if(error.type=='sd')
    error<-sd(na.omit(x))
  return(error)
}

##------end

#############################################################################
##--------------------------------pass word--------------------------------##
#############################################################################

Logged <- FALSE;
LoginPass <- 0; #0: not attempted, -1: failed, 1: passed

login <- box(title = h4("Login to enable plot"),
             textInput(inputId = "userName",value = 'user',label = "Username (user)"),
             # passwordInput("passwd", "Password (test)"),
             passwordInput("passwd", "Password",value = 'AtRTD2profiles'),
             br(),
             actionButton("Login", "Login"),
             width = 3)

loginfail <- box(title = "Login",
                 textInput(inputId = "userName",value = 'user',label = "Username (user)"),
                 passwordInput("passwd", "Password"),
                 p("Username or password incorrect"),
                 br(),
                 actionButton("Login", "Log in"),
                 width = 3)

#############################################################################
##--------------------------------header-----------------------------------##
#############################################################################
header <- dashboardHeader(title = "AtRTD2 Dashboard")


#############################################################################
##--------------------------------sidebar----------------------------------##
#############################################################################
sidebar <- dashboardSidebar(
  sidebarMenu(  
    id = "tabs",
    menuItem("Introduction", tabName = "introduction", icon = icon("book")),
    menuItem("Expression profiles", icon = icon("line-chart"), tabName = "profiles"),
    menuItem("Contact us", icon = icon("address-card-o"), tabName = "map")
    # uiOutput("sidebarpanel"))
  )
)

#############################################################################
##--------------------------------body-------------------------------------##
#############################################################################

mainbody <- dashboardBody(
  tabItems(
    tabItem("introduction",
            # div(p("Dashboard tab content"))
            box(width = NULL, status = "primary",
                uiOutput("page1")
            )
    ),
    tabItem("profiles",
            fluidRow(
              column(width = 10,
                     box(title='Profile plot',width = NULL, solidHeader = F,status = 'primary',
                         plotOutput("profile_plot", height = '550px' )
                     )
              ),
              column(width = 2,
                     box(width = NULL, status = "primary",
                         
                         ##select or input a gene name
                         
                         textInput("genename", label = "Input a gene",
                                   value = ''
                         ),
                         
                         # uiOutput("gene_ui"),#for dropdown list
                         ##action button to make the plot
                         actionButton('plotprofiles','Click to plot',icon("send outline icon")),
                         div(p('')),
                         div(p('Note: Input gene ID of TAIR/Araport11, e.g. AT1G01060'))
                     ),
                     box(width = NULL, status = "primary",
                         radioButtons('plottype',label='Select a type',
                                      choices =list('png'='png','pdf'='pdf'),
                                      selected = 'png',inline = F),
                         div(align='left',downloadButton('profile_download', 
                                                         'Download'))
                         
                     ),
                     box(width=NULL,status = "primary",
                         textOutput('loading_data'),
                         p('Notes: It may take a while to load the RNA-seq 
                           transcript expression dataset.
                           Please wait until the gene name selection box activated.')
                         )
                         )
                     ),
            fluidRow(
              column(width = 12,
                     box(title='Profile in TPM (transcript per million)',
                         width = NULL,status = 'primary', solidHeader = F,
                         shiny::dataTableOutput('profile_table')
                     )
              )
            ) 
              ),
    tabItem("map",
            box(width = 12, status = "primary",
                leafletOutput("mymap")
            ),
            fluidRow(
              column(6,
                     box(width = 12, status = "primary",title = 'Contact us',
                         # div(class="section level4",h4('Address:')),
                         div('John W. S. Brown'),
                         div('j.w.s.brown@dundee.ac.uk'),
                         div('Plant Sciences Division, College of Life Sciences, University of Dundee'),
                         div('Cell and Molecular Sciences, The James Hutton Institute'),
                         div('Invergowrie, Dundee DD2 5DA, Scotland UK')              
                     )
              ),
              column(6,
                     box(width = 12, status = "primary",title = 'Contact us',
                         # div(class="section level4",h4('Address:')),
                         div('Runxuan Zhang'),
                         div('runxuan.zhang@hutton.ac.uk'),
                         div('Information and Computational Sciences'),
                         div('The James Hutton Institute'),
                         div('Invergowrie, Dundee DD2 5DA, Scotland UK')
                     )
              )
            )
            
    )
  )
  )


#############################################################################
##----------------------------------ui-------------------------------------##
#############################################################################

dashboardPage(
  header,
  sidebar,
  body = mainbody
)


