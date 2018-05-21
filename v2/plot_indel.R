## Slimmed down single plotting function
library(ggplot2)

## Set the general formatting criteria
indel_plot_formatting <- theme(axis.text.x=element_text(angle=45, vjust=0.5, size=10,colour = "black"),
        axis.text.y=element_text(size=10,colour = "black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size=10),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA))

feature_barplot <- function(df,...){
  args <- list(...)
  params <- list(x='Subtype',y='aggregate',fill='Type',x_labels=df$Subtype,x_order=df$Subtype,xlab='Indel Channels',ylab='Counts',group_color='',
                ylim=c(0,1))
  params <- modifyList(params,args)
  p <- ggplot(df,aes_string(x=params$x,y=params$y,fill=params$fill)) + geom_bar(stat='identity',position = 'dodge') +
    xlab(label = params$xlab) + ylab(label= params$ylab)
  
  if(all(df[[params$y]] < 1)){
    p <- p+ scale_y_continuous(limits=params$ylim,breaks=(seq(0,1,0.2)),labels = scales::percent) +
      ylab('Percentage')
  }
  
  ## Relabel by x_labels
  p <- p+scale_x_discrete(labels = params$x_labels)
  
  ## Color processing 
  if(nchar(params$group_color) > 0){
    p <- p+scale_fill_manual(values=params$group_color)
  }

  ## Add the theme
  p <- p + indel_plot_formatting
  return(p)
}



