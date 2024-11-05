plot_venn_3sets <- function(data_list, 
                       param, 
                       labels, 
                       title,
                       colors=c("green", "red", "blue"))
{
  if (param == 'All SNPs') {
    set_list <- lapply(data_list, function(x) names(x)) 
  }
  else {
    if (param == "Overlapped SNPs") {
      set_list <- lapply(data_list, function(x) unique(x$snp_id))
    }
    else {
      stop("Error: param can be either 'All SNPs' or 'Overlapped SNPs'")
    }
  }
  
  catnames <- labels
  dists <- c(0.05, 0.05, 0.05)
  catposs <- c(0, 0, 180)
  titlepos <- c(0.5, 0.95) # Adjusted title position
  
  # Create the Venn diagram with the title as part of the diagram
  venn_plot <- venn.diagram(
    x = set_list,
    category.names = catnames,  # Use empty strings to hide labels
    filename = NULL,  # To plot directly instead of saving to a file
    output = TRUE,
    fill = colors , # Colors for each set
    alpha = 0.4,
    cat.col = colors ,
    cat.cex = 1,
    margin = 0.1,
    main = title, # Add title within the Venn diagram
    main.cex = 1, # Adjust the size of the title
    main.pos = titlepos, # Adjust the position of the title
    main.fontfamily = "sans",
    main.fontface = "bold",
    cat.fontfamily = "sans",
    fontfamily = "sans",
    cat.pos = catposs,
    cat.dist = dists
  )
  
  # Draw the Venn Diagram
  grid.newpage() # Opens a new graphics page
  grid.draw(venn_plot) # Draws the Venn Diagram
}



plot_venn_4sets <- function(data_list, 
                       param, 
                       labels, 
                       title,
                       colors=c("green", "red", "blue", "yellow"))
{
  if (param == 'All SNPs') {
    set_list <- lapply(data_list, function(x) names(x))
  }
  else {
    if (param == "Overlapped SNPs") {
      set_list <- lapply(data_list, function(x) unique(x$snp_id))
    }
    else {
      stop("Error: param can be either 'All SNPs' or 'Overlapped SNPs'")
    }
  }
  
  catnames <- rep("", 4)
  dists <- c(0.2, 0.2, 0.2, 0.1)
  catposs <- c(0, 0, 0, 0)
  titlepos <- c(0.5, 0.85)
  
  # Create the Venn diagram with the title as part of the diagram
  venn_plot <- venn.diagram(
    x = set_list,
    category.names = catnames,  # Use empty strings to hide labels
    filename = NULL,  # To plot directly instead of saving to a file
    output = TRUE,
    fill = colors , # Colors for each set
    alpha = 0.4,
    cat.col = colors ,
    cat.cex = 1,
    margin = 0.1,
    main = title, # Add title within the Venn diagram
    main.cex = 1, # Adjust the size of the title
    main.pos = titlepos, # Adjust the position of the title
    main.fontfamily = "sans",
    main.fontface = "bold",
    cat.fontfamily = "sans",
    fontfamily = "sans",
    cat.pos = catposs,
    cat.dist = dists
  )
  
  
  # Create a new graphics page with an adjusted layout for Venn and legend
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(6, 2), "null")))) # Adjust layout heights
  
  # Draw the Venn Diagram in the first row
  pushViewport(viewport(layout.pos.row = 1))
  grid.draw(venn_plot)
  popViewport()
  
  # Create a legend grob with two columns
  legend_grob <- legendGrob(
    labels = labels,
    pch = 15, # Square symbols
    gp = gpar(col = colors , fill = colors ),
    ncol = 2 # Arrange the legend items in 2 columns
  )
  
  # Draw the legend on the second row
  pushViewport(viewport(layout.pos.row = 2))
  grid.draw(legend_grob)
  popViewport()
}