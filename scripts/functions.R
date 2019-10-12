require(tidyverse)
require(annotables)
require(grImport)
require(BuenColors)
require(cowplot)
require(edgeR)

# 
blood_cell_plot <- function(genes, palette = "solar_rojos", output_dir = ".",
                            type = "log2_cpm" ){
   myshape<- readPicture("./data/blank_cells.eps.xml")
   genes = genes
   n <- length(genes)
   df <- import_data(type = type)
   df <- anno_grch38(df)
   df <- filter_genes_long(df, genes)
   df <- generate_plot_info(df, palette, n, type = type)[[1]] %>% distinct()
   legend <- generate_plot_info(df, palette, n, type = type)[[2]]
   for(symb in unique(df$symbol)){
      curr <- df %>% filter(symbol == symb)
      my_shape <- color_cells(myshape, curr) %>% pictureGrob(.) %>% ggdraw(.)
      p <- my_shape + annotation_custom(legend, xmin =.75, ymin = .65)
      ggsave(file.path(output_dir,
                       paste0(symb, "_bloodexp.pdf")),
             p, width = 8, height = 6)
   }
}


# generate color graphic information using ggplot
generate_plot_info <- function(df, palette, n_genes, type){
   palette = jdb_palette(palette)
   df <- data.frame(df)
   if(n_genes == 1){
      p1 <- ggplot(df, aes(x = cell_type, y = expression, fill = expression)) +
         geom_bar(stat = "identity") + 
         coord_flip() +
         pretty_plot() + 
         scale_fill_gradientn(colors = palette)
   }else{
      p1 <- ggplot(df, aes(x = cell_type, y = expression, fill = expression)) +
         geom_bar(stat = "identity") + 
         coord_flip() +
         facet_wrap(~factor(ensg), ncol = n_genes) +
         pretty_plot() + 
         scale_fill_gradientn(colors = palette)
   }
   if(type == "cpm"){
      p1 <- p1 + labs(fill = "cpm")
   }else if(type == "log2_cpm"){
      p1 <- p1 + labs(fill = "log2(cpm + 1)")
   }else if(type == "raw"){
      p1 <- p1 + labs(fill = "counts")
   }
   
   df <- ggplot_build(p1)[[1]] %>%
      data.frame(.) %>% 
      dplyr::select(1,3) %>%
      left_join(df, ., by = c("expression" = "y"))
   legend <- get_legend(p1)
   list(df, legend)
}


# import data
import_data <- function(type){
   df <- read_tsv("https://raw.githubusercontent.com/jeffverboon/blood_expressions_plots/master/data/counts.tsv")
   if(type == "log2_cpm" || type == "cpm"){
      df[-1] <- lapply(df[-1], cpm)
      if(type == "log2_cpm"){
         df[-1] <- lapply(df[-1], function(x) log2(x + 1))
      }
   }
   df
}

# annotate the dataframe
anno_grch38 <- function(df){
   annotables::grch38 %>%
      mutate(ensg = ensgene) %>%
      dplyr::select(ensg, symbol) %>%
      right_join(., df)
}

# subset df by genes
filter_genes_long <- function(df, genes){
      df %>% 
      filter(symbol %in% genes | ensg %in% genes) %>%
      gather(cell_type, expression, -ensg, -symbol)
   if()
}


# function to change colors of cells
color_cells <- function(picture, df){
   my_shape = picture
   B <- df %>% filter(cell_type == "B") %>% pull(fill)
   CD4 <- df %>% filter(cell_type == "CD4") %>% pull(fill)
   CD8 <- df %>% filter(cell_type == "CD8") %>% pull(fill)
   CLP <- df %>% filter(cell_type == "CLP") %>% pull(fill)
   CMP <- df %>% filter(cell_type == "CMP") %>% pull(fill)
   ERY <- df %>% filter(cell_type == "ERY") %>% pull(fill)
   GMPA <- df %>% filter(cell_type == "GMPA") %>% pull(fill)
   GMPB <- df %>% filter(cell_type == "GMPB") %>% pull(fill)
   GMPC <- df %>% filter(cell_type == "GMPC") %>% pull(fill)
   HSC <- df %>% filter(cell_type == "HSC") %>% pull(fill)
   LMPP <- df %>% filter(cell_type == "LMPP") %>% pull(fill)
   MEGA <- df %>% filter(cell_type == "MEGA") %>% pull(fill)
   MEP <- df %>% filter(cell_type == "MEP") %>% pull(fill)
   MONO <- df %>% filter(cell_type == "MONO") %>% pull(fill)
   MPP <- df %>% filter(cell_type == "MPP") %>% pull(fill)
   NK <- df %>% filter(cell_type == "NK") %>% pull(fill)
   GRAN <- df %>% filter(cell_type == "GRAN") %>% pull(fill)
   MDC <- df %>% filter(cell_type == "MDC") %>% pull(fill)
   PDC <- df %>% filter(cell_type == "PDC") %>% pull(fill)
   PLT <- df %>% filter(cell_type == "PLT") %>% pull(fill)
   
   # B cell
   my_shape@paths[59]$path@rgb <- B
   my_shape@paths[181]$path@rgb <- B
   
   #CD4  53, 177
   my_shape@paths[53]$path@rgb <- CD4
   my_shape@paths[177]$path@rgb <- CD4
   
   #CD8 56, 179 
   my_shape@paths[56]$path@rgb <- CD8
   my_shape@paths[179]$path@rgb <- CD8
   
   #CLP  44, 169
   my_shape@paths[44]$path@rgb <- CLP
   my_shape@paths[169]$path@rgb <- CLP
   
   #CMP   46, 171
   my_shape@paths[46]$path@rgb <- CMP
   my_shape@paths[171]$path@rgb <- CMP
   
   #Ery   68, 70, 72, 187, 205, 207
   my_shape@paths[68]$path@rgb <- ERY
   my_shape@paths[70]$path@rgb <- ERY
   my_shape@paths[72]$path@rgb <- ERY
   my_shape@paths[187]$path@rgb <- ERY
   my_shape@paths[205]$path@rgb <- ERY
   my_shape@paths[207]$path@rgb <- ERY
   
   #GMP-A 97, 173
   my_shape@paths[97]$path@rgb <- GMPA
   my_shape@paths[173]$path@rgb <- GMPA
   
   #GMP-B 84, 197
   my_shape@paths[84]$path@rgb <- GMPB
   my_shape@paths[197]$path@rgb <- GMPB
   
   # GMP-C 48, 199
   my_shape@paths[48]$path@rgb <- GMPC
   my_shape@paths[199]$path@rgb <- GMPC
   
   # GRAN 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141,
   #       143, 145, 147, 149, 151, 153
   my_shape@paths[121]$path@rgb <- GRAN
   my_shape@paths[123]$path@rgb <- GRAN
   my_shape@paths[125]$path@rgb <- GRAN
   my_shape@paths[127]$path@rgb <- GRAN
   my_shape@paths[129]$path@rgb <- GRAN
   my_shape@paths[131]$path@rgb <- GRAN
   my_shape@paths[133]$path@rgb <- GRAN
   my_shape@paths[135]$path@rgb <- GRAN
   my_shape@paths[137]$path@rgb <- GRAN
   my_shape@paths[139]$path@rgb <- GRAN
   my_shape@paths[141]$path@rgb <- GRAN
   my_shape@paths[145]$path@rgb <- GRAN
   my_shape@paths[147]$path@rgb <- GRAN
   my_shape@paths[149]$path@rgb <- GRAN
   my_shape@paths[143]$path@rgb <- GRAN
   my_shape@paths[151]$path@rgb <- GRAN
   my_shape@paths[153]$path@rgb <- GRAN
   
   # HSC   25, 163
   my_shape@paths[25]$path@rgb <- HSC
   my_shape@paths[163]$path@rgb <- HSC
   
   # LMPP  31, 167
   my_shape@paths[31]$path@rgb <- LMPP
   my_shape@paths[167]$path@rgb <- LMPP
   
   # mDC   115, 201
   my_shape@paths[115]$path@rgb <- MDC
   my_shape@paths[201]$path@rgb <- MDC
   
   # mega  1, 155, 157, 159, 161
   my_shape@paths[1]$path@rgb <- MEGA
   my_shape@paths[157]$path@rgb <- MEGA
   my_shape@paths[159]$path@rgb <- MEGA
   my_shape@paths[161]$path@rgb <- MEGA
   my_shape@paths[155]$path@rgb <- MEGA
   
   # MEP   50, 175
   my_shape@paths[50]$path@rgb <- MEP
   my_shape@paths[175]$path@rgb <- MEP
   
   # Mono  65, 185
   my_shape@paths[65]$path@rgb <- MONO
   my_shape@paths[185]$path@rgb <- MONO
   
   # MPP   27, 165
   my_shape@paths[27]$path@rgb <- MPP
   my_shape@paths[165]$path@rgb <- MPP
   
   # NK    62, 183
   my_shape@paths[62]$path@rgb <- NK
   my_shape@paths[183]$path@rgb <- NK
   
   # pDC   111, 203
   my_shape@paths[111]$path@rgb <- PDC
   my_shape@paths[203]$path@rgb <- PDC
   
   # Plt   190, 192, 194, 217, 219, 221, 223, 225
   my_shape@paths[190]$path@rgb <- PLT
   my_shape@paths[192]$path@rgb <- PLT
   my_shape@paths[194]$path@rgb <- PLT
   my_shape@paths[217]$path@rgb <- PLT
   my_shape@paths[219]$path@rgb <- PLT
   my_shape@paths[221]$path@rgb <- PLT
   my_shape@paths[223]$path@rgb <- PLT
   my_shape@paths[225]$path@rgb <- PLT
   return(my_shape)
}   
   
