suppressMessages(library(stringr))

# download required packages if not already installed
get_reqs <- function(reqs){
  reqs_not_present <- reqs[!(reqs %in% installed.packages()[,'Package'])]
  if(length(reqs_not_present) > 0){
    reqs_to_install <- paste(reqs_not_present, sep = ', ')
    writeLines(paste(reqs_to_install, 'not found, installing now.'))
    install.packages(reqs_not_present, repos = 'http://cran.us.r-project.org')
  }
}

# download Bioc packages if not already installed
get_bioc_reqs <- function(reqs){
  
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  
  reqs_not_present <- reqs[!(reqs %in% installed.packages()[,'Package'])]
  if(length(reqs_not_present) > 0){
    reqs_to_install <- paste(reqs_not_present, sep = ', ')
    writeLines(paste(reqs_to_install, 'not found, installing now.'))
    BiocManager::install(reqs_not_present)
  }
}

parseArgs <- function(args, keys){
  
  args <- paste(args, collapse = ' ')
  
  parsed_args <- unlist(str_split(args, '--[:alpha:]+_{0,1}[:alpha:]*'))
  parsed_args <- str_trim(parsed_args[-1])
  argnames <- unlist(str_extract_all(args, '(?<=--)[:alpha:]+_{0,1}[:alpha:]*'))
  
  if(length(parsed_args) != length(argnames)){
    stop('Number of argument markers does not match number of arguments provided. Try again.')
  }
  
  names(parsed_args) <- argnames
  
  return(parsed_args)
}

checkArgs <- function(parsed_args, req_args, file_args){

  # check if any of the required arguments have been omitted
  if(any(!(req_args %in% names(parsed_args)))){
    omitted_args <- paste(req_args[!(req_args %in% names(parsed_args))], collapse = ', ')
    stop(paste('The following required arguments were omitted:', omitted_args))
  }
  
  if(length(file_args) > 0){
    # check if all the files can be found
    if(any(!file.exists(parsed_args[file_args]))){
      files_not_found <- paste(parsed_args[file_args][!file.exists(parsed_args[file_args])], collapse = ', ')
      stop(paste0('The following files arguments were not found: ', files_not_found, 
                  '. Please check you have provided the correct filepath(s).'))
    }
  }
  
  writeLines('\nScript arguments parsed, running script...\n')
  
}

# function to make p-value histogram from DESeq results
make_pvalue_histogram <- function(DESeq_results, title, save = FALSE){
  
  df <- as.data.frame(DESeq_results)
  
  plot <- ggplot(df, aes(x = padj)) + geom_histogram() + 
    theme_minimal() +
    labs(x = '\nAdjusted p-value', title = gsub('_', ' ', title)) +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'), text = element_text(size = 14))
  
  if(save){
    filepath = paste('figures/pvalue_histogram_', title, '.png', sep = '')
    ggsave(filepath, width = 33.8, height = 19.5, units = 'cm')
  }
  
  return(plot)
}

# function to make volcano plot from DESeq results
make_volcano_plot <- function(DESeq_results, title, transcript_type = c('Human genes', 'Retroelements'), 
                              save = FALSE){
  
  volcano_plot <- EnhancedVolcano::EnhancedVolcano(DESeq_results, lab = rownames(DESeq_results), 
                                                   x = 'log2FoldChange', y = 'padj', caption = '', 
                                                   title = gsub('_', ' ', title), labSize = 7, 
                                                   subtitle = transcript_type, titleLabSize = 24,
                                                   subtitleLabSize = 22) + xlim(-20, 20)
  
  if(save){
    filepath = paste('figures/volcano_plot_', title, '.png', sep = '')
    ggsave(filepath, width = 33.8, height = 19.5, units = 'cm')
  }
  
  return(volcano_plot)
  
}
