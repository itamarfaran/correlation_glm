dir <- 'tex/simulations'

all_files <- list.files(dir, '*.R$')
all_files <- all_files[!all_files %in% c('aux.R', 'aux_.R', 'run_all.R')]

for(file in all_files){
  tryCatch({
    message(paste0('source ', file, '...'))
    source(paste0(dir, '/', file))
    message('success.')
  },
    warning = function(w) warning(paste0('in ', file, ': ', w)),
    error = function(e) warning(paste0('in ', file, ': ', e)),
    finally = rm(list = setdiff(ls(), c('dir', 'all_files', 'file')))
  )
}
