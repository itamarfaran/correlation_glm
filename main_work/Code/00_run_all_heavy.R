tryCatch(source("main_work/Code/05_main_code_noPlots.R"), error = function(e) print(e))
remove(list = ls()); gc()

tryCatch(source("main_work/Code/06_bootstraps.R"), error = function(e) print(e))
remove(list = ls()); gc()

tryCatch(source("main_work/Code/07_bootstrapsARMA.R"), error = function(e) print(e))
remove(list = ls()); gc()

tryCatch(source("main_work/Code/08_power_comparisions.R"), error = function(e) print(e))
remove(list = ls()); gc()

#tryCatch(source("main_work/Code/09_ndma.R"), error = function(e) print(e))
# remove(list = ls()); gc()

tryCatch(source("main_work/Code/10_adni.R"), error = function(e) print(e))
remove(list = ls()); gc()

tryCatch(source("main_work/Code/11_tga.R"), error = function(e) print(e))
remove(list = ls()); gc()
