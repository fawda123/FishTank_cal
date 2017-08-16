######
# create rdata files in data

source('R/funcs.R')
library(FME)
library(tidyr)
library(dplyr)
library(lubridate)
library(ncdf4)
library(tibble)

# ##
# # sensitivity of all state variables (most)
# 
# # state variables to evaluate
# sts <- c('Chla_mg_tot', 'O2', 'OM1', 'OM2', 'PO4', 'NH4', 'NO3', 'irradiance')
# sts <- sort(sts)
# 
# # parameters to eval
# pars <- par_tst()
# 
# # sens eval for all state variables and relevant pars
# sens_ests_all <- sensfun(pars = pars, out_var = sts, p1z1 = TRUE)
# 
# # save
# save(sens_ests_all, file = 'rdata/sens_ests_all.RData', compress = 'xz')
# 
# ##
# # for all state variables, combine individual sensitivity estimates with categories
# # retain those that are not zero or not NA (NA means the orig value was set to zero)
# 
# load(file = 'rdata/sens_ests_all.RData')
# 
# # get parameter, cat labels
# cats <- reshape2::melt(parcats()) %>% 
#   dplyr::select(L1, L2) %>% 
#   rename(Category = L1, Parameter = L2) 
# 
# # sens_ests_all <- sens_ests_all[!names(sens_ests_all) %in% 'PO4']
# # get sens, merge with cats, arrange by cat
# sens_ests_cat_all <- lapply(sens_ests_all, function(x){
#   
#   out <- dplyr::select(x$sens, Parameter, L1) %>% 
#     rename(
#       error = L1
#     ) %>% 
#     left_join(cats, by = 'Parameter') %>% 
#     dplyr::filter(error > 0 & !is.na(error)) %>% 
#     group_by(Category) %>% 
#     arrange(error) %>% 
#     ungroup %>% 
#     mutate(Parameter = factor(Parameter, levels = Parameter))
#   
#   return(out)
#   
# })
# 
# save(sens_ests_cat_all, file = 'rdata/sens_ests_cat_all.RData', compress = 'xz')

##
# calibrating to parameter subsets

load(file = 'rdata/tocal_all.RData')

idpars <- get_cmbs(tocal_all, 'O2', coll = FALSE, pert = NULL)

cals <- vector('list', length = length(idpars))
names(cals) <- names(idpars)
for(idpar in 'Optics'){#:names(idpars)){

 # output
 cat(idpar, '\n')

 # get heurist and parameter subsets
 parsin <- idpars[[idpar]]
 vals <- parsin$vals
 minv <- parsin$minv
 maxv <- parsin$maxv

 # calibrate
 res <- fishopt(vals, minv, maxv)

 # append to output, save as you
 cals[[idpar]] <- res
 save(cals, file = 'rdata/cals.RData', compress = 'xz')

}

save(cals, file = 'rdata/cals.RData', compress = 'xz')

#######
# reset RData input file from ignore folder when done
GEM_InputFile <- formpars('ignore/GEM_InputFile')
save(GEM_InputFile, file = 'input/GEM_InputFile.RData')  
