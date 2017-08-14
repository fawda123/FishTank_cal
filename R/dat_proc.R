######
# create rdata files in data

source('R/funcs.R')
library(FME)
library(tidyr)
library(dplyr)
library(lubridate)
library(ncdf4)
library(tibble)

##
# calibrating to parameter subsets

rm(list = ls())
source('R/funcs.R')

load(file = 'rdata/tocal_all.RData')

idpars <- get_cmbs(tocal_all, 'O2', coll = FALSE)

cals <- vector('list', length = length(idpars))
names(cals) <- names(idpars)
for(idpar in names(idpars)){

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
