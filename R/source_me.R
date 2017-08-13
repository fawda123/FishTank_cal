source('R/funcs.R')
library(FME)
library(SWMPr)
library(tidyverse)
library(WtRegDO)
library(lubridate)
library(oce)

######

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
  res <- fishopt(vals, minv, maxv, factr = 1e9)

  # append to output, save as you
  cals[[idpar]] <- res
  save(cals, file = 'rdata/cals_srcme.RData', compress = 'xz')

}
