# a caller for state variable labels
labels_fun <- function(){
  
  if(!file.exists('NETCDF/output.000000.nc'))
    stop('Run model to get labels')

  # get model output from netcdf
  nc <- nc_open('NETCDF/output.000000.nc')
  
  # shrt names, all other vectors will be ordered using shrt
  shrt <- sort(names(nc$var))

  # long names
  lngs <- sapply(shrt, function(x) ncatt_get(nc, x, 'description')$value)
  lngs <- lngs[shrt]
  
  # units for each variable
  vals <- sapply(shrt, function(x) ncatt_get(nc, x, 'units')$value)
  vals <- vals[shrt]
  
  # remove names
  names(lngs) <- NULL
  names(vals) <- NULL
  
  # close connection
  nc_close(nc)
  
  out <- list(shrt = shrt, lngs = lngs, vals = vals)
  return(out)

}

# convert the initial conditions between formats
#
# parsin can be a chr string of file location of original InitialConditions.txt or input data frame of input conditions to be converted to ASCII format, the data frame is the output from the ASCII text file, but there's an additional function that replaces parameters in the data frame for conversion to the standard GEM format
forminps <- function(parsin){
  
  # data frame to input file  
  if(inherits(parsin, 'data.frame')){
    
    # combine for easy write with writelines
    out <- with(parsin, paste(value, parm, sep = '\t!'))
    
  }
  
  # input file to data frame
  if(inherits(parsin, 'character')){
    
    # sanity checks
    if(!file.exists(parsin)) stop('Input file does not exist')

    # convert to data frame
    out <- readLines(parsin) %>% 
      gsub('\t+!|\\s+!', '\t', .) %>% # replace strings inbetween value and variable name with tab
      gsub('\\s+$', '', .) %>% # remove trailing spaces
      strsplit('\t') %>% 
      do.call('rbind', .) %>% 
      data.frame(., stringsAsFactors = FALSE)
    names(out) <- c('value', 'parm')
  
  }
  
  return(out)
  
}

# convert the input parameter info between formats
#
# parsin can be a chr string of file location of original GEM_InputFile or input data frame of parameters to be converted to ASCII format, the data frame is the ouput from the ASCII text file, but there's an additional function that replaces parameters in the data frame for conversion to the standard GEM format
formpars <- function(parsin){
  
  # data frame to input file  
  if(inherits(parsin, 'data.frame')){

    # split data frame by parameters
    # reorganize duplicates as single chr vector
    # back to data frame
    tmp <- mutate(parsin, 
        parm = gsub('_[0-9]*$', '', parm),
        parm = factor(parm, levels =  unique(parm))
      ) %>% 
      split(., .$parm) %>% 
      lapply(., function(x) paste(x[1][, 1, drop = TRUE], collapse = ' ')) %>% 
      reshape2::melt(.) %>% 
      rename(parm = L1) %>% 
      mutate(
        ord = as.numeric(row.names(.)), 
        parm = paste0('!', parm), 
        value = paste0(value, '\t\t')
        )

    # add category labels and NA rows above
    # added by index but with decimal change for ordering
    cats <- c('Simulation Specifics', 'Switches in GEM', 'Optics', 'Temperature', 'Phytoplankton, up to 6 types', 'Zooplankton, up to 2 types', 'Organic Matter', 'Other including Boundary Conditions') 
    cats <- paste0('!', cats)
    locs <- grep('starting time|^!Which_fluxes|^!Kw$|^!Tref|^!ediblevector\\(Z1|^!Zeffic|^!KG1$|^!rcNO3|^!Which_VMix', tmp$parm, ignore.case = F)
    locs <- locs - 0.01
    lab <- data.frame(value = cats, parm = '', ord = locs)
    labfill <- data.frame(value = '', parm = '', ord = locs - 0.01)
    
    tmp <- rbind(tmp, lab, labfill) %>% 
      arrange(ord) %>% 
      select(-ord)
      
    # combine columns to vector
    out <- paste(tmp$value, tmp$parm, sep = '')
      
  }
  
  # input file to data frame
  if(inherits(parsin, 'character')){
    
    # sanity checks
    if(!file.exists(parsin)) stop('Input file does not exist')

    # readby lines, one list element per parameter 
    tmp <- readLines(parsin) %>% 
      strsplit('!') %>% 
      lapply(., function(x){
        tospl <- x[1] %>% 
          gsub('\t', '', .) %>%
          strsplit(' ') %>% 
          .[[1]] %>% 
          .[nchar(.) > 0] %>% 
          list(., x[2])
        tospl
      })

    # get parameter names, everything left of colon
    nms <- lapply(tmp, function(x) x[[2]]) %>% 
      unlist %>% 
      gsub(':.*$', '', .)
    
    # get parameter values, rename list elements as parameter names
    tmp <- lapply(tmp, function(x) x[[1]]) 
    names(tmp) <- nms

    # list elements that contain parameters (inverse of notparms)
    cats <- c('Simulation Specifics', 'Switches in GEM', 'Optics', 'Temperature', 'Phytoplankton, up to 6 types', 'Zooplankton, up to 2 types', 'Organic Matter', 'Other including Boundary Conditions') 
    notparms <- which(nms %in% cats)
    if(length(notparms) != length(cats)) stop('Cannot find all parameter categories in input file')
    notparms <- sort(c(as.numeric(which(is.na(tmp))), notparms))
    parms <- seq(length(tmp))
    parms <- parms[!parms %in% notparms]
    tmp <- tmp[parms] 

    # empties are important, need to convert to blank character
    tmp[unlist(lapply(tmp, length)) == 0] <- ' '
    
    # suffix for parameters with more than one value
    suff <- lapply(tmp, function(x) seq(length(x))) %>% 
      unlist #%>% 
      # gsub('0', '1', .) # the blank rows are important
    
    # melt list, add suffix to id more than one parameter value
    out <- reshape2::melt(tmp) %>% 
      mutate(
        L1 = paste0(L1, '_', suff),
        value = as.character(value)
        ) %>% 
      rename(parm = L1)
    
  }
  
  return(out)
  
}

# create file for initial conditions, defaults to existing RData object is inps = NULL, otherwise values are replaced
# 
# inps is a named list where each element is one to many parameter values for each input condition
# partial string matching is used for the names to replace values in a default input list
# passing NULL to inps will return the existing input conditions file
setinps <- function(inps = NULL){
  
  # load default parameter file
  load('input/InitialConditions.RData')
  
  # copy default parameter file in CGEM format if no new parameters found
  if(is.null(inps)){
    
    out <- forminps(InitialConditions)
   
  # replace row values in input file with input list
  } else {
    
    # format parm names in input list for matching with new parm names
    # must remove regex metacharacters 
    inp_nms <- gsub('\\+|\\(|\\)', '', InitialConditions$parm)
    
    # replace each parameter with new
    for(inp in names(inps)){
      
      # index in defaults to replace
      sel <- gsub('\\+|\\(|\\)', '', inp) %>% 
        paste0('^', .) %>% 
        grep(., inp_nms)
      
      # stop if name is not matched in defaults
      if(length(sel) == 0) stop('No matches for initial conditions in ', inp)
      if(length(sel) > 1){
        mtchs <- paste(GEM_InputFile[sel, 'parm'], collapse = ', ')
        stop('Multiple matches for inital conditions in ', inp, ': ', mtchs)
      }
      
      InitialConditions[sel, 'value'] <- as.character(inps[[inp]])
      
    }
    
    # format new file for export
    out <- forminps(InitialConditions)
    
  }
  
  # save output to file
  writeLines(out, 'input/InitialConditions.txt')
  
}

# create parameter file for input with new parameters, uses formpars to change defaults
# 
# pars is a named list where each element is one to many parameter values for each parameter
# partial string matching is used for the element names to replace values in a default parameter list
# passing NULL to pars will return the default parameter list for Weeks Bay
setpars <- function(pars = NULL){

  # load default parameter file
  load('input/GEM_InputFile.RData')
  
  # copy default parameter file in CGEM format if no new parameters found
  if(is.null(pars)){

    out <- formpars(GEM_InputFile)
   
  # replace row values in input file with input list
  } else {

    # format parm names in input list for matching with new parm names
    # must remove regex metacharacters 
    par_nms <- gsub('\\+|\\(|\\)', '', GEM_InputFile$parm)
    
    # replace each parameter with new
    for(par in names(pars)){
      
      # index in defaults to replace
      sel <- gsub('\\+|\\(|\\)', '', par) %>% 
        paste0('^', .) %>% 
        grep(., par_nms)
      
      # stop if name is not matched in defaults
      if(length(sel) == 0) stop('No parameter matches for ', par)
      if(length(sel) > 1){
        mtchs <- paste(GEM_InputFile[sel, 'parm'], collapse = ', ')
        stop('Multiple parameter matches for ', par, ': ', mtchs)
      }
      
      GEM_InputFile[sel, 'value'] <- as.character(pars[[par]])
      
    }
    
    # format new file for export
    out <- formpars(GEM_InputFile)
    
  }
  
  # save output to file
  writeLines(out, 'input/GEM_InputFile')
  
}

# get chr string of parameter names that can be replaced
showpars <- function(){
  
  fl <- 'input/GEM_InputFile.RData'
  if(!file.exists(fl))stop(fl, ' does not exist')
  
  load(file = fl)
  
  out <- GEM_InputFile$parm
  
  return(out)
  
}
  
# import the initial conditions in a useful format
# only for display, not running the model
import_input <- function(){
  
  # import the file, rename columns for merge
  dat <- read.table(paste0('input/InitialConditions.txt'))
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[, 2] <- gsub('^!', '', dat[, 2])
  names(dat) <- c('strt', 'shrt')
  
  # get names from labels_fun, merge with input conditions
  out <- data.frame(Variable = labels_fun()$lng, Value = labels_fun()$val, shrt = labels_fun()$shrt)
  out <- merge(out, dat, by = 'shrt', all.x = TRUE)

  # return
  return(out)

}

# change files in tree to run FishTank with one phyto and one zoop group
# changes the files data/Model_dim.txt, InitialConditions.txt, and GEM_InputFile 
# last two files are in root, created on the fly with run_mod from input file so no need to convert back if to = F
#
# to logical indicating if changes are made from default (six phyto, two zoop) to one each, set to F to go from one group each back to default
p1z1_swtch <- function(to = TRUE){
 
  # forward change
  if(to){
    
    # model dimensions file
    mod_dim <- readLines('data/FishTank/Model_dim.txt')
    phytsel <- grep('nospA', mod_dim)
    zoopsel <- grep('nospZ', mod_dim)
    mod_dim[phytsel] <- gsub('^[6]', '1', mod_dim[phytsel])
    mod_dim[zoopsel] <- gsub('^[2]', '1', mod_dim[zoopsel])
    writeLines(mod_dim, 'data/FishTank/Model_dim.txt')
    
    # InitialConditions, remove 2-6 phyto and 2 zoop
    inits <- readLines('data/FishTank/InitialConditions.txt')
    rm <- grep('A[2-6]$|n[2-6]$|p[2-6]$|G[2]$', inits)
    inits <- inits[-rm]
    writeLines(inits, 'data/FishTank/InitialConditions.txt')
    
    # GEM_InputFile
    gem_inp <- readLines('GEM_InputFile')
    sel <- grep('ediblevector\\(Z[2]\\)$', gem_inp)
    gem_inp <- gem_inp[-sel]
    writeLines(gem_inp, 'GEM_InputFile')
  
  # backward change
  # GEM_InputFile in root does not need to be changed back, overwritten every time model is run
  } else {

    # model dimensions file
    mod_dim <- readLines('data/FishTank/Model_dim.txt')
    phytsel <- grep('nospA', mod_dim)
    zoopsel <- grep('nospZ', mod_dim)
    mod_dim[phytsel] <- gsub('^[1]', '6', mod_dim[phytsel])
    mod_dim[zoopsel] <- gsub('^[1]', '2', mod_dim[zoopsel])
    writeLines(mod_dim, 'data/FishTank/Model_dim.txt')
    
  }
 
  return(NULL)
  
}

# run the model
# copies intial conditions and parameter values from input
# executes model
# formats output to return data frame
#
# out_var is chr string of variable to return
# p1z1 logical if only one phyto and one zoop group are used, passed to p1z1_swtch
#
# any existing netcdf files are deleted before running the model
# the netcdf file is also deleted on function exit
# function will exit with error if model fails
run_mod <- function(pars = NULL, inps = NULL, out_var = 'O2', p1z1 = FALSE){

  # create parameter file based on inputs
  setpars(pars)
  
  # create initial conditions file based on inputs
  setinps(inps)

  # move the input files from input to root
  fls <- c('input/InitialConditions.txt', 'input/GEM_InputFile')
  file.copy(fls[1], 'data/FishTank', overwrite = TRUE)
  file.copy(fls[2], getwd(), overwrite = TRUE)
  
  # set to one phyto, one zoop group if TRUE
  if(p1z1) p1z1_swtch(to = TRUE)

  # delete any previous netcdf files
  fl <- 'NETCDF/output.000000.nc'
  if(file.exists(fl))
    file.remove('NETCDF/output.000000.nc')
  
  # run model, suppress output messages
  system('./FishTank.exe > /dev/null')

  # exist if it didn't work
  if(!file.exists(fl)){
    stop('FishTank broke')  
  }

  # back to default six phyt, two zoop if TRUE
  if(p1z1) p1z1_swtch(to = FALSE)
  
  # remove temp files
  file.remove(c('input/GEM_InputFile', 'input/InitialConditions.txt'))
  
  # get model output from netcdf
  nc <- nc_open('NETCDF/output.000000.nc')
  
  # return all if out_var is NULL
  if(is.null(out_var)) return(nc)

  # get time, add output variable
  time <- ncvar_get(nc, "time")
  
  # output as list
  out_ls <- vector('list', length = length(out_var))
  names(out_ls) <- out_var
  for(ov in out_var){

    # get all relevant organic matter variables for OM1 or OM2
    if(ov %in% c('OM1', 'OM2')){
  
      # need to sum the OM variables
      out_var_OM <- paste(ov, c('A', 'BC', 'R', 'Z'), sep = '_')
      sel <- names(nc$var) %in% out_var_OM
      if(sum(sel) != 4) stop('Check hard-coded organic variables')
      
      # get each variable, take rowsums
      var <- sapply(out_var_OM, function(x) ncvar_get(nc, x)) %>% 
        rowSums
  
    } else {
  
      sel <- names(nc$var) %in% ov
      if(sum(sel) == 0) stop(paste(ov, 'not found in model output'))
      
      # get variable
      var <- ncvar_get(nc, ov)
      
    }
  
    # combine with time
    out <- data.frame(time, var)
    names(out) <- c('time', ov)
    
    # append to list
    out_ls[[ov]] <- out
    
  }
  
  # make sure to close the connection, delete file
  nc_close(nc)
  if(file.exists(fl))
    file.remove(fl)

  return(out_ls)
  
}

# formatting the labels from labels_fun as expressions for plots
expr_fun <- function(lab_in){
 
  sel <- which(labels_fun()$shrt == lab_in)
  val <- labels_fun()$vals[sel]
  
#   if(grepl('-', val)){
#     val <- strsplit(val, '-')[[1]]
#     val <- bquote(.(val[1]) ^ .(paste0('-', val[2])))
#   }
    
  return(val)
   
}

# sensitivity analysis for individual parameters, different approach from FEM::sensFun
#
# pars is chr string of parameter name(s)
# per is perturbation factor, as 1 + per multiplied by original default value
# out_var chr string of output variable to evaluate
# p1z1 logical passed to run_mod for one phyto, one zoop group
# inps list of initial conditions, passed to run_mod
sensfun <- function(pars, per = 0.5, out_var = 'O2', trace = TRUE, p1z1 = FALSE, inps = NULL){
  
  # load lookup table
  load('input/GEM_InputFile.RData')
  
  # format parm names in inputs, doesn't matter which file
  # must remove regex metacharacters 
  par_nms <- gsub('\\+|\\(|\\)', '', GEM_InputFile$parm)
    
  # get default for each 
  pardef <- parper <- vector('numeric', length = length(pars))
  names(pardef) <- names(parper) <- pars

  # find default parameter values for each par in the lookup file
  for(par in pars){
    
    # index in defaults to replace
    sel <- gsub('\\+|\\(|\\)', '', par) %>% 
      paste0('^', .) %>% 
      grep(., par_nms)
    
    # stop if name is not matched in defaults
    if(length(sel) == 0) stop('No parameter matches for ', par)
    if(length(sel) > 1){
      mtchs <- paste(GEM_InputFile[sel, 'parm'], collapse = ', ')
      stop('Multiple parameter matches for ', par, ': ', mtchs)
    }
    
    pardef[par] <- GEM_InputFile[sel, 'value'] %>% 
      as.numeric
    
    # if the default parameter value is zero, no perturbation
    if(pardef[par] == 0) {
      
      parper[par] <- pardef[par]
      
    # otherwise apply perturbation      
    } else {
      
      parper[par] <- pardef[par] * (1 + per)
    
    }
    
  }
  
  # reference model  
  resref <- run_mod(out_var = out_var, p1z1 = p1z1, inps = inps)

  # iterate through each perturbed parameter and get results
  parper <- as.list(parper)
  ressij <- resper <- vector('list', length = length(pars))
  names(ressij) <- names(resper) <- pars
  for(par in names(parper)){
    
    if(trace) cat(par, '\t')

    # run the model with changed parameter
    ests <- run_mod(pars = parper[par], out_var = out_var, p1z1 = p1z1, inps = inps)
    
    # get new estimates from the perturbation
    ests <- lapply(ests, function(x){
        x <- select(x, -time) %>% 
          unlist
        
        names(x) <- NULL
        
        return(x)
        
      })
  
    # get sij for each variables from ests
    sij <- sapply(names(ests), function(x){
      
      ests_tmp <- ests[[x]]
      ref_tmp <- resref[[x]][, x]
      
      out <- (ests_tmp - ref_tmp)/(parper[[par]] - pardef[[par]]) * (pardef[[par]])/ests_tmp
      out[is.infinite(out)] <- NaN
      out
      
    }, simplify = F)

    resper[[par]] <- ests
    ressij[[par]] <- sij
    
  }
  
  # go through each state variable, summarize sensitivity data
  out_ls <- vector('list', length = length(out_var))
  names(out_ls) <- out_var
  for(ov in out_var){

    # raw output from perturbation 
    resper_tmp <- lapply(resper, function(x) x[[ov]])
    raw <- c(list('time' = resref[[ov]]$time, 'default' = resref[[ov]][, ov]), resper_tmp) %>% 
      do.call('cbind', .) %>% 
      as.data.frame %>% 
      tidyr::gather('Parameter', 'Estimate', -time, -default)
  
    # sij values
    ressij_tmp <- lapply(ressij, function(x) x[[ov]])
    sij <- do.call('cbind', ressij_tmp)
    nms <- dimnames(sij)[[2]]
    sij <- data.frame(time = resref[[ov]]$time, sij)
    names(sij) <- c('time', nms)
    sij <- tidyr::gather(sij, 'Parameter', 'sij', -time)

    # sensitivity from sij
    sens <- group_by(sij, Parameter) %>% 
      summarize(
        L1 = sum(abs(sij), na.rm = TRUE)/length(sij), 
        L2 = sqrt(sum(sij^2, na.rm = TRUE)/length(sij)), 
        Min = min(sij, na.rm = TRUE), 
        Max = max(sij, na.rm = TRUE), 
        N = length(sij)
      ) %>%  
      data.frame(., stringsAsFactors = F)
      
    # output
    out_ls[[ov]] <- list('raw' = raw, 'sij' = sij, 'sens' = sens)
    
  }
  
  return(out_ls)
  
}

######
# format parameter inputs for shiny app, input is reactive object from shiny (user inputs from ui), output is list to send to setpars
# react_ls is created from reactiveValuesToList(input)
form_parinps <- function(react_ls){
  
  react_ls <- reactiveValuesToList(react_ls)
  
  # format times argument separately
  times <- react_ls$times
  times <- list(
    '- starting time.*1' = as.numeric(format(times[1], '%Y')),
    '- starting time.*2' = as.numeric(format(times[1], '%m')),
    '- starting time.*3' = as.numeric(format(times[1], '%d')),
    '- ending   time.*1' = as.numeric(format(times[2], '%Y')),
    '- ending   time.*2' = as.numeric(format(times[2], '%m')),
    '- ending   time.*3' = as.numeric(format(times[2], '%d'))
    )
   
  # remove inputs that are not parameters, make sure this works
 
  # format parm names in input list for matching with new parm names
  # must remove regex metacharacters 
  inps <- gsub('\\+|\\(|\\)', '', names(react_ls))
  load('input/GEM_InputFile.RData')
  torm <- gsub('\\+|\\(|\\)*', '', GEM_InputFile$parm) %>% 
    paste0('^', .) %>% 
    paste(., collapse = '|') %>% 
    grep(., inps, value = T, invert = T)
  
  out_ls <- react_ls[!names(react_ls) %in% torm]
  out_ls <- c(out_ls, times)
  
  return(out_ls)

}

######
# format initial condition inputs for shiny app, input is reactive object from shiny (user inputs from ui), output is list to send to setpars
# react_ls is created from reactiveValuesToList(input)
form_iniinps <- function(react_ls){
  
  react_ls <- reactiveValuesToList(react_ls)

  # select the initial condition values from the list
  load('input/InitialConditions.RData')
  tosel <- InitialConditions$parm
  out_ls <- react_ls[tosel]
  
  return(out_ls)

}

######
# return parameter values by categories
parcats <- function(){
  
  load('input/GEM_InputFile.RData')

  raw <- formpars(GEM_InputFile)
  
  cats <- c('Optics', 'Temperature', 'Phytoplankton', 'Zooplankton', 'Organic Matter')
  
  out_ls <- vector('list', length = length(cats))
  names(out_ls) <- cats
  for(cat in cats){
   
    # get parameter block for the category 
    strt <- grep(cat, raw)
    blck <- raw[strt:length(raw)]
    ends <- grep('^$', blck)[1] - 1
    blck <- blck[2:ends]
    
    # get unique parameter names from blck
    parms <- gsub('^.*\\t!', '', blck)
    
    # get parameter values from GEM_InputFile
    sels <- gsub('_[0-9]$', '', GEM_InputFile$parm)
    sels <- sels %in% parms
    out <- as.numeric(GEM_InputFile[sels, 'value'])
    names(out) <- GEM_InputFile[sels, 'parm']
    out <- as.list(out)
   
    # append to output list
    out_ls[[cat]] <- out
    
  }

  return(out_ls)
  
}

######
# parameter identifiability based on residuals
# similar to FME::collin, a lot of hacked from https://github.com/cran/FME/blob/master/R/collin.R
# requires sens object from sensfun
#
# sens_in input sens object with sensitivity matrices
# pars chr string of parameters for evaluating identifiability
# N maximum number of parameters to select from 2:N from parameter str, defaults to all, large N is floored at the number of parameters in pars
# maxN logical indicating if only the maximum N is evaluated for identifiability, otherwise all combinations from 2:N will be evaluated
# maxcomb numeric for upper limit of combinations to evaluate for a given subset of parameters, testing all combinations is prohibitive for large parameter lists, maxcomb is ignored if larger than the number of unique combinations for a given parameter subset
# trace logical for printed progress on console
# shrt logical if output is truncated
par_ident <- function(sens_in, pars = NULL, N = NULL, maxN = FALSE, maxcomb = 5000, trace = FALSE, shrt = FALSE, ...){

  # get pars
  if(is.null(pars)){
    pars <- unique(sens_in$sij$Parameter)
  } else {
    if(length(pars) < 2) 
      stop('Not enough parameters to evaluate')
  }

  # subset parameters, get residuals
  toeval <- tidyr::spread(sens_in$sij, 'Parameter', 'sij') %>% 
    as.matrix

  # remove first row if all NaN
  if(sum(is.na(toeval[1, ])) == length(toeval[1, -1]))
    toeval <- toeval[-1, ]

  # get max number of parameters to eval
  if(is.null(N)){
   
     N <- length(pars)
  
  # use N if provided, but floor at max number of parameters if less than N
  } else {
 
     N <- pmin(N, length(pars))
     
  }
    
  # normalize the residual matrix
  sensest <- sqrt(colSums(toeval^2, na.rm = TRUE))
  toeval_nrm <- t(t(toeval) / sensest)
  
  # get parameter combinations to eval 
  # get all combos in increasing selections if maxN is T
  if(!maxN){
    combs <- combs_fun(pars = pars, N = N, maxcomb = maxcomb, trace = trace)
  } else {
    combs <- list(pars)
  }
  
  # output matrix to allocate
    
  # dims
  np <- length(pars)
  nc <- length(combs)
  collout <- matrix(0, ncol = np + 2, nrow = nc) 
  
  # counter
  if(trace){
    counts <- round(seq(1, length(combs), length = 20))
    cat('Percent estimated..\n\n')
  }

  # iterate through different combination levels
  for(i in seq_along(combs)){
    
    if(trace){
      perc <- 5 * which(i == counts)
      if(length(perc) != 0) cat(perc, '\t') 
    }
    
    # subset matrix with comb, get collinearity
    comb <- as.character(combs[[i]])
    mat <- toeval_nrm[, comb]
    
    # get collinearity, taken from collFun in FME
    # symmetric matrix and prevents errors if negative eigenv (floors at zero)
    coll <- 1 / sqrt(max(0, tail(eigen(crossprod(mat), only.values = TRUE, symmetric = TRUE)$value, 1)))
    # coll <- 1 / sqrt(min(eigen(t(mat) %*% mat)$value))
      
    # prep output for collout
    n  <- length(comb)
    whichp <- rep(0, np)
    whichp[pars %in% comb] <- 1
    collout[i, ] <- c(whichp, n, coll)
    
  }

  if(trace) cat('\n\n')
  
  # prep final output
  collout <- as.data.frame(collout) 
  names(collout) <- c(pars, 'N', 'coll')
  
  if(shrt) collout <- collout[, c('N', 'coll')]
  
  return(collout)
  
}

######
# retrieve list of parameter combinations for evalauting identifiability, used in par_ident
# the function retrieves every unique combintation from sel in 2 to N parameters in the complete parameter set
# the argument maxcomb places an upper limit for any i to prevent memory allocation issues with large parameter sets
# function returns a list with many elements, further processed in par_ident
# 
# pars chr string of parameters for evaluating identifiability
# N maximum number of parameters to select from 2:N from parameter str, defaults to all
# maxcomb numeric for upper limit of combinations to evaluate for a given subset of parameters, set to NULL to prevent this behavior
combs_fun <- function(pars, N, maxcomb, trace = TRUE){
  
  if(trace) cat('\nSelecting parameter combinations...\n\n')
  
  # get parameter combinations to eval
  combs <- vector('list', length = N - 1)
  for(sel in 2:N){

    if(trace) cat(sel, '\t')
    
    # get all unique combos of parameters if no upper limit
    if(is.null(maxcomb)){
      
      # get unique combinations to test, from min to max combs  
      combtmp <- combn(pars, sel, simplify = FALSE)
  
    # otherwise subsample for very large datasets
    } else {

      # number of unique combinations given sel
      unicomb <- choose(length(pars), sel)        
      
      # take all combs if unicomb < maxcomb
      if(unicomb < maxcomb){
        
        # get unique combinations to test, from min to max combs  
        combtmp <- combn(pars, sel, simplify = FALSE)
        
      } else {
      
        # floor maxcomb if > unicomb
        maxcomb <- pmin(maxcomb, unicomb)
        
        # initialize while
        unis <- 1
  
        # get random permutations until number of uniques are greater than or equal to maxcomb
        grw <- list()
        while(unis < maxcomb){
  
          # take sample size of N, repeat maxcomb times, sort, get unique
          smps <- sapply(1:maxcomb, function(x){ sample(x = pars, size = sel, replace = FALSE)}, simplify = FALSE)
          smps <- lapply(smps, sort) %>%
            unique
          
          # append to output, get unique again, check length
          grw <- c(grw, smps)
          grw <- unique(grw)
          unis <- length(grw)
          
        }
   
        # truncate unique by maxcomb, unless number of uniques are less than maxcomb
        combtmp <- grw[1:maxcomb] 
        
      }
  
    }
    
    # add to output
    combs[[sel - 1]] <- combtmp
    
  }
    
  if(trace) cat('\n\n')
  
  combs <- do.call('c', combs)
  return(combs)
  
}
######
# get collinearity matrix (pairwise) for all parameters in a category
# input for ggcorr in GGally
#
# cat chr string of input category (optics, organic, phytoplankton, zooplankton)
# sens_fl chr string of path to file with sensitivity estimates
# sens_cat_fl chr string of path to file with sensitivity estimates by parameter category
# st chr string of state variable to evaluate
collmat <- function(cat, sens_fl = 'rdata/sens_ests_all.RData', sens_cat_fl = 'rdata/sens_ests_cat_all.RData', st = 'O2'){
  
  # load data, assign to arbitrary object
  load(file = sens_fl)
  load(file = sens_cat_fl)
  sens_fl <- get(gsub('^rdata/|\\.RData$', '', sens_fl))
  sens_cat_fl <- get(gsub('^rdata/|\\.RData$', '', sens_cat_fl))
  
  # subset by state variables
  sens_fl <- sens_fl[[st]]
  sens_cat_fl <- sens_cat_fl[[st]]
  
  # subset parameters by category
  if(length(cat) == 1){
    
    pars <- filter(sens_cat_fl, Category == cat) %>% 
      .$Parameter %>% 
      as.character
    
  # otherwise user-supplied    
  } else {
  
    pars <- cat
     
  }

  # get identifiability
  ident <- par_ident(sens_fl, pars, N = 2, maxcomb = NULL)
  
  # format a vector as symmetric matrix for ggcorr
  ident_mat <- ident$coll
  class(ident_mat) <- 'dist'
  attr(ident_mat,'Size') <- length(pars)
  ident_mat <- as.matrix(ident_mat)
  diag(ident_mat) <- 1
  dimnames(ident_mat) <- list(pars, pars)
  
  return(ident_mat)
   
}

######
# get legend from an existing ggplot object
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# alternative list structure for parameter categories
parcats2 <- function(as_df = FALSE){
  
  Optics = list(
    cats = 'Optics',
    shrt = c('Kw_1', 
      'Kcdom_1', 
      'Kspm_1', 
      'Kchla_1',
      'astar490_1',
      'aw490_1',
      'astarOMA_1',
      'astarOMZ_1',
      'astarOMR_1',
      'astarOMbC_1',
      'PARfac_1',
      'sink CDOM_1'
      ),
    lngs = c(
      'Kw: AOP, light attenuation due to water',
      'Kcdom: AOP, light attenuation due to CDOM',
      'Kspm: AOP, light attenuation due to SPM', 
      'Kchla: AOP, light attenuation due to chla ',
      'astar490: Chla specific absorption at 490 nm',
      'aw490: seawater absorption at 490 nm',
      'astarOMA: OM_A specific absorption at 490 nm',
      'astarOMZ: OM_Z specific absorption at 490 nm',
      'astarOMR: OM_R specific absorption at 490 nm',
      'astarOMBC: OM_BC specific absorption at 490 nm',
      'PARfac: Multiplies surface PAR', 
      'sink CDOM: sinking rate'
      ),
    vals = c(
      0.146,
      0.001,
      0.029,
      0.024,
      0.0375,
      0.015,
      0.1,
      0.1,
      0.1,
      0.1,
      1,
      0
    )
  )
  
  Temperature = list(
    cats = 'Temperature',
    shrt = c(
      'Tref(nospA+nospZ)_1',
      'Tref(nospA+nospZ)_2',
      'Tref(nospA+nospZ)_3',
      'Tref(nospA+nospZ)_4',
      'Tref(nospA+nospZ)_5',
      'Tref(nospA+nospZ)_6',
      'Tref(nospA+nospZ)_7',
      'Tref(nospA+nospZ)_8',
      'KTg1(nospA+nospZ)_1',
      'KTg1(nospA+nospZ)_2',
      'KTg1(nospA+nospZ)_3',
      'KTg1(nospA+nospZ)_4',
      'KTg1(nospA+nospZ)_5',
      'KTg1(nospA+nospZ)_6',
      'KTg1(nospA+nospZ)_7',
      'KTg1(nospA+nospZ)_8',
      'KTg2(nospA+nospZ)_1',
      'KTg2(nospA+nospZ)_2',
      'KTg2(nospA+nospZ)_3',
      'KTg2(nospA+nospZ)_4',
      'KTg2(nospA+nospZ)_5',
      'KTg2(nospA+nospZ)_6',
      'KTg2(nospA+nospZ)_7',
      'KTg2(nospA+nospZ)_8',
      'Ea_R(nospA+nospZ)_1',
      'Ea_R(nospA+nospZ)_2',
      'Ea_R(nospA+nospZ)_3',
      'Ea_R(nospA+nospZ)_4',
      'Ea_R(nospA+nospZ)_5',
      'Ea_R(nospA+nospZ)_6',
      'Ea_R(nospA+nospZ)_7',
      'Ea_R(nospA+nospZ)_8'
      ),
    lngs = c(
      'Tref(nospA+nospZ): Optimum temperature for growth(C)',
      'Tref(nospA+nospZ): Optimum temperature for growth(C)',
      'Tref(nospA+nospZ): Optimum temperature for growth(C)',
      'Tref(nospA+nospZ): Optimum temperature for growth(C)',
      'Tref(nospA+nospZ): Optimum temperature for growth(C)',
      'Tref(nospA+nospZ): Optimum temperature for growth(C)',
      'Tref(nospA+nospZ): Optimum temperature for growth(C)',
      'Tref(nospA+nospZ): Optimum temperature for growth(C)',
      'KTg1(nospA+nospZ): Effect of T below Topt(C^2)',
      'KTg1(nospA+nospZ): Effect of T below Topt(C^2)',
      'KTg1(nospA+nospZ): Effect of T below Topt(C^2)',
      'KTg1(nospA+nospZ): Effect of T below Topt(C^2)',
      'KTg1(nospA+nospZ): Effect of T below Topt(C^2)',
      'KTg1(nospA+nospZ): Effect of T below Topt(C^2)',
      'KTg1(nospA+nospZ): Effect of T below Topt(C^2)',
      'KTg1(nospA+nospZ): Effect of T below Topt(C^2)',
      'KTg2(nospA+nospZ): Effect of T above Topt(C^2)',
      'KTg2(nospA+nospZ): Effect of T above Topt(C^2)',
      'KTg2(nospA+nospZ): Effect of T above Topt(C^2)',
      'KTg2(nospA+nospZ): Effect of T above Topt(C^2)',
      'KTg2(nospA+nospZ): Effect of T above Topt(C^2)',
      'KTg2(nospA+nospZ): Effect of T above Topt(C^2)',
      'KTg2(nospA+nospZ): Effect of T above Topt(C^2)',
      'KTg2(nospA+nospZ): Effect of T above Topt(C^2)',
      'Ea_R(nospA+nospZ): Slope of Arrhenius plot',
      'Ea_R(nospA+nospZ): Slope of Arrhenius plot',
      'Ea_R(nospA+nospZ): Slope of Arrhenius plot',
      'Ea_R(nospA+nospZ): Slope of Arrhenius plot',
      'Ea_R(nospA+nospZ): Slope of Arrhenius plot',
      'Ea_R(nospA+nospZ): Slope of Arrhenius plot',
      'Ea_R(nospA+nospZ): Slope of Arrhenius plot',
      'Ea_R(nospA+nospZ): Slope of Arrhenius plot'
    ),
    vals = c(
      22, 
      22, 
      22, 
      22, 
      22, 
      22, 
      22,
      22,
      0.0035,
      0.0035,
      0.0035,
      0.0035,
      0.0035,
      0.0035,
      0.0035,
      0.0035,
      0.001,
      0.001,
      0.001,
      0.001,
      0.001,
      0.001,
      0.001,
      0.001,
      10000,
      10000,
      10000,
      10000,
      10000,
      10000,
      10000, 
      10000
      )
  )
  
  Phytoplankton = list(
    cats = 'Phytoplankton',
    shrt = c(
      'ediblevector(Z1)_1',
      'ediblevector(Z1)_2',
      'ediblevector(Z1)_3',
      'ediblevector(Z1)_4',
      'ediblevector(Z1)_5',
      'ediblevector(Z1)_6',
      'ediblevector(Z2)_1',
      'ediblevector(Z2)_2',
      'ediblevector(Z2)_3',
      'ediblevector(Z2)_4',
      'ediblevector(Z2)_5',
      'ediblevector(Z2)_6',
      'umax_1',
      'umax_2',
      'umax_3',
      'umax_4',
      'umax_5',
      'umax_6',
      'alpha_1',
      'alpha_2',
      'alpha_3',
      'alpha_4',
      'alpha_5',
      'alpha_6',
      'beta_1',
      'beta_2',
      'beta_3',
      'beta_4',
      'beta_5',
      'beta_6',
      'respg_1',
      'respg_2',
      'respg_3',
      'respg_4',
      'respg_5',
      'respg_6',
      'respb_1',
      'respb_2',
      'respb_3',
      'respb_4',
      'respb_5',
      'respb_6',
      'QminN_1',
      'QminN_2',
      'QminN_3',
      'QminN_4',
      'QminN_5',
      'QminN_6',
      'QminP_1',
      'QminP_2',
      'QminP_3',
      'QminP_4',
      'QminP_5',
      'QminP_6',
      'QmaxN_1',
      'QmaxN_2',
      'QmaxN_3',
      'QmaxN_4',
      'QmaxN_5',
      'QmaxN_6',
      'QmaxP_1',
      'QmaxP_2',
      'QmaxP_3',
      'QmaxP_4',
      'QmaxP_5',
      'QmaxP_6',
      'Kn_1',
      'Kn_2',
      'Kn_3',
      'Kn_4',
      'Kn_5',
      'Kn_6',
      'Kp_1',
      'Kp_2',
      'Kp_3',
      'Kp_4',
      'Kp_5',
      'Kp_6',
      'Ksi_1',
      'Ksi_2',
      'Ksi_3',
      'Ksi_4',
      'Ksi_5',
      'Ksi_6',
      'KQn_1',
      'KQn_2',
      'KQn_3',
      'KQn_4',
      'KQn_5',
      'KQn_6',
      'KQp_1',
      'KQp_2',
      'KQp_3',
      'KQp_4',
      'KQp_5',
      'KQp_6',
      'nfQs_1',
      'nfQs_2',
      'nfQs_3',
      'nfQs_4',
      'nfQs_5',
      'nfQs_6',
      'vmaxN_1',
      'vmaxN_2',
      'vmaxN_3',
      'vmaxN_4',
      'vmaxN_5',
      'vmaxN_6',
      'vmaxP_1',
      'vmaxP_2',
      'vmaxP_3',
      'vmaxP_4',
      'vmaxP_5',
      'vmaxP_6',
      'vmaxSi_1',
      'vmaxSi_2',
      'vmaxSi_3',
      'vmaxSi_4',
      'vmaxSi_5',
      'vmaxSi_6',
      'aN_1',
      'aN_2',
      'aN_3',
      'aN_4',
      'aN_5',
      'aN_6',
      'volcell_1',
      'volcell_2',
      'volcell_3',
      'volcell_4',
      'volcell_5',
      'volcell_6',
      'Qc_1',
      'Qc_2',
      'Qc_3',
      'Qc_4',
      'Qc_5',
      'Qc_6',
      'Athresh_1',
      'Athresh_2',
      'Athresh_3',
      'Athresh_4',
      'Athresh_5',
      'Athresh_6',
      'sink A_1',
      'sink A_2',
      'sink A_3',
      'sink A_4',
      'sink A_5',
      'sink A_6',
      'mA_1',
      'mA_2',
      'mA_3',
      'mA_4',
      'mA_5',
      'mA_6', 
      'A_wt_1',
      'A_wt_2',
      'A_wt_3',
      'A_wt_4', 
      'A_wt_5', 
      'A_wt_6'
      ),
    lngs = c(
      "ediblevector(Z1): edibility vector for Z1", 
      "ediblevector(Z1): edibility vector for Z1", 
      "ediblevector(Z1): edibility vector for Z1", 
      "ediblevector(Z1): edibility vector for Z1", 
      "ediblevector(Z1): edibility vector for Z1", 
      "ediblevector(Z1): edibility vector for Z1", 
      "ediblevector(Z2): edibility vector for Z2", 
      "ediblevector(Z2): edibility vector for Z2", 
      "ediblevector(Z2): edibility vector for Z2", 
      "ediblevector(Z2): edibility vector for Z2", 
      "ediblevector(Z2): edibility vector for Z2", 
      "ediblevector(Z2): edibility vector for Z2", 
      "umax: maximum growth rate", 
      "umax: maximum growth rate",
      "umax: maximum growth rate", 
      "umax: maximum growth rate", 
      "umax: maximum growth rate", 
      "umax: maximum growth rate", 
      "alpha: initial slope of the photosynthesis-irradiance relationship", 
      "alpha: initial slope of the photosynthesis-irradiance relationship", 
      "alpha: initial slope of the photosynthesis-irradiance relationship", 
      "alpha: initial slope of the photosynthesis-irradiance relationship", 
      "alpha: initial slope of the photosynthesis-irradiance relationship", 
      "alpha: initial slope of the photosynthesis-irradiance relationship", 
      "beta: photoinhibition constant", 
      "beta: photoinhibition constant", 
      "beta: photoinhibition constant",
      "beta: photoinhibition constant", 
      "beta: photoinhibition constant",
      "beta: photoinhibition constant", 
      "respg: phytoplankton growth respiration coefficient",
      "respg: phytoplankton growth respiration coefficient", 
      "respg: phytoplankton growth respiration coefficient",
      "respg: phytoplankton growth respiration coefficient", 
      "respg: phytoplankton growth respiration coefficient",
      "respg: phytoplankton growth respiration coefficient", 
      "respb: phytoplankton basal respiration coefficient",
      "respb: phytoplankton basal respiration coefficient", 
      "respb: phytoplankton basal respiration coefficient",
      "respb: phytoplankton basal respiration coefficient", 
      "respb: phytoplankton basal respiration coefficient",
      "respb: phytoplankton basal respiration coefficient", 
      "QminN: minimum N cell-quota",
      "QminN: minimum N cell-quota", 
      "QminN: minimum N cell-quota",
      "QminN: minimum N cell-quota", 
      "QminN: minimum N cell-quota",
      "QminN: minimum N cell-quota", 
      "QminP: minimum P cell-quota",
      "QminP: minimum P cell-quota", 
      "QminP: minimum P cell-quota",
      "QminP: minimum P cell-quota", 
      "QminP: minimum P cell-quota",
      "QminP: minimum P cell-quota", 
      "QmaxN: maximum N cell-quota",
      "QmaxN: maximum N cell-quota", 
      "QmaxN: maximum N cell-quota", 
      "QmaxN: maximum N cell-quota", 
      "QmaxN: maximum N cell-quota", 
      "QmaxN: maximum N cell-quota", 
      "QmaxP: maximum P cell-quota",
      "QmaxP: maximum P cell-quota", 
      "QmaxP: maximum P cell-quota", 
      "QmaxP: maximum P cell-quota", 
      "QmaxP: maximum P cell-quota",
      "QmaxP: maximum P cell-quota", 
      "Kn: half-saturation constant for N",
      "Kn: half-saturation constant for N", 
      "Kn: half-saturation constant for N",
      "Kn: half-saturation constant for N", 
      "Kn: half-saturation constant for N",
      "Kn: half-saturation constant for N", 
      "Kp: half-saturation constant for P", 
      "Kp: half-saturation constant for P", 
      "Kp: half-saturation constant for P",
      "Kp: half-saturation constant for P", 
      "Kp: half-saturation constant for P", 
      "Kp: half-saturation constant for P", 
      "Ksi: half-saturation constant for Si uptake", 
      "Ksi: half-saturation constant for Si uptake", 
      "Ksi: half-saturation constant for Si uptake",
      "Ksi: half-saturation constant for Si uptake", 
      "Ksi: half-saturation constant for Si uptake", 
      "Ksi: half-saturation constant for Si uptake", 
      "KQn: Qn constant for Flynn nutrient dependent growth model", 
      "KQn: Qn constant for Flynn nutrient dependent growth model", 
      "KQn: Qn constant for Flynn nutrient dependent growth model", 
      "KQn: Qn constant for Flynn nutrient dependent growth model", 
      "KQn: Qn constant for Flynn nutrient dependent growth model", 
      "KQn: Qn constant for Flynn nutrient dependent growth model", 
      "KQp: Qp constant for Flynn nutrient dependent growth model", 
      "KQp: Qp constant for Flynn nutrient dependent growth model", 
      "KQp: Qp constant for Flynn nutrient dependent growth model", 
      "KQp: Qp constant for Flynn nutrient dependent growth model", 
      "KQp: Qp constant for Flynn nutrient dependent growth model", 
      "KQp: Qp constant for Flynn nutrient dependent growth model", 
      "nfQs: exponent for Geider nutrient uptake model",
      "nfQs: exponent for Geider nutrient uptake model", 
      "nfQs: exponent for Geider nutrient uptake model",
      "nfQs: exponent for Geider nutrient uptake model", 
      "nfQs: exponent for Geider nutrient uptake model",
      "nfQs: exponent for Geider nutrient uptake model", 
      "vmaxN: N-uptake rate measured at umax", 
      "vmaxN: N-uptake rate measured at umax", 
      "vmaxN: N-uptake rate measured at umax", 
      "vmaxN: N-uptake rate measured at umax", 
      "vmaxN: N-uptake rate measured at umax",
      "vmaxN: N-uptake rate measured at umax", 
      "vmaxP: P-uptake rate measured at umax", 
      "vmaxP: P-uptake rate measured at umax", 
      "vmaxP: P-uptake rate measured at umax",
      "vmaxP: P-uptake rate measured at umax", 
      "vmaxP: P-uptake rate measured at umax",
      "vmaxP: P-uptake rate measured at umax", 
      "vmaxSi: Si-uptake rate measured at umax",
      "vmaxSi: Si-uptake rate measured at umax", 
      "vmaxSi: Si-uptake rate measured at umax",
      "vmaxSi: Si-uptake rate measured at umax", 
      "vmaxSi: Si-uptake rate measured at umax",
      "vmaxSi: Si-uptake rate measured at umax", 
      "aN: coefficient for non-limiting nutrient",
      "aN: coefficient for non-limiting nutrient", 
      "aN: coefficient for non-limiting nutrient",
      "aN: coefficient for non-limiting nutrient", 
      "aN: coefficient for non-limiting nutrient",
      "aN: coefficient for non-limiting nutrient", 
      "volcell: phytoplankton volume/cell", 
      "volcell: phytoplankton volume/cell", 
      "volcell: phytoplankton volume/cell", 
      "volcell: phytoplankton volume/cell", 
      "volcell: phytoplankton volume/cell",
      "volcell: phytoplankton volume/cell", 
      "Qc: phytoplankton carbon/cell",
      "Qc: phytoplankton carbon/cell", 
      "Qc: phytoplankton carbon/cell",
      "Qc: phytoplankton carbon/cell", 
      "Qc: phytoplankton carbon/cell", 
      "Qc: phytoplankton carbon/cell", 
      "Athresh: Phytoplankton threshold for grazing, is multiplied by VOLcell", 
      "Athresh: Phytoplankton threshold for grazing, is multiplied by VOLcell", 
      "Athresh: Phytoplankton threshold for grazing, is multiplied by VOLcell", 
      "Athresh: Phytoplankton threshold for grazing, is multiplied by VOLcell", 
      "Athresh: Phytoplankton threshold for grazing, is multiplied by VOLcell", 
      "Athresh: Phytoplankton threshold for grazing, is multiplied by VOLcell", 
      "sink A: sinking rate of phytoplankton cells",
      "sink A: sinking rate of phytoplankton cells", 
      "sink A: sinking rate of phytoplankton cells",
      "sink A: sinking rate of phytoplankton cells", 
      "sink A: sinking rate of phytoplankton cells",
      "sink A: sinking rate of phytoplankton cells", 
      "mA: mortality coefficient", 
      "mA: mortality coefficient",
      "mA: mortality coefficient", 
      "mA: mortality coefficient", 
      "mA: mortality coefficient", 
      "mA: mortality coefficient",
      "A_wt: Relative proportion of total chlA for initializing phytoplankton",
      "A_wt: Relative proportion of total chlA for initializing phytoplankton",
      "A_wt: Relative proportion of total chlA for initializing phytoplankton",
      "A_wt: Relative proportion of total chlA for initializing phytoplankton",
      "A_wt: Relative proportion of total chlA for initializing phytoplankton",
      "A_wt: Relative proportion of total chlA for initializing phytoplankton"
      ),
    vals = c(
      0.25,
      0.25,
      0.25, 
      0.25,
      0.25,
      0.25, 
      0.25,
      0.25, 
      0.25,
      0.25, 
      0.25, 
      0.25, 
      0.41,
      0.41,
      0.41,
      0.41,
      0.41, 
      0.41,
      8.42e-17, 
      8.42e-17, 
      8.42e-17, 
      8.42e-17, 
      8.42e-17, 
      8.42e-17, 
      1.1e-18, 
      1.1e-18, 
      1.1e-18, 
      1.1e-18, 
      1.1e-18, 
      1.1e-18, 
      0.1, 
      0.1, 
      0.1, 
      0.1, 
      0.1, 
      0.1, 
      0.02,
      0.02,
      0.02,
      0.02,
      0.02,
      0.02,
      6.08e-09,
      6.08e-09,
      6.08e-09,
      6.08e-09,
      6.08e-09,
      6.08e-09,
      6.19e-10,
      6.19e-10, 
      6.19e-10,
      6.19e-10,
      6.19e-10,
      6.19e-10,
      2.04e-07,
      2.04e-07,
      2.04e-07,
      2.04e-07,
      2.04e-07,
      2.04e-07,
      1.28e-08,
      1.28e-08,
      1.28e-08,
      1.28e-08,
      1.28e-08,
      1.28e-08,
      4.51,
      4.51,
      4.51,
      4.51,
      4.51,
      4.51,
      2.86,
      2.86,
      2.86,
      2.86,
      2.86,
      2.86,
      4.51,
      4.51,
      4.51,
      4.51,
      4.51,
      4.51,
      5,
      5,
      5,
      5,
      5,
      5,
      0.2,
      0.2,
      0.2,
      0.2,
      0.2,
      0.2,
      1,
      1,
      1,
      1,
      1,
      1,
      4.10e-08,
      4.10e-08,
      4.10e-08,
      4.10e-08,
      4.10e-08,
      4.10e-08,
      2.68e-08,
      2.68e-08,
      2.68e-08,
      2.68e-08,
      2.68e-08,
      2.68e-08,
      4.10e-08,
      4.10e-08,
      4.10e-08,
      4.10e-08,
      4.10e-08,
      4.10e-08,
      1,
      1,
      1,
      1,
      1,
      1,
      33693,
      33693,
      33693,
      33693,
      33693,
      33693,
      1.35e-06,
      1.35e-06,
      1.35e-06,
      1.35e-06,
      1.35e-06,
      1.35e-06,
      1.721e8,
      1.721e8,
      1.721e8, 
      1.721e8,
      1.721e8,
      1.721e8,
      1.49,
      1.49,
      1.49,
      1.49,
      1.49,
      1.49,
      0.1,
      0.1,
      0.1,
      0.1,
      0.1,
      0.1,
      1,
      1,
      1,
      1,
      1,
      1
      )
  )
  
  Zooplankton = list(
    cats = 'Zooplankton',
    shrt = c(
      'Zeffic_1',
      'Zeffic_2',
      'Zslop_1',
      'Zslop_2',
      'Zvolcell_1',
      'Zvolcell_2',
      'ZQc_1', 
      'ZQc_2',
      'ZQn_1',
      'ZQn_2',
      'ZQp_1',
      'ZQp_2',
      'ZKa_1',
      'ZKa_2',
      'Zrespg_1',
      'Zrespg_2',
      'Zrespb_1',
      'Zrespb_2',
      'Zumax_1',
      'Zumax_2',
      'Zm_1',
      'Zm_2'
      ), 
    lngs = c(
      'Zeffic: assimilation efficiency as a fraction of ingestion',
      'Zeffic: assimilation efficiency as a fraction of ingestion',
      'Zslop: proportion of grazed phytoplankton lost to sloppy feeding',
      'Zslop: proportion of grazed phytoplankton lost to sloppy feeding',
      'Zvolcell: zooplankton volume/individual',
      'Zvolcell: zooplankton volume/individual',
      'ZQc: zooplankton carbon/individual',
      'ZQc: zooplankton carbon/individual',
      'ZQn: zooplankton nitrogen/individual',
      'ZQn: zooplankton nitrogen/individual',
      'ZQp: zooplankton phosphorus/individual',
      'ZQp: zooplankton phosphorus/individual',
      'ZKa: half saturation coefficient for grazing',
      'ZKa: half saturation coefficient for grazing',
      'Zrespg: Zooplankton growth-dependent respiration factor',
      'Zrespg: Zooplankton growth-dependent respiration factor',
      'Zrespb: Zooplankton biomass-dependent respiration factor',
      'Zrespb: Zooplankton biomass-dependent respiration factor',
      'Zumax: maximum growth rate of zooplankton',
      'Zumax: maximum growth rate of zooplankton',
      'Zm: Zooplankton mortality constant for quadratic mortality',
      'Zm: Zooplankton mortality constant for quadratic mortality'
      ), 
    vals = c(
      0.4,
      0.4,
      0.25,
      0.25,
      2.98e7,
      2.98e7,
      3.13e-4,
      3.13e-4,
      6.95e-05,
      6.95e-05,
      3.77e-06,
      3.77e-06,
      1.12e12,
      1.12e12,
      0.2,
      0.2,
      0.1,
      0.1,
      9.45e7,
      9.45e7,
      0.00072,
      0.00072
      )
  )
  
  `Organic Matter` = list(
    cats = 'Organic Matter',
    shrt = c(
      'KG1_1',
      'KG2_1',
      'KG1_R_1',
      'KG2_R_1',
      'KG1_BC_1',
      'KG2_BC_1',
      'KNH4_1',
      'nitmax_1',
      'KO2_1',
      'KstarO2_1',
      'KNO3_1',
      'pCO2_1',
      'stoich_x1R_1',
      'stoich_y1R_1',
      'stoich_x2R_1',
      'stoich_y2R_1',
      'stoich_x1BC_1',
      'stoich_y1BC_1',
      'stoich_x2BC_1',
      'stoich_y2BC_1',
      'sink OM1_A_1',
      'sink OM2_A_1',
      'sink OM1_Z_1',
      'sink OM2_Z_1',
      'sink OM1_R_1',
      'sink OM2_R_1',
      'sink OM1_BC_1',
      'sink OM2_BC_1',
      'KGcdom_1',
      'CF_SPM_1'
      ), 
    lngs = c(
      'KG1: turnover rate for OM1_A and OM1_G',	
      'KG2: turnover rate for OM2_A and OM2_G',
      'KG1_R: OM1 turnover rate for riverine',
      'KG2_R: OM2 turnover rate for riverine',
      'KG1_BC: OM1 turnover rate for initial and bc',
      'KG2_BC: OM2 turnover rate for initial and bc',
      'KNH4: NH4 rate constant for nitrification',
      'nitmax: maximum rate of nitrification per day',
      'KO2: half-saturation concentration for O2 utilization',
      'KstarO2: O2 concentration that inhibits denitrification',
      'KNO3: half-saturation concentration for NO3 used in denitrification',
      'pCO2: atmospheric CO2',
      'stoich_x1R:  C:P stoichiometry of OM1_R',
      'stoich_y1R:  N:P stoichiometry of OM1_R',
      'stoich_x2R:  C:P stoichiometry of OM2_R',
      'stoich_y2R:  N:P stoichiometry of OM2_R',
      'stoich_x1BC: C:P stoichiometry of OM1_BC',
      'stoich_y1BC: N:P stoichiometry of OM1_BC',
      'stoich_x2BC: C:P stoichiometry of OM2_BC',
      'stoich_y2BC: N:P stoichiometry of OM2_BC',
      'sink OM1_A:  sinking rate',
      'sink OM2_A:  sinking rate',
      'sink OM1_Z:  sinking rate',
      'sink OM2_Z:  sinking rate',
      'sink OM1_R:  sinking rate',
      'sink OM2_R:  sinking rate',
      'sink OM1_BC: sinking rate',
      'sink OM2_BC: sinking rate',
      'KGcdom: decay rate of CDOM, 1/day',
      'CF_SPM: conversion factor for river OM to river SPM'
      ), 
    vals = c(
      50,
      50,
      11,
      3.7,
      1,
      1,
      1,
      0.52,
      10,
      10,
      10,
      380,
      51,
      4.5,
      700,
      50,
      106,
      16,
      106,
      16,
      10,
      0,
      10,
      0,
      10,
      0,
      10,
      0,
      0.01,
      0.018
      )
  )
  
  # list
  out <- list(Optics = Optics, Temperature = Temperature, Phytoplankton = Phytoplankton, Zooplankton = Zooplankton, `Organic Matter` = `Organic Matter`)
  
  # return as data frame if T
  if(as_df){
    
    out <- lapply(out, data.frame)
    out <- do.call('rbind', out)
    row.names(out) <- 1:nrow(out)
    
    return(out)
  }
  
  return(out)
  
}

######
# convert  parameter names to latex format
#
# parin chr vector of short names to convert
# frm chr string indicating format of output, tex or exp for latex or expression
par_txt <- function(parin, frm = 'tex'){

  # sanity check
  if(!frm %in% c('tex', 'exp'))
    stop('frm argument must be "tex" or "exp"')

  parin <- as.character(parin) 
  
  # all parameters and names
  cats <- parcats2(as_df = TRUE)[, c('cats', 'shrt')]

  # get which row the parameter is in
  sels <- which(cats$shrt %in% parin)
  if(length(sels) != length(unique(parin))) stop('parin not completely matched in shrt')
  
  # split the names by category
  splits <- cats[sels, ] %>% 
    .[match(parin, .$shrt), ] # this is important to make sure the output order matches with input
  
  # gsub the shrt names differnet by category
  subs <- apply(splits, 1,  function(x){

    # 1-6 are phytos, 7-8 are zoops (changed to 1-2)
    if('Temperature' %in% x['cats']){

      x['shrt'] <- gsub('_([1-6])$', '_p\\1', x['shrt'])
      x['shrt'] <- gsub('_[7]$', '_z1', x['shrt']) 
      x['shrt'] <- gsub('_[8]$', '_z2', x['shrt'])  
      
    }
      
    # add p to subscript
    if('Phytoplankton' %in% x['cats']){
      
      x['shrt'] <- gsub('_([1-9])$', '_p\\1', x['shrt']) 
      
    }

    # add z to subscript
    if('Zooplankton' %in% x['cats']){
     
      x['shrt'] <- gsub('_([1-9])$', '_z\\1', x['shrt']) 
      
    }
    
    # remove subscript
    if(any(c('Optics', 'Organic Matter') %in% x['cats'])){
     
      x['shrt'] <- gsub('_[1-9]$', '', x['shrt'])
       
    }
        
    return(x)
    
  }) %>% 
  t %>% 
  data.frame(., stringsAsFactors = FALSE) %>% 
  .$shrt

  # convert output format for tex
  if(frm == 'tex'){
    
    out <- gsub('_([pz][1-9])$', '$_{\\1}$', subs)
    out <- paste0('\\textit{', out, '}')
      
  }
  
  # convert output format as expressions for R
  if(frm == 'exp'){
  
    out <- gsub('_([pz][1-9])$', '[italic(\\1)]', subs)
    out <- paste0('italic(', out, ')')
    out <- gsub('\\((.*)\\s(.*)\\)', '("\\1 \\2")', out)
    out <- parse(text = as.expression(out))
    
  }
    
  return(out)
  
}

######
# this is a convenience function to return the relevant parameters for sensitivity analysis
par_tst <- function(){
  
  # parameters to eval
  pars <- lapply(parcats(), names)
  pars <- unlist(pars)
  
  # remove irrelevant parameters, taken from 'ignore/noparams.csv'
  rmpars <- c('Kw_1', 'Kcdom_1', 'Kspm_1', 'Kchla_1', 'aw490_1', 'astarOMR_1', 'astarOMBC_1', 'PARfac_1', 'KTg1(nospA+nospZ)_1', 'KTg1(nospA+nospZ)_2', 'KTg1(nospA+nospZ)_3', 'KTg1(nospA+nospZ)_4', 'KTg1(nospA+nospZ)_5', 'KTg1(nospA+nospZ)_6', 'KTg1(nospA+nospZ)_7', 'KTg1(nospA+nospZ)_8', 'KTg2(nospA+nospZ)_1', 'KTg2(nospA+nospZ)_2', 'KTg2(nospA+nospZ)_3', 'KTg2(nospA+nospZ)_4', 'KTg2(nospA+nospZ)_5', 'KTg2(nospA+nospZ)_6', 'KTg2(nospA+nospZ)_7', 'KTg2(nospA+nospZ)_8', 'Ea_R(nospA+nospZ)_1', 'Ea_R(nospA+nospZ)_2', 'Ea_R(nospA+nospZ)_3', 'Ea_R(nospA+nospZ)_4', 'Ea_R(nospA+nospZ)_5', 'Ea_R(nospA+nospZ)_6', 'Ea_R(nospA+nospZ)_7', 'Ea_R(nospA+nospZ)_8', 'beta_1', 'beta_2', 'beta_3', 'beta_4', 'beta_5', 'beta_6', 'QmaxN_1', 'QmaxN_2', 'QmaxN_3', 'QmaxN_4', 'QmaxN_5', 'QmaxN_6', 'QmaxP_1', 'QmaxP_2', 'QmaxP_3', 'QmaxP_4', 'QmaxP_5', 'QmaxP_6', 'nfQs_1', 'nfQs_2', 'nfQs_3', 'nfQs_4', 'nfQs_5', 'nfQs_6', 'A_wt_1', 'A_wt_2', 'A_wt_3', 'A_wt_4', 'A_wt_5', 'A_wt_6', 'KG1_R_1', 'KG2_R_1', 'KG1_BC_1', 'KG2_BC_1', 'pCO2_1', 'stoich_x1R_1', 'stoich_y1R_1', 'stoich_x2R_1', 'stoich_y2R_1', 'stoich_x1BC_1', 'stoich_y1BC_1', 'stoich_x2BC_1', 'stoich_y2BC_1', 'sink OM1_R_1', 'sink OM2_R_1', 'sink OM1_BC_1', 'sink OM2_BC_1', 'CF_SPM_1')
  rmpars <- c(rmpars, 'respb_7') # there's an extra of this one
  pars <- pars[!pars %in% rmpars]
  pars <- grep('_1$|_7$', pars, value = TRUE) # get parameters for one phyto, one zoop (_7 is for temp) 
  
  return(pars)

}

######
# get list of parameter combinations to calibrate from tocal_all data object
#
# tocal_all input data object created in dat_proc.R
# out_var chr string of state variable to eval
# thrsh numeric for collinearity thresholds
# coll logical if gamma is returned as an additional element for each heuristic
# errs logical if errors for each parameter in the set are returned for each heuristic
# pert numeric optional value to fix min/max ranges of parameters for optimization as plus/minus the parameter times the proportion in pert, minv/maxv values from rngs are used if not supplied
# 
# temperature is manually added as a category in the first heuristic
# only one parameter so wasn't included in tocal_all output
#
# min/max values for ediblevector and tref were added, 
# former form Eldridge and Roelke, latter from reasonable temp ranges
get_cmbs <- function(tocal_all, out_var = 'O2', thrsh = 15, coll = TRUE, errs = TRUE, pert = 0.5){
  
  # expected parameter values and ranges
  allp <- get(load(file = 'input/GEM_InputFile.RData'))
  rngs <- get(load(file = 'rdata/rngs.RData'))
  
  # get out_var results, melt, create heurist column, filter by thresh, 
  # add temp category (one parameter), combine with input file to get values
  idpars <- tocal_all[[out_var]] %>% 
    reshape2::melt(id.vars = names(.[[1]])) %>% 
    mutate(
      Category = ifelse(!L1 %in% 'cat', '', Category),
      L1 = ifelse(Category %in% '', L1, ''), 
      parm = as.character(Parameter)
      ) %>% 
    unite('heurist', L1, Category, sep = '') %>% 
    filter(coll < thrsh) %>%
    group_by(heurist) %>% 
    mutate(coll = max(coll, na.rm = TRUE)) %>% 
    ungroup %>% 
    select(parm, error, heurist, coll) %>% 
    rbind(
      data.frame(parm = 'Tref(nospA+nospZ)_1', error = 0.02, heurist = 'Temperature', coll = NA)
      ) %>% 
    left_join(., allp, by = 'parm') %>% 
    rename(vals = value)
  
  # get expected ranges
  # add missing values
  rngs <- filter(rngs, parm %in% idpars$parm) %>% 
    separate(value, c('est', 'minv', 'maxv'), sep = '/') %>% 
    mutate(
      minv = ifelse(parm %in% 'Tref(nospA+nospZ)_1', 0, 
        ifelse(parm %in% 'ediblevector(Z1)_1', 0.05, 
          minv)),
      maxv = ifelse(parm %in% 'Tref(nospA+nospZ)_1', 40, 
        ifelse(parm %in% 'ediblevector(Z1)_1', 1, 
          maxv))
      ) %>% 
    gather('var', 'val', -parm) %>% 
    mutate(val = as.numeric(val)) %>% 
    spread(var, val) %>% 
    select(-est) 
  
  # combine both, create list of lists by heuristic and val/min/max
  out <- left_join(idpars, rngs, by = 'parm') %>% 
    split(.$heurist) %>% 
    lapply(., function(x){
      
      # collinearity value
      collv <- unique(x$coll)
      
      # error valus
      errsv <- x$error
      names(errsv) <- x$parm
      
      # change minv maxv to +- pert if not null
      if(!is.null(pert)){
        x$minv <- as.numeric(x$vals) - pert * as.numeric(x$vals)
        x$maxv <- as.numeric(x$vals) + pert * as.numeric(x$vals)
      }
      
      # get value, min, max of each parameter
      x <- select(x, -heurist) %>% 
        gather('var', 'val', -parm, -coll) %>% 
        split(.$var) %>% 
        lapply(., function(y){
          
          nms <- y$parm
          out <- as.list(as.numeric(y$val))
          names(out) <- nms
          
          return(out)
          
        })
      
      # add gamma to list for heurist
      if(coll)
        x$coll <- collv
      
      # add individual errors to list
      if(errs)
        x$errs <- errsv
      
      return(x)

    })
  
  # make sure parameter values in range of min/max
  chks <- sapply(names(out), function(x){
      
      prms <- out[[x]]
      vals <- unlist(prms$vals)
      minv <- unlist(prms$minv)
      maxv <- unlist(prms$maxv)
      
      chk <- minv > vals | maxv < vals
      if(any(chk)) warning('parameter values for ', x , ' not in expected range')
      
    })
  
  return(out)
  
}

######
# optimization function
#
# ... additional arguments to optim
#
# pars list of parameters to optimize
# minv list of minimum expected values for parameters
# maxv list of maximum expected values for parameters
#
# see get_cmbs
fishopt <- function(pars, minv, maxv, ...){
  
  ##
  # parameter scales
  # parscl <- apply(do.call('cbind', pars), 1, function(x) diff(range(x))) 		
  control <- c(list(parscale = unlist(pars)), list(...))
  
  ##
  # observed data for calibration
  load(file = 'rdata/flaskdo.RData')
  obs <- flaskdo
  
  ##
  # matrix for estimates, list for model run results
  mat <- matrix(nrow = 0, ncol = length(pars) + 1)
  res <- list()
  
  ##
  # error function used in fishopt
  #
  # pars is list of parameters passed to fisherr from optim
  # obs is observed values for comparison of predictions
  #
  fisherr <- function(pars = NULL, obs){
    
    cat(pars, '\n')
    
    # change start, end times
    # make system completely closed
    parsdts <- c(
      pars,
      '^- starting.*_1$' = 2013,
      '^- starting.*_2$' = 9,
      '^- starting.*_3$' = 25,
      '^- starting.*_4$' = 8,
      '^- starting.*_4$' = 40,
      '^- ending.*_1$' = 2013,
      '^- ending.*_2$' = 9,
      '^- ending.*_3$' = 25,
      '^- ending.*_4$' = 19,
      'ReadVars_1' = 1,
      'ReadVars_3' = 1,
      '- dT (timestep, seconds); dT_out (output interval, seconds)_1' = 300, 
      '- dT (timestep, seconds); dT_out (output interval, seconds)_2' = 600, 
      'Which_fluxes_1' = 0, # O2 exchange off
      'Which_fluxes_2' = 0  # CO2 exchange off
      )
    
    # add column in obs to add modelled data
    obs$modO2 <- NA
    
    # iterate through treatments
    # restmp is for intermediate model predictions
    tmt_rep <- unique(obs$tmt_rep)
    restmp <- vector('list', length = length(tmt_rep))
    names(restmp) <- tmt_rep
    for(i in tmt_rep){
      
      ##
      # initial conditions, spam 9/25/13
      
      # same between treatments, reps
      inps <- list(
        'CDOM' = 30,
        'Si' = 38.4
      )

      # these change depending on treatment
      
      # O2
      inps$O2 <- obs[obs$tmt_rep %in% i, 'sttO2'] %>% 
        as.numeric
      
      # nh4
      inps$NH4 <- obs[obs$tmt_rep %in% i, 'NH4'] %>% 
        as.numeric
      
      # no3
      inps$NO3 <- obs[obs$tmt_rep %in% i, 'NO3'] %>% 
        as.numeric
      
      # po4
      inps$PO4 <- obs[obs$tmt_rep %in% i, 'PO4'] %>% 
        as.numeric
      
      ## set input PAR and temp files    
      load(file = 'rdata/flaskhobo.RData')
      exp_inp(flaskhobo, getvar = 'par', gettmt = gsub('_.*$', '', i))
      exp_inp(flaskhobo, getvar = 'temp', gettmt = gsub('_.*$', '', i))
      
      ## run the model with parameters
      # sometimes fishtank crashes, so try again until it works
      p <- try({run_mod(parsdts, inps = inps, p1z1 = T)[[1]]})
      while(class(p) == 'try-error'){
        cat('mod error, running again\n')
        p <- try({run_mod(parsdts, inps = inps, p1z1 = T)[[1]]})
      }
    
      # convert time to posix
      tocmp <- mutate(p, 
         time = as.POSIXct(time + (60 * 60 * 6), origin = c('2002-01-01'), tz = 'America/Regina')
        )
 
      # find end time in obs to grab in tocmp
      endtm <- obs %>% 
        filter(tmt_rep %in% i) %>% 
        .$end
      
      # get modeled O2 at end time
      modO2 <- tocmp %>% 
        filter(time %in% endtm) %>% 
        .$O2
  
      # append to output
      obs[obs$tmt_rep %in% i, 'modO2'] <- modO2
      restmp[[i]] <- tocmp
      
    }

    # get errors, append to output
    err <- with(obs, sqrt(mean((endO2 - modO2)^2)))
    
    cat('\n\t', err, '\n\n')
    mat <<- rbind(mat, c(unlist(pars), err))
    
    # get model results for each treatment, append to output
    restmp <- restmp %>% 
      enframe %>% 
      unnest
    res <<- c(res, list(restmp))
    
    return(err)
    
  }
    
  # format output
  opt <- optim(
    pars, 
    fisherr, 
    obs = obs, 
    method = 'L-BFGS-B', 
    lower = as.numeric(minv), 
    upper = as.numeric(maxv), 
    control = control
    )
  opt$mat <- mat
  opt$res <- res
  
  return(opt)
  
}

######
#
# setup temperature and solar inputs if provided
#
# flaskhobo data frame with temperature, solar data for three treatments on a given date
# var variable to create, temp or lt
exp_inp <- function(flaskhobo, getvar = c('temp', 'par'), gettmt = c('lt', 'dk', 'ltnt')){
  
  # get arguments
  getvar <- match.arg(getvar)
  gettmt <- match.arg(gettmt)
  
  # pull var and tmt data from flaskhobo
  # format times
  fl <- flaskhobo %>% 
    filter(var %in% getvar & tmt %in% gettmt) %>% 
    mutate(
      iYr = year(datetime),
      iMon = month(datetime),
      iDay = day(datetime), 
      iHour = hour(datetime),
      iMin = minute(datetime), 
      iSec = second(datetime), 
      Var1 = val
    ) %>% 
    select(iYr, iMon, iDay, iHour, iMin, iSec, Var1)

  # switch getvar for correct name
  getvar <- switch(getvar,
                  temp = 'Temp', 
                  par = 'Solar'
    )
  
  # output directory and write
  outdir <- paste0('data/FishTank/INPUT/', getvar, '.dat')
  write.table(fl, outdir, sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
  
}
