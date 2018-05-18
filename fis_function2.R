FIS <- function(layers) {
  
  
  ##### Gather parameters and layers #####
  ### * ram_b_bmsy, ram_f_fmsy, ram_catch

  status_yr_span <- c(2001:2016)
  
  ram_b_bmsy      <- read_csv("output/ram_b_bmsy.csv") %>%
    rename(b_bmsy = value) %>%
    select(-param)
  ram_f_fmsy      <- read_csv("output/ram_f_fmsy.csv") %>%
    rename(f_fmsy = value) %>%
    select(stock_id, year, f_fmsy)
  dfo_catch        <- read_csv("output/dfo_catch.csv")
  # rgn_stock_area  <- read_csv("output/ram_stock_area.csv") %>%
  #   rename(region_id = rgn_id, stock_id = stockid)
  # #get the prop of catch in each rgn for groundfish too
  # rgn_gfish_prop  <- read_csv("data/2_dfo_groundfish_catch_prop.csv") %>%
  #   rename(region_id = rgn_id)
  
  ### These parameters are based on conversation with Ian Perry, Karen Hunter,
  ### and Karin Bodtker on May 24 2017.
  b_bmsy_underexploit_penalty <- 0.25
  b_bmsy_underexploit_thresh  <- 3.00
  f_fmsy_underfishing_penalty <- 0.25
  f_fmsy_overfishing_thresh   <- 2.00
  
  ### Apply rolling mean to F/Fmsy
  ram_f_fmsy <- ram_f_fmsy %>%
    mutate(f_fmsy_raw = f_fmsy) %>%
    arrange(stock_id, year) %>%
    group_by(stock_id) %>%
    mutate(f_fmsy = zoo::rollmean(f_fmsy_raw, k = 4, align = 'right', fill = NA)) %>%
    ungroup()
  
  stock_status_layers <- ram_b_bmsy %>%
    full_join(ram_f_fmsy, by = c('year', 'stock_id'))
  
  ########################################################.
  ##### run each fishery through the Kobe plot calcs #####
  ########################################################.
  ### * ram_b_bmsy, ram_f_fmsy
  
  
  ### Function for converting B/Bmsy values into a 0 - 1 score
  rescale_bprime_crit <- function(fish_stat_df,
                                  bmax, bmax_val) {
    
    ### parameter from DFO harvest control rule:
    overfished_th  <- 0.8
    ### parameter from OHI California Current:
    underfished_th <- 1.5
    
    bmax_adj <- (bmax - underfished_th) / (1 - bmax_val) + underfished_th
    ### this is used to create a "virtual" B/Bmsy max where score drops
    ### to zero.  If bmax_val == 0, this is bmax; if bmax_val > 0, bmax_adj
    ### extends beyond bmax, to create a gradient where bmax_val occurs at bmax.
    
    fish_stat_df <- fish_stat_df %>%
      # group_by(stock) %>% ### grouping by stock will set b_max by max per stock, instead of max overall
      mutate(b_max     = max(b_bmsy, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(bPrime = NA,
             bPrime = ifelse(b_bmsy < overfished_th,  ### overfished stock
                             b_bmsy / overfished_th,
                             bPrime),
             bPrime = ifelse(b_bmsy >= overfished_th & b_bmsy < underfished_th,
                             1,                       ### appropriately fished stock
                             bPrime),
             bPrime = ifelse(b_bmsy >= underfished_th,
                             (bmax_adj - b_bmsy) / (bmax_adj - underfished_th), ### underfished stock
                             bPrime),
             bPrime = ifelse(bPrime < 0, 0, bPrime))
    return(fish_stat_df)
  }
  
  
  ### Function to create vertical gradient based on distance from
  ### ideal F/Fmsy value to actual F/Fmsy value
  f_gradient <- function(f, over_f, under_f, fmax, fmin_val) {
    x <- ifelse(f < over_f & f > under_f, 1, NA)
    x <- ifelse(f <= under_f, (f * (1 - fmin_val) / under_f + fmin_val), x)
    x <- ifelse(f >= over_f,  (fmax - f) / (fmax - over_f), x)
    x <- ifelse(f > fmax, NA, x)
    return(x)
  }
  
  ### Function to convert F/Fmsy values into 0 - 1 score
  rescale_fprime_crit <- function(fish_stat_df,
                                  fmax, fmin_val) {
    
    ### params from DFO harvest control rule:
    Bcrit <- 0.4; overfished_th <- 0.8
    ### params from OHI California Current:
    underfishing_th <- 0.8; overfishing_th  <- 1.2
    
    bcritslope = 1 / (overfished_th - Bcrit)
    ### connecting from (Bcrit, 0) to (overfished_th, 1)
    
    fish_stat_df <- fish_stat_df %>%
      mutate(fPrime = ifelse(b_bmsy < overfished_th & f_fmsy < fmax,
                             f_gradient(f_fmsy + (overfished_th - b_bmsy) * bcritslope,
                                        over_f = overfishing_th,
                                        under_f = underfishing_th,
                                        fmax = fmax,
                                        fmin_val = fmin_val),
                             NA),
             fPrime = ifelse(b_bmsy >= overfished_th & f_fmsy < fmax,
                             f_gradient(f_fmsy,
                                        over_f = overfishing_th,
                                        under_f = underfishing_th,
                                        fmax = fmax,
                                        fmin_val = fmin_val),
                             fPrime),
             fPrime = ifelse(is.na(fPrime), 0, fPrime), ### fill zeros everywhere unscored
             fPrime = ifelse(is.na(f_fmsy), NA, fPrime) ### but if no f_fmsy, reset to NA
      )
    return(fish_stat_df)
  }
  
  stock_status_df <- stock_status_layers %>%
    rescale_bprime_crit(bmax     = b_bmsy_underexploit_thresh,
                        bmax_val = b_bmsy_underexploit_penalty) %>%
    rescale_fprime_crit(fmax     = f_fmsy_overfishing_thresh,
                        fmin_val = f_fmsy_underfishing_penalty) %>%
    mutate(x_prod = ifelse(!is.na(fPrime), (fPrime * bPrime), bPrime),
           basis  = ifelse(!is.na(fPrime), 'F_Fmsy, B_Bmsy', 'B_Bmsy only')) %>%
    dplyr::select(year, stock_id,
                  score = x_prod,
                  basis,
                  bPrime, fPrime,
                  b_bmsy, f_fmsy)  %>%
    group_by(stock_id) %>%
    complete_years(status_yr_span, method = 'carry', dir = 'forward') %>%
    ungroup()
  
  ## joing stock status data with dfo catch data
  
  stock_status_catch <- stock_status_df %>%
    full_join()
  
  #################################################################
  ##### identify regions where each stock is fished each year #####
  #################################################################
  
  ### Take the df that lists area proportion of each stock in each region (not including groundfish) and
  ### add all years to this df. Since we don't have year by year spatial distributions of catch, we have to
  ### assume that each stock is fished in overlapping regions for all years.
  
  rgns_stocks <- rgn_stock_area %>%
    select(region_id, stock_id)
  
  #groundfish data
  gfish <- rgn_gfish_prop %>%
    mutate(stock_id = case_when(
      species == "Hake" ~ "PHAKEPCOAST",
      species == "Halibut" ~ "PHALNPAC",
      species == "Sablefish" ~ "SABLEFPCAN",
      species == "Pacific Ocean Perch" & rgn_name == "West Vancouver Island" ~ "PERCHWCVANI",
      species == "Pacific Ocean Perch" & rgn_name == "Haida Gwaii" ~ "PERCHQCI"
    )) %>%
    filter(!is.na(stock_id)) %>%  ## This removes some rows of pacific ocean perch outside of wc van island and haida gwaii. B/c we dont have assessments for those areas...
    select(year, region_id, stock_id) %>%
    group_by(region_id, stock_id) %>%
    complete_rgn_years(status_yr_span, method = 'carry', dir = 'forward') %>%
    ungroup() %>%
    mutate(rgn = region_id) %>%
    select(-region_id)
  
  #dataframe of all stocks, years and their stock scores and regions where they are found
  stock_score_df <- stock_status_df %>%
    filter(!is.na(score)) %>%
    left_join(rgns_stocks) %>%
    left_join(gfish) %>%
    mutate(region_id = ifelse(is.na(rgn), region_id, rgn)) %>%
    select(-rgn)
  
  write.csv(stock_score_df, "output/stock_scores.csv")
  

  #the status is 
  fis_status <- stock_score_df %>%
    group_by(year, region_id) %>%
    summarize(status = mean(score)*100) %>%
    filter(year >= min(status_yr_span)) %>%
    mutate(goal = 'FIS',
           score = status,
           dimension = 'status')
  
  #save for plotting
  write_csv(fis_status, 'output/fis_status.csv')
  
  return(fis_status)
  
}