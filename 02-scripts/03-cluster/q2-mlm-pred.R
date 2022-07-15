##### Question 2: MLM ######

setwd("/projects/p31385/ps-oscillators")

library(rstan)
library(plyr)
library(tidyverse)
library(psych)
library(brms)
# library(ggquiver)

sessionInfo()

jobid = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))+1
print(jobid)
args <- read.table("scripts/cluster/args/q2-mlm.txt"
                   , header = T, stringsAsFactors = F)[jobid,]
print(args)


q2_pred_fun <- function(study, pname, pitem, sitem){
  study2 <- mapvalues(
    study
    , c("01-pairs", "02-ipcs", "03-sherman", "04-horstmann")
    , c("pairs", "ipcs", "sherman", "horstmann")
    , warn_missing = F
  )
  
  load(sprintf("05-results/02-q2/02-mlm/01-models/%s-%s-%s-%s.RData", study2, pname, pitem, sitem))
  fx <- broom.mixed::tidy(m)
  rx <- coef(m, probs = c(0.025, 0.975))[[1]] %>% array_tree(3) %>% 
    tibble(term = names(.), data = .) %>% 
    mutate(data = map(data, ~(.) %>% data.frame %>% 
                        rownames_to_column("SID"))) %>% 
    unnest(data) %>% 
    select(term, SID, estimate = Estimate, conf.low = Q2.5, conf.high = Q97.5)
  save(fx, rx, file = sprintf("05-results/02-q2/02-mlm/02-summary/%s-%s-%s-%s.RData", study2, pname, pitem, sitem))
  
  draws <- as_draws_df(m) %>% as_tibble()
  fx_draws <- draws %>% select(-contains("r_SID"), -contains("sd_"), -contains("cor_"), -contains("sigma"))
  rx_draws <- draws %>% 
    mutate(n = 1:n()) %>%
    select(starts_with("r_SID"), n) %>% 
    pivot_longer(
      names_to = "term"
      , values_to = "estimate"
      , cols = -n
    ) %>%
    mutate(term = str_remove_all(term, "r_SID__")) %>%
    separate("term", c("outcome", "term"), sep = "\\[") %>%
    separate("term", c("SID", "term"), sep = ",") %>%
    mutate(term = str_remove_all(term, "\\]")) %>%
    pivot_wider(
      names_from = c("outcome", "term")
      , names_sep = "__"
      , values_from = "estimate"
    ) %>%
    full_join(fx_draws %>% mutate(n = 1:n()) %>% select(n, contains("b_"))) %>%
    mutate(pdeltavalue__p_value_c = pdeltavalue__p_value_c + b_pdeltavalue_p_value_c
           , pdeltavalue__s_value_c = pdeltavalue__s_value_c + b_pdeltavalue_s_value_c
           , `pdeltavalue__p_value_c:s_value_c` = `pdeltavalue__p_value_c:s_value_c` + `b_pdeltavalue_p_value_c:s_value_c`
           , sdeltavalue__p_value_c = sdeltavalue__p_value_c + b_sdeltavalue_p_value_c
           , sdeltavalue__s_value_c = sdeltavalue__s_value_c + b_sdeltavalue_s_value_c
           , `sdeltavalue__p_value_c:s_value_c` = `sdeltavalue__p_value_c:s_value_c` + `b_sdeltavalue_p_value_c:s_value_c`) %>%
    select(-starts_with("b_"))
  save(fx_draws, rx_draws
       , file = sprintf("05-results/02-q2/02-mlm/05-draws/%s-%s-%s-%s.RData", study2, pname, pitem, sitem))
  
  # find attractor location by solving for x
  # linear               # quadratic
  # y = b0 + b1*x        x = (-b +/ sqrt(b2 - 4ac))/2a
  # 0 = b0 + b1*x        
  # -b0 = b1*x         
  # x = -b0/b1
  gr_sd <- m$data %>% 
    group_by(SID) %>% 
    summarize_at(vars(p_value_c, s_value_c), sd) %>% 
    ungroup() %>% 
    summarize_at(vars(p_value_c, s_value_c), mean) 
  fx_equil <- fx_draws %>% 
    as.data.frame() %>%
    # fixed estimates
    mutate(p__attr_loc_0 = -1*b_pdeltavalue_Intercept/b_pdeltavalue_p_value_c
           , p__attr_loc_p = -1*(b_pdeltavalue_Intercept + gr_sd$s_value_c*b_pdeltavalue_s_value_c)/(b_pdeltavalue_p_value_c + gr_sd$s_value_c * `b_pdeltavalue_p_value_c:s_value_c`)
           , p__attr_loc_m = -1*(b_pdeltavalue_Intercept + -1*gr_sd$s_value_c*b_pdeltavalue_s_value_c)/(b_pdeltavalue_p_value_c + -1*gr_sd$s_value_c* `b_pdeltavalue_p_value_c:s_value_c`)
           # delta situations at 0, 5, 10
           , s__attr_loc_0 = -1*b_sdeltavalue_Intercept/b_sdeltavalue_s_value_c
           , s__attr_loc_p = -1*(b_sdeltavalue_Intercept + gr_sd$p_value_c*b_sdeltavalue_p_value_c)/(b_sdeltavalue_s_value_c + gr_sd$p_value_c * `b_pdeltavalue_p_value_c:s_value_c`)
           , s__attr_loc_m = -1*(b_sdeltavalue_Intercept + -1*gr_sd$p_value_c*b_sdeltavalue_p_value_c)/(b_sdeltavalue_s_value_c + -1*gr_sd$p_value_c * `b_sdeltavalue_p_value_c:s_value_c`))
  
  # participant-level estimates
  ind_sd <- m$data %>% 
    group_by(SID) %>% 
    summarize_at(vars(p_value_c, s_value_c), sd) %>%
    ungroup() %>%
    mutate(SID = as.character(SID))
  
  rx_equil <- rx_draws %>%
    full_join(ind_sd) %>%
    # delta personality at 0, 5, 10
    mutate(p__attr_loc_0 = -1*0/pdeltavalue__p_value_c
           , p__attr_loc_p = -1*(0 + s_value_c*pdeltavalue__s_value_c)/(pdeltavalue__p_value_c + s_value_c * `pdeltavalue__p_value_c:s_value_c`)
           , p__attr_loc_m = -1*(0 + -1*s_value_c*pdeltavalue__s_value_c)/(pdeltavalue__p_value_c + -1*s_value_c * `pdeltavalue__p_value_c:s_value_c`)
           # delta situations at 0, 5, 10
           , s__attr_loc_0 = -1*0/sdeltavalue__s_value_c
           , s__attr_loc_p = -1*(0 + p_value_c*sdeltavalue__p_value_c)/(sdeltavalue__s_value_c + p_value_c * `pdeltavalue__p_value_c:s_value_c`)
           , s__attr_loc_m = -1*(0 + -1*p_value_c*sdeltavalue__p_value_c)/(sdeltavalue__s_value_c + -1*p_value_c * `sdeltavalue__p_value_c:s_value_c`))
  
  # find attractor strength
  equil <- fx_equil %>% 
    # fixed estimates
    mutate(
      n = 1:n()
      , SID = "overall"
      , p__attr_str_0 = b_pdeltavalue_p_value_c
      , p__attr_str_p = b_pdeltavalue_p_value_c + gr_sd$s_value_c*`b_pdeltavalue_p_value_c:s_value_c`
      , p__attr_str_m = b_pdeltavalue_p_value_c + -1*gr_sd$s_value_c*`b_pdeltavalue_p_value_c:s_value_c`
      , s__attr_str_0 = b_sdeltavalue_s_value_c
      , s__attr_str_p = b_sdeltavalue_s_value_c + gr_sd$p_value_c*`b_sdeltavalue_p_value_c:s_value_c`
      , s__attr_str_m = b_sdeltavalue_s_value_c + -1*gr_sd$p_value_c*`b_sdeltavalue_p_value_c:s_value_c`
    ) %>% 
    select(n, SID, contains("attr_loc"), contains("attr_str")) %>%
    full_join(
      # participant-level estimates
      rx_equil %>% 
        mutate(
          SID = as.character(SID)
          , p__attr_str0 = pdeltavalue__p_value_c
          , p__attr_str_p = pdeltavalue__p_value_c + s_value_c*`pdeltavalue__p_value_c:s_value_c`
          , p__attr_str_m = pdeltavalue__p_value_c + -1*s_value_c*`pdeltavalue__p_value_c:s_value_c`
          , s__attr_str0 = sdeltavalue__s_value_c
          , s__attr_str_p = sdeltavalue__s_value_c + p_value_c*`sdeltavalue__p_value_c:s_value_c`
          , s__attr_str_m = sdeltavalue__s_value_c + -1*p_value_c*`sdeltavalue__p_value_c:s_value_c`
        ) #%>%
      # select(n, SID, contains("attr_loc"), contains("attr_str"))
    )
  
  # summarize the posterior draws
  equil <- equil %>% 
    select(SID, contains("attr_loc"), contains("attr_str")) %>%
    group_by(SID) %>%
    nest() %>% 
    ungroup() %>%
    mutate(
      data = map(data
                 , ~(.) %>% 
                   posterior_summary %>%
                   as.data.frame() %>% 
                   rownames_to_column("term")
      )
    ) %>%
    unnest(data)
  save(equil, file = sprintf("05-results/02-q2/02-mlm/03-equil/%s-%s-%s-%s.RData", study, pname, pitem, sitem))
  
  # setup predictions for plotting
  # fixed estimates
  gr_sd <- m$data %>% 
    group_by(SID) %>% 
    summarize_at(vars(p_value_c, s_value_c), sd) %>% 
    ungroup() %>% 
    summarize_at(vars(p_value_c, s_value_c), mean) 
  
  frame <- crossing(p_value_c = seq(-2*gr_sd$p_value_c, 2*gr_sd$p_value_c, length.out = 25), s_value_c = seq(-2*gr_sd$s_value_c, 2*gr_sd$s_value_c, length.out = 25))
  fx_pred <- bind_cols(
    bind_rows(frame, frame)
    , fitted(m, newdata = frame, re_formula = NA) %>% array_tree(3) %>% ldply(., .id = "outcome")
  ) %>%
    pivot_wider(
      names_from = "outcome"
      , values_from = c("Estimate", "Est.Error", "Q2.5", "Q97.5")
    )
  
  crossing_fun <- function(d){
    crossing(
      p_value_c = seq(-2*d$p_value_c, 2*d$p_value_c, length.out = 25)
      , s_value_c = seq(-2*d$s_value_c, 2*d$s_value_c, length.out = 25)
    )
  }
  # participant-level estimates
  frame <- m$data %>% 
    group_by(SID) %>% 
    summarize_at(vars(p_value_c, s_value_c), sd) %>%
    ungroup() %>%
    mutate(SID = as.character(SID)) %>%
    group_by(SID) %>%
    nest() %>%
    ungroup() %>%
    mutate(data = map(data, crossing_fun)) %>%
    unnest(data)
  
  # frame <- crossing(SID = unique(m$data$SID), p_value_c = seq(-5,5,.25), s_value_c = seq(-5,5,.25))
  rx_pred <- bind_cols(
    bind_rows(frame, frame)
    , fitted(m, newdata = frame) %>% array_tree(3) %>% ldply(., .id = "outcome")
  ) %>%
    pivot_wider(
      names_from = "outcome"
      , values_from = c("Estimate", "Est.Error", "Q2.5", "Q97.5")
    )
  save(fx_pred, rx_pred
       , file = sprintf("05-results/02-q2/02-mlm/04-predictions/%s-%s-%s-%s.RData", study, pname, pitem, sitem))
  
  # clean up and return nothing
  rm(list = c("frame", "fx_pred", "m1", "m", "fx", "rx", "fx_draws", "rx_draws", "equil", "draws"))
  gc()
  return(T)
}

q2_pred_fun(args[,1], args[,2], args[,3], args[,4])