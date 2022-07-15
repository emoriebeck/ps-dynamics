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
args <- read.table("scripts/cluster/args/q2-mlm-uncentered.txt"
                   , header = T, stringsAsFactors = F)[jobid,]
print(args)


q2_mlm_model_fun <- function(study, pname, pitem, sitem){
  study2 <- mapvalues(
    study
    , c("01-pairs", "02-ipcs", "03-sherman", "04-horstmann")
    , c("pairs", "ipcs", "sherman", "horstmann")
    , warn_missing = F
  )
  
  # load the data 
  load(sprintf("03-data/%s/06-q2-mlm/%s-%s-%s.RData", study, pname, pitem, sitem))
  
  # setup the formula 
  f1 <- bf(p_delta_value ~ 1 + p_value + s_value + p_value:s_value + (1 + p_value + s_value + p_value:s_value | p | SID))
  f2 <- bf(s_delta_value ~ 1 + p_value + s_value + p_value:s_value + (1 + s_value + p_value + p_value:s_value | p | SID))
  
  # load the sample model 
  load("05-results/02-q2/03-mlm-uncentered/sample-mlm-linear.RData")
  
  # run the model using update for faster compilation 
  # m <- update(
  #   m1
  #   , newdata = d
  #   , iter = 2000
  #   , warmup = 1000
  #   , cores = 4
  # )
  m <- brm(
    mvbrmsformula(f1, f2) + set_rescor(TRUE)
    , data = d
    , iter = 2000
    , warmup = 1000
    , cores = 12
    , threads = threading(3)
    , backend = "cmdstanr"
    , chains = 4
  )
  save(m, file = sprintf("05-results/02-q2/03-mlm-uncentered/01-models/%s-%s-%s-%s.RData", study2, pname, pitem, sitem))
  
  fx <- broom.mixed::tidy(m)
  rx <- coef(m, probs = c(0.025, 0.975))[[1]] %>% array_tree(3) %>% 
    tibble(term = names(.), data = .) %>% 
    mutate(data = map(data, ~(.) %>% data.frame %>% 
                        rownames_to_column("SID"))) %>% 
    unnest(data) %>% 
    select(term, SID, estimate = Estimate, conf.low = Q2.5, conf.high = Q97.5)
  save(fx, rx, file = sprintf("05-results/02-q2/03-mlm-uncentered/02-summary/%s-%s-%s-%s.RData", study2, pname, pitem, sitem))
  
  draws <- as_draws_df(m) %>% as_tibble()
  fx_draws <- draws %>% select(-contains("r_SID"), -contains("sd_"), -contains("cor_"))
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
    mutate(term = str_replace_all(term, "s_value:p_value", "p_value:s_value")) %>%
    pivot_wider(
      names_from = c("outcome", "term")
      , names_sep = "__"
      , values_from = "estimate"
    ) %>%
    full_join(fx_draws %>% mutate(n = 1:n()) %>% select(n, contains("b_"))) %>%
    mutate(pdeltavalue__Intercept = pdeltavalue__Intercept + b_pdeltavalue_Intercept
           , pdeltavalue__p_value = pdeltavalue__p_value + b_pdeltavalue_p_value
           , pdeltavalue__s_value = pdeltavalue__s_value + b_pdeltavalue_s_value
           , `pdeltavalue__p_value:s_value` = `pdeltavalue__p_value:s_value` + `b_pdeltavalue_p_value:s_value`
           , sdeltavalue__Intercept = sdeltavalue__Intercept + b_sdeltavalue_Intercept
           , sdeltavalue__p_value = sdeltavalue__p_value + b_sdeltavalue_p_value
           , sdeltavalue__s_value = sdeltavalue__s_value + b_sdeltavalue_s_value
           , `sdeltavalue__p_value:s_value` = `sdeltavalue__p_value:s_value` + `b_sdeltavalue_p_value:s_value`) %>%
    select(-starts_with("b_"))
  save(fx_draws, rx_draws
       , file = sprintf("05-results/02-q2/03-mlm-uncentered/05-draws/%s-%s-%s-%s.RData", study2, pname, pitem, sitem))
  
  # find attractor location by solving for x
  # linear               # quadratic
  # y = b0 + b1*x        x = (-b +/ sqrt(b2 - 4ac))/2a
  # 0 = b0 + b1*x        
  # -b0 = b1*x         
  # x = -b0/b1
  fx_equil <- fx_draws %>% 
    as.data.frame() %>%
    # fixed estimates
    mutate(p__attr_loc_0 = -1*b_pdeltavalue_Intercept/b_pdeltavalue_p_value
           , p__attr_loc5 = -1*(b_pdeltavalue_Intercept + 5*b_pdeltavalue_s_value)/(b_pdeltavalue_p_value + 5 * `b_pdeltavalue_p_value:s_value`)
           , p__attr_loc10 = -1*(b_pdeltavalue_Intercept + 10*b_pdeltavalue_s_value)/(b_pdeltavalue_p_value + 10 * `b_pdeltavalue_p_value:s_value`)
           # delta situations at 0, 5, 10
           , s__attr_loc_0 = -1*b_sdeltavalue_Intercept/b_sdeltavalue_s_value
           , s__attr_loc5 = -1*(b_sdeltavalue_Intercept + 5*b_sdeltavalue_p_value)/(b_sdeltavalue_s_value + 5 * `b_pdeltavalue_p_value:s_value`)
           , s__attr_loc10 = -1*(b_sdeltavalue_Intercept + 10*b_sdeltavalue_p_value)/(b_sdeltavalue_s_value + 10 * `b_sdeltavalue_p_value:s_value`))
  # participant-level estimates
  rx_equil <- rx_draws %>%
    # delta personality at 0, 5, 10
    mutate(p__attr_loc_0 = -1*pdeltavalue__Intercept/pdeltavalue__p_value
           , p__attr_loc5 = -1*(pdeltavalue__Intercept + 5*pdeltavalue__s_value)/(pdeltavalue__p_value + 5 * `pdeltavalue__p_value:s_value`)
           , p__attr_loc10 = -1*(pdeltavalue__Intercept + 10*pdeltavalue__s_value)/(pdeltavalue__p_value + 10 * `pdeltavalue__p_value:s_value`)
           # delta situations at 0, 5, 10
           , s__attr_loc_0 = -1*sdeltavalue__Intercept/sdeltavalue__s_value
           , s__attr_loc5 = -1*(sdeltavalue__Intercept + 5*sdeltavalue__p_value)/(sdeltavalue__s_value + 5 * `pdeltavalue__p_value:s_value`)
           , s__attr_loc10 = -1*(sdeltavalue__Intercept + 10*sdeltavalue__p_value)/(sdeltavalue__s_value + 10 * `sdeltavalue__p_value:s_value`))
  
  # find attractor strength
  equil <- fx_equil %>% 
    # fixed estimates
    mutate(
      n = 1:n()
      , SID = "overall"
      , p__attr_str0 = b_pdeltavalue_p_value
      , p__attr_str5 = b_pdeltavalue_p_value + 5*`b_pdeltavalue_p_value:s_value`
      , p__attr_str10 = b_pdeltavalue_p_value + 10*`b_pdeltavalue_p_value:s_value`
      , s__attr_str0 = b_sdeltavalue_s_value
      , s__attr_str5 = b_sdeltavalue_s_value + 5*`b_sdeltavalue_p_value:s_value`
      , s__attr_str10 = b_sdeltavalue_s_value + 10*`b_sdeltavalue_p_value:s_value`
    ) %>% 
    select(n, SID, contains("attr_loc"), contains("attr_str")) %>%
    full_join(
      # participant-level estimates
      rx_equil %>% 
        mutate(
          p__attr_str0 = pdeltavalue__p_value
          , p__attr_str5 = pdeltavalue__p_value + 5*`pdeltavalue__p_value:s_value`
          , p__attr_str10 = pdeltavalue__p_value + 10*`pdeltavalue__p_value:s_value`
          , s__attr_str0 = sdeltavalue__s_value
          , s__attr_str5 = sdeltavalue__s_value + 5*`sdeltavalue__p_value:s_value`
          , s__attr_str10 = sdeltavalue__s_value + 10*`sdeltavalue__p_value:s_value`
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
  save(equil, file = sprintf("05-results/02-q2/03-mlm-uncentered/03-equil/%s-%s-%s-%s.RData", study2, pname, pitem, sitem))
  
  # setup predictions for plotting
  # fixed estimates
  frame <- crossing(p_value = seq(0,10,.5), s_value = seq(0,10,.5))
  fx_pred <- bind_cols(
    bind_rows(frame, frame)
    , fitted(m, newdata = frame, re_formula = NA) %>% array_tree(3) %>% ldply(., .id = "outcome")
  ) %>%
    pivot_wider(
      names_from = "outcome"
      , values_from = c("Estimate", "Est.Error", "Q2.5", "Q97.5")
    )
  # participant-level estimates
  frame <- crossing(SID = unique(m$data$SID), p_value = seq(0,10,.5), s_value = seq(0,10,.5))
  rx_pred <- bind_cols(
    bind_rows(frame, frame)
    , fitted(m, newdata = frame) %>% array_tree(3) %>% ldply(., .id = "outcome")
  ) %>%
    pivot_wider(
      names_from = "outcome"
      , values_from = c("Estimate", "Est.Error", "Q2.5", "Q97.5")
    )
  save(fx_pred, rx_pred
       , file = sprintf("05-results/02-q2/03-mlm-uncentered/04-predictions/%s-%s-%s-%s.RData", study2, pname, pitem, sitem))
  
  # clean up and return nothing
  rm(list = c("frame", "fx_pred", "m1", "m", "fx", "rx", "fx_draws", "rx_draws", "equil", "draws"))
  gc()
  return(T)
}
q2_mlm_model_fun(args[,1], args[,2], args[,3], args[,4])