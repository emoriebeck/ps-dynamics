##### Question 1: MLM ######

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
args <- read.table("scripts/cluster/args/q1-mlm.txt"
                   , header = T, stringsAsFactors = F)[jobid,]
print(args)

q1_mlm_model_fun <- function(study, cat, name, item){
  study2 <- mapvalues(
    study
    , c("01-pairs", "02-ipcs", "03-sherman", "04-horstmann")
    , c("pairs", "ipcs", "sherman", "horstmann")
    , warn_missing = F
  )
  
  # load the data 
  load(sprintf("03-data/%s/04-q1-mlm/%s-%s-%s.RData", study, cat, name, item))
  
  # setup the formula 
  # actually stays consistent, so here more for convenience
  f <- bf(delta_value ~ 1 + value + (1 + value | SID))
  
  # load the sample model 
  load("05-results/01-q1/02-mlm/sample-mlm-linear.RData")
  
  # run the model using update for faster compilation 
  m <- update(
    m1
    , f = f
    , newdata = d
    , iter = 2000
    , warmup = 1000
    , cores = 4
  )
  save(m, file = sprintf("05-results/01-q1/02-mlm/01-models/%s-%s-%s-%s.RData", study2, cat, name, item))
  
  # get fixed and participant-specific term estimates
  fx <- broom.mixed::tidy(m)
  rx <- coef(m, probs = c(0.025, 0.975))[[1]] %>% array_tree(3) %>% 
    tibble(term = names(.), data = .) %>% 
    mutate(data = map(data, ~(.) %>% data.frame %>% 
                        rownames_to_column("SID"))) %>% 
    unnest(data) %>% 
    select(term, SID, estimate = Estimate, conf.low = Q2.5, conf.high = Q97.5)
  save(fx, rx, file = sprintf("05-results/01-q1/02-mlm/02-summary/%s-%s-%s-%s.RData", study2, cat, name, item))
  
  # get and wrangle fixed and random draws
  draws <- as_draws_df(m)
  fx_draws <- draws %>% select(-contains("r_SID["))
  rx_draws <- draws %>% 
    mutate(n = 1:n()) %>%
    select(contains("r_SID["), n) %>% 
    pivot_longer(
      names_to = "term"
      , values_to = "estimate"
      , cols = -n
    ) %>%
    mutate(
      term = str_remove_all(term, "r_SID\\[")
      , term = str_remove_all(term, "[, ]")
      , term = str_remove_all(term, "\\]")
      , SID = str_remove_all(term, "[A-Z a-z]")
      , term = str_remove_all(term, "[0-9]")
    ) %>%
    pivot_wider(names_from = "term", values_from = "estimate") %>%
    full_join(fx_draws %>% mutate(n = 1:n()) %>% select(b_Intercept, b_value, n)) %>%
    mutate(Intercept = Intercept + b_Intercept
           , value = value + b_value) %>%
    select(-b_Intercept, -b_value)
  save(fx_draws, rx_draws
       , file = sprintf("05-results/01-q1/02-mlm/05-draws/%s-%s-%s-%s.RData", study2, cat, name, item))
  
  # find attractor location by solving for x
  # could have used hypothesis() because we don't have to solve for quadratics but leaving it 
  # using posterior draws and then later summarizing for consistency with idiographic Q1
  # linear               # quadratic
  # y = b0 + b1*x        x = (-b +/ sqrt(b2 - 4ac))/2a
  # 0 = b0 + b1*x        
  # -b0 = b1*x         
  # x = -b0/b1
  fx_equil <- fx_draws %>% 
    as.data.frame() %>%
    mutate(attr_loc = -1*b_Intercept/b_value)
  rx_equil <- rx_draws %>%
    mutate(attr_loc = -1*Intercept/value)
  
  # find attractor strength
  equil <- fx_equil %>% 
    mutate(n = 1:n(), SID = "overall") %>% 
    select(SID, n, attr_loc, attr_str = b_value) %>%
    full_join(
      rx_equil %>% 
        select(n, SID, attr_loc, attr_str = value)
    )
  
  # summarize the posterior draws
  equil <- equil %>% 
    select(SID, starts_with("b_"), starts_with("attr")) %>%
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
  
  save(equil, file = sprintf("05-results/01-q1/02-mlm/03-equil/%s-%s-%s-%s.RData", study2, cat, name, item))
  
  # setup predictions for plotting
  frame <- tibble(value = seq(0,10,.2))
  fx_pred <- bind_cols(frame, fitted(m, newdata = frame, re_formula = NA))
  frame <- crossing(SID = unique(rx_equil$SID), value = seq(0,10,.2))
  rx_pred <- bind_cols(frame, fitted(m, newdata = frame))
  save(fx_pred, rx_pred
       , file = sprintf("05-results/01-q1/02-mlm/04-predictions/%s-%s-%s-%s.RData", study2, cat, name, item))
  
  # clean up and return nothing
  rm(list = c("frame", "fx_pred", "m1", "m", "fx", "rx", "fx_draws", "rx_draws", "equil", "draws"))
  gc()
  return(T)
}

q1_mlm_model_fun(args[,1], args[,2], args[,3], args[,4])