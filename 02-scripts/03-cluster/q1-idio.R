##### Question 1: Idiographic ######

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
args <- read.table("scripts/cluster/args/q1-idio.txt"
                   , header = T, stringsAsFactors = F)[jobid,]
print(args)

q1_idio_model_fun <- function(study, sid, cat, name, item){
  study2 <- mapvalues(
    study
    , c("01-pairs", "02-ipcs", "03-sherman", "04-horstmann")
    , c("pairs", "ipcs", "sherman", "horstmann")
    , warn_missing = F
  )
  load(sprintf("03-data/%s/03-q1-idio/%s-%s-%s-%s.RData", study, sid, cat, name, item))
  
  # fit the loess model
  ml <- loess(delta_value ~ value, data = d)
  
  # get loess predictions
  pl <- predict(ml, newdata = tibble(value = seq(0, 10, .2)))
  
  # count sign changes = number of equilibria
  num_equil <- sum(sign(pl) != sign(lag(pl)), na.rm = T)
  
  # setup rhs of formula 
  rhs <- "value"
  if(num_equil > 1) rhs <- c(rhs, "value2")
  rhs <- paste0(rhs, collapse = " + ")
  
  # add in the lhs
  f <- paste(c("delta_value", rhs), collapse = " ~ ")
  
  # load sample model 
  if(num_equil <= 1) load("05-results/01-q1/01-idio/sample-idio-linear.RData") else load("05-results/01-q1/01-idio/sample-idio-quad.RData")
  
  # updated from sample run to prevent compilation and speed runtime
  m <- update(m1
              , f = f
              , newdata = d
              , iter = 2000
              , warmup = 1000
              # , cores = 4
  ) 
  save(m, file = sprintf("05-results/01-q1/01-idio/01-models/%s-%s-%s-%s-%s.RData"
                         , study2, sid, cat, name, item))
  
  fx <- broom.mixed::tidy(m)
  save(fx, file = sprintf("05-results/01-q1/01-idio/02-summary/%s-%s-%s-%s-%s.RData"
                          , study2, sid, cat, name, item))
  
  draws <- as_draws_df(m) %>% as_tibble()
  save(draws, file = sprintf("05-results/01-q1/01-idio/05-draws/%s-%s-%s-%s-%s.RData"
                             , study2, sid, cat, name, item))
  
  # find attractor location by solving for x
  # linear               # quadratic
  # y = b0 + b1*x        x = (-b +/ sqrt(b2 - 4ac))/2a
  # 0 = b0 + b1*x        
  # -b0 = b1*x         
  # x = -b0/b1
  equil <- if(num_equil == 1){
    bind_cols(draws, tibble(attr_loc1 = -1*draws$b_Intercept/draws$b_value))
  } else {
    bind_cols(draws, pmap_df(list(draws$b_value2, draws$b_value, draws$b_Intercept), quad_form))
  }
  
  # find attractor strength
  equil <- if(num_equil == 1){
    bind_cols(equil, tibble(attr_str1 = draws$b_value))
  } else {
    # b1*x + b2*x2
    # dx = b1 + 2*b2*x, where x is the location of an equilibrium
    bind_cols(
      equil
      , tibble(
        attr_str1 = draws$b_value + 2*draws$b_value2*draws$attr_loc1
        , attr_str2 = draws$b_value + 2*draws$b_value2*draws$attr_loc2
      )
    )
  }
  
  equil <- equil %>% 
    select(starts_with("b_"), starts_with("attr")) %>%
    posterior_summary() %>%
    as.data.frame() %>% 
    rownames_to_column("term")
  save(equil, file = sprintf("05-results/01-q1/01-idio/03-equil/%s-%s-%s-%s-%s.RData"
                             , study2, sid, cat, name, item))
  
  frame <- tibble(value = seq(0,10,.2))
  if(num_equil > 1) frame <- frame %>% mutate(value = value^2)
  
  pred <- bind_cols(frame, predict(m, newdata = frame))
  save(pred, file = sprintf("05-results/01-q1/01-idio/04-predictions/%s-%s-%s-%s-%s.RData"
                            , study2, sid, cat, name, item))
  
  rm(list = c("draws", "equil", "m", "m1", "d", "frame", "pred", "fx", "rhs", "f"))
}

quad_form <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    x_1    <- (-b+sqrt(delta(a,b,c)))/(2*a)
    x_2    <- (-b-sqrt(delta(a,b,c)))/(2*a)
    result <- tibble(attr_loc1 = x_1, attr_loc2 = x_2)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    result <- tibble(attr_loc1 = -b/(2*a), attr_loc2 = NA_real_)
  }
  else { # third case D<0
    result <- tibble(attr_loc1 = NA_real_, attr_loc2 = NA_real_)
  } 
}

# Constructing delta
delta<-function(a,b,c){
  b^2-4*a*c
}
q1_idio_model_fun(args[,1], args[,2], args[,3], args[,4], args[,5])