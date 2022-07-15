##### Question 2: Idiographic ######

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
args <- read.table("scripts/cluster/args/q2-idio.txt"
                   , header = T, stringsAsFactors = F)[jobid,]
print(args)



q2_idio_model_fun <- function(study, sid, pname, pitem, sname){
  study2 <- mapvalues(
    study
    , c("01-pairs", "02-ipcs", "03-sherman", "04-horstmann")
    , c("pairs", "ipcs", "sherman", "horstmann")
    , warn_missing = F
  )
  load(sprintf("03-data/%s/05-q2-idio/%s-%s-%s-%s.RData", study, sid, pname, pitem, sname))
  
  
  # load sample model 
  load("05-results/02-q2/01-idio/sample-idio-linear.RData")
  
  # formula (not added to update because update doesn't support new multivariate formulas)
  # f1 <- bf(p_delta_value ~ 1 + p_value + s_value + p_value:s_value)
  # f2 <- bf(s_delta_value ~ 1 + s_value + p_value + s_value:p_value)
  # f <- mvbrmsformula(f1, f1) + set_rescor(TRUE)
  
  # updated from sample run to prevent compilation and speed runtime
  m <- update(m1
              , newdata = d
              , iter = 2000
              , warmup = 1000
              # , cores = 4
  ) 
  save(m, file = sprintf("05-results/02-q2/01-idio/01-models/%s-%s-%s-%s-%s.RData", study2, sid, pname, pitem, sname))
  
  fx <- broom.mixed::tidy(m)
  save(fx, file = sprintf("05-results/02-q2/01-idio/02-summary/%s-%s-%s-%s-%s.RData", study2, sid, pname, pitem, sname))
  
  draws <- as_draws_df(m)
  save(draws, file = sprintf("05-results/02-q2/01-idio/05-draws/%s-%s-%s-%s-%s.RData", study2, sid, pname, pitem, sname))
  
  # find attractor location by solving for x
  # linear               # quadratic
  # y = b0 + b1*x        x = (-b +/ sqrt(b2 - 4ac))/2a
  # 0 = b0 + b1*x        
  # -b0 = b1*x         
  # x = -b0/b1
  equil <- draws %>%
    mutate(
      p__attr_loc0 = -1 * b_pdeltavalue_Intercept / b_pdeltavalue_p_value
      , p_attr_loc5 = -1*(b_pdeltavalue_Intercept + 5*b_pdeltavalue_s_value)/(b_pdeltavalue_p_value + 5*`b_pdeltavalue_p_value:s_value`)
      , p_attr_loc10 = -1*(b_pdeltavalue_Intercept + 10*b_pdeltavalue_s_value)/(b_pdeltavalue_p_value + 10*`b_pdeltavalue_p_value:s_value`)
      , s__attr_loc0 = -1 * b_sdeltavalue_Intercept / b_sdeltavalue_s_value
      , s_attr_loc5 = -1*(b_sdeltavalue_Intercept + 5*b_sdeltavalue_p_value)/(b_sdeltavalue_s_value + 5*`b_sdeltavalue_p_value:s_value`)
      , s_attr_loc10 = -1*(b_sdeltavalue_Intercept + 10*b_sdeltavalue_p_value)/(b_sdeltavalue_s_value + 10*`b_sdeltavalue_p_value:s_value`)
    )
  
  # find attractor strength
  equil <- equil %>%
    mutate(
      p__attr_str0    = b_pdeltavalue_p_value
      , p__attr_str5  = b_pdeltavalue_p_value + 5 *`b_pdeltavalue_p_value:s_value`
      , p__attr_str10 = b_pdeltavalue_p_value + 10*`b_pdeltavalue_p_value:s_value`
      , s__attr_str0  = b_sdeltavalue_s_value
      , s__attr_str5  = b_sdeltavalue_s_value + 5 *`b_sdeltavalue_p_value:s_value`
      , s__attr_str10 = b_sdeltavalue_s_value + 10*`b_sdeltavalue_p_value:s_value`
    ) #%>%
  # select(n = .draw, contains("attr_"))
  
  equil <- equil %>% 
    select(contains("attr")) %>%
    posterior_summary() %>%
    as.data.frame() %>% 
    rownames_to_column("term")
  save(equil, file = sprintf("05-results/02-q2/01-idio/03-equil/%s-%s-%s-%s-%s.RData", study2, sid, pname, pitem, sname))
  
  frame <- crossing(p_value = seq(0,10,.5), s_value = seq(0,10,.5))
  
  pred <- bind_cols(
    bind_rows(frame, frame)
    , fitted(m, newdata = frame) %>% array_tree(3) %>% ldply(., .id = "outcome")
  ) %>%
    pivot_wider(
      names_from = "outcome"
      , values_from = c("Estimate", "Est.Error", "Q2.5", "Q97.5")
    )
  save(pred, file = sprintf("05-results/02-q2/01-idio/04-predictions/%s-%s-%s-%s-%s.RData", study2, sid, pname, pitem, sname))
  
  # pred %>% 
  #   mutate(col = sqrt(Estimate_pdeltavalue^2 + Estimate_sdeltavalue^2)) %>%
  #   ggplot(aes(x = p_value, y = s_value, u = Estimate_pdeltavalue, v = Estimate_sdeltavalue)) +
  #     geom_quiver(aes(color = col))   +
  #     # geom_hline(aes(yintercept = m_h)) +
  #     # geom_vline(aes(xintercept = m_p)) +
  #     scale_color_gradient(low = "gray80", high = "black") +
  #     labs(x = "Personality State Level", y = "Situation Characteristic Level", title = "") +
  #     theme_classic() +
  #     theme(legend.position = "none",
  #           strip.background = element_rect(fill = "black"),
  #           strip.text = element_text(face = "bold", color = "white", size = rel(1.2)),
  #           axis.text = element_text(face = "bold", color = "black"),
  #           axis.title = element_text(face = "bold", size = rel(1.1)),
  #           plot.title = element_text(face = "bold", hjust = .5))
  
  rm(list = c("draws", "equil", "m", "m1", "d", "frame", "pred", "fx", "rhs", "f"))
}
q2_idio_model_fun(args[,1], args[,2], args[,3], args[,4], args[,5])