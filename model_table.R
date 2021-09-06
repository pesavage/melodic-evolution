## This script aggregates the model output for each data subset
## into a single table

library(brms)
library(dplyr)


build_table = function(suffix){
  model_files = list.files('results/', 
                           suffix, 
                           full.names = TRUE)
  
  model_names = model_files %>% 
    basename() %>%
    tools::file_path_sans_ext()
  
  models = lapply(model_files, function(f) brm(file = f))
  names(models) = model_names
  
  coefficients = lapply(models, fixef)
  coefficients_formatted = lapply(coefficients, function(cf){
    cf = signif(cf, 2)
    cf = data.frame(cf)
    formatted = 
      paste0(cf$Estimate, " (", cf$Q2.5, ", ", cf$Q97.5, ")")
    names(formatted) = rownames(cf)
    formatted
  })
  
  
  results_table = 
    matrix(NA, 
           nrow = 8, 
           ncol = 8, 
           dimnames = 
             list(model_names, 
                  c("N", 
                    "Intercept",
                    "functional_changew",
                    "semitonal_distance",
                    "functional_changew:semitonal_distance",
                    "frequency1:frequency2",
                    "frequency1:frequency2:functional_total",
                    "ELPD"
                  )))
  
  
  ## Fill table 
  for(i in seq_along(coefficients_formatted)){
    m_name = names(coefficients_formatted)[i]
    v_names = names(coefficients_formatted[[i]])
    
    results_table[m_name, v_names] = coefficients_formatted[[i]]
    
    N = nrow(models[[m_name]]$data)
    results_table[m_name, "N"] = N
  }
  results_table = data.frame(results_table)
  
  # Get model comparisons
  loos = lapply(models, function(l){
    values = l$criteria$loo$estimates[1,]
    
    if(is.null(values)){
      out = "-"
    } else {
      values = round(values, 2)
      out = paste0(values[1], " (", values[2], ")")
    }
  })
  results_table$ELPD = unlist(loos)
  
  results_table
}

english_results = build_table("*_english.rds")
japanese_results = build_table("*_japanese.rds")

english_results$Sample = "English"
japanese_results$Sample = "Japanese"

english_results$Model = factor(c("Function", 
                                 "Function + Distance",
                                 "Function * Distance",
                                 "Note frequency",
                                 "Baseline: Note Function & Frequency",
                                 "Null",
                                 "Functional model - strong",
                                 "Functional model - weak"),
                               levels = 
                                 c("Null",
                                   "Note frequency",
                                   "Baseline: Note Function & Frequency",
                                   "Functional model - strong",
                                   "Functional model - weak",
                                   "Function", 
                                   "Function + Distance",
                                   "Function * Distance"
                                    ))

english_results = english_results %>% 
  arrange(Model)

japanese_results$Model = factor(c("Function", 
                                 "Function + Distance",
                                 "Function * Distance",
                                 "Note frequency",
                                 "Baseline: Note Function & Frequency",
                                 "Null",
                                 "Functional model - strong",
                                 "Functional model - weak"),
                               levels = 
                                 c("Null",
                                   "Note frequency",
                                   "Baseline: Note Function & Frequency",
                                   "Functional model - strong",
                                   "Functional model - weak",
                                   "Function", 
                                   "Function + Distance",
                                   "Function * Distance"
                                 ))
japanese_results = japanese_results %>% 
  arrange(Model)

result_table = rbind(english_results, japanese_results)
result_table = result_table %>% 
  relocate(Model, Sample, .before = N)

write.csv(result_table, 'results/modeloutput_table.csv')