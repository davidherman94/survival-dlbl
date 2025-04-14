rm(list=ls()) 

#Libraries
library(survival)
library(dplyr)
library(car)
library(writexl)
library(gtsummary)
library(finalfit)
library(data.table)

## Working directory
setwd("C:/Users/david/Desktop/EPH/Lymphoma/Analysis")

## Load data
datos_sup <- read.csv("datos_sup_rec.csv", sep =  ";", header = TRUE, stringsAsFactors = FALSE)

## Select variables without ROC (only medians)
datos_sup <- datos_sup %>%
  filter(!is.na(Evento)) %>%
  select(-contains("ROC"))

# Ensuring integer values for time and age
datos_sup$tiempo_evento_meses <- as.numeric(datos_sup$tiempo_evento_meses)
datos_sup$edad <- as.numeric(datos_sup$edad)

## Transform to factor variables of interest variables since column event
datos_sup <- mutate_at(datos_sup, vars(4:ncol(datos_sup)), as.factor)

## Cox model Univariate
# function to conduct univariable Cox regression for each variable of interest
function_cox_univariable <- function(data, tiempo_evento_meses, Evento, variables) {
  models <- list()
  results <- list()
  
  for (variable in variables) {
    formula <- as.formula(paste("Surv(", tiempo_evento_meses, ",", Evento, ") ~", variable))
    model <- coxph(formula, data = data, method = "efron", id =  pacientes)
    models[[variable]] <- model

    # Create a summary dataframe for each variable
    summary_df <- fit2df(model, condense = FALSE)
    p <- summary(model)$coefficients[,"Pr(>|z|)"] # Extract p-value
    HR <- exp(coef(model)) # Extract exp(coef)
    conf_int <- exp(confint(model)) # Extract confidence intervals
    CI <- paste(sprintf("%.2f", conf_int[1]), sprintf("%.2f", conf_int[2]), sep = " - ") # Formatting CI and storing
    
    # Store results
    summary_df$HR <- format(HR, digits = 4, scientific = FALSE) # Format HR
    summary_df$p <- sprintf("%.4f", as.numeric(p))
    summary_df$CI <- CI
    results[[variable]] <- summary_df
  }
  
  return(results)
}

# Assuming datos_sup is your dataframe
variables <- names(datos_sup)[-c(1, 2, 4)]  # Exclude "patients", "time in months", and "event" columns
cox_models <- function_cox_univariable(datos_sup, "tiempo_evento_meses", "Evento", variables)

# Combine results into a single data frame
combined_results <- rbindlist(cox_models, idcol = "Variable")
result_table <- knitr::kable(combined_results, digits = c(0, 2, 3, 4, 5))
result_df <- as.data.frame(combined_results) %>%
  select(-explanatory ,-L95, -U95)

# Save as excel output
writexl::write_xlsx(result_df, "resultados_univariable.xlsx")


## Cox model multivariable for significative variables in the univariable analysis
CoxModel_NT5EMEDIANA_mult <- coxph(Surv(tiempo_evento_meses, Evento) ~ NT5EMEDIANA + log(edad) + Fem + EBV + Subtipo, id=pacientes, 
                               method="efron", data=datos_sup)

summary(CoxModel_NT5EMEDIANA_mult)
cox.zph(CoxModel_NT5EMEDIANA_mult)


CoxModel_CD34MEDIANA_mult <- coxph(Surv(tiempo_evento_meses, Evento) ~ CD34MEDIANA + log(edad) + Fem + EBV + Subtipo, id=pacientes, 
                               method="efron", data=datos_sup)

summary(CoxModel_CD34MEDIANA_mult)
cox.zph(CoxModel_CD34MEDIANA_mult)

CoxModel_MYCMEDIANA_mult <- coxph(Surv(tiempo_evento_meses, Evento) ~ MYCMEDIANA + log(edad) + Fem + EBV + Subtipo, id=pacientes, 
                                   method="efron", data=datos_sup)

summary(CoxModel_MYCMEDIANA_mult)
cox.zph(CoxModel_MYCMEDIANA_mult)
