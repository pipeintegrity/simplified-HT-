
library(readxl)
library(tidyverse)
library(tidymodels)
library(here)
library(minpack.lm)
library(patchwork)
library(latex2exp)

theme_set(theme_light(14))
pal <- RColorBrewer::brewer.pal(3,"Set1")

## Read data in------------------------------------
tim_cvn <-
  read_excel(here("timp_met_charpy_clean_rev02.xlsx"),
             sheet = "Charpy") %>%
  janitor::clean_names() %>%
  rename(
    location = specimen_location,
    temp_f = test_temperature_f,
    charpy_thickness_mm = specimen_thickness_mm,
    cvn_ss_ft_lbs = measured_energy_absorbed_ft_lbs,
    sa_ss = shear_area_percent
  ) %>%
  mutate(across(temp_f:sa_ss, as.numeric),
         sa_ss = as.numeric(str_remove(
           string = sa_ss, pattern = c("<", ">")
         ))) %>%
  filter(
    str_detect(location, "Base Metal"),
    #removing weld tests
    str_detect(single_or_multiple_temperature_test, "Multiple"),
    #remove single temp tests
    str_detect(full_transition_curve, "Yes"),
    sa_ss <= 100
  ) %>% #removing shear areas greater than 100
  select(id,
         bar_orientation,
         charpy_thickness_mm,
         temp_f,
         cvn_ss_ft_lbs,
         sa_ss) %>% 
  rename(cvn = cvn_ss_ft_lbs)

#list of id that max is less than 90%
idlist <- tim_cvn %>% 
  group_by(id) %>% 
  summarise(mx_sa=max(sa_ss)) %>% 
  filter(mx_sa<90)

temp_id_list <- tim_cvn %>% 
  group_by(id, temp_f) %>% 
  summarise(unique(temp_f)) %>% 
  mutate(count = n()) %>% 
  filter(count<4)

`%notin%` <- Negate(`%in%`) # not in special function

tim_cvn <- tim_cvn %>% 
  filter(id %notin% idlist$id,
         id %notin% temp_id_list$id)

ids <- unique(tim_cvn$id) #create vector of ids

## Select which id to use ----------------------------------------
idx <- 3 #start index at 1

# initial estimates for A, B, C, D for SA solution

cvn_data <- tim_cvn %>% 
  filter(id==ids[idx])

plot(sa_ss ~ temp_f, cvn_data)
plot(cvn ~ temp_f, cvn_data)

Ai_sa <- 51 #initial estimates A + B <= 100
Bi_sa <- 49 #initial estimates
Ci <-  IQR(cvn_data$temp_f) / 2 #inner quartile range/2 = Ci
Di <- (cvn_data$temp_f[which.min(abs(50 - cvn_data$sa_ss))]) 
# Di = what temp is closest to 50% SA

# linearize the data -----------------------------------------
cvn_data <- cvn_data %>%
  mutate(k1 = atanh((sa_ss - Ai_sa) / Bi_sa))


# hyperbolic tangent function
func <- function(t, A, B, C, D) {
  A + B * tanh((t - D) / C)
}


# SA solution ----------------------------------------------
fit_sa <-
  nlsLM(
    sa_ss ~ func(temp_f, A, B, C, D),
    data = cvn_data ,
    start = list(
      A = Ai_sa,
      B = Bi_sa,
      C = Ci,
      D = Di
    ),
    trace = F
  );tidy(fit_sa)

# simplified SA HT method ------------
fit_sa_simp <-
  nlsLM(
    sa_ss ~ func(temp_f,Ai_sa, Bi_sa, C, D),
    data = cvn_data ,
    start = list(
      C = Ci,
      D = Di
    ),
    trace = F
  );tidy(fit_sa_simp)

fit_sa_lm <- lm(temp_f ~ k1, data = cvn_data) 


SA_coefs_gen <- tidy(fit_sa,conf.int = FALSE) %>% 
  mutate(model="General Solution")

SA_coefs_simp <- tidy(fit_sa_simp,conf.int = FALSE) %>% 
  mutate(model="Simplified HT")

sa_lm_coefs <- tidy(fit_sa_lm, conf.int = FALSE) %>%
  mutate(term = ifelse(term == "k1", "C", "D"))

sa_coefs_all <-  sa_lm_coefs %>%
  bind_rows(tibble(term = c("A", "B"),
                   estimate = c(51, 49))) %>%
  mutate(model = "Simplified OLS") %>%
  bind_rows(SA_coefs_simp) %>% 
  bind_rows(tibble(term = c("A", "B"),
                   estimate = c(51, 49), model = "Simplified HT")) %>%
  mutate(term = case_when(term == "(Intercept)" ~ "D",
                          term == "k1" ~ "C",
                          TRUE ~ term)) %>%
  bind_rows(SA_coefs_gen) %>%
  mutate(across(.cols = everything(), ~ ifelse(is.na(.), 0, .))) %>%
  arrange(term, model)

satt_summary <- sa_coefs_all %>% 
  select(term, model, estimate) %>% 
  pivot_wider(id_cols = model, names_from = term, values_from = estimate) %>% 
  group_by(model) %>% 
  summarise(satt = C*atanh((85 - A)/B)+D)

func_ht <- function(temp_f, A, B, C, D) {
  A + B * tanh((temp_f - D) / C)
}

cols <-
  c(
    "General" = "red",
    "Simplified OLS" = "blue",
    "Simplified HT" = "green"
  )

cvn_data %>%
  ggplot(aes(temp_f, sa_ss)) +
  geom_point() +
  stat_function(
    fun = func_ht,
    args = list(
      sa_coefs_all$estimate[1],
      sa_coefs_all$estimate[4],
      sa_coefs_all$estimate[7],
      sa_coefs_all$estimate[10]
    ),
    aes(col = "General"),
    lwd = 1, 
    n = 301
  ) +
  stat_function(
    fun = func_ht,
    args = list(
      sa_coefs_all$estimate[2],
      sa_coefs_all$estimate[5],
      sa_coefs_all$estimate[8],
      sa_coefs_all$estimate[11]
    ),
    aes(col = "Simplified HT"),
    lwd =1, 
    n = 301
  ) +
  stat_function(
    fun = func_ht,
    args = list(
      sa_coefs_all$estimate[3],
      sa_coefs_all$estimate[6],
      sa_coefs_all$estimate[9],
      sa_coefs_all$estimate[12]
    ),
    aes(col = "Simplified OLS"),
    lwd =1, 
    n = 301
  ) +
  geom_hline(yintercept = 85, 
             col='grey50', 
             lty =2, 
             lwd = 0.65)+
  scale_y_continuous(breaks = scales::pretty_breaks())+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  labs(x = "Temperature (\u00B0F)",
       y = "Shear Area (%)") +
  annotate(
    "text",
    x = 20,
    y = 85,
    label = "SATT",
    vjust = -0.5
  ) +
  scale_color_manual(
    values = cols,
    breaks = c("General", "Simplified HT", "Simplified OLS"),
    name = "Solution"
  )


# AE solution -------------------------------------------------

# initial values for AE solution

ae_init <- cvn_data %>% 
  summarise(Ai_ae = max(cvn)/2,
            Bi_ae = (max(cvn) - min(cvn))/2)


fit_cvn <-
  nlsLM(
    cvn ~ func(temp_f, A, B, C, D),
    data = cvn_data ,
    start = list(
      A = ae_init$Ai_ae,
      B = ae_init$Bi_ae,
      C = Ci,
      D = Di
    ),
    trace = F
  )

tdy_cvn <- tidy(fit_cvn, conf.int = TRUE);tdy_cvn

sa_simp_sols <- tibble(C = SA_coefs_simp$estimate[1], D =SA_coefs_simp$estimate[2])

# AE solutions ------------------------------------------------------------

#Simplified HT Solution

fit_cvn_simp <-
  nlsLM(
    cvn ~ func(temp_f, A, B, sa_simp_sols$C, sa_simp_sols$D),
    data = cvn_data ,
    start = list(
      A = ae_init$Ai_ae,
      B = ae_init$Bi_ae
    ),
    trace = F
  );tidy(fit_cvn_simp)

AE_gen_coefs <- tidy(fit_cvn, conf.int = FALSE) %>% 
  mutate(model = "General Solution")

AE_coefs_simp <- tidy(fit_cvn_simp, conf.int = FALSE) %>% 
  bind_rows(SA_coefs_simp) %>% 
  mutate(model = "Simplified HT")

# linearize AE Data with SA Solution
cvn_data <- cvn_data %>% 
  mutate(k2 = tanh((temp_f - sa_lm_coefs$estimate[1]) / sa_lm_coefs$estimate[2]))

ae_lm_mod <- lm(cvn ~ k2, data = cvn_data) # linear model

ae_lm_coefs <- tidy(ae_lm_mod, conf.int = FALSE) %>% 
  mutate(term = ifelse(term=="k2", "B", "A")) %>% 
  bind_rows(sa_lm_coefs) %>% 
  arrange(term) %>% 
  mutate(model = "Simplified OLS");ae_lm_coefs

ae_coefs_all <- AE_gen_coefs %>% 
  bind_rows(ae_lm_coefs,AE_coefs_simp) %>% 
  arrange(term, model);ae_coefs_all

ae_us_all <- ae_coefs_all %>% 
  pivot_wider(id_cols = model,names_from = term, values_from = estimate) %>% 
  group_by(model) %>% 
  mutate(us = A + B);ae_us_all


# plot AE solutions -------------------------------------------------------

cvn_data %>%
  ggplot(aes(temp_f, cvn)) +
  geom_point() +
  stat_function(
    fun = func_ht,
    args = list(
      ae_coefs_all$estimate[1],
      ae_coefs_all$estimate[4],
      ae_coefs_all$estimate[7],
      ae_coefs_all$estimate[10]
    ),
    aes(col = "General"),
    lwd = 1, 
    n = 301
  ) +
  stat_function(
    fun = func_ht,
    args = list(
      ae_coefs_all$estimate[2],
      ae_coefs_all$estimate[5],
      ae_coefs_all$estimate[8],
      ae_coefs_all$estimate[11]
    ),
    aes(col = "Simplified HT"),
    lwd =1, 
    n = 301
  ) +
  stat_function(
    fun = func_ht,
    args = list(
      ae_coefs_all$estimate[3],
      ae_coefs_all$estimate[6],
      ae_coefs_all$estimate[9],
      ae_coefs_all$estimate[12]
    ),
    aes(col = "Simplified OLS"),
    lwd =1, 
    n = 301
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(x = "Temperature (\u00B0F)",
       y = "Absorbed Energy") +
  scale_color_manual(
    values = cols,
    breaks = c("General", "Simplified HT", "Simplified OLS"),
    name = "Solution"
  )

# NLS bootstrap -----------------------------------------------------------

nsim <- 1e3
# cvn_data <- cvn_data %>% rename(cvn = ae_ss_ft_lbs)
# AE Bootstrap ------------------------------------------------------------
ae_init<- cvn_data %>% 
  summarise(Ai_ae = max(cvn)/2,
            Bi_ae = (max(cvn) - min(cvn))/2)
# initial estimates for A, B, C, D for SA solution
Ai_sa <- 51 #initial estimates A + B <= 100
Bi_sa <- 49 #initial estimates
Ci <-  IQR(cvn_data$temp_f) / 2 #inner quartile range/2 = Ci
Di <- (cvn_data$temp_f[which.min(abs(50 - cvn_data$sa_ss))]) 



fit_ae_on_bootstrap2 <- function(split) {
  nlsLM(
    cvn ~  A + B * tanh((temp_f - D) / C),
    data = analysis(split),
    start = list(
      A = ae_init$Ai_ae,
      B = ae_init$Bi_ae,
      C = Ci,
      D = Di
    ),
    trace = F,
    control = nls.lm.control(maxiter = 100)
  )
}

fit_sa_on_bootstrap2 <- function(split) {
  nlsLM(
    sa_ss ~  A + B * tanh((temp_f - D) / C),
    data = analysis(split),
    start = list(
      A = Ai_sa,
      B = Bi_sa,
      C = Ci,
      D = Di
    ),
    trace = F
  )
}

#Create Bootstrap data frame
boots <- cvn_data %>%
  bootstraps(times = nsim,
             apparent = FALSE) 

# Create nested models for AE
ae_models <-
  boots %>%
  mutate(
    model = map(splits, fit_ae_on_bootstrap2),
    coef_info = map(model, tidy),
    glanced = map(model, glance)
  )


# Bootstrap coefficients for AE
ae_boot_coefs <-
  ae_models %>%
  unnest(coef_info)

ae_coefs_clean <- ae_boot_coefs %>%
  select(term, estimate, id) %>%
  pivot_wider(id_cols = id,
              names_from = term,
              values_from = estimate) %>%
  filter(A < 50, B < 50)

# Histogram of Coefficients 
ae_coefs_clean %>%
  pivot_longer(A:D,
               names_to = "term",
               values_to = "estimate") %>%
  ggplot(aes(estimate)) +
  geom_histogram(aes(fill = term),
                 col = 'black',
                 show.legend = FALSE) +
  facet_wrap(~ term, scales = "free")

# plotting temperatures
temps <-
  tibble(temp_f = seq(min(cvn_data$temp_f),
                      max(cvn_data$temp_f),
                      length.out = 200)) 

# grab sample of bootstraps
ae_coefs_samp <- ae_coefs_clean %>% 
  slice_sample(n=100)

#Generate HT predictions for each temp and bootstrap combo
ae_cross <- crossing(temp_f = temps$temp_f, ae_coefs_samp) %>%
  arrange(id) %>%
  mutate(f = func(temp_f, A, B, C, D))

# Upper shelf confidence interval (general)
use_CI <- tdy_cvn %>%
  summarise(
    use_lci = conf.low[1] + conf.low[2],
    use_uci = conf.high[1] + conf.high[2],
    use = estimate[1] + estimate[2]
  )
  
# plot bootstraps for AE
ae_boot <- ggplot(cvn_data, aes(temp_f, cvn)) +
  geom_point() +
  geom_line(
    data = ae_cross,
    aes(temp_f, f, group = id),
    col = 'blue',
    lwd = 0.1,
    alpha = 0.1
  ) +
  geom_errorbar(
    aes(
      x = 170,
      ymax = use_CI$use_uci,
      ymin = use_CI$use_lci
    ),
    lwd = 0.75,
    width = 5
  ) +
  geom_point(aes(x = 170, y = use_CI$use),
             col = 'red') +
  labs(
    title = "Absorbed Energy HT Bootstrap Solutions",
    y = "Absorbed Energy (ft-lbs)",
    x = "Temperature (\u00B0F)",
    caption = "100 Random Bootstrapped Solutions"
  )

# Create nested models for SA
sa_models <-
  boots %>%
  mutate(
    model = map(splits, fit_sa_on_bootstrap2),
    coef_info = map(model, tidy),
    glanced = map(model, glance)
  )

# grab sample of bootstraps
sa_boot_coefs <-
  sa_models %>%
  unnest(coef_info)

sa_coefs_clean <- sa_boot_coefs %>%
  select(term, estimate, id) %>%
  pivot_wider(id_cols = id,
              names_from = term,
              values_from = estimate) %>% 
  mutate(AB = A + B) #%>% 
  # filter(AB<=100)

sa_coefs_samp <- sa_coefs_clean %>% 
  slice_sample(n=100)

#Generate HT predictions for each temp and bootstrap combo
sa_cross <- crossing(temp_f = temps$temp_f, sa_coefs_samp) %>%
  arrange(id) %>%
  mutate(f = func(temp_f, A, B, C, D))


saboot <- ggplot(cvn_data, aes(temp_f, sa_ss)) +
  geom_point() +
  geom_line(
    data = sa_cross,
    aes(temp_f, f, group = id),
    col = 'blue',
    lwd = 0.1,
    alpha = 0.1
  ) +
  geom_errorbarh(
    aes(
      y = 85,
      xmax = satt_gen_CI$UCI,
      xmin = satt_gen_CI$LCI
    ),
    lwd = 0.75,
    height = 5
  ) +
  geom_point(aes(x =satt_gen_CI$estimate , y = 85),
             col = 'red') +
  labs(
    title = "Shear Area HT Bootstrap Solutions",
    y = "Shear Area (%)",
    x = "Temperature (\u00B0F)",
    caption = "100 Random Bootstrapped Solutions"
  )

# Histogram SA of Coefficients 
sa_coefs_clean %>%
  pivot_longer(A:D,
               names_to = "term",
               values_to = "estimate") %>%
  ggplot(aes(estimate)) +
  geom_histogram(aes(fill = term),
                 col = 'black',
                 show.legend = FALSE,
                 bins = 45) +
  facet_wrap(~ term, scales = "free")+
  labs(title = "Shear Area Coefficients")

ae_boot + saboot
