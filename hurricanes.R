###################################################################################################
#' Replicability crisis in Science?
#' Summer School 2023 - Padova, Italy
#' 
#' Teamwork 6
#'  Patric Dolmeta, Bocconi University
#'  Anna Gerna, Roma3 University 
#'  Daniele Girolimetto, University of Padova
#'
#' Original code: https://mucollective.github.io/multiverse/articles/hurricane.html
###################################################################################################

####### Library ----
library(tidyverse)
library(tidybayes)
library(multiverse)
library(gridExtra)

####### Functions ----
pred2expectation_2 <- function(model, mu, sigma) {
  ifelse(model == "linear", exp(mu + sigma^2/2) - 1,  mu)
}
set_panel_heights <- function(g, x, y){
  id_panels <- unique(g$layout[g$layout$name=="panel", "t"])
  for(i in 1:length(id_panels)){
    g$heights[id_panels[i]] <- unit(y[i], x[i])
  }
  g
}
extract_ICmean <- function(x, y){
  methods <- c("norm","basic", "perc")
  out <- t(sapply(setNames(methods, methods),
                  function(z){
                    DescTools::MeanDiffCI(x, y,
                                          method = z)
                  }
  ))
  colnames(out) <- c("mean", "ICd", "ICu")
  list(as_tibble(out, rownames = "type"))
}

####### Multiverse analysis ----
data("hurricane")

# read and process data
hurricane_data <- hurricane |>
  # rename some variables
  rename(
    year = Year,
    name = Name,
    dam = NDAM,
    death = alldeaths,
    female = Gender_MF,
    masfem = MasFem,
    category = Category,
    pressure = Minpressure_Updated_2014,
    wind = HighestWindSpeed
  ) |>
  # create new variables
  mutate(
    post = ifelse(year>1979, 1, 0),
    zcat = as.numeric(scale(category)),
    zpressure = -scale(pressure),
    zwind = as.numeric(scale(wind)),
    z3 = as.numeric((zpressure + zcat + zwind) / 3)
  )

# Original model
df <- hurricane_data |>
  filter(name != "Katrina" & name != "Audrey")
fit <- glm(death ~ masfem * dam + masfem * zpressure, data = df, family = "poisson")

# Multiverse approach
M <- multiverse()
inside(M, {
  df <- hurricane_data |>
    filter(branch(death_outliers, 
                  "no_exclusion" ~ TRUE,
                  "most_extreme_deaths" ~ name != "Katrina",
                  "most_extreme_two_deaths" ~ ! (name %in% c("Katrina", "Audrey"))
    )) |>
    filter(branch(damage_outliers,
                  "no_exclusion" ~ TRUE,
                  "most_extreme_one_damage" ~ ! (name %in% c("Sandy")),
                  "most_extreme_two_damage" ~ ! (name %in% c("Sandy", "Andrew")),
                  "most_extreme_three_damage" ~ ! (name %in% c("Sandy", "Andrew", "Donna"))
    ))
  
  df <- df |>
    mutate(femininity = branch(femininity_calculation,
                          "masfem" ~ masfem,
                          "female" ~ female),
      damage = branch(damage_transform,
                      "no_transform" ~ identity(dam),
                      "log_transform" ~ log(dam)))
  
  fit <- glm(branch(model, "linear" ~ log(death + 1), "poisson" ~ death) ~ 
               branch(main_interaction,
                      "no" ~ femininity + damage,
                      "yes" ~ femininity * damage) + 
               branch(other_predictors,
                          "none" ~ NULL,
                          "pressure" %when% (main_interaction == "yes") ~ femininity * zpressure,
                          "wind" %when% (main_interaction == "yes") ~ femininity * zwind,
                          "category" %when% (main_interaction == "yes") ~ femininity * zcat,
                          "all" %when% (main_interaction == "yes") ~ femininity * z3,
                          "all_no_interaction" %when% (main_interaction == "no") ~ z3) + 
               branch(covariates, "1" ~ NULL, "2" ~ year:damage, "3" ~ post:damage), 
             family = branch(model, "linear" ~ "gaussian", "poisson" ~ "poisson"),  
             data = df)
  
  pred <- predict(fit, se.fit = TRUE, type = "response")
  
  pred2expectation <- function(mu, sigma) {
    branch(model, "linear" ~ exp(mu + sigma^2/2) - 1, "poisson" ~ mu)
  }
  
  disagg_fit <- df  |>
    mutate(fitted = pred$fit, # add fitted predictions and standard errors to dataframe
      se.fit = pred$se.fit,
      deg_f = df.residual(fit), # get degrees of freedom
      sigma = sigma(fit), # get residual standard deviation
      se.residual = sqrt(sum(residuals(fit)^2) / deg_f)) # get residual standard errors
  
  # aggregate fitted effect of female storm name
  expectation <- disagg_fit |>
    mutate(expected_deaths = pred2expectation(fitted, sigma)) |> 
    group_by(female) |>
    summarise(mean_deaths = mean(expected_deaths), .groups = "drop_last") |> 
    compare_levels(mean_deaths, by = female)
})

execute_multiverse(M, progress = TRUE)

# Extract multiverse results
multiverse_object <- multiverse::expand(M)

# Universe information dataset (original .universe = 1375)
info <- multiverse_object |>
  select(.universe, death_outliers, damage_outliers, femininity_calculation, 
         damage_transform, model, main_interaction, other_predictors, covariates)



# Plot dataset
dataplot <- multiverse_object |>
  extract_variables(disagg_fit) |>
  select(.universe, disagg_fit, model) |>
  unnest(disagg_fit) |>
  mutate(expected_deaths = pred2expectation_2(model, fitted, sigma)) |> 
  select(.universe, female, expected_deaths, se.residual, deg_f, se.fit) |>
  group_by(.universe) |>
  summarise(
    data = extract_ICmean(expected_deaths[female == 1], 
                                  expected_deaths[female == 0])) |> 
  unnest(data) |>
  left_join(info, by = ".universe")  |>
  arrange(mean) %>%
  mutate(.id = 1:nrow(.),
         check = ifelse(ICd>0 | ICu<0, "outside the confidence interval", "inside"),
         original = .universe == 1375)
  
# Plots
g1 <- ggplotGrob(dataplot |>
                   filter(type == "norm") |>
                   arrange(mean) %>%
                   mutate(.id = 1:nrow(.)) |>
                   ggplot(aes(x = .id)) +
                   geom_line(aes(y = mean)) +
                   theme_minimal() +
                   labs(x = "universe", y = NULL,
                        title = "Mean difference in expected deaths") +
                   theme(axis.title.x = element_blank(),
                         plot.margin = margin(b = 0, t = 4, r = 4, l = 4),
                         plot.background = element_rect(fill = "white", colour = "white"),
                         axis.text.x=element_blank()))

g1_norm <- ggplotGrob(dataplot |>
                        filter(type == "norm") |>
                        arrange(mean) %>%
                        mutate(.id = 1:nrow(.)) |>
                        ggplot(aes(x = .id)) +
                        geom_line(aes(y = mean))+
                        geom_errorbar(aes(ymin = ICd, ymax = ICu), linewidth = 0.05)+
                        #geom_vline(aes(xintercept = .id, col = check), data = function(x) filter(x, .universe == 1375)) +
                        theme_minimal() +
                        labs(x = "universe", y = NULL,
                             title = "Mean difference in expected deaths", 
                             subtitle = "95% CI norm approach") +
                        theme(axis.title.x = element_blank(),
                              legend.position = "none",
                              plot.margin = margin(b = 0, t = 4, r = 4, l = 4),
                              plot.background = element_rect(fill = "white", colour = "white"),
                              axis.text.x=element_blank()))

g1_basic <- ggplotGrob(dataplot |>
                        filter(type == "basic") |>
                         arrange(mean) %>%
                         mutate(.id = 1:nrow(.)) |>
                        ggplot(aes(x = .id)) +
                        geom_line(aes(y = mean))+
                        geom_errorbar(aes(ymin = ICd, ymax = ICu), linewidth = 0.05)+
                        #geom_vline(aes(xintercept = .id, col = check), data = function(x) filter(x, .universe == 1375)) +
                        theme_minimal() +
                        labs(x = "universe", y = NULL,
                             title = "Mean difference in expected deaths", 
                             subtitle = "95% CI boot approach") +
                        theme(axis.title.x = element_blank(),
                              legend.position = "none",
                              plot.margin = margin(b = 0, t = 4, r = 4, l = 4),
                              plot.background = element_rect(fill = "white", colour = "white"),
                              axis.text.x=element_blank()))

g1_perc <- ggplotGrob(dataplot |>
                        filter(type == "perc") |>
                        arrange(mean) %>%
                        mutate(.id = 1:nrow(.)) |>
                        ggplot(aes(x = .id)) +
                        geom_line(aes(y = mean))+
                        geom_errorbar(aes(ymin = ICd, ymax = ICu), linewidth = 0.05)+
                        #geom_vline(aes(xintercept = .id, col = check), data = function(x) filter(x, .universe == 1375)) +
                        theme_minimal() +
                        labs(x = "universe", y = NULL,
                             title = "Mean difference in expected deaths", 
                             subtitle = "95% CI boot (perc) approach") +
                        theme(axis.title.x = element_blank(),
                              legend.position = "none",
                              plot.margin = margin(b = 0, t = 4, r = 4, l = 4),
                              plot.background = element_rect(fill = "white", colour = "white"),
                              axis.text.x=element_blank()))

legend_tmp <- theme(legend.position = "bottom",
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 10),
                    plot.background = element_rect(fill = "white", colour = "white"),
                    plot.margin = margin(b = 4, t = 0, r = 4, l = 4),
                    panel.grid = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y=element_text(color = "white"), 
                    axis.ticks.y=element_blank())

g2 <- ggplotGrob(dataplot |>
                   filter(type == "norm") |>
                   arrange(mean) %>%
                   mutate(.id = 1:nrow(.)) |>
                   ggplot(aes(x = .id)) +
                   geom_tile(aes(y = 1, fill = death_outliers)) +
                   labs(x = "universe") +
                   theme_minimal() + legend_tmp)

g3 <- ggplotGrob(dataplot |>
                   filter(type == "norm") |>
                   arrange(mean) %>%
                   mutate(.id = 1:nrow(.)) |>
                   ggplot(aes(x = .id)) +
                   geom_tile(aes(y = 1, fill = model)) +
                   #facet_grid("death_outliers"~.)  +
                   labs(x = "universe") +
                   theme_minimal() + legend_tmp)

g4_norm <- ggplotGrob(dataplot |>
                   filter(type == "norm") |>
                     arrange(mean) %>%
                     mutate(.id = 1:nrow(.)) |>
                   ggplot(aes(x = .id)) +
                   geom_tile(aes(y = 1, fill = check)) +
                   #facet_grid("death_outliers"~.)  +
                   labs(x = "universe", fill = "Zero is ") +
                   theme_minimal() + legend_tmp)

g4_basic <- ggplotGrob(dataplot |>
                   filter(type == "basic") |>
                     arrange(mean) %>%
                     mutate(.id = 1:nrow(.)) |>
                   ggplot(aes(x = .id)) +
                   geom_tile(aes(y = 1, fill = check)) +
                   #facet_grid("death_outliers"~.)  +
                   labs(x = "universe", fill = "Zero is ") +
                   theme_minimal() + legend_tmp)

g4_perc <- ggplotGrob(dataplot |>
                         filter(type == "perc") |>
                        arrange(mean) %>%
                        mutate(.id = 1:nrow(.)) |>
                         ggplot(aes(x = .id)) +
                         geom_tile(aes(y = 1, fill = check)) +
                         #facet_grid("death_outliers"~.)  +
                         labs(x = "universe", fill = "Zero is ") +
                         theme_minimal() + legend_tmp)


g <- set_panel_heights(rbind(g1, g2), y = c(1,1), x = c("null", "line"))
plot(g) 

gmod <- set_panel_heights(rbind(g1, g3), y = c(1,1), x = c("null", "line"))
plot(gmod) 

gcheck_norm <- set_panel_heights(rbind(g1_norm, g4_norm), 
                                 y = c(1,1), x = c("null", "line"))
plot(gcheck_norm) 

gcheck_basic <- set_panel_heights(rbind(g1_basic, g4_basic), 
                                 y = c(1,1), x = c("null", "line"))
plot(gcheck_basic) 

gcheck_perc <- set_panel_heights(rbind(g1_perc, g4_perc), 
                                  y = c(1,1), x = c("null", "line"))
plot(gcheck_perc) 
save(g, gmod, gcheck_norm, gcheck_basic, gcheck_perc, file = "slideplot.RData")


