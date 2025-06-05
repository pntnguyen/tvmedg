library(tidyverse)

#----- Theme
#===============================================================================
mytheme <- function(...) {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14,color = "grey10",  face = "bold", hjust = 0.5),
      axis.line = element_line(linetype = "solid"),
      axis.text = element_text(color = "gray10", size = 10),
      axis.title = element_text(color = "gray10", size = 12),
      # plot.background = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction = "horizontal",
      legend.position = "top",
      legend.background = element_rect(fill = NA, color = NA),
      legend.text = element_text(size = 12),
      legend.key.width = unit(1, "line"),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = NA, color = NA)
    )
}


#----- Create dataset
#===============================================================================
# df : original dataset
# datMC: output from g-formula

# Proportion for outcome in the original data
tv_cumsum <- dat |>
  group_by(mm) |>
  summarise(
    y_prop = mean(Yp),
    y_count = sum(Yp)
    ) |> 
  ungroup() |>
  mutate(y_prob_cum = cumsum(y_count)/nrow(df))


# Proportion for outcome in the output data
dat1M <- datMC |> 
  mutate(group = if_else(Ay == 1 & Am == 1, "Q(1,1)",
                         if_else(Ay == 1 & Am ==0, "Q(1,0)", "Q(0,0)")),
         group = factor(group, labels = c("Q(0,0)", "Q(1,0)", "Q(1,1)")),
         groupM = factor(Am, label = c("No", "Yes"))) |> 
  group_by(group, mm) |>
  summarise(Y = mean(Yp)) |>
  ungroup() |>
  group_by(group) |>
  mutate(Ysum = cumsum(Y))|>
  ungroup() 


# Plot for cumulative Y
(f_cumY <- dat1M |> 
    left_join(tv_cumsum, by = "mm") |>
    ggplot() +
    geom_line(aes(x = mm, y = Ysum, color = group), linewidth = 1) +
    geom_line(aes(x = mm, y = y_prob_cum), color = "gray20", 
              linewidth = 1, linetype = 2) +
    scale_color_brewer(palette = "Set1", direction = -1) +
    scale_y_continuous(limits = c(0, 0.6)) +
    scale_x_continuous(breaks = seq(0, 60, by = 12)) +
    mytheme() +
    labs(x = "Month",
       y = "%",
       title = "Cumulative Y (%)",
       caption = "The dashed line represents the observed %",
       color = NULL))

# Plot for time-varying Y
(f_tvY <- dat1M |> 
    left_join(tv_cumsum, by = "mm") |>
    ggplot() +
    geom_line(aes(x = mm, y = Y, color = group), linewidth = 1) +
    geom_line(aes(x = mm, y = y_prop), color = "gray20", 
              linewidth = 1, linetype = 2) +
    scale_color_brewer(palette = "Set1", direction = -1) +
    scale_y_continuous(limits = c(0, 0.1)) +
    scale_x_continuous(breaks = seq(0, 60, by = 12)) +
    mytheme() +
    labs(x = "Month",
         y = "%",
         title = "Y by month (%)",
         caption = "The dashed line represents the observed %",
         color = NULL))










