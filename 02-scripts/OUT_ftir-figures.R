## Setup

grind_sample_order <- c( # make sure the grind samples are ordered correctly in plots
  "Archaeological calculus",
  "Artificial calculus day 16",
  "Artificial calculus day 20",
  "Artificial calculus day 24",
  "Archaeological bone",
  "Bone-Dentine",
  "Bone-Dentine_2",
  "Enamel",
  "Enamel_2",
  "Enamel_3")

## Spectra

# day 7
ftir_day7 <- ftir_data %>%
  filter(
    analysis_id == "F7.1A6_b"
  ) %>%
  ggplot(aes(x = wavenumber, y = abs, col = analysis_id)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse()

# day 12
ftir_day12 <- ftir_data %>%
  filter(
    analysis_id == "F12.1A5+F12.B1_B" |
      analysis_id == "F12.1D1+F12.1D2"
  ) %>%
  ggplot(aes(x = wavenumber, y = abs, col = analysis_id)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse()

# day 16
ftir_day16 <- ftir_data %>%
  filter(
    #day == 16,
    #grind == FALSE,
    analysis_id == "F16.1A2" |
      analysis_id == "F16.1C6" |
      analysis_id == "F16.2D2_b"
  ) %>%
  ggplot(aes(x = wavenumber, y = abs, col = analysis_id)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse()

# day 24 (final product)
ftir_day24 <- ftir_data %>%
  filter(
    analysis_id == "F24.1A3" |
      analysis_id == "F24.1C2" |
      analysis_id == "F24.1D3"
  ) %>%
  ggplot(aes(x = wavenumber, y = abs, col = analysis_id)) +
  geom_line() +
  theme_classic() +
  scale_x_reverse()

calc_compar <- ftir_data %>%
  filter(
    analysis_id == "ArchDC_MB11_grind_c" | 
      analysis_id == "F24.1A3" | 
      analysis_id == "modern-ref_1"
  ) %>%
  mutate(
    analysis_id = factor(analysis_id, levels = c("F24.1A3", "ArchDC_MB11_grind_c", "modern-ref_1")),
    abs = case_when(analysis_id == "F24.1A3" ~ abs + 0.45,
      TRUE ~ abs)
  ) %>% 
  ggplot(aes(x = wavenumber, y = abs, col = analysis_id)) +
  geom_line() +
  theme_classic() +
  # scale_x_continuous(
  #   n.breaks = length(as.character(seq(500, 4000, by = 500)))
  # ) +
  scale_x_reverse() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

## Grind curves

asscher_enamel_1 <- function(x) {-0.018 * x + 6.0} #y = -0.018x + 6.0 #Qesem   
asscher_enamel_2 <- function(x) {-0.021 * x + 6.1} #y = -0.021x + 6.1 #Ateret 
asscher_enamel_3 <- function(x) {-0.017 * x + 5.5} #y = -0.017x + 5.5 #Neve-Yarak
asscher_enamel_4 <- function(x) {-0.015 * x + 5.2} #y = -0.015x + 5.2 #Modern

# produce plot with grind curve for all samples
grind_all_plot <- ftir_grind_data %>%
  mutate(
    Sample_day = factor(Sample_day, levels = grind_sample_order)) %>%
  group_by(day, Sample) %>%
  ggplot(aes(x = FWHM, y = IRSF, col = Sample_day, shape = Sample_day)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = F) +
  geom_function(inherit.aes = F, fun = asscher_enamel_1, xlim = c(100,200), aes(col = "grey20", linetype = "dotted"), alpha = 0.3) + #Qesem   
  geom_function(inherit.aes = F, fun = asscher_enamel_2, xlim = c(100,200), aes(col = "grey20", linetype = "dotted"), alpha = 0.3) + #Ateret
  geom_function(inherit.aes = F, fun = asscher_enamel_3, xlim = c(100,200), aes(col = "grey20", linetype = "dotted"), alpha = 0.3) + #Neve-Yarak
  geom_function(inherit.aes = F, fun = asscher_enamel_4, xlim = c(100,200), aes(col = "grey20", linetype = "dashed"), alpha = 0.3) + #Modern
  theme_minimal()

# isolate calculus samples to see diffs between days and the 'real deal'
grind_calc_plot <- ftir_grind_data %>%
  filter(Sample == "Artificial calculus" | Sample == "Archaeological calculus") %>%
  ggplot(aes(x = FWHM, y = IRSF, col = Sample_day, shape = Sample_day)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = F) +
  labs(col = "Sample", shape = "Sample",
    x = "FWHM of the 1035 peak",
    y = "Splitting factor") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(
      colour = "grey", 
      fill = "transparent", linewidth = 1),
    panel.background = element_rect(fill = "white")
  ) +
  scale_colour_viridis_d()
