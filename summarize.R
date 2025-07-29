library(openxlsx)
library(dplyr)
library(readr)
### Settings ================================
pT.true <- rbind(
  c(0.1, 0.18, 0.35, 0.40, 0.50),
  c(0.05, 0.15, 0.25, 0.35, 0.50),
  c(0.02, 0.06, 0.10, 0.20, 0.35),
  c(0.01, 0.03, 0.05, 0.12, 0.22),
  c(0.10, 0.20, 0.35, 0.43, 0.50),
  c(0.02, 0.06, 0.10, 0.20, 0.35),
  c(0.05, 0.10, 0.20, 0.35, 0.40),
  c(0.01, 0.05, 0.15, 0.18, 0.35)
)
pE.true <- rbind(
  c(0.35, 0.35, 0.37, 0.39, 0.39),
  c(0.10, 0.35, 0.35, 0.38, 0.39),
  c(0.05, 0.10, 0.35, 0.35, 0.40),
  c(0.05, 0.10, 0.15, 0.35, 0.36),
  c(0.10, 0.36, 0.37, 0.40, 0.41),
  c(0.05, 0.10, 0.15, 0.35, 0.37),
  c(0.35, 0.36, 0.37, 0.40, 0.41),
  c(0.05, 0.35, 0.36, 0.37, 0.38)
)

all_config <- expand.grid(
  Scenarrio = 1:nrow(pT.true),
  DLT_window = c(1, 2, 3),
  eff_window = c(1, 2, 3), 
  accural_rate = c(1, 2, 3)
)
all_config <- cbind(setting.idx = 1:nrow(all_config), all_config)

tab_sel <- read_csv("./results/output_selection.csv") %>% 
  arrange(setting.idx) %>% left_join(all_config, ., by = "setting.idx")
tab_EN <- read_csv("./results/output_EN.csv") %>% 
  arrange(setting.idx) %>% left_join(all_config, ., by = "setting.idx")
tab_settings <- read_csv("./results/output_settings.csv") %>% 
  arrange(setting.idx) %>% left_join(all_config, ., by = "setting.idx")

file_output <- list(
  settings = tab_settings, selection = tab_sel, EN = tab_EN
)
write.xlsx(file_output, file = "./results/0729summary.xlsx", asTable = TRUE)
