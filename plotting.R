library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

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
w0 <- 2/3
true.utility <- (pE.true - w0 * pT.true)
OBD <- apply(true.utility, 1, which.max)
OBD.df <- data.frame(Scenario = seq_along(OBD), OBD = OBD)
MTD <- apply(abs(pT.true - 0.25), 1, which.min)
MTD.df <- data.frame(Scenario = seq_along(MTD), MTD = MTD)

all_config <- expand.grid(
  Scenario = 1:nrow(pT.true),
  DLT_window = c(1, 2, 3),
  eff_window = c(1, 2, 3), 
  accural_rate = c(1, 2, 3),
  eff_complete = c(0.21, 0.31, 0.41)
)
all_config <- cbind(setting.idx = 1:nrow(all_config), all_config)

scenarios <- (all_config %>%
  filter(Scenario %in% c(1:4) & DLT_window == 2 & eff_window == 3 & accural_rate == 2 & eff_complete == 0.31))$setting.idx



tab_setting <-  read_csv("./results/output_settings.csv") %>% 
  arrange(setting.idx) %>% left_join(., all_config[, c("setting.idx", "Scenario")])

### Setting Curves =============================================

p_setting <- list()
for (ii in 1:8) {
  
  df_tmp <- data.frame(
    Toxicity = pT.true[ii, ], Efficacy = pE.true[ii, ], Utility = true.utility[ii, ]
  )
  
  df_tmp <- df_tmp %>%
    mutate(dose = 1:n()) %>%  
    pivot_longer(cols = c("Toxicity", "Efficacy", "Utility"), names_to = "group", values_to = "value") %>% mutate(group = factor(group))
  
  p_setting[[ii]] <- ggplot(df_tmp, aes(x = dose, y = value, color = group, shape = group, linetype = group)) +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1)) + 
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.25, color = "gray", linetype = "dashed", size = 1) + 
    labs(x = "Dose Level", y = "Probability", title = paste("Scenario", ii)) +
    scale_color_manual(name = "", values = c("Toxicity" = "red", "Efficacy" = "blue", "Utility" = "orange")) +
    scale_linetype_manual(name = "", values = c("Toxicity" = "dashed", "Efficacy" = "dashed", "Utility" = "solid")) +
    scale_shape_manual(name = "", values = c("Toxicity" = 17, "Efficacy" = 16, "Utility" = 15)) +  
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank())+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
          axis.line = element_line(color = "black"),   # add axis lines
          axis.line.y.right = element_line(color = "black"),
          axis.line.x.top = element_line(color = "black"))
  
}


plot_list_clean <- lapply(seq_along(p_setting), function(rr) {
  p <- p_setting[[rr]]
  
  if (rr == 1) {
    p <- p + theme(legend.position = "left")
  }
  
  if (rr == 4) {
    p <- p + 
      scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1), position = "right", name = "Utility")
  } 
  
  if ((rr - 1) %% 4 != 0 & (rr - 1) %% 4 != 3) {
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  p
})

p_setting_all <- wrap_plots(plot_list_clean[1:4], ncol = 4)
# p_setting_all <- wrap_plots(plot_list_clean[1:8], ncol = 4)
p_setting_all

ggsave(file="~/Intern/03slides/figures/settings.jpg",plot = p_setting_all,dpi = 400,width = 12,height = 3)



### OBD Selection ==============================================
methods <- c("OBD.bf", "OBD.ours1", "OBD.titeet", "OBD.bfet")
tab_sel <- read_csv("~/Intern/01project/01_code/backfilling/results/output_selection.csv") %>% 
  arrange(setting.idx) %>% filter(setting.idx %in% scenarios & method %in% methods) %>% 
  left_join(., all_config[, c("setting.idx", "Scenario")], by = "setting.idx") %>% 
  left_join(., OBD.df, by = "Scenario")

df_tmp <- matrix(NA, nrow = nrow(tab_sel), ncol = 4)
for(ii in 1:nrow(tab_sel)) {
  df_tmp[ii, ] <- tab_sel[ii, c(1, 2, unlist(2 + tab_sel[ii,"OBD"]), 9)] %>% unlist()
}
df_tmp <- data.frame(df_tmp) %>%  setNames(c("setting.idx", "method", "obd_sel", "Scenario"))

# ### add BF-BOIN OBD selection
# new_rows <- data.frame(
#   setting.idx = unique(df_tmp$setting.idx), 
#   method = rep("bfboin", nrow(df_tmp)/length(methods)),
#   obd_sel = rep(0, nrow(df_tmp)/length(methods)),
#   Scenario = unique(df_tmp$Scenario)
#   )
# 
# df_tmp <- df_tmp %>% rbind(., new_rows)

df_tmp$obd_sel <- as.numeric(df_tmp$obd_sel)
df_tmp$Scenario<-as.factor(df_tmp$Scenario)
df_tmp$method <- factor(df_tmp$method, levels = c("OBD.bf", "OBD.titeet", "OBD.bfet", "OBD.ours1"))


p1<-ggplot(data=df_tmp, aes(x = Scenario, y = obd_sel, fill = method)) +
  geom_bar(stat="identity",position=position_dodge(width = 0.8),color="black",width = 0.8)+
  geom_text(
    aes(label = round(obd_sel, 1)),  # You can customize the number format
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 65)) + 
  ylab("OBD Selection (%)")+
  scale_fill_manual(
    name = " ", 
    values = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
    labels = c("BF-BOIN", "TITE-BOIN-ET", "BF-BOIN-ET", "Proposed")
  ) +
  theme_classic() + 
  theme(legend.position = "top",
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank()
        )
p1
ggsave(file="~/Intern/03slides/figures/OBD_sel.jpg",plot = p1,dpi = 400,width = 6.5,height = 2.5)

### MTD Selection ==============================================
methods <- c("MTD.bf", "MTD.ours1")
tab_sel <- read_csv("~/Intern/01project/01_code/backfilling/results/output_selection.csv") %>% 
  arrange(setting.idx) %>% filter(setting.idx %in% scenarios & method %in% methods) %>% 
  left_join(., all_config[, c("setting.idx", "Scenario")], by = "setting.idx") %>% 
  left_join(., MTD.df, by = "Scenario")

df_tmp <- matrix(NA, nrow = nrow(tab_sel), ncol = 4)
for(ii in 1:nrow(tab_sel)) {
  df_tmp[ii, ] <- tab_sel[ii, c(1, 2, unlist(2 + tab_sel[ii,"MTD"]), 9)] %>% unlist()
}
df_tmp <- data.frame(df_tmp) %>%  setNames(c("setting.idx", "method", "mtd_sel", "Scenario"))

### add BF-BOIN OBD selection
new_rows1 <- data.frame(
  setting.idx = unique(df_tmp$setting.idx), 
  method = rep("titeet", nrow(df_tmp)/length(methods)),
  mtd_sel = rep(0, nrow(df_tmp)/length(methods)),
  Scenario = unique(df_tmp$Scenario)
)
new_rows2 <- data.frame(
  setting.idx = unique(df_tmp$setting.idx), 
  method = rep("bfet", nrow(df_tmp)/length(methods)),
  mtd_sel = rep(0, nrow(df_tmp)/length(methods)),
  Scenario = unique(df_tmp$Scenario)
)

df_tmp <- df_tmp %>% rbind(., new_rows1) %>% rbind(., new_rows2)

df_tmp$mtd_sel <- as.numeric(df_tmp$mtd_sel)
df_tmp$Scenario<-as.factor(df_tmp$Scenario)
df_tmp$method <- factor(df_tmp$method, levels = c("MTD.bf", "titeet", "bfet", "MTD.ours1"))


p6<-ggplot(data=df_tmp, aes(x = Scenario, y = mtd_sel, fill = method)) +
  geom_bar(stat="identity",position=position_dodge(),color="black",width = 0.8)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  ylab("MTD Selection (%)")+
  scale_fill_manual(
    name = " ", 
    values = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
    labels = c("BF-BOIN", "TITE-BOIN-ET", "BF-BOIN-ET", "Proposed")
  ) +
  theme_classic() + 
  theme(legend.position = "top",
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank()
  )
p6
ggsave(file="~/Intern/03slides/figures/MTD_sel.jpg",plot = p6,dpi = 400,width = 6.5,height = 2.5)

### Average Trial Sample Sizes ==============================================
methods <- c("bf.total", "ours1.total", "titeet", "bfet")
tab_EN <- read_csv("./results/output_EN.csv") %>% 
  arrange(setting.idx) %>% filter(setting.idx %in% scenarios & method %in% methods) %>% 
  left_join(., all_config[, c("setting.idx", "Scenario")], by = "setting.idx") %>% 
  left_join(., OBD.df, by = "Scenario") %>% 
  select(c(method, EN, Scenario))

df_tmp <- data.frame(tab_EN)

df_tmp$EN <- as.numeric(df_tmp$EN)
df_tmp$method <-factor(df_tmp$method, levels = c("bf.total", "titeet", "bfet", "ours1.total"))
df_tmp$Scenario<-as.factor(df_tmp$Scenario)

p2<-ggplot(data=df_tmp, aes(x = Scenario, y = EN, fill = method)) +
  geom_bar(stat="identity",position=position_dodge(width = 0.8),color="black",width = 0.8)+
  geom_text(
    aes(label = round(EN, 1)),  # You can customize the number format
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 45)) + 
  ylab("EN")+
  scale_fill_manual(
    name = " ", 
    values = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
    labels = c("BF-BOIN", "TITE-BOIN-ET", "BF-BOIN-ET", "Proposed")
  ) +
  theme_classic() + 
  theme(legend.position = "top",
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank()
  )
p2
ggsave(file="~/Intern/03slides/figures/EN.jpg",plot = p2,dpi = 400,width = 6.5,height = 2.5)

### Number of Patients at Overdoses =======================================
methods <- c("bf.total", "ours1.total", "titeet", "bfet")
tab_overdose <- read_csv("./results/output_EN.csv") %>% 
  arrange(setting.idx) %>% filter(setting.idx %in% scenarios & method %in% methods) %>% 
  left_join(., all_config[, c("setting.idx", "Scenario")], by = "setting.idx") %>% 
  left_join(., OBD.df, by = "Scenario") %>% 
  select(c(method, EN, EN.Overdose, Scenario)) %>% 
  mutate(perc = EN.Overdose / EN * 100)

df_tmp <- data.frame(tab_overdose)
df_tmp$method <- factor(df_tmp$method, levels = c("bf.total", "titeet", "bfet", "ours1.total"))
df_tmp$Scenario<-as.factor(df_tmp$Scenario)
df_tmp$perc[df_tmp$perc < 0] <- 0

p3<-ggplot(data=df_tmp, aes(x = Scenario, y = perc, fill = method)) +
  geom_bar(stat="identity",position=position_dodge(width = 0.8),color="black",width = 0.8)+
  geom_text(
    aes(label = round(perc, 1)),  # You can customize the number format
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 20)) + 
  ylab("Overdose Pts %")+
  scale_fill_manual(
    name = " ", 
    values = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
    labels = c("BF-BOIN", "TITE-BOIN-ET", "BF-BOIN-ET", "Proposed")
  ) +
  theme_classic() + 
  theme(legend.position = "top",
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank()
  )
p3
ggsave(file="~/Intern/03slides/figures/overdose.jpg",plot = p3,dpi = 400,width = 6.5,height = 2.5)

### Trial Duration =======================================
methods <- c("bf.total", "ours1.total", "titeet", "bfet")
tab_duration <- read_csv("./results/output_EN.csv") %>% 
  arrange(setting.idx) %>% filter(setting.idx %in% scenarios & method %in% methods) %>% 
  left_join(., all_config[, c("setting.idx", "Scenario")], by = "setting.idx") %>% 
  left_join(., OBD.df, by = "Scenario") %>% 
  select(c(method, duration, Scenario)) 

df_tmp <- data.frame(tab_duration)
df_tmp$method <-factor(df_tmp$method, levels = c("bf.total", "titeet", "bfet", "ours1.total"))
df_tmp$Scenario<-as.factor(df_tmp$Scenario)

p4<-ggplot(data=df_tmp, aes(x = Scenario, y = duration, fill = method)) +
  geom_bar(stat="identity",position=position_dodge(width = 0.8),color="black",width = 0.8)+
  geom_text(
    aes(label = round(duration, 1)),  # You can customize the number format
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 40)) + 
  ylab("Duration (month)")+
  scale_fill_manual(
    name = " ", 
    values = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
    labels = c("BF-BOIN", "TITE-BOIN-ET", "BF-BOIN-ET", "Proposed")
  ) +
  theme_classic() + 
  theme(legend.position = "top",
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank()
  )
p4
ggsave(file="~/Intern/03slides/figures/duration.jpg",plot = p4,dpi = 400,width = 6.5,height = 2.5)

### Trial Enrollment Efficiency =======================================
methods <- c("bf.total", "ours1.total", "titeet", "bfet")
tab_enrol <- read_csv("./results/output_EN.csv") %>% 
  arrange(setting.idx) %>% filter(setting.idx %in% scenarios & method %in% methods) %>% 
  left_join(., all_config[, c("setting.idx", "Scenario")], by = "setting.idx") %>% 
  left_join(., OBD.df, by = "Scenario") %>% 
  select(c(method, EN, duration, Scenario)) %>% 
  mutate(efficiency = EN / duration)

df_tmp <- data.frame(tab_enrol)
df_tmp$method <- factor(df_tmp$method, levels = c("bf.total", "titeet", "bfet", "ours1.total"))
df_tmp$Scenario<- factor(df_tmp$Scenario)

p5<-ggplot(data=df_tmp, aes(x = Scenario, y = efficiency, fill = method)) +
  geom_bar(stat="identity",position=position_dodge(),color="black",width = 0.8)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  ylab("Relative Efficiency")+
  scale_fill_manual(
    name = " ", 
    values = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
    labels = c("BF-BOIN", "TITE-BOIN-ET", "BF-BOIN-ET", "Proposed")
  ) +
  theme_classic() + 
  theme(legend.position = "right",
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank()
  )
p5
ggsave(file="~/Intern/03slides/figures/trial_efficiency.jpg",plot = p5,dpi = 400,width = 6.5,height = 2.5)
