library(tidyverse)

d4 <- read_csv("./plots/fig4.csv")

d4 %>%
    pivot_longer(-time) %>%
    mutate(var = ifelse(grepl("mu", name), "mu", "sd")) %>%
    mutate(Genotype = ifelse(grepl("SS", name), "SS/RS", "RR")) %>%
    select(-name) %>%
    pivot_wider(values_from = value, names_from = var) %>%
    ggplot(aes(x = time)) +
    geom_line(aes(y = mu, linetype = Genotype)) +
    geom_ribbon(aes(ymin = mu - sd, ymax = mu + sd, fill = Genotype), alpha = 0.1) +
    theme_classic() +
    scale_x_continuous(breaks = 0:10, expand = c(0.005, 0)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10), expand = c(0.005, 0)) +
    scale_fill_manual(values = c("black", "black")) +
    theme(
        axis.title = element_text(size = 14, face = "bold", color = "#595959"),
        axis.text = element_text(size = 12, color = "#595959")
    ) +
    xlab("Time (h)") +
    ylab("Cumulative mortality (%)")
ggsave("./Rplot/Figure4.png", width = 8, height = 5)


d5 <- read_csv("./plots/fig5.csv")

d5 %>%
    pivot_longer(-time) %>%
    filter(name != "n1000") %>%
    mutate(Probability = 100 * value) %>%
    mutate(`Number tested` = factor(str_extract(name, "\\d+"), levels = (c(30, 100, 300, 1000)))) %>%
    ggplot(aes(x = time)) +
    geom_line(aes(y = Probability, alpha = `Number tested`, linetype = `Number tested`)) +
    geom_point(aes(y = Probability, alpha = `Number tested`, shape = `Number tested`), size = 3) +
    theme_classic() +
    scale_x_continuous(breaks = 0:10, expand = c(0.005, 0)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10), expand = c(0.005, 0)) +
    theme(
        axis.title = element_text(size = 14, face = "bold", color = "#595959"),
        axis.text = element_text(size = 12, color = "#595959")
    ) +
    xlab("N survivors") +
    ylab("Probability of observing N survivors (%)")
ggsave("./Rplot/Figure5.png", width = 8, height = 5)
