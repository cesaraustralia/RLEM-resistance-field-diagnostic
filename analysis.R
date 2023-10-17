library(tidyverse)

# make trials smaller for faster run time (but less smooth plots)
trials <- 10000 # number of sims for each settings

# vars estimated from data
p_dead_SS <- 0.996 # 98/101 genotyped
p_dead_SR <- 0.996 # 56/59 genotyped
p_dead_RR <- 0.14 # 14/32

# vars for sensitivity analysis
n_mites_screened <- 100 # number of mites bioassayed
p_R <- 0.10 # 10% resistant allele

simulate_bioassay <- function(trials, p_dead_SS, p_dead_SR, p_dead_RR, n_mites_screened, p_R) {
    n_RR <- rbinom(trials, n_mites_screened, p_R^2)
    n_SS <- rbinom(trials, n_mites_screened, (1 - p_R)^2)
    n_SR <- n_mites_screened - n_RR - n_SS

    n_dead_SS <- rbinom(trials, n_SS, p_dead_SS)
    n_dead_SR <- rbinom(trials, n_SR, p_dead_SR)
    n_dead_RR <- rbinom(trials, n_RR, p_dead_RR)

    n_alive_SS <- n_SS - n_dead_SS
    n_alive_SR <- n_SR - n_dead_SR
    n_alive_RR <- n_RR - n_dead_RR

    outcome <- rep("", trials)
    outcome[(n_alive_RR + n_alive_SR + n_alive_SS) > 0 & p_R > 0] <- "true_positive"
    outcome[n_alive_SS > 0 & p_R == 0] <- "false_positive"
    outcome[(n_alive_SS + n_alive_SR + n_alive_RR) == 0 & p_R > 0] <- "false_negative"
    outcome[(n_alive_SS + n_alive_SR + n_alive_RR) == 0 & p_R == 0] <- "true_negative"

    return(tibble(n_mites_screened, n_alive_SS, n_alive_SR, n_alive_RR, p_R, outcome))
}

# run sensitivity analysis
p_R_sensitivity <- c(0, 0.01, 0.05, 0.10, 0.50, 1.0)
n_mites_screened_sensitivity <- c(1:50, seq(50, 500, by = 10))
sims <- list() # init
i <- 0 # init
for (p_R in p_R_sensitivity) {
    for (n_mites_screened in n_mites_screened_sensitivity) {
        i <- i + 1
        sims[[i]] <- simulate_bioassay(trials, p_dead_SS, p_dead_SR, p_dead_RR, n_mites_screened, p_R)
    }
}
d <- bind_rows(sims) %>%
    group_by(n_mites_screened, p_R) %>%
    summarise(
        false_negatives = mean(outcome == "false_negative"),
        false_positives = mean(outcome == "false_positive"),
        true_negatives  = mean(outcome == "true_negative"),
        true_positives  = mean(outcome == "true_positive")
    ) %>%
    pivot_longer(cols = c(false_negatives, false_positives, true_negatives, true_positives), names_to = "outcome", values_to = "mean_rate") %>%
    mutate(resistance = if_else(p_R == 0, "susceptible", "resistant")) %>%
    mutate(`R (%)` = factor(p_R * 100))

ggplot(d, aes(n_mites_screened, mean_rate,
    colour = `R (%)`, linetype = `R (%)`
)) +
    geom_line() +
    facet_wrap(~outcome) +
    theme_bw() +
    scale_color_viridis_d(option = "B", direction = -1, begin = 0.1, end = 0.8)

ggsave(sprintf("plots/sim.p_dead_SS_%s.p_dead_SR_%s.p_dead_RR_%s.png", p_dead_SS, p_dead_SR, p_dead_RR))
