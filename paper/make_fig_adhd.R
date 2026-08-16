library(dplyr); library(ggplot2); library(brms); library(cogmod); library(patchwork); library(ggdist)
m <- readRDS("m_ln.rds"); d <- readRDS("d.rds"); sc <- readRDS("sc.rds")
cols <- c(CONTROL = "#3B7DD8", ADHD = "#E4572E")
p <- posterior::as_draws_df(m); sp <- function(x) log1p(exp(x)); minrt <- min(d$RT)

pA <- ggplot(sc, aes(x = 1000*cost, fill = Diagnosis)) +
  geom_density(alpha = .45, colour = NA) +
  geom_vline(data = sc |> group_by(Diagnosis) |> summarise(m = 1000*mean(cost)),
             aes(xintercept = m, colour = Diagnosis), linetype = "dashed", linewidth = .6) +
  scale_fill_manual(values = cols) + scale_colour_manual(values = cols) +
  labs(x = "Switch cost (ms)", y = NULL, subtitle = "A. The conventional index") +
  theme_classic(base_size = 8) +
  theme(legend.position = "top", legend.title = element_blank(), legend.key.size = unit(.3,"cm"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

pB <- ggplot(d, aes(x = RT, fill = Diagnosis)) +
  geom_density(alpha = .45, colour = NA) + scale_fill_manual(values = cols) +
  labs(x = "Reaction time (s)", y = NULL, subtitle = "B. The distributions") +
  theme_classic(base_size = 8) +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

ndt_c <- plogis(p$b_tau_Intercept)*minrt; ndt_a <- plogis(p$b_tau_Intercept+p$b_tau_DiagnosisADHD)*minrt
s_c <- sp(p$b_sigma_Intercept); s_a <- sp(p$b_sigma_Intercept+p$b_sigma_DiagnosisADHD)
dec_c <- exp(p$b_Intercept+s_c^2/2); dec_a <- exp(p$b_Intercept+p$b_DiagnosisADHD+s_a^2/2)
eff <- bind_rows(
  tibble(par = "Total\nmean RT",      v = 1000*((ndt_a+dec_a)-(ndt_c+dec_c))),
  tibble(par = "Non-decision\ntime",  v = 1000*(ndt_a-ndt_c)),
  tibble(par = "Decision\ncomponent", v = 1000*(dec_a-dec_c))
) |> mutate(par = factor(par, levels = c("Total\nmean RT","Non-decision\ntime","Decision\ncomponent")))

pC <- ggplot(eff, aes(x = v, y = par)) +
  stat_halfeye(aes(fill = after_stat(x > 0)), .width = .95, normalize = "groups",
               slab_alpha = .75, point_size = 1.2, interval_size = .6) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35") +
  scale_fill_manual(values = c(`TRUE`="#E4572E", `FALSE`="#3B7DD8"), guide = "none") +
  labs(x = "ADHD \u2212 control (ms)", y = NULL,
       subtitle = "C. What the model separates") +
  theme_classic(base_size = 8)

g <- (pA | pB) / pC + plot_layout(heights = c(1, 1.15))
ggsave("fig_adhd.png", g, width = 6.4, height = 4.3, dpi = 300, bg = "white")
cat("saved\n")
