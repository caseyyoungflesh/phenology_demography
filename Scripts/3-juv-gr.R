####################
# 4 - juv ~ gr (sensitivity)
#  
####################


# set dir -----------------------------------------------------------------

dir <- 'XXXX'
run_date <- '2022-04-08'
pc_run_date <- '2022-04-08'


# load packages -----------------------------------------------------------

library(tidyverse)
library(rstan)
library(MCMCvis)


# load data ---------------------------------------------------------------

# #load data used in PC model
fit_data <- readRDS(paste0(dir, 'Results/pro-PC-',
                           pc_run_date, '/pro-PC-data-', pc_run_date, '.rds'))

tdata <- fit_data$pro_data

#sp_id for each st_id
sp_id_df <- dplyr::distinct(tdata, st_id, .keep_all = TRUE) %>%
  dplyr::select(st_id, sp_id)

#center gr each species
tdata2 <- tdata %>%
  dplyr::group_by(sci_name) %>%
  mutate(sc_gr_sp = scale(gr_mid, scale = FALSE)[,1],
         mn_gr_sp = mean(gr_mid)) %>%
  ungroup()


# Run Stan model --------------------------------------------------------------

DATA <- list(N = NROW(tdata2),
             Nst = length(unique(tdata2$st_id)),
             Nsp = length(unique(tdata2$sp_id)),
             juv = tdata2$juv_meanday,
             st_id = tdata2$st_id, #id for each species/station
             sp_id = tdata2$sp_id,
             sp_id2 = sp_id_df$sp_id,
             gr = tdata2$sc_gr_sp,
             pro_data = tdata2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.98
TREE_DEPTH <- 12
STEP_SIZE <- 0.05
CHAINS <- 4
ITER <- 5000

fit <- rstan::stan(paste0(dir, 'Scripts/Model_files/3-juv-gr.stan'),
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'mu_alpha',
                            'mu_beta',
                            'sigma_alpha',
                            'sigma_beta',
                            'sigma',
                            'mma',
                            'mmb',
                            'sma',
                            'smb',
                            'nu',
                            'lambda',
                            'kappa',
                            'y_rep'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))


# save summary ------------------------------------------------------------

setwd(paste0(dir, 'Scripts'))
#save out summary, model fit, data
MCMCvis::MCMCdiag(fit, 
                  round = 4,
                  file_name = paste0('juv-gr-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('juv-gr-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('juv-gr-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('juv-gr-data-', run_date),
                  cp_file = c('Model_files/3-juv-gr.stan', 
                              '3-juv-gr.R'),
                  cp_file_names = c(paste0('3-juv-gr-', run_date, '.stan'),
                                    paste0('3-juv-gr-', run_date, '.R')))

# library(shinystan)
# shinystan::launch_shinystan(fit)


# PPC ---------------------------------------------------------------------

#figure dir
fig_dir <- paste0(dir, 'Results/juv-gr-', run_date, '/')

#posterior predictive check
y_val <- DATA$juv
y_rep <- MCMCvis::MCMCchains(fit, params = 'y_rep')

pdf(paste0(fig_dir, 'juv-gr-PPC-', run_date, '.pdf'), height = 5, width = 5)
bayesplot::ppc_dens_overlay(y_val, y_rep[1:100,])
dev.off()


# cat plots ---------------------------------------------------------------

pdf(paste0(fig_dir, 'juv-sens-', run_date, '.pdf'),
    height = 10, width = 5)
MCMCvis::MCMCplot(fit,
                  params = 'mu_beta',
                  labels = distinct(tdata2, sci_name, sp_id)$sci_name,
                  sz_labels = 0.75,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = 'Fledge sensitivity to green-up',
                  guide_lines = TRUE)
dev.off()


# plot fledge ~ gr sens --------------------------------------------------------------------

mu_alpha_mn <- MCMCvis::MCMCpstr(fit, params = 'mu_alpha')[[1]]
mu_beta_mn <- MCMCvis::MCMCpstr(fit, params = 'mu_beta')[[1]]

gr_sim <- seq(-20, 20, by = 0.5)

#mu_alpha_k + mu_beta_k * gr_sim
#SAME AS
#X %*% beta_vec
xmat <- cbind(rep(1, length(gr_sim)), gr_sim)
bb <- rbind(mu_alpha_mn, mu_beta_mn)
res_mat <- xmat %*% bb

#long format
tplt <- data.frame(sc_gr = xmat[,2], res_mat) %>%
  pivot_longer(cols = -sc_gr) %>%
  arrange(name, sc_gr) %>%
  group_by(name) %>%
  mutate(sc_val = scale(value, scale = FALSE)[,1]) %>%
  ungroup() %>%
  mutate(mu_beta_mn = rep(mu_beta_mn, each = length(gr_sim)))

#lines for background (for reference of slope = 1)
sim_int <- seq(100, 250, by = 3)
bb2 <- rbind(sim_int, 1)
res_mat2 <- xmat %*% bb2
tplt2 <- data.frame(sc_gr = xmat[,2], res_mat2) %>%
  pivot_longer(cols = -sc_gr) %>%
  arrange(name, sc_gr)

plt2 <- ggplot(tplt2, aes(sc_gr, value, group = name)) + 
  #background
  geom_line(size = 0.9, alpha = 0.12, col = 'black') +
  xlim(c(-20, 20)) +
  ylim(range(tplt$value)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18)) +
  geom_line(data = tplt, aes(sc_gr, value, group = name),
            size = 1.3, alpha = 0.7) +
  xlab('Green-up anomaly (days)') +
  ylab('Fledge date (ordinal date)')

#save image
ggsave(plt2, filename = paste0(fig_dir, 'juv-gr-', run_date, '.pdf'),
       height = 5, width = 5)
