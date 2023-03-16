####################################
# 4 - juv_sens ~ pc (dis and fledge date) + phylo
#
####################################


# set dir -----------------------------------------------------------------

dir <- 'XXXX'
run_date <- '2023-03-16'
fl_sens_run_date <- '2023-03-16'
pc_run_date <- '2022-04-08'
phylo_date <- '2022-09-15'


# load packages -----------------------------------------------------------

library(tidyverse)
library(rstan)
library(MCMCvis)
library(factoextra)
library(ape)
library(phytools)


# load data and process ------------------------------------------------------------

#read in pro PC model
tfit <- readRDS(paste0(dir, '/Results/pro-PC-', pc_run_date, 
                       '/pro-PC-fit-', pc_run_date, '.rds'))

#read in pro PC data
tdata <- readRDS(paste0(dir, '/Results/pro-PC-', pc_run_date, 
                        '/pro-PC-data-', pc_run_date, '.rds'))$pro_data

#read in migration distances (from 4a-process-range-maps.R)
mig_dis <- readRDS(paste0(dir, '/Data/mig_distance.rds'))

#read in juv sens
t2 <- readRDS(paste0(dir, '/Results/juv-gr-unc-', 
                     fl_sens_run_date, '/juv-gr-fit-', 
                     fl_sens_run_date, '.rds'))

#sensitivity to gr
mb_mn <- MCMCvis::MCMCpstr(t2, params = 'mu_beta')[[1]]
mb_sd <- MCMCvis::MCMCpstr(t2, params = 'mu_beta', fun = sd)[[1]]

#mean fledge date for each species
mn_juv <- tdata %>%
  group_by(sci_name, station) %>%
  summarize(mn_juv = mean(juv_meanday)) %>%
  ungroup() %>%
  group_by(sci_name) %>%
  summarize(mjd = mean(mn_juv)) %>%
  ungroup()

#merge
tmrg <- dplyr::distinct(tdata, sci_name, sp_id) %>%
  dplyr::left_join(mig_dis, by = 'sci_name') %>%
  dplyr::mutate(juv_sens = mb_mn, juv_sens_sd = mb_sd) %>%
  dplyr::left_join(mn_juv, by = 'sci_name') %>%
  dplyr::select(sci_name, sp_id, juv_sens, juv_sens_sd, dis, mjd)

#sens of mig species
mean(dplyr::filter(tmrg, dis > 0)$juv_sens)
mean(dplyr::filter(tmrg, dis == 0)$juv_sens)


# PCA ---------------------------------------------------------------------

covs <- dplyr::select(tmrg, dis, mjd) * -1
colnames(covs) <- c('Mig Distance', 'Mean Fledge Date')
cor(covs)

#PCA on variables
covs_pca <- prcomp(covs, center = TRUE, scale. = TRUE)
summary(covs_pca)


# phylo data --------------------------------------------------------------

#read in phylo species names
setwd(paste0(dir, '/Data/bird_phylo'))
bn <- read.table(paste0('phylo_names-', phylo_date, '.txt'), sep = ',')[,1]
bn2 <- gsub(' ', '_', bn)
bn3 <- data.frame(name = bn2, num = 1:length(bn2))

setwd(paste0(dir, "/Data/bird_phylo/"))
pnk <- read.csv(paste0('phylo_names_key-', phylo_date, '.csv'))

#read in phylo tree
setwd(paste0(dir, '/Data/bird_phylo/tree-pruner-9cd43878-f76b-48ec-aa78-a75a03fae18a'))
ptree <- ape::read.nexus('output.nex')

#consensus tree
tree <- phytools::consensus.edges(ptree)

#calculate covariance matrix
V <- ape::vcv.phylo(tree)

#scale by max variance -> correlation matrix
R <- V[bn2, bn2] / max(V)


# stan data ---------------------------------------------------------------

#create data list for Stan
DATA <- list(juv_sens = tmrg$juv_sens,
             juv_sens_sd = tmrg$juv_sens_sd,
             cov = covs_pca$x[,1] * -1,
             R = R,
             I = diag(1, nrow = NROW(tmrg), ncol = NROW(tmrg)),
             N = NROW(tmrg))


# Call Stan model space --------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.93
TREE_DEPTH <- 12
STEP_SIZE <- 0.1
CHAINS <- 4
ITER <- 5000

fit1 <- rstan::stan(paste0(dir, 'Scripts/Model_files/4-sens-juv-pc.stan'),
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'lambda',
                            'sigma',
                            'mu_js',
                            'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))


# save summary space ------------------------------------------------------------

setwd(paste0(dir, '/Scripts'))
#save out summary, model fit, data
MCMCvis::MCMCdiag(fit1, 
                  round = 4,
                  file_name = paste0('sens-juv-pc-results-', run_date),
                  dir = paste0(dir, '/Results'),
                  mkdir = paste0('sens-juv-pc-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('sens-juv-pc-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('sens-juv-pc-data-', run_date),
                  cp_file = c('Model_files/4-sens-juv-pc.stan', 
                              '4-sens-juv-pc.R'),
                  cp_file_names = c(paste0('4-sens-juv-pc-', run_date, '.stan'),
                                    paste0('4-sens-juv-pc-', run_date, '.R')))

# library(shinystan)
# shinystan::launch_shinystan(fit1)


# PPC ---------------------------------------------------------------------

y_val <- DATA$juv_sens
y_rep <- MCMCvis::MCMCchains(fit1, params = 'y_rep')

pdf(paste0(dir, 'Results/sens-juv-pc-', run_date, '/PPC.pdf'), height = 5, width = 5)
bayesplot::ppc_dens_overlay(y_val, y_rep[1:100,])
dev.off()


# plots --------------------------------------------------------------------

#which species are resident (migratory distance = 0)
tmrg$PC <- DATA$cov
mig_tmrg <- dplyr::filter(tmrg, dis > 0)
nm_tmrg <- dplyr::filter(tmrg, dis == 0)

#model fit
post_ch <- MCMCvis::MCMCchains(fit1, params = c('alpha', 'beta'))

#latent juv sens
post_js_mn <- MCMCvis::MCMCpstr(fit1, params = c('mu_js'), fun = mean)[[1]]
post_js_sd <- MCMCvis::MCMCpstr(fit1, params = c('mu_js'), fun = sd)[[1]]

rr <- range(tmrg$PC)

xsim <- seq(rr[1] - 0.25, rr[2] + 0.25, by = 0.01)
ysim <- matrix(NA, nrow = NROW(post_ch), ncol = length(xsim))
for (i in 1:NROW(post_ch))
{
  #i <- 1
  ysim[i,] <- post_ch[i,1] + post_ch[i,2] * xsim
}

FIT_PLOT<- data.frame(xsim = xsim,
                      LCI = apply(ysim, 2, function(x) quantile(x, probs = 0.055)),
                      UCI = apply(ysim, 2, function(x) quantile(x, probs = 0.945)),
                      MED = apply(ysim, 2, function(x) quantile(x, probs = 0.5)))

#point with
plt1 <- ggplot() +
  geom_ribbon(data = FIT_PLOT,
              aes(x = xsim, ymin = LCI, ymax = UCI),
              fill = 'grey', alpha = 0.7,
              inherit.aes = FALSE) +
  geom_errorbar(data = mig_tmrg,
                aes(x = PC, 
                    ymin = juv_sens - juv_sens_sd,
                    ymax = juv_sens + juv_sens_sd),
                width = 0.1, size = 0.8,
                color = 'black', alpha = 0.2) +
  geom_errorbar(data = nm_tmrg,
                aes(x = PC, 
                    ymin = juv_sens - juv_sens_sd,
                    ymax = juv_sens + juv_sens_sd),
                width = 0.1, size = 0.8,
                color = 'black', alpha = 0.2) +
  geom_point(data = mig_tmrg, aes(PC, juv_sens),
               color = 'black', alpha = 0.6, size = 3) +
  geom_point(data = nm_tmrg, aes(PC, juv_sens),
             color = 'black', alpha = 0.6, size = 3.25) +
  geom_point(data = nm_tmrg, aes(PC, juv_sens),
             color = 'white', alpha = 1, size = 2) +
  geom_line(data = FIT_PLOT, aes(xsim, MED), color = 'red',
            alpha = 0.9,
            inherit.aes = FALSE,
            size = 1.4) +
  theme_bw() +
  xlab('PC (Migration Distance and Fledge Date)') +
  ylab('Breeding sensitivity (days / day)') +
  ylim(c(0.1, 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(size = 1.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
        axis.ticks.length = unit(0.2, 'cm')) #length of axis tick

#save to pdf
ggsave(plt1, filename = paste0(dir, 'Results/sens-juv-pc-', 
                                          run_date, '/sens-juv-pc-', 
                                          run_date, '.pdf'),
       height = 5, width = 5.5)


# PCA plot ----------------------------------------------------------------

pdf(paste0(dir, 'Results/sens-juv-pc-', run_date,
           '/pca_biplot-', run_date, '.pdf'), height = 5.75, width = 5.75)
factoextra::fviz_pca_var(covs_pca,
                         repel = FALSE)
dev.off()
