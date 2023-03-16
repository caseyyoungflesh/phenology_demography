####################
# 2 - pro ~ PC
# 
####################


# set dir -----------------------------------------------------------------

dir <- 'XXXX'
run_date <- '2022-04-08'
gr_date <- '2022-03-10'


# load packages -----------------------------------------------------------

library(tidyverse)
library(rstan)
library(MCMCvis)
library(factoextra)
library(viridis)


# load data ---------------------------------------------------------------

#MAPS data
#only species/stations/years with at least 15 total captures, 5 juvs
#only sites <= 50 degrees lat
#only species/stations with at least 5 years data
#only species with at least 15 station/years data
data <- readRDS(paste0(dir, 'Data/MAPS_data/MAPS-data.rds'))

#gr data
gr_mid_master <- readRDS(paste0(dir, 'Data/environment/processed/',
                                gr_date, '/MidGreenup-10km-', gr_date, '-forest.rds'))


# process data ------------------------------------------------------------

#calculate mean gr and mean n_pix at each site
gr_df <- gr_mid_master %>%
  dplyr::group_by(station) %>%
  dplyr::summarize(mn_site_gr = mean(gr_mn),
                   mn_site_np = median(n_pix)) %>%
  dplyr::filter(mn_site_gr >= 60, #mean day March 1 or later
                mn_site_np >= 150) #at least 150 pixels on avg

#circle area = pi * 10km ^2 = 314.2 km^2
#threshold = 0.5km * 0.5km * 150 = 37.5km^2
#37.5/314.2 ~ 12 % forest cover

#merge
gr_mid_master2 <- gr_mid_master %>%
  dplyr::inner_join(gr_df, by = 'station')

#remove extreme outliers green-up  
gr_mid_master3 <- gr_mid_master2 %>%
  dplyr::group_by(station) %>%
  dplyr::summarize(MAD = mad(gr_mn),
            med = median(gr_mn)) %>%
  dplyr::mutate(UL = med + MAD * 6,
         LL = med - MAD * 6) %>%
  dplyr::left_join(gr_mid_master2, by = 'station') %>%
  dplyr::filter(gr_mn < UL & gr_mn > LL)

#join MAPS and gr data
tdata <- dplyr::full_join(data, gr_mid_master3, by = c('year', 'station')) %>%
  dplyr::filter(!is.na(sci_name)) %>%
  dplyr::group_by(sci_name, station) %>%
  dplyr::mutate(sc_gr_mid = scale(gr_mn, scale = FALSE)[,1],
                sc_juv = scale(juv_meanday, scale = FALSE)[,1]) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(sci_name) %>%
  #scale proportion hours
  dplyr::mutate(sc_ph = scale(prop_hours, scale = FALSE)[,1]) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sci_name) %>%
  dplyr::mutate(sp_id = cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sp_id, station) %>%
  dplyr::mutate(st_id = cur_group_id()) %>%
  dplyr::ungroup()


#179 stations
NROW(dplyr::count(tdata, station))
#2001-2018
range(tdata$year)
#41 species
NROW(dplyr::distinct(tdata, sci_name))
#species/station/years
NROW(tdata)


# PCA ---------------------------------------------------------------------

#cor.test(tdata$sc_gr_mid, tdata$sc_juv)
tcov <- tdata[,c('sc_gr_mid', 'sc_juv')]
colnames(tcov) <- c('Green-up (stdized)', 'Fledge (stdized)')
covs_pca <- prcomp(tcov, center = TRUE, scale. = TRUE)

#multiply by -1 to make + PC1 later date and - PC1 earlier date
tdata$PC1 <- covs_pca$x[,1] * -1
tdata$PC2 <- covs_pca$x[,2]


# Run Stan model --------------------------------------------------------------

#sp_id for each st_id
sp_id_df <- dplyr::distinct(tdata, st_id, .keep_all = TRUE) %>%
  dplyr::select(st_id, sp_id)

#poly data from PCs (orthogonal polynomials)
poly_PC1 <- poly(tdata$PC1, 2, raw = FALSE)
poly_PC2 <- poly(tdata$PC2, 2, raw = FALSE)

#data for model
DATA <- list(N = NROW(tdata),
             Nsp = length(unique(tdata$sp_id)),
             Nst = length(unique(tdata$st_id)),
             y = tdata$njuv,
             ntot = tdata$njuv + tdata$nad,
             sp_id = tdata$sp_id, #species id for each data point
             sp_id2 = sp_id_df$sp_id, #species id for each site
             st_id = tdata$st_id, #species/station id for each data point
             PC1_1 = poly_PC1[,1],
             PC1_2 = poly_PC1[,2],
             PC2_1 = poly_PC2[,1],
             PC2_2 = poly_PC2[,2],
             EF = tdata$sc_ph,
             poly_PC1 = poly_PC1,
             poly_PC2 = poly_PC2,
             covs_pca = covs_pca,
             pro_data = tdata,
             gr_mid_master = gr_mid_master3)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.99
TREE_DEPTH <- 12
STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 5000

t2 <- rstan::stan(paste0(dir, 'Scripts/Model_files/2-pro-pc.stan'),
                  data = DATA,
                  chains = CHAINS,
                  iter = ITER,
                  cores = CHAINS,
                  pars = c('alpha',
                           'mu_alpha',
                           'sigma_alpha',
                           'mma',
                           'sma',
                           'beta1',
                           'mu_beta1',
                           'beta2',
                           'mu_beta2',
                           'gamma1',
                           'mu_gamma1',
                           'gamma2',
                           'mu_gamma2',
                           'kappa',
                           'mu_kappa',
                           'sigma_kappa',
                           'sigma_bg',
                           'Rho_bg',
                           'nu',
                           'epsilon',
                           'mu_se',
                           'sd_se',
                           'sigma_eps',
                           'p',
                           'lpnre'),
                  control = list(adapt_delta = DELTA,
                                 max_treedepth = TREE_DEPTH,
                                 stepsize = STEP_SIZE))


# save summary space ------------------------------------------------------------

setwd(paste0(dir, 'Scripts'))
#save out summary, model fit, data
MCMCvis::MCMCdiag(t2, 
                  round = 4,
                  file_name = paste0('pro-PC-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('pro-PC-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('pro-PC-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('pro-PC-data-', run_date),
                  cp_file = c('Model_files/2-pro-pc.stan', 
                              '2-pro-PC.R'),
                  cp_file_names = c(paste0('2-pro-pc-', run_date, '.stan'),
                                    paste0('2-pro-PC-', run_date, '.R')))

# library(shinystan)
# shinystan::launch_shinystan(t2)


# PCA plot -----------------------------------------------------------------

fig_dir <- paste0(dir, 'Results/pro-PC-', run_date, '/')

covs_pca$rotation[,1] <- covs_pca$rotation[,1] * -1

pdf(paste0(fig_dir,'pca_biplot-', run_date, '.pdf'),
    height = 5.75, width = 5.75)
factoextra::fviz_pca_var(covs_pca,
                         repel = FALSE)
dev.off()


# Summary -----------------------------------------------------------------

#model summary - effect of EL index, effect of AI index
MCMCvis::MCMCsummary(t2, round = 3, 
                     params = c('mu_beta1', 'mu_beta2', 'mu_gamma1', 'mu_gamma2'))


# PPC ---------------------------------------------------------------------

#posterior predictive check (proportion successes)
y_val <- DATA$y
p_rep <- MCMCvis::MCMCchains(t2, params = 'p')
prop_val <- DATA$y / DATA$ntot

pdf(paste0(fig_dir, 'pro_PC_PPC-', run_date, '.pdf'), height = 5, width = 5)
bayesplot::ppc_dens_overlay(prop_val, p_rep[1:100,])
dev.off()


# R^2 ---------------------------------------------------------------------

#Gelman et al. 2019 - American Statistician
#var predicted / (var predicted + var residuals)
#in other words:
#explained var / (explained var + resid var)

# For R^2, we take variance of predicted values, considering intercept and covariates of interest (pheno indices) -> lpnre
# We then take the variance of epsilon. We do not consider kappa, the effect of effort on productivity, because it isn't a qty of interest -> eps_ch
# We calculate these in probability space (take inv.logit transform)
# We then take: y_prep / (y_pred + y_resid)

#for each st_id
lpnre_ch <- MCMCvis::MCMCchains(t2, params = 'lpnre')
eps_ch <- MCMCvis::MCMCchains(t2, params = 'epsilon')
ust_id <- unique(tdata$st_id)

#prob space
out_r2 <- matrix(NA, nrow = NROW(lpnre_ch), ncol = length(ust_id))
for (i in 1:length(ust_id))
{
  #i <- 1
  print(paste0('st_id: ', i, ' of ', length(ust_id)))
  idx <- which(tdata$st_id == ust_id[i])
  sp_idx <- tdata$sp_id[which(tdata$st_id == ust_id[i])][1]
  
  #variation in data - not considering effort and epsilon
  cov_var <- apply(boot::inv.logit(lpnre_ch[, idx]), 1, var)
  #variation contributed by epsilon
  eps_var <- apply(boot::inv.logit(eps_ch[, idx]), 1, var)
  
  out_r2[,i] <- cov_var / (cov_var + eps_var)
}

#species-specific r2
u_sp_id <- unique(tdata$sp_id)
sp_r2_out <- data.frame(sp_id = rep(NA, length(u_sp_id)),
                        sci_name = NA,
                        med_r2 = NA)
for (k in 1:length(u_sp_id))
{
  #k <- 1
  #st_id for each species
  idx <- unique(dplyr::filter(tdata, sp_id == u_sp_id[k])$st_id)
  
  sp_r2_out$sp_id[k] <- u_sp_id[k]
  sp_r2_out$sci_name[k] <- dplyr::filter(tdata, sp_id == u_sp_id[k])$sci_name[1]
  sp_r2_out$med_r2[k] <- median(apply(out_r2[,idx], 2, median))
}

#overall r2 - median of species-level medians
median(sp_r2_out$med_r2)
#species-specific range
range(sp_r2_out$med_r2)



# cat plots ---------------------------------------------------------------

pdf(paste0(fig_dir, 'pro_PC1_linear-', run_date, '.pdf'),
    height = 10, width = 5)
MCMCvis::MCMCplot(t2,
                  params = 'beta1',
                  labels = distinct(tdata, sci_name, sp_id)$sci_name,
                  sz_labels = 0.75,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = 'PC1 - Linear effect',
                  guide_lines = TRUE)
dev.off()
pdf(paste0(fig_dir, '/pro_PC1_quad-', run_date, '.pdf'),
    height = 10, width = 5)
MCMCvis::MCMCplot(t2,
                  params = 'beta2',
                  labels = distinct(tdata, sci_name, sp_id)$sci_name,
                  sz_labels = 0.75,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = 'PC1 - Quadratic effect',
                  guide_lines = TRUE)
dev.off()
pdf(paste0(fig_dir, '/pro_PC2_linear-', run_date, '.pdf'),
    height = 10, width = 5)
MCMCvis::MCMCplot(t2,
                  params = 'gamma1',
                  labels = distinct(tdata, sci_name, sp_id)$sci_name,
                  sz_labels = 0.75,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = 'PC2 - Linear effect',
                  guide_lines = TRUE)
dev.off()
pdf(paste0(fig_dir, '/pro_PC2_quad-', run_date, '.pdf'),
    height = 10, width = 5)
MCMCvis::MCMCplot(t2,
                  params = 'gamma2',
                  labels = distinct(tdata, sci_name, sp_id)$sci_name,
                  sz_labels = 0.75,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = 'PC2 - Quadratic effect',
                  guide_lines = TRUE)
dev.off()


# plot overall effect -----------------------------------------------------

#get coefficients from poly to use in prediction
#https://stats.stackexchange.com/questions/295816/how-to-correctly-predict-from-orthogonalized-covariates-on-new-data-in-polynomia
coefs_PC1 <- attr(poly_PC1, 'coefs')
coefs_PC2 <- attr(poly_PC2, 'coefs')

#simulate data for overall fit
out_all <- array(NA, c(2500, NROW(p_rep)))
out_all_PC1 <- array(NA, c(50, NROW(p_rep)))
out_all_PC2 <- array(NA, c(50, NROW(p_rep)))

opt_all_PC1 <- rep(NA, NROW(p_rep))
opt_all_PC2 <- rep(NA, NROW(p_rep))

new_PC1_all <- seq(range(tdata$PC1)[1], range(tdata$PC1)[2], length.out = 50)
new_PC2_all <- seq(range(tdata$PC2)[1], range(tdata$PC2)[2], length.out = 50)
surface_PC_all <- as.matrix(expand.grid(PC1 = seq(range(tdata$PC1)[1], 
                                                  range(tdata$PC1)[2], length.out = 50), 
                                        PC2 = seq(range(tdata$PC2)[1], 
                                                  range(tdata$PC2)[2], length.out = 50), 
                                        KEEP.OUT.ATTRS = FALSE))

#convert to poly for prediction
poly_PC1_all <- poly(surface_PC_all[, 1], 2, raw = FALSE, coefs = coefs_PC1)
poly_PC2_all <- poly(surface_PC_all[, 2], 2, raw = FALSE, coefs = coefs_PC2)

#just PC1 and just PC2
poly_just_PC1_all <- poly(new_PC1_all, 2, raw = FALSE, coefs = coefs_PC1)
poly_just_PC2_all <- poly(new_PC2_all, 2, raw = FALSE, coefs = coefs_PC2)

mma_ch <- MCMCvis::MCMCchains(t2, params = 'mma')
mu_beta1_ch <- MCMCvis::MCMCchains(t2, params = 'mu_beta1')
mu_beta2_ch <- MCMCvis::MCMCchains(t2, params = 'mu_beta2')
mu_gamma1_ch <- MCMCvis::MCMCchains(t2, params = 'mu_gamma1')
mu_gamma2_ch <- MCMCvis::MCMCchains(t2, params = 'mu_gamma2')

#loop through iter
for (i in 1:NROW(mma_ch))
{
  #i <- 1
  print(paste0('iter : ', i, ' of ', NROW(mma_ch)))
  
  #fill prediction vec
  out_all[,i] <- mma_ch[i,1] +
    mu_beta1_ch[i,1] * poly_PC1_all[,1] + 
    mu_beta2_ch[i,1] * poly_PC1_all[,2] +
    mu_gamma1_ch[i,1] * poly_PC2_all[,1] + 
    mu_gamma2_ch[i,1] * poly_PC2_all[,2]
  
  # #value with PC2 = 0
  out_all_PC1[,i] <- mma_ch[i,1] +
    mu_beta1_ch[i,1] * poly_just_PC1_all[,1] + 
    mu_beta2_ch[i,1] * poly_just_PC1_all[,2]
  
  # #value with PC1 = 0
  out_all_PC2[,i] <- mma_ch[i,1] +
    mu_gamma1_ch[i,1] * poly_just_PC2_all[,1] + 
    mu_gamma2_ch[i,1] * poly_just_PC2_all[,2]
  
  #get 'optimal' productivity for each species
  opt_all_PC1[i] <- new_PC1_all[which.max(out_all_PC1[,i])]
  opt_all_PC2[i] <- new_PC2_all[which.max(out_all_PC2[,i])]
}

#inverse logit to convert to probability space
out_all_il <- boot::inv.logit(out_all)

#just PC1 and just PC2
out_all_PC1_il <- boot::inv.logit(out_all_PC1)
out_all_PC2_il <- boot::inv.logit(out_all_PC2)

#calculate mean and CI for each species/obs
pred_mean_all <- apply(out_all_il, 1, mean)
pred_LCI_all <- apply(out_all_il, 1, function(x) quantile(x, probs = c(0.025)))
pred_UCI_all <- apply(out_all_il, 1, function(x) quantile(x, probs = c(0.975)))

pred_mean_all_PC1 <- apply(out_all_PC1_il, 1, mean)
pred_LCI_all_PC1 <- apply(out_all_PC1_il, 1, function(x) quantile(x, probs = c(0.025)))
pred_UCI_all_PC1 <- apply(out_all_PC1_il, 1, function(x) quantile(x, probs = c(0.975)))

pred_mean_all_PC2 <- apply(out_all_PC2_il, 1, mean)
pred_LCI_all_PC2 <- apply(out_all_PC2_il, 1, function(x) quantile(x, probs = c(0.025)))
pred_UCI_all_PC2 <- apply(out_all_PC2_il, 1, function(x) quantile(x, probs = c(0.975)))

pred_tplt_all <- data.frame(NA, NA, pred_mean_all,
                            pred_LCI_all, 
                            pred_UCI_all)
colnames(pred_tplt_all) <- c('PC1', 'PC2', 'pro_index_mean',
                             'pro_index_LCI', 'pro_index_UCI')
pred_tplt_all$PC1 <- surface_PC_all[,1]
pred_tplt_all$PC2 <- surface_PC_all[,2]

#raster/contour plot
rc_plt <- ggplot(pred_tplt_all, aes(PC1, PC2)) +
  geom_raster(aes(fill = pro_index_mean)) +
  scale_fill_viridis() +
  geom_contour(aes(PC1, PC2, z = pro_index_mean),
               inherit.aes = FALSE,
               binwidth = 0.03, color = 'white',
               alpha = 0.6, size = 0.8) +
  theme_bw() +
  xlab('PC1') +
  ylab('PC2') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

ggsave(filename = paste0(fig_dir, '/pro_PC1_PC2_fit-', run_date, '.pdf'),
       rc_plt, width = 6, height = 5)


# plot per-species --------------------------------------------------------

#simulate data to plot fit per species - sim data only over ranges of actual data for each species
mu_alpha_ch <- MCMCvis::MCMCchains(t2, params = 'mu_alpha')
beta1_ch <- MCMCvis::MCMCchains(t2, params = 'beta1')
beta2_ch <- MCMCvis::MCMCchains(t2, params = 'beta2')
gamma1_ch <- MCMCvis::MCMCchains(t2, params = 'gamma1')
gamma2_ch <- MCMCvis::MCMCchains(t2, params = 'gamma2')

#num species
usp <- unique(tdata$sp_id)
#dims: obs, species, iteration
new_PC1 <- rep(NA, (50 * length(usp)))
new_PC2 <- rep(NA, (50 * length(usp)))
surface_PC <- matrix(NA, nrow = (50 * 50) * length(usp), ncol = 2)

out <- array(NA, c(2500, length(usp), NROW(mu_alpha_ch)))
out_PC1 <- array(NA, c(50, length(usp), NROW(mu_alpha_ch)))
out_PC2 <- array(NA, c(50, length(usp), NROW(mu_alpha_ch)))
#dims: species, iteration
opt_PC1 <- array(NA, c(length(usp), NROW(mu_alpha_ch)))
opt_PC2 <- array(NA, c(length(usp), NROW(mu_alpha_ch)))
counter <- 1
counter2 <- 1

#loop through each species
for (k in 1:length(usp))
{
  #k <- 1
  print(paste0('species: ', k, ' of ', length(usp)))
  
  #async index to predict on
  df <- dplyr::filter(tdata, sp_id == usp[k])

  new_PC1[counter:(counter+49)] <- seq(range(df$PC1)[1], 
                                       range(df$PC1)[2], length.out = 50)
  new_PC2[counter:(counter+49)] <- seq(range(df$PC2)[1], 
                                       range(df$PC2)[2], length.out = 50)
  surface_PC[counter2:(counter2+2499), ] <- as.matrix(expand.grid(PC1 = seq(range(df$PC1)[1], 
                                                                            range(df$PC1)[2],
                                                                            length.out = 50), 
                                                                  PC2 = seq(range(df$PC2)[1], 
                                                                            range(df$PC2)[2], 
                                                                            length.out = 50), 
                                                                  KEEP.OUT.ATTRS = FALSE))
  
  #convert to poly for prediction
  poly_PC1 <- poly(surface_PC[counter2:(counter2+2499), 1], 2, raw = FALSE, coefs = coefs_PC1)
  poly_PC2 <- poly(surface_PC[counter2:(counter2+2499), 2], 2, raw = FALSE, coefs = coefs_PC2)
  
  #just PC1 and just PC2
  poly_just_PC1 <- poly(new_PC1[counter:(counter+49)], 2, raw = FALSE, coefs = coefs_PC1)
  poly_just_PC2 <- poly(new_PC2[counter:(counter+49)], 2, raw = FALSE, coefs = coefs_PC2)
  
  #loop through iter
  for (i in 1:NROW(mma_ch))
  {
    #i <- 1
    #print(paste0('iter : ', i, ' of ', NROW(sp_int_ch)))
    
    #fill prediction vec
    out[,k,i] <- mu_alpha_ch[i,k] +
      beta1_ch[i,k] * poly_PC1[,1] + 
      beta2_ch[i,k] * poly_PC1[,2] + 
      gamma1_ch[i,k] * poly_PC2[,1] + 
      gamma2_ch[i,k] * poly_PC2[,2]
    
    #value with PC2 = 0
    out_PC1[,k,i] <- mu_alpha_ch[i,k] +
      beta1_ch[i,k] * poly_just_PC1[,1] + 
      beta2_ch[i,k] * poly_just_PC1[,2]
    
    #value with PC1 = 0
    out_PC2[,k,i] <- mu_alpha_ch[i,k] +
      gamma1_ch[i,k] * poly_just_PC2[,1] + 
      gamma2_ch[i,k] * poly_just_PC2[,2]
    
    #get 'optimal' productivity for each species
    opt_PC1[k,i] <- new_PC1[counter:(counter+49)][which.max(out_PC1[,k,i])]
    opt_PC2[k,i] <- new_PC2[counter:(counter+49)][which.max(out_PC2[,k,i])]
  }
  counter <- counter + 50
  counter2 <- counter2 + 2500
}

#inverse logit to convert to probability space
out_il <- boot::inv.logit(out)
#just PC1 and just PC2
out_PC1_il <- boot::inv.logit(out_PC1)
out_PC2_il <- boot::inv.logit(out_PC2)

#calculate mean and CI for each species/obs
pred_mean <- apply(out_il, c(1,2), mean)
pred_LCI <- apply(out_il, c(1,2), function(x) quantile(x, probs = c(0.025)))
pred_UCI <- apply(out_il, c(1,2), function(x) quantile(x, probs = c(0.975)))

pred_mean_PC1 <- apply(out_PC1_il, c(1,2), mean)
pred_LCI_PC1 <- apply(out_PC1_il, c(1,2), function(x) quantile(x, probs = c(0.025)))
pred_UCI_PC1 <- apply(out_PC1_il, c(1,2), function(x) quantile(x, probs = c(0.975)))

pred_mean_PC2 <- apply(out_PC2_il, c(1,2), mean)
pred_LCI_PC2 <- apply(out_PC2_il, c(1,2), function(x) quantile(x, probs = c(0.025)))
pred_UCI_PC2 <- apply(out_PC2_il, c(1,2), function(x) quantile(x, probs = c(0.975)))

#wide to long format
pred_mean_long <- reshape2::melt(pred_mean)
pred_mean_PC1_long <- reshape2::melt(pred_mean_PC1)
pred_mean_PC2_long <- reshape2::melt(pred_mean_PC2)

pred_tplt <- cbind(NA, pred_mean_long,
                   reshape2::melt(pred_LCI)[,3], 
                   reshape2::melt(pred_UCI)[,3])
colnames(pred_tplt) <- c('PC1', 'PC2', 'sci_name', 'pro_index_mean',
                         'pro_index_LCI', 'pro_index_UCI')
pred_tplt$PC1 <- surface_PC[,1]
pred_tplt$PC2 <- surface_PC[,2]
pred_tplt$sci_name <- rep(unique(tdata$sci_name), each = 2500)
#just PC1
pred_tplt_PC1 <- cbind(pred_mean_PC1_long,
                       reshape2::melt(pred_LCI_PC1)[,3], 
                       reshape2::melt(pred_UCI_PC1)[,3])
colnames(pred_tplt_PC1) <- c('PC1', 'sci_name', 'pro_index_mean',
                             'pro_index_LCI', 'pro_index_UCI')
pred_tplt_PC1$PC1 <- new_PC1
pred_tplt_PC1$sci_name <- rep(unique(tdata$sci_name), each = 50)
#just PC2
pred_tplt_PC2 <- cbind(pred_mean_PC2_long,
                       reshape2::melt(pred_LCI_PC2)[,3], 
                       reshape2::melt(pred_UCI_PC2)[,3])
colnames(pred_tplt_PC2) <- c('PC2', 'sci_name', 'pro_index_mean',
                             'pro_index_LCI', 'pro_index_UCI')
pred_tplt_PC2$PC2 <- new_PC2
pred_tplt_PC2$sci_name <- rep(unique(tdata$sci_name), each = 50)

pred_tplt_PC1 %>%
  group_by(sci_name) %>%
  summarize(dd = round(sum(range(diff(pro_index_mean))), 4))

p1 <- ggplot(pred_tplt_PC1, aes(PC1, pro_index_mean, by = factor(sci_name))) +
  geom_line(alpha = 0.3, size = 1.5) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('Breeding productivity index') +
  xlab('PC1: early (-) -> late (+)') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18)) +
  ylim(c(0.1, 0.9)) +
  xlim(c(-5, 5))

p2 <- ggplot(pred_tplt_PC2, aes(PC2, pro_index_mean, by = factor(sci_name))) +
  geom_line(alpha = 0.3, size = 1.5) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylab('Breeding productivity index') +
  xlab('PC2: early juv/late gr (-) -> late juv/early gr (+)') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18)) +
  ylim(c(0.1, 0.9)) +
  xlim(c(-5, 5))

ggsave(p1, filename = paste0(fig_dir, 'pro_juv_gr_PC1_sp_fits-', run_date, '.pdf'),
       width = 5, height = 5)
ggsave(p2, filename = paste0(fig_dir, 'pro_juv_gr_PC2_sp_fits-', run_date, '.pdf'),
       width = 5, height = 5)
