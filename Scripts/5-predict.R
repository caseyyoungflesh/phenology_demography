####################################
# 5 - predict productivity in response to gr change
#
####################################


# set dir -----------------------------------------------------------------

dir <- 'XXXX'
run_date <- '2022-04-08'
fl_sens_run_date <- '2022-04-08'
pc_run_date <- '2022-04-08'

#create output dir if it doesn't exist
fig_dir <- paste0(dir, 'Results/predict-', run_date, '/')
ifelse(!dir.exists(fig_dir),
       dir.create(fig_dir),
       FALSE)


# load packages -----------------------------------------------------------

library(boot)
library(MCMCvis)
library(RColorBrewer)


# load and process data ------------------------------------------------------------

# #load fl sens model fit
fl_sens_fit <- readRDS(paste0(dir, 'Results/juv-gr-',
                              fl_sens_run_date, '/juv-gr-fit-', 
                              fl_sens_run_date, '.rds'))

#posterior fledge sens - cross-species
mmb_ch <- MCMCvis::MCMCchains(fl_sens_fit, params = 'mmb')
mu_beta_ch <- MCMCvis::MCMCchains(fl_sens_fit, params = 'mu_beta')

# #load PC model fit
fit <- readRDS(paste0(dir, 'Results/pro-PC-',
                      pc_run_date, '/pro-PC-fit-', pc_run_date, '.rds'))

# #load data used in PC model
fit_data <- readRDS(paste0(dir, 'Results/pro-PC-',
                           pc_run_date, '/pro-PC-data-', pc_run_date, '.rds'))

#get posteriors
mma_ch <- MCMCvis::MCMCchains(fit, 'mma')
mb1_ch <- MCMCvis::MCMCchains(fit, 'mu_beta1')
mb2_ch <- MCMCvis::MCMCchains(fit, 'mu_beta2')
mg1_ch <- MCMCvis::MCMCchains(fit, 'mu_gamma1')
mg2_ch <- MCMCvis::MCMCchains(fit, 'mu_gamma2')
mu_alpha_ch <- MCMCvis::MCMCchains(fit, 'mu_alpha')
beta1_ch <- MCMCvis::MCMCchains(fit, 'beta1')
beta2_ch <- MCMCvis::MCMCchains(fit, 'beta2')
gamma1_ch <- MCMCvis::MCMCchains(fit, 'gamma1')
gamma2_ch <- MCMCvis::MCMCchains(fit, 'gamma2')

#read in migration distance data
mig_dis <- readRDS(paste0(dir, '/Data/mig_distance.rds'))
#resident species
res_sp <- dplyr::filter(mig_dis, dis == 0)$sci_name

#convert predicted PC to poly to feed into model coefs
#use then when next model run has the poly objects saved in DATA
coefs_PC1 <- attr(fit_data$poly_PC1, 'coefs')
coefs_PC2 <- attr(fit_data$poly_PC2, 'coefs')

#prcomp object from 2-pro-PC.R
covs_pca <- fit_data$covs_pca

usp <- unique(fit_data$pro_data$sci_name)

#temp at 2100 - CMIP6 ensemble
CMIP_45_2100 <- read.csv(paste0(dir, 'Data/environment/climate/Ens_45_2100.csv'))
CMIP_85_2100 <- read.csv(paste0(dir, 'Data/environment/climate/Ens_85_2100.csv'))

#historical temp data
hist_2001_2010 <- read.csv(paste0(dir, 'Data/environment/climate/hist_2001_2010.csv'))

#change in temp between 2001-2010 and 2100
(delta_T_45 <- mean(CMIP_45_2100$Tave_sp - hist_2001_2010$Tave_sp))
(delta_T_85 <- mean(CMIP_85_2100$Tave_sp - hist_2001_2010$Tave_sp))


# functions ----------------------------------------------------------------

#convert gr and fl to PC
ph_to_PC_fun <- function(gr_delta, 
                         fl_delta, 
                         pro_pca = covs_pca)
{
  #amount of change in gr and fledge
  pp <- data.frame(sc_gr_mid = gr_delta,
                   sc_juv = fl_delta)
  colnames(pp) <- c('Green-up (stdized)', 'Fledge (stdized)')
  
  #predicted PC values following specified period of time (and assuming rate of pheno change remains constant)
  pred_pc <- predict(pro_pca, pp)
  pred_PC1 <- pred_pc[,1] * -1
  pred_PC2 <- pred_pc[,2]
  
  pp <- list(pred_PC1 = pred_PC1, pred_PC2 = pred_PC2)
  
  return(pp)
}

#output is posterior estimate of overall productivity across species
pro_pred_fun <- function(gr_delta,
                         fl_delta, #distribution or single value
                         c_PC1,
                         c_PC2,
                         pro_pca = covs_pca, #from 2-pro-PC.R
                         mma,
                         mb1,
                         mb2,
                         mg1,
                         mg2)
{
  pp <- ph_to_PC_fun(gr_delta = gr_delta,
                     fl_delta = fl_delta,
                     pro_pca = covs_pca)
  
  pred_PC1_poly <- poly(pp$pred_PC1, 2, raw = FALSE, coefs = c_PC1)
  pred_PC2_poly <- poly(pp$pred_PC2, 2, raw = FALSE, coefs = c_PC1)
  
  #i <- 2
  beta_vec <- cbind(mma, mb1, mb2, mg1, mg2)
  xmat <- rbind(1, pred_PC1_poly[,1],
                pred_PC1_poly[,2], pred_PC2_poly[,1],
                pred_PC2_poly[,2])
  if (length(fl_delta) == 1)
  {
    res_mat2 <- beta_vec %*% xmat
  } else {
    res_mat <- matrix(NA, nrow = NROW(mma), ncol = ncol(xmat))
    for (j in 1:NCOL(xmat))
    {
      #j <- 3
      res_mat[,j] <- beta_vec %*% xmat[,j]
    }
    #sample one col for each row (realization of parameters)
    res_mat2 <- apply(res_mat, 1, function(x) sample(x, 1))
  }
  
  res_iv <- boot::inv.logit(res_mat2)
  
  out <- list(post = res_iv, 
              PC1 = as.numeric(pp$pred_PC1), 
              PC2 = as.numeric(pp$pred_PC2))
  
  return(out)
}


#output is posterior estimate of overall productivity for each species
pro_pred_ss_fun <- function(gr_delta,
                         fl_delta, #distribution or single value
                         c_PC1,
                         c_PC2,
                         pro_pca = covs_pca, #from 2-pro-PC.R
                         mu_alpha,
                         mu_beta1,
                         mu_beta2,
                         mu_gamma1,
                         mu_gamma2)
{
  pp <- ph_to_PC_fun(gr_delta = gr_delta,
                     fl_delta = fl_delta,
                     pro_pca = covs_pca)
  
  pred_PC1_poly <- poly(pp$pred_PC1, 2, raw = FALSE, coefs = c_PC1)
  pred_PC2_poly <- poly(pp$pred_PC2, 2, raw = FALSE, coefs = c_PC1)
  
  #i <- 2
  beta_vec <- cbind(mu_alpha, mu_beta1, mu_beta2, mu_gamma1, mu_gamma2)
  xmat <- rbind(1, pred_PC1_poly[,1],
                pred_PC1_poly[,2], pred_PC2_poly[,1],
                pred_PC2_poly[,2])
  if (length(fl_delta) == 1)
  {
    res_mat2 <- beta_vec %*% xmat
  } else {
    res_mat <- matrix(NA, nrow = NROW(mu_alpha), ncol = ncol(xmat))
    for (j in 1:NCOL(xmat))
    {
      #j <- 3
      res_mat[,j] <- beta_vec %*% xmat[,j]
    }
    #sample one col for each row (realization of parameters)
    res_mat2 <- apply(res_mat, 1, function(x) sample(x, 1))
  }

    res_iv <- boot::inv.logit(res_mat2)
  
  out <- list(post = res_iv, 
              PC1 = as.numeric(pp$pred_PC1), 
              PC2 = as.numeric(pp$pred_PC2))
  
  return(out)
}


# all species - run function ---------------------------------------------------

#pro at 0 years
pred0 <- pro_pred_fun(gr_delta = 0,
                      fl_delta = 0,
                      c_PC1 = coefs_PC1, #coefs from poly fun
                      c_PC2 = coefs_PC2,
                      pro_pca = covs_pca, 
                      mma = mma_ch, #post means for each species
                      mb1 = mb1_ch,
                      mb2 = mb2_ch,
                      mg1 = mg1_ch,
                      mg2 = mg2_ch)

# PREDICTED RANGE OF GREEN-UP CHANGE CONSIDERING DELTA TEMP 
#temp change in 100 years = (SSP2 4.5 = 2.83 degrees C, SSP5 8.5 = 5.22 degrees C)
#associated pheno change (given Wolkovich et al. 2012 pheno response) in 100 years = -18 (-10, -26), -33 (-19, -48)

grs <- 0:-25
tplt <- data.frame(gr = grs,
                   fl = NA,
                   mn = NA,
                   med = NA,
                   lci = NA,
                   uci = NA,
                   PC1 = NA,
                   PC2 = NA,
                   med_delta = NA)
for (i in 1:length(grs))
{
  print(paste0('processing ', i, ' of ', length(grs)))
  
  #i <- 1
  FL_DELTA <- grs[i] * mmb_ch
  predi <- pro_pred_fun(gr_delta = grs[i], #change in gr over period
                        fl_delta = FL_DELTA, #change in fl over period
                        c_PC1 = coefs_PC1, #coefs from poly fun
                        c_PC2 = coefs_PC2,
                        pro_pca = covs_pca, 
                        mma = mma_ch, #post means for each species
                        mb1 = mb1_ch,
                        mb2 = mb2_ch,
                        mg1 = mg1_ch,
                        mg2 = mg2_ch)
  
  #posterior difference vs. year 0
  delta <- predi$post - pred0$post
  prop <- delta / pred0$post
  
  #store in df
  tplt$gr[i] <- grs[i]
  tplt$fl[i] <- median(FL_DELTA)
  tplt$mn[i] <- mean(prop)
  tplt$med[i] <- quantile(prop, probs = 0.5)
  tplt$lci[i] <- quantile(prop, probs = 0.055)
  tplt$uci[i] <- quantile(prop, probs = 0.945)
  tplt$PC1[i] <- median(predi$PC1)
  tplt$PC2[i] <- median(predi$PC2)
  tplt$med_delta[i] <- median(delta)
}

saveRDS(tplt, paste0(fig_dir, 'pred_pro_time_mean_data_', run_date, '.rds'))

#stats at 25-day advance
predi <- pro_pred_fun(gr_delta = grs[26],
                      fl_delta = grs[i] * mmb_ch, #change in fl over period
                      c_PC1 = coefs_PC1, #coefs from poly fun
                      c_PC2 = coefs_PC2,
                      pro_pca = covs_pca, 
                      mma = mma_ch, #post means for each species
                      mb1 = mb1_ch,
                      mb2 = mb2_ch,
                      mg1 = mg1_ch,
                      mg2 = mg2_ch)
delta <- predi$post - pred0$post
prop <- delta / pred0$post
mean(prop)
quantile(prop, probs = c(0.055, 0.945))
sum(prop < 0) / NROW(prop)


# per species - run function ----------------------------------------------

ps_fun <- function(idx)
{
  #pro at 0 years
  pred0 <- pro_pred_ss_fun(gr_delta = 0,
                        fl_delta = 0,
                        c_PC1 = coefs_PC1, #coefs from poly fun
                        c_PC2 = coefs_PC2,
                        pro_pca = covs_pca, 
                        mu_alpha = mu_alpha_ch[,idx],
                        mu_beta1 = beta1_ch[,idx],
                        mu_beta2 = beta2_ch[,idx],
                        mu_gamma1 = gamma1_ch[,idx],
                        mu_gamma2 = gamma2_ch[,idx])

  grs <- 0:-25
  tplt <- data.frame(sci_name = NA,
                     gr = grs,
                     fl = NA,
                     mn = NA,
                     med = NA,
                     lci = NA,
                     uci = NA,
                     PC1 = NA,
                     PC2 = NA,
                     med_delta = NA)
  for (i in 1:length(grs))
  {
    print(paste0('processing ', i, ' of ', length(grs)))
    
    FL_DELTA <- grs[i] * mu_beta_ch[,idx]
    predi <- pro_pred_ss_fun(gr_delta = grs[i], #change in gr over period
                          fl_delta = FL_DELTA, #change in fl over period
                          c_PC1 = coefs_PC1, #coefs from poly fun
                          c_PC2 = coefs_PC2,
                          pro_pca = covs_pca, 
                          mu_alpha = mu_alpha_ch[,idx],
                          mu_beta1 = beta1_ch[,idx],
                          mu_beta2 = beta2_ch[,idx],
                          mu_gamma1 = gamma1_ch[,idx],
                          mu_gamma2 = gamma2_ch[,idx])
    
    #posterior difference vs. year 0
    delta <- predi$post - pred0$post
    prop <- delta / pred0$post
    
    #store in df
    tplt$sci_name[i] <- unique(fit_data$pro_data$sci_name)[idx]
    tplt$gr[i] <- grs[i]
    tplt$fl[i] <- median(FL_DELTA)
    tplt$mn[i] <- mean(prop)
    tplt$med[i] <- quantile(prop, probs = 0.5)
    tplt$lci[i] <- quantile(prop, probs = 0.055)
    tplt$uci[i] <- quantile(prop, probs = 0.945)
    tplt$PC1[i] <- median(predi$PC1)
    tplt$PC2[i] <- median(predi$PC2)
    tplt$med_delta[i] <- median(delta)
  }
  
  return(tplt)
}

tplt_ss <- data.frame(sci_name = rep(NA, 26 * length(usp)),
                      gr = NA,
                      fl = NA,
                      mn = NA,
                      med = NA,
                      lci = NA,
                      uci = NA,
                      PC1 = NA,
                      PC2 = NA,
                      med_delta = NA)
counter <- 1
for (i in 1:length(usp))
{
  #i <- 1
  print(paste0('species: ', i, ' of ', length(usp)))
  tplt_ss[counter:(counter+25),] <- ps_fun(idx = i)
  counter <- counter + 26
}

saveRDS(tplt_ss, paste0(fig_dir, 'pred_pro_time_ss_data_', run_date, '.rds'))


# plot over gr advance --------------------------------------------------------------------

#line transparency
L_ALPHA <- 0.8
#ribbon transparency
R_ALPHA <- 0.3
#line color
med_col <- 'black'

pred_pro_year <- ggplot(tplt, aes(-gr, med)) +
  geom_hline(yintercept = 0,
             linetype = 'dashed',
             size = 1,
             alpha = 0.4) +
  geom_ribbon(aes(x = -gr, ymin = lci, ymax = uci),
              alpha = R_ALPHA, fill = med_col) +
  geom_line(size = 2, alpha = L_ALPHA, col = med_col) +
  xlab('Days green-up advancement') +
  ylab('Projected proportion change in pro') +
  ylim(c(-0.28, 0.1)) + #to match species-specific
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
        axis.ticks.length = unit(0.2, 'cm')) #length of axis tick


ggsave(pred_pro_year, filename = paste0(fig_dir, 'pred_pro_time_', run_date, '.pdf'),
       height = 5,
       width = 5)


# plot species specific ---------------------------------------------------

#get resident species
'%ni%' <- Negate('%in%')
res_tplt_ss <- dplyr::filter(tplt_ss, sci_name %in% res_sp)
mig_tplt_ss <- dplyr::filter(tplt_ss, sci_name %ni% res_sp)

pred_pro_year_ss <- ggplot(mig_tplt_ss, aes(-gr, med, group = sci_name)) +
  geom_line(size = 1.5, alpha = 0.2) +
  geom_line(data = res_tplt_ss,
            alpha = 0.2, size = 1.5, col = 'red') +
  xlab('Days green-up advancement') +
  ylab('Projected proportion change in pro') +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
        axis.ticks.length = unit(0.2, 'cm')) #length of axis tick

ggsave(pred_pro_year_ss, filename = paste0(fig_dir, 'pred_pro_time_ss_', run_date, '.pdf'),
       height = 5,
       width = 5)

#magnitude of responses from species - which has largest and smallest decline
tplt_ss_end <- tplt_ss %>%
  filter(gr == min(grs)) %>%
  arrange(med)

tplt_end <- tplt %>%
  filter(grs == min(grs))

#posterior medians for final change in pro
pnt_plt <- ggplot(tplt_ss_end, aes(0, med)) +
  geom_jitter(shape = 1, alpha = 0.5, size = 4) +
  geom_hline(yintercept = 0,
             linetype = 'dashed',
             size = 1,
             alpha = 0.4) +
  geom_segment(aes(y = tplt_end$lci, yend = tplt_end$uci,
                   x = 0, xend = 0), 
               lineend = "round",
               size = 3, alpha = 0.025) +
    geom_point(aes(0, tplt_end$med), 
             size = 8, alpha = 0.025) +
  xlim(-2, 2) +
  ylab('Proportion change productivity') +
  ylim(c(-0.28, 0.1)) + #to match community-level
  xlab('') +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 18),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
        axis.ticks.length = unit(0.2, 'cm')) #length of axis tick

pnt_plt

ggsave(pnt_plt,
       filename = paste0(fig_dir, 'pnt_plt_ss_', run_date, '.pdf'),
       height = 5,
       width = 5)


# grid of pro (and PC) for different gr and juv ----------------------------------------

grs2 <- seq(-25, 25, by = 0.5)
fls <- seq(-25, 25, by = 0.5) #fledge anomaly
#matrix of proportion increase in prod (to get uncertainty in optimal values)
pr_mat <- array(NA, c(length(fls), length(grs2), 10000))
tplt_grid <- data.frame(gr = rep(NA, length(fls) * length(grs2)),
                   fl = NA,
                   mn = NA,
                   med = NA,
                   lci = NA,
                   uci = NA,
                   PC1 = NA,
                   PC2 = NA,
                   med_delta = NA)
counter <- 1
for (i in 1:length(grs2))
{
  for (j in 1:length(fls))
  {
    #i <- 19
    #j <- 1
    print(paste0('processing grs: ', i, '; fl: ', j))

    FL_DELTA <- fls[j]
    predi <- pro_pred_fun(gr_delta = grs2[i], #change in gr over period
                          fl_delta = FL_DELTA, #change in fl over period
                          c_PC1 = coefs_PC1, #coefs from poly fun
                          c_PC2 = coefs_PC2,
                          pro_pca = covs_pca,
                          mma = mma_ch, #post means for each species
                          mb1 = mb1_ch,
                          mb2 = mb2_ch,
                          mg1 = mg1_ch,
                          mg2 = mg2_ch)

    #posterior difference vs. year 0
    delta <- predi$post - pred0$post
    prop <- delta / pred0$post
  
    #fill posterior matrix
    pr_mat[j,i,] <- prop[, 1]
    
    #store in df
    tplt_grid$gr[counter] <- grs2[i]
    tplt_grid$fl[counter] <- FL_DELTA
    tplt_grid$mn[counter] <- mean(prop)
    tplt_grid$med[counter] <- quantile(prop, probs = 0.5)
    tplt_grid$lci[counter] <- quantile(prop, probs = 0.055)
    tplt_grid$uci[counter] <- quantile(prop, probs = 0.945)
    tplt_grid$PC1[counter] <- predi$PC1
    tplt_grid$PC2[counter] <- predi$PC2
    tplt_grid$med_delta[counter] <- median(delta)

    counter <- counter + 1
  }
}

#optimal values
tplt_grid[which(tplt_grid$mn > 0.011),]
opt <- tplt_grid[which.max(tplt_grid$mn),]

#difference between long-term average and optimum ~ 1%
opt$mn

#get uncertainty in maximized productivity - calc gr and fl values at each posterior iter
gr_opt <- rep(NA, dim(pr_mat)[3])
fl_opt <- rep(NA, dim(pr_mat)[3])
PC1_opt <- rep(NA, dim(pr_mat)[3])
PC2_opt <- rep(NA, dim(pr_mat)[3])
for (i in 1:dim(pr_mat)[3])
{
  #i <- 3485
  tidx <- which(pr_mat[,,i] == max(pr_mat[,,i]), arr.ind = TRUE)
  
  fl_opt[i] <- fls[tidx[1]] #fledge is in rows above
  gr_opt[i] <- grs2[tidx[2]]
  
  #convert to PC
  tpp <- ph_to_PC_fun(gr_delta = gr_opt[i], 
               fl_delta = fl_opt[i], 
               pro_pca = covs_pca)
  PC1_opt[i] <- tpp$pred_PC1
  PC2_opt[i] <- tpp$pred_PC2
}


#plot kde contour of max productivity at each iter
kde_PC <- MASS::kde2d(PC1_opt, PC2_opt)
kde_PH <- MASS::kde2d(gr_opt, fl_opt)

#turn kde into df to plot on ggplot figs
probs <- c(0.89, 0.5)
dx <- diff(kde_PC$x[1:2])
dy <- diff(kde_PC$y[1:2])
sz <- sort(kde_PC$z)
c1 <- cumsum(sz) * dx * dy
dimnames(kde_PC$z) <- list(kde_PC$x, kde_PC$y)
PC_dc <- reshape2::melt(kde_PC$z)
PC_dc$prob <- approx(sz, 1-c1, PC_dc$value)$y

dx <- diff(kde_PH$x[1:2])
dy <- diff(kde_PH$y[1:2])
sz <- sort(kde_PH$z)
c1 <- cumsum(sz) * dx * dy
dimnames(kde_PH$z) <- list(kde_PH$x, kde_PH$y)
PH_dc <- reshape2::melt(kde_PH$z)
PH_dc$prob <- approx(sz, 1-c1, PH_dc$value)$y


#plot pro ~ gr and fl
pro_gr_fl_plt <- ggplot(tplt_grid, aes(gr, fl, fill = mn)) +
  geom_tile() +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(7, 'YlOrRd'))) + 
  geom_contour(data = PH_dc, aes(x = Var1, y = Var2, z = prob),
               color = 'black',
               breaks = probs,
               alpha = 0.8, size = 1,
               inherit.aes = FALSE) +
  geom_segment(data = tplt,
               aes(x = 0,
                   y = 0,
                   xend = gr[NROW(tplt)-1],
                   yend = fl[NROW(tplt)-1]),
               color = 'grey70',
               arrow = arrow(), size = 1.1, 
               inherit.aes = FALSE) +
  geom_point(data = opt, aes(gr, fl),
             col = 'black',
             size = 2, shape = 19) +
  theme_bw() +
  ggtitle('pro ~ gr and fl') +
  xlab('green-up anomaly') +
  ylab('fledge anomaly') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18)) +
  geom_vline(xintercept = 0, size = 1.5, 
             color = 'black', alpha = 0.4,
             linetype = 'dashed') +
  geom_hline(yintercept = 0, size = 1.5, 
             color = 'black', alpha = 0.4,
             linetype = 'dashed')

#plot pro ~ PC1 and PC2
pro_pc_plt <- ggplot(tplt_grid, aes(PC1, PC2, fill = mn)) +
  geom_tile(width = 0.24, height = 0.24) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(7, 'YlOrRd'))) + 
  geom_contour(data = PC_dc, aes(x = Var1, y = Var2, z = prob),
               color = 'black',
               breaks = probs,
               alpha = 0.8, size = 1,
               inherit.aes = FALSE) +
  geom_segment(data = tplt,
               aes(x = 0,
                   y = 0,
                   xend = PC1[NROW(tplt)-1],
                   yend = PC2[NROW(tplt)-1]),
               color = 'grey70',
               arrow = arrow(), size = 1.1, 
               inherit.aes = FALSE) +
  geom_point(data = opt, aes(PC1, PC2),
             col = 'black', 
             size = 2, shape = 19) +
  theme_bw() +
  ggtitle('pro ~ PC1 and PC2') +
  xlab('Early/Late Index') +
  ylab('Asynchrony Index') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18)) +
  geom_vline(xintercept = 0, size = 1.5, 
             color = 'black', alpha = 0.4,
             linetype = 'dashed') +
  geom_hline(yintercept = 0, size = 1.5, 
             color = 'black', alpha = 0.4,
             linetype = 'dashed')


#gr and fl anomaly from observed data
data_avail_gr_fl <- ggplot(fit_data$pro_data, aes(sc_gr_mid, sc_juv)) +
  stat_bin_hex() +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(7, 'Reds')) + 
  theme_bw() +
  xlab('Green-up anomaly') +
  ylab('Fledge anomaly') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

#plot PC1 and PC2 from observed data
data_avail_PC <- ggplot(fit_data$pro_data, aes(PC1, PC2)) +
  stat_bin_hex() +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(7, 'Greens')) + 
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
        axis.title.y = element_text(size = 18)) +
  xlim(c(-5, 5)) +
  ylim(c(-4, 4))


#save plots
ggsave(pro_gr_fl_plt, filename = paste0(fig_dir, 'pro_gr_fl_plt_', run_date, '.pdf'),
       height = 5,
       width = 6)

ggsave(pro_pc_plt, filename = paste0(fig_dir, 'pro_PC_plt_', run_date, '.pdf'),
       height = 5,
       width = 6)

ggsave(data_avail_PC, filename = paste0(fig_dir, 'data_avail_PC_', run_date, '.pdf'),
       height = 5,
       width = 5)

ggsave(data_avail_gr_fl, filename = paste0(fig_dir, 'data_avail_gr_fl_', run_date, '.pdf'),
       height = 5,
       width = 5)


# translation of juv/gr to PC ---------------------------------------------

#amount of change in gr and fledge
pp <- data.frame(sc_gr_mid = grs,
                 sc_juv = grs*median(mmb_ch))
colnames(pp) <- c('Green-up (stdized)', 'Fledge (stdized)')

#predicted PC values following specified period of time (and assuming rate of pheno change remains constant)
PC_df <- data.frame(gr = pp$`Green-up (stdized)`,
                    juv = pp$`Green-up (stdized)`,
                    PC1 = predict(covs_pca, pp)[,1] * -1,
                    PC2 = predict(covs_pca, pp)[,2])


# Observed PC values compared to predictions --------------------------------

#predicted PC at median gr change

pc1_hist <- ggplot(fit_data$pro_data, aes(PC1)) +
  geom_histogram(color = rgb(0,0,0,0.3),
                 alpha = 0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18)) +
  geom_vline(xintercept = min(tplt$PC1), col = 'red',
             linetype = 'dashed', size = 1.5,
             alpha = 0.6)

pc2_hist <- ggplot(fit_data$pro_data, aes(PC2)) +
  geom_histogram(color = rgb(0,0,0,0.3),
                 alpha = 0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18)) +
  geom_vline(xintercept = max(tplt$PC2), col = 'red',
             linetype = 'dashed', size = 1.5,
             alpha = 0.6)

ggsave(pc1_hist, 
       filename = paste0(fig_dir, 'pc1_hist_', run_date, '.pdf'),
       height = 5,
       width = 5)

ggsave(pc2_hist, 
       filename = paste0(fig_dir, 'pc2_hist_', run_date, '.pdf'),
       height = 5,
       width = 5)


