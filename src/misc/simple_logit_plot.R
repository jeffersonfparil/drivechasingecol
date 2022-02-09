### Accompanying function for lostic regression plot:
### calculate the area under the receiver operating curve, i.e. FPR vs TPR curve
#' @export
simple_auc = function(y_binary, y_continuous, seed=123){
  ### to prevent insufficient memory leading to R being killed, limit the size to a maximum of 30,000
  if (length(y_binary)>30000){
    set.seed(seed)
    idx_rand = sample(x=c(1:length(y_binary)), size=30000)
    y_binary = y_binary[idx_rand]
    y_continuous = y_continuous[idx_rand]
  }
  idx_sort = order(y_continuous, decreasing=TRUE)
  y_binary = y_binary[idx_sort]
  y_continuous = y_continuous[idx_sort]
  positives = y_continuous[y_binary]
  negatives = y_continuous[!y_binary]
  ### if positives[i] == negatives[j] then 1/2
  ### if positives[i] >  negatives[j] then 2/2
  ### if positives[i] <  negatives[j] then 0/2
  M = outer(length(positives):1, 1:length(negatives),
            FUN=function(i, j) (1 + sign(positives[i]-negatives[j]))/2)
  ### above equivalent to below but above is way faster:
  # M = matrix(NA, nrow=sum(y_binary), ncol=sum(!y_binary))
  # for (i in 1:sum(y_binary)){
  #   for (j in 1:sum(!y_binary)){
  #     M[i,j] = 1 + sign(positives[i]-negatives[j])/2
  #   }
  # }
  auc = mean(M)
  return(auc)
}
### plot the logistic regression curve with the distribution of the binary response variable along the continuous explanatiry variable
#' @export
simple_logit_plot = function(x, y, level=0.95, col_fit="black", lwd_fit=1, col_ci_ribbon="#fc9272", lwd_ci_ribbon=1, col_hist_top="#af8dc3", col_hist_bottom="#7fbf7b", col_hist_bord=NA, main=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, zlim=NULL, x_ticks_len=6, y_ticks_len=6, z_ticks_len=6, x_digits=2, y_digits=2, z_digits=2, show_auc=TRUE, show_legend=FALSE){
  ### test ##########################
  # v_mu = apply(cbind(dat$v_x_Drive_mu, dat$v_y_Drive_mu), MAR=1, FUN=mean)
  # drive_loss = dat$drive_loss
  # x = v_mu
  # y = drive_loss
  # level=0.95
  # col_fit="black"
  # lwd_fit=1
  # col_ci_ribbon="#fc9272"
  # lwd_ci_ribbon=1
  # col_hist="#bdbdbd"
  # col_hist_bord=NA
  # main=NULL
  # xlab=NULL
  # ylab=NULL
  # xlim=NULL
  # ylim=NULL
  # zlim=NULL
  # x_ticks_len=6
  # y_ticks_len=6
  # z_ticks_len=6
  # x_digits=2
  # y_digits=2
  # z_digits=2
  # show_auc=TRUE
  # show_legend=FALSE
  # ####################################
  if (is.null(main)){
    par(mar=c(5,5,1,5))
  } else {
    par(mar=c(5,5,5,5))
  }
  ### prepare plot area
  plot(x, y, type="n", yaxt="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
  grid()
  # axis(side=1, at=seq(min(c(0,x),na.rm=TRUE),max(x,na.rm=TRUE),length=x_ticks_len), lab=formatC(seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=x_ticks_len), format="f", flag="0", digits=x_digits), las=1)
  axis(side=2, at=seq(0,1,length=y_ticks_len), lab=formatC(seq(0,1,length=y_ticks_len), format="f", flag="0", digits=y_digits), las=2)
  ### sort x and y values for ease of plotting lines
  idx = order(x)
  x = x[idx]
  y = y[idx]
  ### remove missing data
  idx = (!is.na(x)) & (!is.na(y))
  x = x[idx]
  y = y[idx]
  ### plot hist
  vec_x_breaks = hist(x, plot=FALSE)$br
  vec_y0_counts = hist(x[y==0], breaks=vec_x_breaks, plot=FALSE)$counts
  vec_y1_counts = hist(x[y==1], breaks=vec_x_breaks, plot=FALSE)$counts
  y_counts_max = max(vec_y0_counts) + max(vec_y1_counts)
  ### bottom hist
  for (i in 1:length(vec_y0_counts)){
    vec_x = c(vec_x_breaks[i], vec_x_breaks[i+1], vec_x_breaks[i+1], vec_x_breaks[i])
    vec_y = c(vec_y0_counts[i], vec_y0_counts[i], 0, 0) / y_counts_max
    polygon(x=vec_x, y=vec_y, col=col_hist_bottom, bord=col_hist_bord)
  }
  ### top hist
  for (i in 1:length(vec_y1_counts)){
    vec_x = c(vec_x_breaks[i], vec_x_breaks[i+1], vec_x_breaks[i+1], vec_x_breaks[i])
    vec_y = 1 - c(0, 0, vec_y1_counts[i], vec_y1_counts[i]) / y_counts_max
    polygon(x=vec_x, y=vec_y, col=col_hist_top, bord=col_hist_bord)
  }
  ### add hist axis (z-axis)
  if (y_counts_max>1000){
    nth = length(strsplit(as.character(y_counts_max), "")[[1]])-1
    nth_factor = as.numeric(paste0(c("1", rep("0", times=nth)), collapse=""))
    hist_axis_lab = seq(0,y_counts_max,length=z_ticks_len)/nth_factor
    hist_axis_lab = formatC(hist_axis_lab, format="f", flag="0", digits=z_digits)
    hist_axis_text = paste0("Count (x ", nth_factor, ")")
  } else {
    hist_axis_lab = seq(0,y_counts_max,length=z_ticks_len)
    hist_axis_lab = formatC(hist_axis_lab, format="f", flag="0", digits=z_digits)
    hist_axis_text = "Count"
  }
  axis(side=4, at=seq(0,1,length=z_ticks_len), lab=hist_axis_lab, las=2)
  mtext(side=4, padj=10*par()$adj, text=hist_axis_text, cex=par()$cex)
  ### model
  model_logit = glm(y ~ x, family=binomial(link="logit"))
  ### predict
  fit = predict(model_logit, level=level, type="resp", se.fit=TRUE)
  # fit_link = predict(model_logit, level=level, type="link", se.fit=TRUE)
  # y_pred = 1/(1 + exp(-fit_link$fit))
  ### confindence interval ribbon
  lower = fit$fit - fit$se.fit
  upper = fit$fit + fit$se.fit
  df_ribbon = data.frame(x=c(x, rev(x)), y=c(upper, rev(lower)))
  polygon(x=df_ribbon$x, y=df_ribbon$y, col=col_ci_ribbon, border=NA)
  lines(x=x, y=upper, lwd=lwd_ci_ribbon, col=col_ci_ribbon)
  lines(x=x, y=lower, lwd=lwd_ci_ribbon, col=col_ci_ribbon)
  ### plot logistic regression line
  lines(x=x, y=fit$fit, lwd=lwd_fit, col=col_fit)
  ### Area under the receiver operating curve, FPR vs TPR
  auc = simple_auc(y_binary=y, y_continuous=fit$fit)
  ### legend
  if (show_auc){
    legend("right", legend=paste0("AUC = ", round(auc*100,2), "%"), bty="n")
  }
  if (show_legend){
    legend("left", legend=c("Logistic regression fit", paste0(round(level*100), "% confidence interval"),
          "Distribution of binary response variable\nalong the continuous explanatory variable"),
          col=c(col_fit, col_ci_ribbon, NA), lty=c(1,1,NA), lwd=c(1,5,NA), bord=c(NA,NA,1), fill=c(NA,NA,col_hist))
  }
}