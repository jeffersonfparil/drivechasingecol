### A simple heatmap plotting function using base R
#' @export
simple_heatmap = function(mat, vec_mat_lim=NULL, vec_col_lim=c("white", "red"), resolution=101, main="", xlab="", ylab="", vec_las=c(1,1), vec_numeric_axes=c(FALSE,FALSE), vec_axes_decplac=c(2,2), show.text=FALSE){
  # ### test:
  # mat = mat_mu
  # vec_mat_lim = c(min(mat), max(mat))
  # vec_col_lim = c("white", "red")
  # resolution = 101
  # main="No chasing\nCrash"
  # xlab="Dispersion"
  # ylab="Rmax"
  # show.text=TRUE
  ### generate colour gradient
  func_col_grad = colorRampPalette(vec_col_lim)
  vec_colours = func_col_grad(resolution)
  ### convert mat into an index of vec_colours
  if (is.null(vec_mat_lim)){
    vec_mat_lim = c(min(mat, na.rm=TRUE), max(mat, na.rm=TRUE))
  }
  mat_min = vec_mat_lim[1]
  mat_max = vec_mat_lim[2]
  mat_idx = (round((mat-mat_min)*100/(mat_max-mat_min)))+1
  ### plot
  vec_x = seq(0,ncol(mat), by=1)
  vec_y = seq(0,nrow(mat), by=1)
  plot(x=c(0, max(vec_x)), y=c(0, max(vec_y)), bty="n", type="n", xaxt="n", yaxt="n", xlab=xlab, ylab=ylab, main=main)
  for (i in 1:ncol(mat)){
    for (j in 1:nrow(mat)){
      ### test: i=2; j=2
      x0 = vec_x[i]
      x1 = vec_x[i+1]
      y0 = vec_y[j]
      y1 = vec_y[j+1]
      col = vec_colours[mat_idx[j,i]]
      polygon(x=c(x0,x1,x1,x0), y=c(y0,y0,y1,y1), col=col, border=NA)
      if (show.text) {
        perc = round(mat[j,i],2)
        if (!((vec_mat_lim[1]==0) & (vec_mat_lim[2]==1) & ((perc==0) | (perc==1)))){
          perc = formatC(perc, format='f', flag='0', digits=2)
        }
        text(x=(x0+x1)/2, y=(y0+y1)/2, lab=perc)
      }
    }
  }
  vec_xtick_pos = (head(vec_x,-1) + tail(vec_x,-1))/2
  vec_ytick_pos = (head(vec_y,-1) + tail(vec_y,-1))/2
  if (vec_numeric_axes[1]==TRUE){
    xtick_lab = formatC(as.numeric(colnames(mat)), format='f', flag='0', digits=vec_axes_decplac[1])
  } else {
    xtick_lab = colnames(mat)
  }
  if (vec_numeric_axes[1]==TRUE){
    ytick_lab = formatC(as.numeric(rownames(mat)), format='f', flag='0', digits=vec_axes_decplac[2])
  } else {
    ytick_lab = rownames(mat)
  }
  axis(side=1, lab=xtick_lab, at=vec_xtick_pos, las=vec_las[1])
  axis(side=2, lab=ytick_lab, at=vec_ytick_pos, las=vec_las[2])
  return(list(mat=mat, vec_colours=vec_colours, mat_idx=mat_idx))
}