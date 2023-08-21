postPlot_v2 = function (post,
                        prior = NULL,
                        postpred,
                        priorpred = NULL,
                        dataMats,
                        classifs,
                        var.names = c("mu_m",
                                      "tau_m",
                                      "tau_norm",
                                      "mu_c",
                                      "tau_c",
                                      "probdiff",
                                      "m",
                                      "c"),
                        hierParams = c("m", "c"),
                        xlabs = NULL,
                        mitoPlot_xlab = "",
                        mitoPlot_ylab = "",
                        ...)
{
  if (!is.null(xlabs)) {
    if (is.null(names(xlabs)))
      names(xlabs) = var.names
  }
  else {
    xlabs = var.names
    names(xlabs) = var.names
  }
  if (is.character(post)) {
    if (file.exists(post)) {
      post = data.table::fread(post, header = TRUE)
      post = as.data.frame(post)
    }
    else {
      stop(paste0("The file, ", post, ", does not exist."))
    }
  }
  if (is.character(postpred)) {
    if (file.exists(postpred)) {
      postpred = data.table::fread(postpred, header = TRUE)
      postpred = as.data.frame(postpred)
    }
    else {
      stop(paste0("The file, ", postpred, ", does not exist."))
    }
  }
  if (!is.null(prior)) {
    if (is.character(prior)) {
      if (file.exists(prior)) {
        prior = data.table::fread(prior, header = TRUE)
        prior = as.data.frame(prior)
      }
      else {
        message(paste0("The file, ", prior, ", does not exist."))
      }
    }
  }
  if (!is.null(priorpred)) {
    if (is.character(priorpred)) {
      if (file.exists(priorpred)) {
        priorpred = data.table::fread(priorpred, header = TRUE)
        priorpred = as.data.frame(priorpred)
      }
      else {
        message(paste0("The file, ", priorpred, ", does not exist."))
      }
    }
  }
  
  nPts = length(unique(dataMats$indexPat))
  nCrl = length(unique(dataMats$indexCtrl))
  
  for( pat_ind in unique(dataMats$indexPat) ) {
    colnames = colnames(post)
    for (var in var.names) {
      if (var %in% hierParams) {
        post_var = post[, grepl(paste0(var, "\\["), colnames)]
        post_dens = list()
        priorParent = prior[, paste0(var, "_pred")]
        priorParent_dens = stats::density(priorParent)
        postParent = post[, paste0(var, "_pred")]
        postParent_dens = stats::density(postParent)
        xlim_varMin = min(postParent_dens$x)
        xlim_varMax = max(postParent_dens$x)
        ylim_var = max(postParent_dens$y)
        for (param in colnames(post_var)) {
          post_dens[[paste(param)]] = stats::density(post_var[[param]])
          ylim_var = max(c(ylim_var, post_dens[[param]]$y))
          xlim_varMax = max(c(xlim_varMax, post_dens[[param]]$x))
          xlim_varMin = min(c(xlim_varMin, post_dens[[param]]$x))
        }
        plot(
          NA,
          main = "",
          xlab = xlabs[var],
          ylab = "",
          xlim = c(xlim_varMin, xlim_varMax),
          ylim = c(0, ylim_var),
          ...
        )
        graphics::lines(priorParent_dens,
                        col = alphaPink(1),
                        lty = 4,
                        ...)
        graphics::lines(postParent_dens, col = alphaGreen(1), lty = 4, ...)
        ind = grep(paste0(var, "\\["), colnames, value = TRUE)
        ind = gsub("\\]", "", gsub(paste0(var, "\\["), "", ind))
        max_ind = max(as.numeric(ind))
        for (param in colnames(post_var)[paste0(var, "[", pat_ind, "]") != colnames(post_var)]) {
          graphics::lines(post_dens[[param]], col = alphaGreen(0.1), ...)
        }
        graphics::lines(post_dens[[paste0(var, "[", pat_ind, "]")]], col = alphaGreen(1), ...)
      }
      else {
        if( var=="probdiff") var = paste0("probdiff[", pat_ind-nCrl, "]")
        post_den = stats::density(post[, paste(var)])
        prior_den = stats::density(prior[, paste(var)])
        xlims = range(c(prior_den$x, post_den$x))
        yMax = max(c(prior_den$y, post_den$y))
        plot(
          NA,
          xlim = xlims,
          ylim = c(0, yMax),
          xlab = xlabs[var],
          ylab = "",
          main = "",
          ...
        )
        graphics::lines(prior_den, col = alphaPink(1), ...)
        graphics::lines(post_den, col = alphaGreen(1), ...)
      }
    }
    xlims = range((c(dataMats$ctrl[, 1], dataMats$pts[, 1])))
    ylims = range((c(dataMats$ctrl[, 2], dataMats$pts[, 2])))
    plot(
      NULL,
      xlab = mitoPlot_xlab,
      ylab = mitoPlot_ylab,
      main = "",
      xlim = xlims,
      ylim = ylims,
      ...
    )
    graphics::points(dataMats$ctrl, pch = 20, col = alphaBlack(0.05), ...)
    graphics::points(dataMats$pts[dataMats$indexPat==pat_ind, ], 
                     pch = 20, 
                     col = classcols( apply(classifs[[pat_ind-nCrl]], 2, mean)),
                     ...)
    graphics::lines(postpred[[pat_ind-nCrl]][, "mitochan"], postpred[[pat_ind-nCrl]][, "lwr_norm"],
                    lty = 2, col = alphaGreen(1), ...)
    graphics::lines(postpred[[pat_ind-nCrl]][, "mitochan"], postpred[[pat_ind-nCrl]][, "med_norm"],
                    lty = 1, col = alphaGreen(1), ...)
    graphics::lines(postpred[[pat_ind-nCrl]][, "mitochan"], postpred[[pat_ind-nCrl]][, "upr_norm"],
                    lty = 2, col = alphaGreen(1), ...)
  }
}
