## ----'preamble', include=FALSE, warning=FALSE, message=FALSE------------------
library(knitr)
library(ggplot2)

theme_set(theme_minimal(base_size=12))
theme_update(legend.position="bottom")

if (require(viridis,quietly=TRUE)) {
  scale_colour_discrete <- function(...) {
    require(viridis,quietly=TRUE)
    scale_color_viridis(discrete=TRUE,end=0.85,option='D',...)
  }
  scale_fill_discrete <- function(...) {
    require(viridis,quietly=TRUE)
    scale_fill_viridis(discrete=TRUE,end=0.85,option='D',...)
  }
  # wait, can I also default for continuous cases?
  scale_colour_continuous <- function(...) {
    require(viridis,quietly=TRUE)
    scale_color_viridis(discrete=FALSE,end=0.85,option='D',...)
  }
  scale_fill_continuous <- function(...) {
    require(viridis,quietly=TRUE)
    scale_fill_viridis(discrete=FALSE,end=0.85,option='D',...)
  }
}

# set the knitr options ... for everyone!
# if you unset this, then vignette build bonks. oh, joy.
#opts_knit$set(progress=TRUE)
opts_knit$set(eval.after='fig.cap')
# for a package vignette, you do want to echo.
# opts_chunk$set(echo=FALSE,warning=FALSE,message=FALSE)
opts_chunk$set(echo=FALSE,warning=FALSE,message=FALSE)
#opts_chunk$set(results="asis")
opts_chunk$set(cache=TRUE,cache.path="cache/nmf_")

#opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
opts_chunk$set(fig.path="figure/nmf_",dev=c("pdf"))
opts_chunk$set(fig.width=5,fig.height=4,dpi=64)
# stolen from my book:
resol <- 7
aspr <- 1.41
opts_chunk$set(fig.width=resol,
               fig.height=resol / aspr,
               out.width="0.975\\textwidth",
               out.height=sprintf('%.3f\\textwidth',0.975/aspr),
               dpi=450)

# a more compact figure use
# fig.width=cfg_width,fig.height=cfg_height,out.width=cot_width,out.height=cot_height,
cfg_width <- resol
cfg_height <- resol / 2
cot_width <- "0.975\\textwidth"
cot_height <- "0.4875\\textwidth"


# for cup plots.
# fig.width=cpg_width,fig.height=cpg_height,out.width=cpt_width,out.height=cpt_height,
cpg_width <- 1.5 * resol
cpg_height <- 1.5* resol / 1.75
cpt_width <- "0.975\\textwidth"
cpt_height <- paste0(signif(0.975 * cpg_height / cpg_width,4),"\\textwidth")

# for tall cup plots.
# fig.width=tcg_width,fig.height=tcg_height,out.width=tct_width,out.height=tct_height,
tcg_width <- 1.5 * resol
tcg_height <- 1.5* resol * 1.4
tct_width <- "0.975\\textwidth"
tct_height <- paste0(signif(0.975 * tcg_height / tcg_width,4),"\\textwidth")


# a wide figure use
# fig.width=wfg_width,fig.height=wfg_height,out.width=wot_width,out.height=wot_height,
wfg_width <- 1.2 * resol
wfg_height <- 0.6 * resol
wot_width <- "0.99\\textwidth"
wot_height <- "0.495\\textwidth"

# a tall figure use
# fig.width=tfg_width,fig.height=wfg_height,out.width=tot_width,out.height=tot_height,
tfg_width <- resol
tfg_height <- aspr * resol
tot_width <- "0.99\\textwidth"
tot_height <- paste0(round(0.99 * aspr,3),"\\textwidth")

opts_chunk$set(fig.pos='h')
# doing this means that png files are made of figures;
# the savings is small, and it looks like shit:
#opts_chunk$set(fig.path="figure/",dev=c("png","pdf","cairo_ps"))
#opts_chunk$set(fig.width=4,fig.height=4)
# for figures? this is sweave-specific?
#opts_knit$set(eps=TRUE)

# this would be for figures:
#opts_chunk$set(out.width='.8\\textwidth')
# for text wrapping:
# CRAN complains if you do this:
# options(width=64,digits=2)
opts_chunk$set(size="small")
opts_chunk$set(tidy=TRUE,tidy.opts=list(width.cutoff=50,keep.blank.line=TRUE))

## ----'mc_sims',eval=TRUE------------------------------------------------------
# ok.
library(dplyr)
library(rnnmf)

frobenius_err <- function(Y, L, R) {
	sqrt(sum(abs(Y - L %*% R)^2))
}

runifmat <- function(nr,nc,...) {
	matrix(pmax(0,runif(nr*nc,...)),nrow=nr)
}

test_a_bunch <- function(Y_t, L_0, R_0, niter=1e4L) {
	iter_hist <- new.env()
	iter_hist[['history']] <- rep(NA_real_, niter)
	on_iteration_end <- function(iteration, Y, L, R, ...) {
		iter_hist[['history']][iteration] <<- frobenius_err(Y,L,R)
	}
	wuz <- aurnmf(Y_t, L_0, R_0, max_iterations=niter, on_iteration_end=on_iteration_end)
	df1 <- tibble(x=seq_along(iter_hist[['history']]),y=iter_hist[['history']]) %>% mutate(method='additive, optimal step')

	#iter_hist[['history']] <- rep(NA_real_, niter)
	#wuz <- aurnmf(Y_t, L_0, R_0, max_iterations=niter, check_optimal_step=FALSE, on_iteration_end=on_iteration_end)
	#df15 <- tibble(x=seq_along(iter_hist[['history']]),y=iter_hist[['history']]) %>% mutate(method='additive, naive step')

	#iter_hist[['history']] <- rep(NA_real_, niter)
	#wuz <- aurnmf(Y_t, L_0, R_0, max_iterations=niter, check_optimal_step=FALSE, tau=0.9, on_iteration_end=on_iteration_end)
	#df17 <- tibble(x=seq_along(iter_hist[['history']]),y=iter_hist[['history']]) %>% mutate(method='additive, naive step, tau=0.9')

	iter_hist[['history']] <- rep(NA_real_, niter)
	wuz <- murnmf(Y_t, L_0, R_0, max_iterations=niter, on_iteration_end=on_iteration_end)
	df2 <- tibble(x=seq_along(iter_hist[['history']]),y=iter_hist[['history']]) %>% mutate(method='multiplicative')

	#retv <- bind_rows(df1,df15,df17,df2) %>%
	retv <- bind_rows(df1,df2) %>%
		mutate(nr=nrow(Y_t),
					 nc=ncol(Y_t),
					 nd=ncol(L_0),
					 max_iter=niter)
	return(retv)
}

nr <- 30
nc <- 8
ynd <- 2
set.seed(1234)
L_t <- runifmat(nr,ynd)
R_t <- runifmat(ynd,nc)
Y_t <- L_t %*% R_t

L_0 <- runifmat(nrow(Y_t),ynd+1)
R_0 <- runifmat(ncol(L_0),ncol(Y_t))

test1 <- test_a_bunch(Y_t, L_0, R_0, niter=1e4L) %>%
	mutate(true_nd=ynd)

nr <- 40
nc <- 10 
ynd <- 3

set.seed(4579)
L_t <- runifmat(nr,ynd,min=-1,max=1)
R_t <- runifmat(ynd,nc,min=-1,max=1)
Y_t <- L_t %*% R_t
L_0 <- runifmat(nrow(Y_t),ynd+1,min=-0.5,max=1)
R_0 <- runifmat(ncol(L_0),ncol(Y_t),min=-0.5,max=1)

test2 <- test_a_bunch(Y_t, L_0, R_0, niter=1e4L) %>%
	mutate(true_nd=ynd)

set.seed(6789)
L_0 <- runifmat(nrow(Y_t),ynd+1,min=1e-4,max=1)
R_0 <- runifmat(ncol(L_0),ncol(Y_t),min=1e-4,max=1)

test2b <- test_a_bunch(Y_t, L_0, R_0, niter=1e4L) %>%
	mutate(true_nd=ynd)

## ----'mc_sims_plot1',dependson=c('mc_sims'),eval=TRUE,fig.cap=paste0("The Frobenius norm is plotted versus step for two methods for a small problem."),eval.after='fig.cap'----
# ell2 <- '\u2113\u2082'
ell2 <- 'Frobenius'

test1 %>%
	ggplot(aes(x,y,color=method)) + 
	geom_line() + 
	scale_x_log10(labels=scales::comma) + scale_y_log10() +
	labs(x='Step',y='Frobenius Norm of Error',
			 title='Frobenius Norm of Error vs Step',
			 color='Method',
			 caption=paste0('Factoring ',test1$nr[1],' x ',test1$nc[1],' matrix down to ',test1$nd[1],' dimensions. Y matrix has rank ',test1$true_nd[1],'.'))

## ----'mc_sims_plot2',dependson=c('mc_sims'),eval=TRUE,fig.cap=paste0("The Frobenius norm is plotted versus step for two methods for a small problem. Starting iterates are taken to be sparse or dense."),eval.after='fig.cap'----

bind_rows(test2 %>% mutate(starting_iterate='sparse'),
	test2b %>% mutate(starting_iterate='dense')) %>%
	ggplot(aes(x,y,color=method)) + 
	geom_line() + 
	scale_x_log10(labels=scales::comma) + scale_y_log10() +
	facet_grid(.~starting_iterate,labeller=label_both) +
	labs(x='Step',y='Frobenius Norm of Error',
			 title='Frobenius Norm of Error vs Step',
			 color='Method',
			 caption=paste0('Factoring ',test2$nr[1],' x ',test2$nc[1],' matrix down to ',test2$nd[1],' dimensions. Y matrix has rank ',test2$true_nd[1],'.'))

