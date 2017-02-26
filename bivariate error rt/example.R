#Load packages ----
library(tidyverse) #for data wrangling & vis
library(rstan) #for Stan stuff
library(ez) #for the ANT data set

#load useful functions
source('helpers.R')

#define a useful contrasts function
best_contrasts_ever = function(n, contrasts = TRUE, sparse = FALSE){
	contr.sum(n,contrasts,sparse)*.5
}

#set options to use the above contrasts by default
options(
	contrasts = c('best_contrasts_ever','contr.poly')
)

# load in the ANT data ----
data(ANT)
head(ANT)

# Prep the data for Stan ----

#generate within-subjects matrix
W = get_contrast_matrix(
	data = ANT
	, formula = ~ cue*flank
)
head(W)

#for the between-subjects contrast matrix, first reduce data to just the subject
# and between-subject predictors
ANT %>%
	dplyr::group_by(subnum,group) %>%
	dplyr::summarize(n=length(rt)) -> b

#generate between-subjects contrast matrix
B = get_contrast_matrix(
	data = b
	, formula = ~ group
)
head(B)

#package in list for Stan
data_for_stan = list(
	#nTrials: num trials
	nTrials = nrow(ANT)
	#RT: RT outcomes
	, RT = scale(ANT$rt)[,1] #scaling for easy priors
	#ER: Error outcomes
	, ER = as.numeric(ANT$error)
	#accInd: Indices of accurate outcomes
	, accInd = (1:nrow(ANT))[!ANT$error]
	#nS: num subjects
	, nS = length(unique(ANT$subnum))
	#S: trial-by-trial subject labels
	, S = as.numeric(factor(ANT$subnum))
	#nWerr: num within predictors on error
	, nWerr = ncol(W)
	#nBerr: num group predictors on error
	, nBerr = ncol(B)
	#nWmrt: num within predictors on mean RT
	, nWmrt = ncol(W)
	#nBmrt: num group predictors on mean RT
	, nBmrt = ncol(B)
	#nWsrt: num within predictors on sd RT
	, nWsrt = ncol(W)
	#nBsrt: num group predictors on sd RT
	, nBsrt = ncol(B)
	#Werr: within predictors for error
	, Werr = W
	#Berr: between predictors for error
	, Berr = B
	#Wmrt: within predictors for mean RT
	, Wmrt = W
	#Bmrt: between predictors for mean RT
	, Bmrt = B
	#Wsrt: within predictors for sd RT
	, Wsrt = W #for just intercepts do: matrix(W[,1],ncol=1)
	#Bsrt: between predictors for sd RT
	, Bsrt = B #for just intercepts do: matrix(B[,1],ncol=1)
)


# Compile & sample the model ----
post = rstan::stan(
	file = 'er_rt.stan'
	, data = data_for_stan
	, seed = 1
	, chains = 4
	, cores = 4 #set this to # of physical cores on your system
	, iter = 2e3
	, init = 0
	, pars = c('dummy','ZcorsW') #don't bother saving these variables
	, include = FALSE
)

# Check the posterior ----

#coefficients relating effects on error rate
stan_summary(
	object = post
	, par = 'coefErr'
	, W = W
	, B = B
)

#coefficients relating effects on mean RT
stan_summary(
	object = post
	, par = 'coefMrt'
	, W = W
	, B = B
)

#coefficients relating effects on SD of RT
stan_summary(
	object = post
	, par = 'coefSrt'
	, W = W
	, B = B
)

#check for correlations amongst the predictors' effects
cors = stan_summary(post,'corsW',return_table=T)
any(((cors[,2]>0)|(cors[,3]<0))&(cors[,1]!=1))
cors[which(((cors[,2]>0)|(cors[,3]<0))&(cors[,1]!=1)),]

#get the names of predictors in cors
# for example, if cor[11,13]:
rep(dimnames(W)[[2]],3)[11]
rep(dimnames(W)[[2]],3)[13]



# Compute & visualize the condition mean RTs

conditions = get_condition_post(
	post = post
	, par = 'coefMrt'
	, W = W
	, B = B
	, data = ANT
)

#convert from z-score to original RT scale (because we passed "scale(ANT$rt)[,1]" as the RT data to the model above)
conditions$value = conditions$value*sd(ANT$rt)+mean(ANT$rt)

#visualize the full flank*cue*group interaction
conditions %>%
	dplyr::group_by(flank,cue,group) %>%
	dplyr::summarize(
		med = quantile(value,.5)
		, lo95 = quantile(value,.25)
		, hi95 = quantile(value,.75)
	) %>%
	ggplot(
		mapping = aes(
			x = flank
			, y = med
			, ymin = lo95
			, ymax = hi95
			, color = cue
			, group = cue
		)
	)+
	facet_grid(
		. ~ group
	)+
	geom_point()+
	geom_line()+
	geom_errorbar(width=0)+
	labs(
		y = 'Mean RT (ms)'
	)

#visualize the main effect of flank interaction
conditions %>%
	dplyr::group_by(flank,sample) %>%
	dplyr::summarize(
		value = mean(value)
	) %>%
	dplyr::group_by(flank) %>%
	dplyr::summarize(
		med = quantile(value,.5)
		, lo95 = quantile(value,.25)
		, hi95 = quantile(value,.75)
	) %>%
	ggplot(
		mapping = aes(
			x = flank
			, y = med
			, ymin = lo95
			, ymax = hi95
			, group = 1
		)
	)+
	geom_point()+
	geom_line()+
	geom_errorbar(width=0)+
	labs(
		y = 'Mean RT (ms)'
	)

#compute the 95%CrI on the I-C flanker effect
conditions %>%
	dplyr::group_by(flank,sample) %>%
	dplyr::summarize(
		value = mean(value)
	) %>%
	dplyr::group_by(sample) %>%
	dplyr::summarize(
		value = value[flank=='Incongruent'] - value[flank=='Congruent']
	) %>%
	dplyr::summarize(
		med = quantile(value,.5)
		, lo95 = quantile(value,.25)
		, hi95 = quantile(value,.75)
	)


#Etc.
