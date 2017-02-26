data {
	#nTrials: num trials
	int<lower=1> nTrials ;
	#RT: RT outcomes
	vector[nTrials] RT ;
	#ER: Error outcomes
	int<lower=0,upper=1> ER[nTrials] ;
	#accInd: Indices of accurate outcomes
	int<lower=1,upper=nTrials> accInd[nTrials-sum(ER)] ;
	#nS: num subjects
	int<lower=1> nS ;
	#S: trial-by-trial subject labels
	int<lower=1,upper=nS> S[nTrials] ;
	#nWerr: num within predictors on error
	int<lower=1> nWerr ;
	#nBerr: num group predictors on error
	int<lower=1> nBerr ;
	#nWmrt: num within predictors on mean RT
	int<lower=1> nWmrt ;
	#nBmrt: num group predictors on mean RT
	int<lower=1> nBmrt ;
	#nWsrt: num within predictors on sd RT
	int<lower=1> nWsrt ;
	#nBsrt: num group predictors on sd RT
	int<lower=1> nBsrt ;
	#Werr: within predictors for error
	matrix[nTrials, nWerr] Werr ;
	#Berr: between predictors for error
	matrix[nS,nBerr] Berr ;
	#Wmrt: within predictors for mean RT
	matrix[nTrials, nWmrt] Wmrt ;
	#Bmrt: between predictors for mean RT
	matrix[nS,nBmrt] Bmrt ;
	#Wsrt: within predictors for sd RT
	matrix[nTrials, nWsrt] Wsrt ;
	#Bsrt: between predictors for sd RT
	matrix[nS,nBsrt] Bsrt ;
}
transformed data{
	#nWtot: total number of within predictors
	int nWtot ;
	#erIntObs: observed intercept for the error data
	real errIntObs ;
	#compute nWtot
	nWtot = nWerr + nWmrt + nWsrt ;
	#compute observed intercept for the error data
	errIntObs = logit(mean(to_vector(ER))) ;
}
parameters {
	#dummy: a helper variable
	matrix[nWtot, nS] dummy;
	#ZsdW: population-level sds (on the z-scale) for each within-subject predictor
	vector<lower=0>[nWtot] sdsW ;
	#ZcorW: population-level correlations (on cholesky factor scale) amongst within-subject predictors
	cholesky_factor_corr[nWtot] ZcorsW ;
	#coefErr: coefficients for between and within subject predictors on error
	matrix[nBerr, nWerr] coefErr ;
	#coefMrt: coefficientsfor between and within subject predictors on mrt
	matrix[nBmrt, nWmrt] coefMrt ;
	#coefSrt: coefficients for between and within subject predictors on srt
	matrix[nBsrt, nWsrt] coefSrt ;
	#u: a dummy variable to implement multivariate cauchy on Sdevs
	real<lower=0> u ;
}
model {
	#u must be chi_square(1) for multivariate cauchy
	u ~ chi_square(1) ;
	#dummy must have normal(0,1) prior for multivariate trick
	to_vector(dummy) ~ normal(0, 1) ;
	#flat prior on correlations
	ZcorsW ~ lkj_corr_cholesky(1) ;
	#weibull prior on subject deviations
	sdsW ~ weibull(2, 1) ;
	#normal(0,4) priors on all error coefs
	to_vector(coefErr) ~ normal(0, 4) ;
	#except the error intercept, which is centered on the observed value
	coefErr[1,1] ~ normal(errIntObs,1) ;
	#normal(0,1) priors on all mrt intercept & coefs
	to_vector(coefMrt) ~ normal(0, 1) ;
	#normal(0,1) priors on all srt intercept & coefs
	to_vector(coefSrt) ~ normal(0, 1) ;
	{
		#Sdevs: subject-by-subject deviations
		matrix[nS, nWtot] Sdevs ;
		#SvalsErr: subject-by-subject err values
		matrix[nS, nWerr] SvalsErr ;
		#SvalsMrt: subject-by-subject mrt values
		matrix[nS, nWmrt] SvalsMrt ;
		#SvalsSrt: subject-by-subject srt values
		matrix[nS, nWsrt] SvalsSrt ;
		#err: trial-by-trial err values
		vector[nTrials] err ;
		#mrt: trial-by-trial mrt values
		vector[nTrials] mrt ;
		#srt: trial-by-trial srt values
		vector[nTrials] srt ;
		#compute
		Sdevs = sqrt(1 / u) * transpose(diag_pre_multiply(sdsW,ZcorsW) * dummy) ;
		SvalsErr = Berr * coefErr + Sdevs[,1:nWerr] ;
		SvalsMrt = Bmrt * coefMrt + Sdevs[,(nWerr+1):(nWerr+nWmrt)] ;
		SvalsSrt = Bsrt * coefSrt + Sdevs[,(nWerr+nWmrt+1):nWtot] ;
		err = rows_dot_product( SvalsErr[S] , Werr ) ;
		mrt = rows_dot_product( SvalsMrt[S] , Wmrt ) ;
		srt = rows_dot_product( SvalsSrt[S] , Wsrt ) ;
		target += cauchy_lpdf( RT[accInd] | mrt[accInd] , exp(srt[accInd]) );
		target += bernoulli_logit_lpmf( ER | err );
	}
}
generated quantities {
	corr_matrix[nWtot] corsW;
	corsW = multiply_lower_tri_self_transpose(ZcorsW);
}
