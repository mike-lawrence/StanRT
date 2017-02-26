get_contrast_matrix = function(data,formula){
	mm = model.matrix(data=data,object=formula)
	attr(mm,'formula') = formula
	return(mm)
}
names_from_mm = function(W,B,reverse = F){
	B_names = dimnames(B)[[2]]
	W_names = dimnames(W)[[2]]
	new_names = rep(NA,length(W_names)*length(B_names))
	k = 1
	if(!reverse){
		for(i in B_names){
			for(j in W_names){
				new_names[k] = paste(j,i,sep=':')
				k = k + 1
			}
		}
	}else{
		for(i in W_names){
			for(j in B_names){
				new_names[k] = paste(i,j,sep=':')
				k = k + 1
			}
		}
	}
	new_names = gsub('(Intercept):','',new_names,fixed=T)
	new_names = gsub(':(Intercept)','',new_names,fixed=T)
	return(new_names)
}

stan_summary = function(object,pars,probs=c(.5,.025,.975),digits=2,W=NULL,B=NULL,return_table=F){

	s = summary(object,pars,probs,digits,use_cache=FALSE)$summary
	s = array(
		s[,4:ncol(s)]
		, dim = c(dim(s)[1],ncol(s)-3)
		, dimnames = list(
			dimnames(s)[[1]]
			, dimnames(s)[[2]][4:ncol(s)]
		)
	)
	if(!is.null(W)){
		dimnames(s)[[1]] = names_from_mm(W,B)
	}
	if(return_table){
		return(s)
	}else{
		print(s,digits=digits)
	}
}

get_condition_post <-
	function(
		post
		, par
		, W
		, B
		, data
		, numeric_res = 0
	){
		if (inherits(data, "tbl_df")) {
			data = as.data.frame(data)
		}
		W_vars = attr(terms(attr(W,'formula')),'term.labels')
		W_vars = W_vars[!stringr::str_detect(W_vars,':')]
		B_vars = attr(terms(attr(B,'formula')),'term.labels')
		B_vars = B_vars[!stringr::str_detect(B_vars,':')]
		data_vars = c(W_vars,B_vars)
		temp = list()
		for(i in 1:length(data_vars)){
			this_var = data_vars[i]
			this_fixed_data = data[,names(data)==this_var]
			if(is.numeric(this_fixed_data)&(numeric_res>0)){
				temp[[i]] = seq(
					min(this_fixed_data)
					, max(this_fixed_data)
					, length.out=numeric_res
				)
			}else{
				temp[[i]] = sort(unique(this_fixed_data))
				if(!is.numeric(this_fixed_data)){
					contrasts(temp[[i]]) = contrasts(this_fixed_data)
				}
			}
		}
		to_return = data.frame(expand.grid(temp))
		names(to_return) = data_vars
		new_names = names_from_mm(W,B,reverse=T)
		both_formula = eval(parse(text=paste('~', paste(data_vars, collapse = '*'))))
		unique_mm = model.matrix(both_formula,to_return)
		unique_mm = unique_mm[,new_names]
		samples = rstan::extract(post,par)[[1]]
		samples = matrix(samples,nrow=dim(samples)[1])
		# samples = samples[,match(new_names,dimnames(unique_mm)[[2]])]
		mat = matrix(NA,nrow=nrow(to_return),ncol=dim(samples)[1])
		for(i in 1:dim(samples)[1]){
			mat[,i] <- unique_mm%*%samples[i,]
		}
		to_return = cbind(to_return,as.data.frame(mat))
		to_return = tidyr::gather(to_return,key='sample',value='value',-one_of(data_vars))
		to_return$sample = as.numeric(gsub('V','',as.character(to_return$sample),fixed=T))
		return(to_return)
	}
