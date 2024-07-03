

*******************************************
**** Chris Walters
**** This program applies empirical Bayes methods
**** to study employment discrimination using data
**** from Kline, Rose, and Walters (QJE 2022)
********************************************


	********************
	*** Basic setup ****
	********************
	
		*Set up stata
		clear
		cap set more off
		cap set trace off
		cap log close
		set seed 1028
		cd "/Users/christopherwalters/Dropbox/teaching/short_courses/ucla/labs" 
		
	**************************************
	*** Lab 1
	***************************************
		
		****Load data
		insheet using "krw_data.csv", names comma clear
		
		****Summarize variables in this data set
		sum	
		
		*****Overall treatment effects -- compare robust and job-clustered SE
		reg callback white, r
		reg callback white, cluster(job)
	
		****Collapse to data set of firm-specific gaps and SEs
		gen callback_white=callback if white==1
		gen callback_black=callback if white==0
		bys job_id: egen ybar_white=mean(callback_white)
		bys job_id: egen ybar_black=mean(callback_black)
		bys job_id: gen n_j=_N
		bys job_id: keep if _n==1
		gen thetahat_jf=ybar_white - ybar_black
		drop if thetahat_jf==.
		bys firm_id: egen thetahat_f = mean(thetahat_jf)
		bys firm_id: egen std=sd(thetahat_jf)
		bys firm_id: gen J=_N
		gen s_f=sqrt((1/J)*(std^2))
		bys firm_id: keep if _n==1
		keep thetahat_f s_f

		****Estimate mean and variance of mixing distribution
		qui sum thetahat_f
		local mu=r(mean)
		local F=r(N)
		gen v=(thetahat_f - (`mu'))^2 - (s_f^2)
		qui sum v
		local sigma=sqrt(r(mean))
		disp "Estimated SD of firm gaps: `sigma'"

		****Linear shrinkage estimates
		gen lambda_f=((`sigma')^2)/(((`sigma')^2)+(s_f^2))
		gen thetastar_f=lambda_f*thetahat_f + (1-lambda_f)*(`mu')
		sum thetahat_f thetastar_f

		****Histogram of estimates vs. linear shrinkage posteriors
		local N=_N
		local width=0.005
		local bwidth=`width'/1.5
		local min = -4*`sigma'
		local max = 4*`sigma'
		graph twoway hist thetahat_f, freq width(`width') color(navy) fintensity(0) ///
					|| hist thetastar_f, width(`width') barwidth(`bwidth') color(maroon) ///
					|| function y = (`N'*`width')*normalden(x,`mu',`sigma'), range(`min' `max') lcolor(black) ///
					legend(lab(1 "Unbiased estimates") lab(2 "Linear shrinkage estimates") lab(3 "Mixing distribution") order(1 3 2)) ///
					xtitle("White/Black contact gap") ytitle("Number of firms") ///
					scheme(s1color)
					
	**************************************
	*** Lab 2
	***************************************	
						
		****Execute R command for nonparametric deconvolution
		gen id_f=_n
		preserve
		keep id_f thetahat_f s_f
		order id_f thetahat_f s_f
		outsheet using "disc_estimates.csv", nonames comma replace
		rsource using "decon_example.R", rpath("/usr/local/bin/R") roptions(`"--vanilla"')
		restore
		
		****Histogram of estimates with nonparametric mixing distribution estimate
		preserve
		gen sample1=1
		tempfile sample
		save "`sample'"
		insheet using "decon_estimates.csv", nonames comma clear
		ren v1 supp_theta
		ren v2 g_theta
		gen sample1=0
		append using "`sample'"
		graph twoway hist thetahat_f if sample1==1,  color(navy) fintensity(0) width(`width') density ///
					|| line g_theta supp_theta if sample1==0,  lcolor(maroon) ///
					legend(lab(1 "Unbiased estimates") lab(2 "Nonparametric prior") order(1 2)) ///
					xtitle("White/Black contact gap") ytitle("Density") ///
					scheme(s1color)	
		restore

			
		*****Plot of linear shrinkage vs. non-parametric posterior means
		
			*Bring in non-parametric posteriors from R
			tempfile mergefile
			save "`mergefile'", replace
			insheet using "postmean_estimates.csv", nonames comma clear
			ren v1 id_f
			ren v2 thetastar_nonpar_f
			merge 1:1 id_f using "`mergefile'"
			tab _merge
			drop _merge
				
			*Plot linear shrinkage vs. nonparametric posteriors
			graph twoway scatter thetastar_f thetastar_nonpar_f, color(navy) ///
				ytitle("Linear shrinkage estimate") xtitle("Nonparametric posterior mean") ///
				scheme(s1color)
				
			preserve
			gen sample1=1
			tempfile sample
			save "`sample'"
			insheet using "decon_estimates.csv", nonames comma clear
			ren v1 supp_theta
			ren v2 g_theta
			gen sample1=0
			append using "`sample'"
			graph twoway hist thetahat if sample1==1,  color(navy) fintensity(0) width(0.005) density ///
						|| line g_theta supp_theta if sample1==0,  lcolor(maroon) ///
						legend(lab(1 "Unbiased estimates") lab(2 "Nonparametric prior") order(1 2)) ///
						xtitle("White/Black contact gap") ytitle("Density") ///
						scheme(s1color)	
			restore

	
	**************************************
	*** Multiple testing: Classify firm discrimination while controlling FDR
	***************************************
		if `multiple_test'==1 {
		
			***Load estimates
			use "`basic_estimates'", clear
		
			*Choose cutoff lambda for calculating bound on share of true nulls (pi0)
			local lambda=0.5
			
			*Choose False Discovery Rate control threshold
			local FDR=0.05
		
			*Compute one-tailed tests of H0: theta=0 vs. HA: theta>0
			gen z_f=thetahat_f/s_f
			gen p_f=1-normal(z_f)
						
			*Bound pi_0
			gen c=(p_f*(p_f>`lambda'))/(1-`lambda')
			qui sum c
			local pi_0=r(mean)
			disp "Bound on pi_0: `pi_0'"
				
			*Compute q-values
			sort p_f
			gen F_p=_n/_N
			gen q_f=(p_f*(`pi_0'))/F_p
				
			*Count firms with q-vals below FDR threshold, and find p-val cutoff
			count if q_f<`FDR'
			local N_firms=r(N)
			sum p_f if q_f<`FDR', detail
			local p_cutoff=r(max)
			
			*Plot histogram of p-vals
			local pi0=round(`pi_0',0.01)
			local lambda_loc=`lambda'+0.06
			local pi0_loc=`pi0'+0.2
			local p_loc=`p_cutoff'+0.175
			graph twoway hist p_f, color(navy) fintensity(0) width(0.05) /// 
					legend(off) ///
					xtitle("P-value from test of no disc. against Black applicants") ytitle("Density") ///
					scheme(s1color) 
					
				
					
			*Plot histogram of p-vals with lambda, pi0, and p-val cutoff for FDR control
			graph twoway hist p_f, color(navy) fintensity(0) width(0.05) /// 
					|| function y=`pi_0', range(0 1) lcolor(maroon) lpattern(dash) ///
					xline(`lambda', lcolor(black) lpattern(dash)) ///
					legend(off) ///
					xtitle("P-value from test of no disc. against Black applicants") ytitle("Density") ///
					scheme(s1color) ///
					text(3.5 `lambda_loc' "{&lambda} = 0`lambda'") ///
					text(`pi0_loc' 0.725 "{&pi}{subscript:0} = 0`pi0'")	
				
			*Plot histogram of p-vals with lambda, pi0, and p-val cutoff for FDR control
			local pi0=round(`pi_0',0.01)
			local lambda_loc=`lambda'+0.06
			local pi0_loc=`pi0'+0.2
			local p_loc=`p_cutoff'+0.175
			graph twoway hist p_f, color(navy) fintensity(0) width(0.05) /// 
					|| function y=`pi_0', range(0 1) lcolor(maroon) lpattern(dash) ///
					xline(`lambda', lcolor(black) lpattern(dash)) ///
					xline(`p_cutoff', lcolor(black) lpattern(dash)) ///
					legend(off) ///
					xtitle("P-value from test of no disc. against Black applicants") ytitle("Density") ///
					scheme(s1color) ///
					text(3.5 `lambda_loc' "{&lambda} = 0`lambda'") ///
					text(`pi0_loc' 0.725 "{&pi}{subscript:0} = 0`pi0'") ///
					text(5.5 `p_loc' "`N_firms' firms with q-vals < 0`FDR'")

		}
	
	
