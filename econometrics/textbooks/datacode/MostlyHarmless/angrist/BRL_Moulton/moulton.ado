*! version 2.3 12/16/08  -- by Steve Pischke -- use at your own risk
*! Moulton correction for regression


program moulton, eclass

    version 9.0


	  syntax varlist(min=2) [if] [in], CLuster(string) [noIC moulton force]
	  marksample touse

	  tempname vc b
        tempvar _r _n _one _n1 _xdev
	  tempfile _temp

	  quietly 	{
	
	  preserve
        keep if `touse'
		
	  reg `varlist'

        predict `_r', residuals
        global S_E_nobs = _result(1)
        global S_E_r2 = _result(7)
        global S_E_ar2 = _result(8)
        global S_E_rmse = _result(9)
        global S_E_df = _result(5)

	  matrix `vc' = e(V)
	  matrix `b' = e(b)

	  save _temp, replace
	
	  collapse (count) `_n'=`_r', by(`cluster')
	  sum `_n'
	  local nvar = r(Var)
	  local nmean = r(mean)
	  local nclus = r(N)

	  /* calculate denominator part in case moulton calculation of ICC is requested */

	  gen `_n1' = `_n'*(`_n' - 1)
	  sum `_n1'
	  local n1sum = r(mean)*r(N)	

	  /* use loneway to ICC of residual */

	  use _temp, clear

	  if ("`moulton'" == "") {
		  loneway `_r' `cluster'
		  local rhor = r(rho)
		 }

	   /* use moulton style ICC calculation */

	   else {
 		  local rvar = $S_E_rmse^2
              if ("`moulton'" ~= "") mvar `_r' `cluster'
		  local rhor = min(max($xv/(`rvar'*`n1sum'),0),1)
		  }

	  tokenize `varlist'
	  local depn = "`1'"

	  /* count rhs variables */

	  local k = 2
	  while "``k''" ~= "" {
	 	local ++k
		}
	  local k = `k' - 2			/* k counter is number of vars plus 1 */

        mac shift


	  /* intraclass correlations of rhs variables */

		
	  matrix rho = J(`k',1,0)
	  local i = 1
	  while "``i''" != "" {
	
	   if "`force'" == "" {
			
	
	      if ("`moulton'" == "") {
	
		loneway ``i'' `cluster'
			if (r(sd_w)==0 | r(sd_w)==.) {
				matrix rho[`i',1] = 1
			}
			else {
				matrix rho[`i',1] = r(rho )
			}
	      }

		else {
			 sum ``i''
			 local xvar = r(Var)*(r(N)-1)/$S_E_df
	             gen `_xdev' = ``i'' - r(mean)
                   if ("`moulton'" ~= "") mvar `_xdev' `cluster'		
     		       matrix rho[`i',1] = min(max($xv/(`xvar'*`n1sum'),0),1)
			 drop `_xdev'
			}		

          } 	/* end no force */
          else {
          matrix rho[`i',1] = 1
	    }		/* end force to 1 */
	
            local moul = 1 + ((`nvar'/`nmean' + `nmean' - 1)*`rhor'*rho[`i',1])

	      matrix `vc'[`i',`i'] = `vc'[`i',`i'] * `moul'

		local ++i
	  }			/* end while */

		

	  /* adjust constant */	
		
        local moul = 1 + (`nvar'/`nmean' + `nmean' - 1)*`rhor'
	  matrix `vc'[`i',`i'] = `vc'[`i',`i'] * `moul'

	  }

        /* write header */

        #delimit ;
        di _n in gr
        "OLS Regression: standard errors " _col(55)
        "Number of obs  =" in yel %8.0f $S_E_nobs _n
        in gr "adjusted for cluster effects using Moulton"
        _col(55) in gr "R-squared      ="
        in yel %8.4f $S_E_r2 _n
        _col(55) in gr "Adj R-squared  ="
        in yel %8.4f $S_E_ar2 _n
        in gr "Number of clusters (`cluster') = " in yel `nclus'
        _col(55) in gr "Root MSE       ="
        in yel %8.0g $S_E_rmse _n `addline'  ;
        #delimit cr

        ereturn post `b' `vc', esample(`touse') depname(`depn') dof($S_E_df) obs($S_E_nobs)
	  ereturn display

	  ereturn local clustvar "`cluster'"
	  ereturn local cmd "moulton"

        /* Footer */

	  if "`ic'" == "" {
		  local i = 1
		  while "``i''" != "" {
	      	dis in gr "Intraclass correlation in " %8s abbrev("``i''",8) " = " in yel %7.4f rho[`i',1]
		      local ++i
		  }
        dis in gr "Intraclass correlation in residual = " in yel %7.4f `rhor'
	  }



    end


capture program drop mvar
program mvar
	args var cluster
      tempvar _one
	
      sort `cluster'
      gen `_one' = 1
      matrix opaccum A = `var', group(`cluster') opvar(`_one') noconst
      matrix accum B = `var', noconst
      matrix A = (A - B)					/* subtract out diagonal */
      global xv = A[1,1]

end
