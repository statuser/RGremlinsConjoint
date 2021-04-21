# Internal functions should not be exported

calculate_segment_memberships <- function(betai, bayesm_data) {
  # The procedure is:
  # 1. Calculate the insample hit rates for each individual
  # 2. Draw a random quantile for each of the information poor segments
  # 3. Find the target hit rate based on the random quantile
  # 4. Split the sample into segment 2 (gremlins) for values
  #    below the target hit rate.  All others are segment 1

  num_respondents <- dim(betai)[1]
  num_iterations <- dim(betai)[3]
  num_concepts <- bayesm_data$p
  num_tasks <- nrow(bayesm_data$lgtdata[[1]]$X)/num_concepts

  hit_rates <- double(num_respondents)
  # Convert to vectorized version to avoid doubel loop and speed things up
  for( resp in 1:num_respondents) {
    x_beta <- array((bayesm_data$lgtdata[[resp]]$X %*% betai[resp,,]), c(num_concepts, num_tasks, num_iterations))
    choices <- apply(x_beta, c(2, 3), which.max)
    hit_rates[resp] <- mean(choices == matrix(bayesm_data$lgtdata[[resp]]$y, nrow=16, ncol=num_iterations))
  }

  # This is actually a difficult problem if slightly arbitrary choice when you
  # have more than 2 lambda segments.  This code seems like a reasonable option
  # since we don't really have any information on the appropriate sizes for each
  # segment beyond 2.  Basically we ensure that the Information rich segment is
  # at least 1/num_lambda_segments of the total sample  We also ensure that the
  # poorest information segment is at least 5% of the total sample.  These are
  # just starting values the algorithm will settle on the final segment sizes as
  # it converges
  target_quantiles <- sort(runif(num_lambda_segments - 1, 0.05, 1 - 1/num_lambda_segments))
  hit_rate_cut_offs <- as.numeric(quantile(hit_rates, prob=target_quantiles))
  hit_rate_cut_offs <- c(hit_rate_cut_offs, 1.0)

  segments <- double(num_respondents)
  segment_numbers <- num_lambda_segments:1
  lower <- 0.0
  upper <- hit_rate_cut_offs[1]
  for(i in 1:num_lambda_segments) {
    segments[hit_rates >= lower & hit_rates < upper] <- segment_numbers[i]
    lower <- upper
    upper <- hit_rate_cut_offs[i+1]
  }


  return(segments)
}

# Generate starting values for slope Atchade parameters for tuning runs
# Adjust tuningFactor until accept rates are between 20-80%
# (larger tuningFactor reduce accept rate)
generate_starting_atchade_slopes <- function(slopeBar, tuningFactor, nParams, nUnits) {
  Atch_MCMC_list_slopes <- rep(
    list(
      list(metstd = rep(tuningFactor, times=nParams),
           VMet = rep(1, times=nParams),
           mu_adapt=slopeBar,
           Gamma_adapt=rep(1, times=nParams)
      )
    ),
    times=nUnits)
  return(Atch_MCMC_list_slopes)
}

# Generate starting values for slope Atchade parameters for tuning runs
# Adjust Atch_tau_tune_slopes up or down to target
# Adjust Atch_tau_tune_slopes and Atch_tau_tune_lambda_in until acceptrates are about 20-80%
# (larger Atch_tau_tune_slopes reduce accept rate)
# ncluster is number of Germlins clusters
# est_mean_lambda      a ncluster by 1 vector, first element has to be one; needs to be good guess of posterior mean lambda
# est_var_lambda       a ncluster by 1 vector, first element has to be one; needs to be good guess of posterior variance lambda
generate_starting_atchade_lambdas <- function(Atch_tau_tune_lambda, ncluster, est_mean_lambda, est_var_lambda){

  est_mean_lambda[1] <- 1
  est_var_lambda[1] <- 1

  Atch_MCMC_list_lambda <-
    list( metstd = c(NA, rep(Atch_tau_tune_lambda, times= (ncluster-1))),
          Vmet = c(NA, est_var_lambda[2:ncluster]),
          mu_adapt=est_mean_lambda,
          Gamma_adapt=est_var_lambda
    )

  return(Atch_MCMC_list_lambda)

}

generateIndividualSlope_ATCH <- function(data, design, currentSlope, lambda, slopeBar, slopeCov, constraints,
                                         ind_in, cur_mcmc_iter, metstd_in,Vmet_in,Gamma_adapt_in,mu_adapt_in  ) {
  proposedSlope <- double(length(currentSlope))

  # Break proposals up by constraint type 0,-1,1
  # How many constraints of each type
  ntypeconstraints <- c(sum(constraints==0), sum(constraints==-1), sum(constraints==1))

  accept_rate_est <- double(length(currentSlope))


  #NO CONSTRAINTS
  if (ntypeconstraints[1]>0){
    which_elements <- which((constraints==0))
    VMH_i <- metstd_in[which_elements] * Vmet_in[which_elements]

    proposedSlope_i <- currentSlope
    proposedSlope_i[which_elements] = currentSlope[which_elements] + rnorm(length(which_elements), 0, 1)*sqrt(VMH_i)

    currentLL <- llmnl(currentSlope/lambda, data, design)
    proposedLL <- llmnl(proposedSlope_i/lambda, data, design)

    currentPrior <- -0.5 * t(currentSlope - slopeBar)%*%solve(slopeCov)%*%(currentSlope - slopeBar)
    proposedPrior <- -0.5 * t(proposedSlope_i - slopeBar)%*%solve(slopeCov)%*%(proposedSlope_i - slopeBar)

    logAcceptProb <- proposedLL + proposedPrior - currentLL - currentPrior

    alpha <- log(runif(1))

    if(alpha < logAcceptProb) {
      currentSlope <- proposedSlope_i
      gremlinsEnv$acceptanceRate_slopes[ind_in,1] <- gremlinsEnv$acceptanceRate_slopes[ind_in,1] + 1

      accept_rate_est[which_elements] <- 1

    }else{    # Reject proposal


      accept_rate_est[which_elements] <- 0

    }



  }

  #NEGATIVE CONSTRAINTS
  if (ntypeconstraints[2]>0){

    which_elements <- which((constraints==-1))
    VMH_i <- metstd_in[which_elements] * Vmet_in[which_elements]
    proposedSlope_i <- currentSlope
    proposedSlope_i[which_elements] = currentSlope[which_elements] + rnorm(length(which_elements), 0, 1)*sqrt(VMH_i)


    #Auto reject sign constraints; do whole batch if there are multiple sign constrains (could be better to loop)
    if (all(proposedSlope_i[which_elements]<0)){
      currentLL <- llmnl(currentSlope/lambda, data, design)
      proposedLL <- llmnl(proposedSlope_i/lambda, data, design)

      currentPrior <- -0.5 * t(currentSlope - slopeBar)%*%solve(slopeCov)%*%(currentSlope - slopeBar)
      proposedPrior <- -0.5 * t(proposedSlope_i - slopeBar)%*%solve(slopeCov)%*%(proposedSlope_i - slopeBar)

      logAcceptProb <- proposedLL + proposedPrior - currentLL - currentPrior

      alpha <- log(runif(1))

      if(alpha < logAcceptProb) {
        currentSlope <- proposedSlope_i
        gremlinsEnv$acceptanceRate_slopes[ind_in,2] <- gremlinsEnv$acceptanceRate_slopes[ind_in,2] + 1
        accept_rate_est[which_elements] <- 1
      }else{    # Reject proposal
        accept_rate_est[which_elements] <- 0
      }
    }else{      # Auto reject (proposal not satisfying constraints)
      accept_rate_est[which_elements] <- 0
    }# END IF statement auto reject if contraint is not satisfied
  }


  #POSITIVE CONSTRAINTS
  if (ntypeconstraints[3]>0){
    which_elements <- which((constraints==1))
    VMH_i <- metstd_in[which_elements] * Vmet_in[which_elements]
    proposedSlope_i <- currentSlope
    proposedSlope_i[which_elements] = currentSlope[which_elements] + rnorm(length(which_elements), 0, 1)*sqrt(VMH_i)

    #Auto reject sign constraints; do whole batch if there are multiple sign constrains (could be better to loop)
    if (all(proposedSlope_i[which_elements]>0)){
      currentLL <- llmnl(currentSlope/lambda, data, design)
      proposedLL <- llmnl(proposedSlope_i/lambda, data, design)

      currentPrior <- -0.5 * t(currentSlope - slopeBar)%*%solve(slopeCov)%*%(currentSlope - slopeBar)
      proposedPrior <- -0.5 * t(proposedSlope_i - slopeBar)%*%solve(slopeCov)%*%(proposedSlope_i - slopeBar)

      logAcceptProb <- proposedLL + proposedPrior - currentLL - currentPrior

      alpha <- log(runif(1))

      if(alpha < logAcceptProb) {
        currentSlope <- proposedSlope_i
        gremlinsEnv$acceptanceRate_slopes[ind_in,3] <- gremlinsEnv$acceptanceRate_slopes[ind_in,3] + 1
        accept_rate_est[which_elements] <- 1
      }else{    # Reject proposal
        accept_rate_est[which_elements] <- 0
      }
    }else{      # Auto reject (proposal not satisfying constraints)
      accept_rate_est[which_elements] <- 0
    }# END IF statement auto reject if contraint is not satisfied
  }

  # Update tuning parameters
  if (cur_mcmc_iter>1000){

    # Fixed constants; set in gremlinsEnv$
    g_adapt <- 10/cur_mcmc_iter;
    Vmet_in <- Gamma_adapt_in + gremlinsEnv$eps2
    mu_adapt_in <- mu_adapt_in + g_adapt*(currentSlope - mu_adapt_in)

    dd_x2 <- sqrt(sum(mu_adapt_in*mu_adapt_in))
    if (dd_x2>gremlinsEnv$A1){
      mu_adapt_in <- (gremlinsEnv$A1/dd_x2)*mu_adapt_in
    }

    Gamma_adapt_in <- Gamma_adapt_in +g_adapt*(((currentSlope - mu_adapt_in)^2) - Gamma_adapt_in)
    dd_x1 <- sqrt(sum(Gamma_adapt_in*Gamma_adapt_in))
    if (dd_x1>gremlinsEnv$A1){
      Gamma_adapt_in <- (gremlinsEnv$A1/dd_x1)*Gamma_adapt_in
    }
    metstd_in <- metstd_in +g_adapt*(accept_rate_est - gremlinsEnv$tau_bar)

    tau_too_small <- which(metstd_in<gremlinsEnv$eps1)
    tau_too_big <- which(metstd_in>gremlinsEnv$A1)

    if(length(tau_too_small)>0){
      metstd_in[tau_too_small] <- gremlinsEnv$eps1
    }else if(length(tau_too_big)>0){
      metstd_in[tau_too_big] <- gremlinsEnv$A1
    }
  }    # End update Atchade constants if statement (need enough iterations)
  Atch_out_respi_list <- list(
    metstd = metstd_in,
    Vmet = Vmet_in,
    mu_adapt = mu_adapt_in,
    Gamma_adapt = Gamma_adapt_in
  )

  return(list(
    currentSlope_out = currentSlope,
    Atch_out_respi_list = Atch_out_respi_list
  ))
}


