# Internal functions should not be exported

#' @importFrom bayesm llmnl
#' @importFrom stats optim
find_starting_values <- function(respondent_level_data, respondent_level_designs) {
  slope_start <- matrix(0, nrow=gremlinsEnv$num_respondents, ncol=gremlinsEnv$num_parameters)
  unconverged_respondents <- c()

  for(ind in 1:gremlinsEnv$num_respondents) {
    mle <- optim(matrix(0, nrow=1, ncol=gremlinsEnv$num_parameters), llmnl, y=respondent_level_data[[ind]], X=respondent_level_designs[[ind]], method = "BFGS", control=list(fnscale=-1))
    if(mle$convergence == 0) {
      slope_start[ind, ] <- mle$par
    } else {
      unconverged_respondents <- append(unconverged_respondents, ind)
    }
  }

  if(length(unconverged_respondents) > 0) {
    ## Calculate the aggregate logit and use the results for starting values for any individuals whose MLE didn't
    overall_design <- do.call(rbind, respondent_level_designs)
    overall_data <- do.call(c, respondent_level_data)
    overall_mle <- optim(matrix(0, nrow=1, ncol=gremlinsEnv$num_parameters), llmnl, y=overall_data, X=overall_design, method = "BFGS", control=list(fnscale=-1))
    slope_start[unconverged_respondents,] <- matrix(rep(overall_mle$par, times=length(unconverged_respondents)), ncol = gremlinsEnv$num_parameters, byrow=TRUE)
  }

  slopeBar_start <- colMeans(slope_start)
  slopeCov_start <- diag(length(slopeBar_start))
  lambda_start <- seq(1, 50, length.out = gremlinsEnv$num_lambda_segments)
  segment_membership_start <- calculate_segment_memberships(slope_start, respondent_level_data, respondent_level_designs)
  phi_lambda <- as.numeric(table(segment_membership_start)/length(segment_membership_start))

  startingValues <- list(slope = slope_start,
                         slopeBar = slopeBar_start,
                         slopeCov = slopeCov_start,
                         lambda = lambda_start,
                         segMembership = segment_membership_start,
                         phi_lambda = phi_lambda)

  return(startingValues)
}


#' @importFrom stats runif quantile

calculate_segment_memberships <- function(betai, respondent_level_data, repondent_level_design) {
  # The procedure is:
  # 1. Calculate the in-sample hit rates for each individual
  # 2. Draw a random quantile for each of the information poor segments
  # 3. Find the target hit rate based on the random quantile
  # 4. Split the sample into segment 2 (gremlins) for values
  #    below the target hit rate.  All others are segment 1

  num_respondents <- gremlinsEnv$num_respondents
  num_concepts <- gremlinsEnv$num_concepts
  num_tasks <- gremlinsEnv$num_tasks


  hit_rates <- double(num_respondents)

  # Convert to vectorized version to avoid double loop and speed things up
  for(resp in 1:num_respondents) {
    x_beta <- array((repondent_level_design[[resp]] %*% betai[resp,]), c(num_concepts, num_tasks))
    choices <- apply(x_beta, 2, which.max)
    hit_rates[resp] <- mean(choices == respondent_level_data[[resp]])
  }

  # This is actually a difficult problem if slightly arbitrary choice when you
  # have more than 2 lambda segments.  This code seems like a reasonable option
  # since we don't really have any information on the appropriate sizes for each
  # segment beyond 2.  Basically we ensure that the Information rich segment is
  # at least 1/num_lambda_segments of the total sample  We also ensure that the
  # poorest information segment is at least 5% of the total sample.  These are
  # just starting values the algorithm will settle on the final segment sizes as
  # it converges
  target_quantiles <- sort(runif(gremlinsEnv$num_lambda_segments - 1, 0.05, 1 - 1/gremlinsEnv$num_lambda_segments))
  hit_rate_cut_offs <- as.numeric(quantile(hit_rates, prob=target_quantiles))
  hit_rate_cut_offs <- c(hit_rate_cut_offs, 1.0)

  segments <- double(num_respondents)
  segment_numbers <- gremlinsEnv$num_lambda_segments:1
  lower <- 0.0
  upper <- hit_rate_cut_offs[1]
  for(i in 1:gremlinsEnv$num_lambda_segments) {
    segments[hit_rates >= lower & hit_rates <= upper] <- segment_numbers[i]
    lower <- upper
    upper <- hit_rate_cut_offs[i+1]
  }


  return(segments)
}

# Generate starting values for slope Atchade parameters for tuning runs
# Adjust tuningFactor until accept rates are between 20-80%
# (larger tuningFactor reduce accept rate)
generate_starting_atchade_slopes <- function(slopeBar, tuningFactor) {
  Atch_MCMC_list_slopes <- rep(
    list(
      list(metstd = rep(tuningFactor, times=gremlinsEnv$num_parameters),
           Vmet = rep(1, times=gremlinsEnv$num_parameters),
           mu_adapt=slopeBar,
           Gamma_adapt=rep(1, times=gremlinsEnv$num_parameters)
      )
    ),
    times=gremlinsEnv$num_respondents)
  return(Atch_MCMC_list_slopes)
}

# Generate starting values for slope Atchade parameters for tuning runs
# Adjust tuning_factor up or down to target
# Adjust tuning_factor until accept rates are about 20-80%
# (larger tuning_factor reduce accept rate)
# est_mean_lambda      a number of lambda segments by 1 vector, first element has to be one; needs to be good guess of posterior mean lambda
# est_var_lambda       a number of lambda segments by 1 vector, first element has to be one; needs to be good guess of posterior variance lambda
generate_starting_atchade_lambdas <- function(tuning_factor, est_mean_lambda){

  est_mean_lambda[1] <- 1
  est_var_lambda <- rep(10, length(est_mean_lambda))
  est_var_lambda[1] <- 1

  Atch_MCMC_list_lambda <-
    list( metstd = c(NA, rep(tuning_factor, times= (gremlinsEnv$num_lambda_segments-1))),
          Vmet = c(NA, est_var_lambda[2:gremlinsEnv$num_lambda_segments]),
          mu_adapt=est_mean_lambda,
          Gamma_adapt=est_var_lambda
    )

  return(Atch_MCMC_list_lambda)

}

# Generate the individual slope draws using the Atchade algorithm.  This function should not be exported
#'
#' @importFrom stats rnorm
#' @importFrom bayesm llmnl
#'
# @param data The individuals data
# @param design The individual design
# @param currentSlope The previous slope for the individual
# @param lambda The current lambda
# @param slopeBar The current draw for slopeBar
# @param slopeCov The current draw for slopeCov
# @param constraints The list of constraints on the parameters (none, positive, or negative)
# @param ind_in The individual number used to update the acceptance rate
# @param cur_mcmc_iter Current MCMC Iteration for deciding whether to apply the Atchade step adjustment
# @param metstd_in ...
# @param Vmet_in ...
# @param Gamma_adapt_in ...
# @param mu_adapt_in ...
# @return A list containing the individuals draws for the slope and the the current atchade parameters
# @example # Do not call directly
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

    currentPrior <- -0.5 * (currentSlope - slopeBar)%*%solve(slopeCov)%*%t(currentSlope - slopeBar)
    proposedPrior <- -0.5 * (proposedSlope_i - slopeBar)%*%solve(slopeCov)%*%t(proposedSlope_i - slopeBar)

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

      currentPrior <- -0.5 * (currentSlope - slopeBar)%*%solve(slopeCov)%*%t(currentSlope - slopeBar)
      proposedPrior <- -0.5 * (proposedSlope_i - slopeBar)%*%solve(slopeCov)%*%t(proposedSlope_i - slopeBar)

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

      currentPrior <- -0.5 * (currentSlope - slopeBar)%*%solve(slopeCov)%*%t(currentSlope - slopeBar)
      proposedPrior <- -0.5 * (proposedSlope_i - slopeBar)%*%solve(slopeCov)%*%t(proposedSlope_i - slopeBar)

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

# Generate the individual slope draws using the Atchade algorithm.  This funciton should not be exported
#
# This is a Metropolis-Hastings step with one little twist.  The lambdas are constrained to be in
# increasing magnitude.  This is accomplished using rejection sampling for the truncated distributions.
# This is equivalent to applying a truncated prior, but is more efficient due to the potentially high rejection
# rate on the draws.  Because of this it is necessary to correct for the non-symetric proposal distribution in
# the acceptance probability step.
#'
#' @importFrom stats rnorm
#' @importFrom bayesm llmnl
#'
# @param currentLambda_in The current values for lambda
# @param data_in The data
# @param design_in The design object
# @param currentSlope_in The current slope values
# @param K_in The current segment assignment for each individual
# @param priors_in The Priors
# @param cur_mcmc_iter The current iteration number
# @param metstd_in ...
# @param Vmet_in ...
# @param Gamma_adapt_in ...
# @param mu_adapt_in ...
#
# @example #This code should not be called directly
# @export
generateLambdaMix_ATCH <- function(currentLambda_in, data_in, design_in, currentSlope_in, K_in, priors_in,
                                   cur_mcmc_iter, metstd_in, Vmet_in, Gamma_adapt_in, mu_adapt_in) {




  # currentLambda_in <- lambda_MCMC
  # data_in <-  data
  # design_in <- design
  # currentSlope_in <- slope_MCMC
  # K_in <- K_MCMC
  # priors_in <- Priors
  # cur_mcmc_iter <- (Atch_mcmc_cnt_in+nreps_total)
  # metstd_in <- Atch_MCMC_list_lambda$metstd        # scale parameter for MH variance proposal; one per lambda, first one irrelevant
  # Vmet_in <- Atch_MCMC_list_lambda$Vmet          # variance for MH proposal; one per lambda, first one irrelevant
  # Gamma_adapt_in <- Atch_MCMC_list_lambda$Gamma_adapt            # Approx post varcov of parameters
  # mu_adapt_in <- Atch_MCMC_list_lambda$mu_adapt





  # NOTE: other proposals may work better with Atchade; e.g. sum of exponentials with draws from normal
  # But this would change the prior structure


  #; Draw the lambdas.  (Lambda[1] is always 1)
  if(gremlinsEnv$num_lambda_segments == 1) {
    return(1)
  }

  #  print(priors_in)

  # Set lambda
  lambda_c <- currentLambda_in

  accept_rate_est <- rep(NA, times=gremlinsEnv$num_lambda_segments)

  # Do MH for each lambda separate (except for first one); construct proposal satisfying the order constraint
  # Automatically reject if it doesnt satisfy the order constraint
  for(k in 2:gremlinsEnv$num_lambda_segments){


    #k <- 2


    VMH_k <- metstd_in[k] * Vmet_in[k]

    if(k < gremlinsEnv$num_lambda_segments) {       # UPDATE any lambda that is not the last lambda

      #proposedLambda <- lambda_c[k] + ((lambda_c[k+1] - lambda_c[k-1])/6) * rnorm(1, 0, 1)

      proposedLambda <- lambda_c[k] + sqrt(VMH_k) * rnorm(1, 0, 1)


      # Auto reject of proposal does not satisfy constraint
      if (proposedLambda > lambda_c[k-1] & proposedLambda < lambda_c[k+1]){

        cIndLL <- double(gremlinsEnv$num_respondents)
        pIndLL <- double(gremlinsEnv$num_respondents)

        for(ind in 1:gremlinsEnv$num_respondents) {
          if(K_in[ind] == k) {

            # cIndLL[ind] <- gremlinsLogLikelihood(data_in[ind, -c(1:2)], design_in[design_in[,1] == data_in[ind, 2], -c(1:3)], currentSlope_in[ind,], lambda_c[k])
            # pIndLL[ind] <- gremlinsLogLikelihood(data_in[ind, -c(1:2)], design_in[design_in[,1] == data_in[ind, 2], -c(1:3)], currentSlope_in[ind,], proposedLambda)
            cIndLL[ind] <- llmnl(currentSlope_in[ind,]/lambda_c[k], data_in[[ind]], design_in[[ind]])
            pIndLL[ind] <- llmnl(currentSlope_in[ind,]/proposedLambda, data_in[[ind]], design_in[[ind]])
          }
        }

        cLL <- sum(cIndLL)
        pLL <- sum(pIndLL)
        cPrior <- (priors_in$lambdaShape[k] - 1) * log(lambda_c[k]) - lambda_c[k]/priors_in$lambdaScale[k] # * sum(K_in[i == k])
        pPrior <- (priors_in$lambdaShape[k] - 1) * log(proposedLambda) - proposedLambda/priors_in$lambdaScale[k] # * sum(K_in[i == k])

        lap <- pLL + pPrior - cLL - cPrior
        alpha <- log(runif(1))


        if(alpha < lap) {             # Accept proposal

          #cat("ACCEPT LAMBDA")

          lambda_c[k] <- proposedLambda


          gremlinsEnv$acceptanceRate_lambda[k] <- gremlinsEnv$acceptanceRate_lambda[k] + 1

          accept_rate_est[k] <- 1

        }else{      # REJECT proposal


          accept_rate_est[k] <- 0


        }



      }else{      # Auto reject (proposal not satisfying constraint)


        accept_rate_est[k] <- 0



      } # END IF statement auto reject if contraint is not satisfied




    } else {        # Repeat for LAST lambda; this is automatic for two Gremlin clusters

      #proposedLambda <- lambda_c[k] + ((lambda_c[k] - lambda_c[1])/(6*gremlinsEnv$num_lambda_segments)) * rnorm(1, 0, 1)

      proposedLambda <- lambda_c[k] + sqrt(VMH_k) * rnorm(1, 0, 1)

      #lambda_c[k-1]
      #proposedLambda

      if (proposedLambda > lambda_c[k-1]){

        cIndLL <- double(gremlinsEnv$num_respondents)
        pIndLL <- double(gremlinsEnv$num_respondents)

        for(ind in 1:gremlinsEnv$num_respondents) {
          if(K_in[ind] == k) {

            # cIndLL[ind] <- gremlinsLogLikelihood(data_in[ind, -c(1:2)], design_in[design_in[,1] == data_in[ind, 2], -c(1:3)], currentSlope_in[ind,], lambda_c[k])
            # pIndLL[ind] <- gremlinsLogLikelihood(data_in[ind, -c(1:2)], design_in[design_in[,1] == data_in[ind, 2], -c(1:3)], currentSlope_in[ind,], proposedLambda)
            cIndLL[ind] <- llmnl(currentSlope_in[ind,]/lambda_c[k], data_in[[ind]], design_in[[ind]])
            pIndLL[ind] <- llmnl(currentSlope_in[ind,]/proposedLambda, data_in[[ind]], design_in[[ind]])
          }
        }

        cLL <- sum(cIndLL)
        pLL <- sum(pIndLL)
        cPrior <- (priors_in$lambdaShape[k] - 1) * log(lambda_c[k]) - lambda_c[k]/priors_in$lambdaScale[k] # * sum(K_in[i == k])
        pPrior <- (priors_in$lambdaShape[k] - 1) * log(proposedLambda) - proposedLambda/priors_in$lambdaScale[k] # * sum(K_in[i == k])

        lap <- pLL + pPrior - cLL - cPrior
        alpha <- log(runif(1))


        if(alpha < lap) {             # Accept proposal

          #cat("ACCEPT LAMBDA")

          lambda_c[k] <- proposedLambda


          gremlinsEnv$acceptanceRate_lambda[k] <- gremlinsEnv$acceptanceRate_lambda[k] + 1

          accept_rate_est[k] <- 1

        }else{      # REJECT proposal


          accept_rate_est[k] <- 0


        }



      }else{      # Auto reject (proposal not satisfying constraint)


        accept_rate_est[k] <- 0



      }# END IF statement auto reject if contraint is not satisfied




    }     # END ifelse for last lambda




  }  # END LOOP OVER NUMBER OF LAMBDAS (FIRST ONE IS SKIPPED)


  #lambda_c
  #lambda[rep - 1,]
  #accept_rate_est

  # Need to update all Atchade parameters
  # Loop again over number of Gremlin clusters


  # Update tuning parameters
  if (cur_mcmc_iter>1000){

    g_adapt <- 10/cur_mcmc_iter;

    #    for(ss in 2:gremlinsEnv$num_lambda_segments){


    #ss <- 2


    #     Vmet_in[ss] <- Gamma_adapt_in[ss] + gremlinsEnv$eps2
    Vmet_in <- Gamma_adapt_in + gremlinsEnv$eps2



    #currentSlope_in
    #Gamma_adapt_in
    #mu_adapt_in
    #g_adapt

    mu_adapt_in <- mu_adapt_in + g_adapt*(lambda_c - mu_adapt_in)

    #      dd_x2 <- mu_adapt_in[ss]
    dd_x2 <- sqrt(sum(mu_adapt_in*mu_adapt_in))
    if (dd_x2>gremlinsEnv$A1){
      mu_adapt_in <- (gremlinsEnv$A1/dd_x2)*mu_adapt_in
    }



    # if (dd_x2>gremlinsEnv$A1){
    #   mu_adapt_in[ss] <- (gremlinsEnv$A1/dd_x2)*mu_adapt_in[ss]
    # }


    #      Gamma_adapt_in[ss] <- Gamma_adapt_in[ss] +g_adapt*(((lambda_c[ss] - mu_adapt_in[ss])^2) - Gamma_adapt_in[ss])
    Gamma_adapt_in <- Gamma_adapt_in +g_adapt*(((lambda_c - mu_adapt_in)^2) - Gamma_adapt_in)


    dd_x1 <- sqrt(sum(Gamma_adapt_in*Gamma_adapt_in))
    #dd_x1
    if (dd_x1>gremlinsEnv$A1){
      Gamma_adapt_in <- (gremlinsEnv$A1/dd_x1)*Gamma_adapt_in
    }


    # dd_x1 <- Gamma_adapt_in[ss]
    # if (dd_x1>gremlinsEnv$A1){
    #   Gamma_adapt_in[ss] <- (gremlinsEnv$A1/dd_x1)*Gamma_adapt_in[ss]
    # }



    #metstd_in
    #accept_rate_est

    #metstd_in[ss] <- metstd_in[ss] +g_adapt*(accept_rate_est - gremlinsEnv$tau_bar)
    metstd_in <- metstd_in +g_adapt*(accept_rate_est - gremlinsEnv$tau_bar)

    tau_too_small <- which(metstd_in<gremlinsEnv$eps1)
    tau_too_big <- which(metstd_in>gremlinsEnv$A1)

    if(length(tau_too_small)>0){
      #      if (metstd_in[ss] < gremlinsEnv$eps1){

      metstd_in[tau_too_small] <- gremlinsEnv$eps1

      #      }else if (metstd_in[ss] > gremlinsEnv$A1){

    }else if(length(tau_too_big)>0){



      metstd_in[tau_too_big] <- gremlinsEnv$A1

    }


    #metstd_in


    # if (metstd_in[ss] < gremlinsEnv$eps1){
    #
    #   metstd_in[ss] <- gremlinsEnv$eps1
    #
    # }else if (metstd_in[ss] > gremlinsEnv$A1){
    #
    #   metstd_in[ss] <- gremlinsEnv$A1
    #
    # }


    #metstd_in






    #    } # End loop update Atchade constants for element lambda's


  }  # End update Atchade constants if statement (need enough iterations)




  # Create a list of the Atchade parameters to return
  Atch_out_lambda_list <- list(
    metstd = metstd_in,
    Vmet = Vmet_in,
    mu_adapt = mu_adapt_in,
    Gamma_adapt = Gamma_adapt_in
  )



  ret_list <- list(
    lambda_out = lambda_c,
    Atch_out_lambda_list = Atch_out_lambda_list
  )

  #ret_list

  #return(lambda_c)
  return(ret_list)


}



#' @importFrom stats rmultinom
generateSegmentMembership <- function(data, design, slope, lambda, phi_lambda, priors) {
  prob <- double(gremlinsEnv$num_lambda_segments)
  for(k in 1:gremlinsEnv$num_lambda_segments) {
    # prob[k] <- gremlinsLogLikelihood(data, design, slope, lambda[k])
    prob[k] <- llmnl(slope/lambda[k], data, design)
    prob[k] <- prob[k] + log(phi_lambda[k])
  }

  prob <- exp(prob)
  return(which.max(rmultinom(1, 1, prob)))
}

# Utility function to compute a multivariate regression.  Based on information from "Bayesian Statistics and Marketing (2005)"
#' @importFrom stats rWishart
multivariateRegression <- function(data, design, priors) {
  U = chol(priors$Ainv)
  R = rbind(design, U)
  Q = rbind(data, U %*% priors$mu_not)
  BTilde = chol2inv(chol(crossprod(R))) %*% crossprod(R, Q)

  S = crossprod(Q - R%*%BTilde)

  Sigma <- chol2inv(chol(drop(rWishart(1, priors$nu_not + nrow(design), chol2inv(chol(priors$V_not + S))))))

  mean <- as.vector(BTilde)
  variance <- Sigma %x% chol2inv(chol(crossprod(design) + priors$Ainv))

  Beta <- mean + rnorm(length(mean), 0, 1) %*% chol(variance)

  return(list(Beta = Beta, Sigma = Sigma))

}


set_atchade_slope_parameters <- function(slopeBar, atchade_tuning_factor, nParameters, nRespondents) {
  slope_parameters <- rep(
    list(
      list(
        metstd = rep(atchade_tuning_factor, times=nParameters),
        Vmet = rep(1, times=nParameters),
        mu_adapt = slopeBar,
        Gamma_adapt = rep(1, times=nParameters)
      )
    ),
  times= nRespondents)

  return(slope_parameters)
}

set_atchade_lambda_parameters <- function(lambda, nClusters, estimated_mean_lambda, estimated_var_lambda) {
  estimated_mean_lambda[1] <- 1
  estimated_var_lambda[1] <- 1

  lambda_parameters <- list(
    metstd = c(NA, rep(lambda, times = (nClusters - 1))),
    Vmet = c(NA, estimated_var_lambda[2:nClusters]),
    mu_adapt = estimated_mean_lambda,
    Gamma_adapt = estimated_var_lambda
  )

  return(lambda_parameters)
}
