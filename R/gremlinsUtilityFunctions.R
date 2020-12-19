# Internal Fucntion should not be exported

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
  num_tasks <- nrow(bayesm_data$lgtdata[[resp]]$X)/num_concepts

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
