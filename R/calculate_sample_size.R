
CalcSWTSampleSize <- function(verbose=FALSE){
  library(ggplot2)
  
  # initialize clusters / individuals 
  n_clust_per_seq = seq(5,15,by=1) 
  n_ind_per_clust = c(1:100)
  power = 0
  epsilon_power = 0.01
  maxIter = 15
  powerTarget = 0.8
  
  # initialize placeholders for sample size/powers 
  sample_size_list = rep(0, times=length(n_clust_per_seq))
  power_list = rep(0, times=length(n_clust_per_seq))
  clust_index = 1
  
  for (n_clust in n_clust_per_seq) { 
    low = 1
    high = length(n_ind_per_clust)
    power = 0
    iterCount = 1
    maxIterCheck = FALSE
    convergenceCheck = FALSE 
    powerPrevValue = 0
    
    # loop to calculate the closest power for specific sample size using binary search 
    while(power > powerTarget + epsilon_power | power < powerTarget-epsilon_power) { 
      
      # stop if can't find close enough power 
      if (iterCount > maxIter) { 
        maxIterCheck = TRUE 
        break
      }
      
      # check midpoint for current bounds  
      mid =  low + (high - low) %/% 2 # floor((low+high) / 2)
      
      
      ##### calculate the power for example
      power = ReplicatePrevSciPaperExample_CalcSampleSize(n_clust_per_seq=n_clust, n_ind_per_clust=n_ind_per_clust[mid])$power
      #####
    
      # stop if new power did not change - meaning the search stopped 
      if (power == powerPrevValue) { 
        convergenceCheck = TRUE 
        break
      }
      
      # reassign previous value 
      powerPrevValue = power
      
      # if verbose = report output for each iteration 
      if (verbose) { 
        print(n_clust)
        print(paste('power', power))
        print(paste('low:', low, 'n_indiv:', mid, 'high:', high))
      }
      
      
      # begin actual binary search 
      if (high - low > 1) { 
        
        # if power too high, check lower half 
        if (power > powerTarget + epsilon_power) {
          high = mid
        }
        
        # if power too low, check upper half 
        else if (power < powerTarget - epsilon_power) {
          low = mid
        }
      }
      
      # case where low and high are just 1 apart, floor wont reach high 
      else { 
          low = high
      }
      # binary search iteration 
      iterCount = iterCount + 1
    }
    
    # if max iterations reached 
    if (maxIterCheck) { 
      print(paste('Warning: max iterations reached for epsilon=', epsilon_power, 'of', powerTarget, 'and sample size=', n_clust))
    }
    # if convergence not within epsilon 
    if (convergenceCheck) { 
      print(paste('Warning: power did not converge within epsilon=', epsilon_power, 'of', powerTarget, 'and for sample size=', n_clust))
    }
    
    # add power / sample to output 
    sample_size_list[clust_index] = mid
    power_list[clust_index] = round(power, 3) 
    clust_index = clust_index + 1 
  }
  
  # create final dataframe for export 
  results = data.frame(num_clust = n_clust_per_seq, sample_size = sample_size_list, power = power_list)
  
  # plot results 
  plot_sample_vs_cluster = ggplot(results, aes(x = num_clust, y = sample_size)) +
    geom_point(alpha = 0.6) +  # plot each point 
    geom_smooth(method = "loess", col = "blue", se = FALSE) +  # add curve 
    geom_text(aes(label = power), vjust = -0.5, hjust = 0.5) +  # add overlay of powers 
    labs(x = "Cluster Size", y = "Estimated Sample Size", title = "Sample Size vs Cluster Size")
  print(plot_sample_vs_cluster)
  return(results)

}

.GuessSampleUpperBound <- function(){
  stop("Not implemented yet")
}

.BisectionSearch <- function(){
  stop("Not implemented yet")
}