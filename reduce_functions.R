############ Get.Stable Function ###############
#This function reduces stable Sr measurements using the derivation described in Krabbenhoft et al. 2009

### Inputs:
# Spiked Raw Ratios - these should be vectors the length of however many cycles were run for the sample
# id_86_84 > raw spiked 86Sr/84Sr
# id_87_84 > raw spiked 87Sr/84Sr
# id_88_84 > raw spiked 88Sr/84Sr

# Unspiked Measurement - this should be a double
# ic_87_86 > average unspiked 87Sr/86Sr

# Spike Composition - this is a tibble with 1 row and 3 columns, sp_86_84, sp_87_84, and sp_88_84 respectively
# spike_comp > the spike composition you're using

get.stable <- function(id_86_84, id_87_84, id_88_84, ic_87_86, spike_comp) {
  
  #Get variables organized
  id_86_84 <- na.omit(id_86_84)
  id_87_84 <- na.omit(id_87_84)
  id_88_84 <- na.omit(id_88_84)
  
  if(length(id_86_84) != length(id_87_84) | length(id_86_84) != length(id_88_84)) {
    stop("spiked ratio columns are not the same length")
  }
  
  id_86_84 <- sapply(id_86_84, as.numeric)
  id_87_84 <- sapply(id_87_84, as.numeric)
  id_88_84 <- sapply(id_88_84, as.numeric)
  
  ic_87_86 <- as.numeric(ic_87_86)
  
  ic_84_86 <- 0.056532661
  ic_88_86 <- 8.375209
  
  sp_86_84 <- as.numeric(spike_comp[["sp_86_84"]])
  sp_87_84 <- as.numeric(spike_comp[["sp_87_84"]])
  sp_88_84 <- as.numeric(spike_comp[["sp_88_84"]])
  
  amu_84 <- 83.91343
  amu_86 <- 85.90926
  amu_87 <- 86.90888
  amu_88 <- 87.90561
  
  ic_86_84 <- 1/ic_84_86
  ic_87_84 <- ic_87_86/ic_84_86
  ic_88_84 <- ic_88_86/ic_84_86
  sp_88_86 <- sp_88_84/sp_86_84
  

  #Iteration 1
  first_iteration <- TRUE
  if(first_iteration == TRUE) {
    Q_86_84 <- (id_86_84 - sp_86_84) / (ic_86_84 - id_86_84) #eq.11
    Q_88_84 <- (id_88_84 - sp_88_84) / (ic_88_84 - id_88_84) #eq.12
    
    Q_84 <- (Q_86_84 + Q_88_84) / 2 #eq.13
    
    calc_87_84 <- ((Q_84*ic_87_84) + sp_87_84)/(1+Q_84) #eq.1 (mix_calc)
    
    beta <- log(calc_87_84/id_87_84)/log(amu_87/amu_84) #eq.2
    
    calc_86_84 <- id_86_84*((amu_86/amu_84)^beta) #eq.3 (mix_calc)
    calc_88_84 <- id_88_84*((amu_88/amu_84)^beta) #eq.4 (mix_calc)
    calc_88_86 <- calc_88_84/calc_86_84 #eq.5 (mix_calc)
    
    Q_86 <- Q_84*(ic_86_84/sp_86_84) #eq.6
    
    new_88_86 <- ((1 + (1/Q_86))*(calc_88_86)) - ((1/Q_86)*sp_88_86) #eq.7 (pureasmple_recalc)
    
    beta_new <- log(new_88_86/ic_88_86)/log(amu_88/amu_86) #eq.8 (alpha)
    
    new_88_84 <- ic_88_84*((amu_88/amu_84)^beta_new) #eq.9 (pureasmple_recalc)
    new_86_84 <- ic_86_84*((amu_86/amu_84)^beta_new) #eq.10 (pureasmple_recalc)
    new_87_84 <- ic_87_84*((amu_87/amu_84)^beta_new) #(pureasmple_recalc)
    
    first_iteration <- FALSE
    
  }
  
  #Iterations 2-200
  if(first_iteration == FALSE) {
    for(i in 1:200) {
      prev_88_86 <- new_88_86 #this is the reassignment of puresample_recalc to puresample in maincalc.txt; initially puresample is the ic measurement
      prev_88_84 <- new_88_84
      prev_86_84 <- new_86_84
      prev_87_84 <- new_87_84
      
      mix_88_84 <- calc_88_84 #this is the reassignment of mix_calc to mix in maincalc.txt; initially mix is the id measurement
      mix_86_84 <- calc_86_84
      mix_87_84 <- calc_87_84
      
      Q_86_84 <- (mix_86_84 - sp_86_84) / (prev_86_84 - mix_86_84) #eq.11
      Q_88_84 <- (mix_88_84 - sp_88_84) / (prev_88_84 - mix_88_84) #eq.12
      
      Q_84 = (Q_86_84 + Q_88_84) / 2 #eq.13
      
      calc_87_84 <- ((Q_84*prev_87_84) + sp_87_84)/(1+Q_84) #eq.1 (mix_calc)
      
      beta <- log(calc_87_84/mix_87_84)/log(amu_87/amu_84) #eq.2
      
      calc_86_84 <- mix_86_84*((amu_86/amu_84)^beta) #eq.3 (mix_calc)
      calc_88_84 <- mix_88_84*((amu_88/amu_84)^beta) #eq.4 (mix_calc)
      calc_88_86 <- calc_88_84/calc_86_84 #eq.5 (mix_calc)
      
      Q_86 <- Q_84*(prev_86_84/sp_86_84) #eq.6
      
      new_88_86 <- ((1 + (1/Q_86))*(calc_88_86)) - ((1/Q_86)*sp_88_86) #eq.7 (puresample_recalc)
      
      beta_new <- log(new_88_86/prev_88_86)/log(amu_88/amu_86) #eq.8 (alpha)
      
      new_88_84 <- prev_88_84*((amu_88/amu_84)^beta_new) #eq.9 (puresample_recalc)
      new_86_84 <- prev_86_84*((amu_86/amu_84)^beta_new) #eq.10 (puresample_recalc)
      new_87_84 <- prev_87_84*((amu_87/amu_84)^beta_new) # (puresample-recalc)
      
      #Check if we need more iterations...
      diff <- abs(Q_86_84 - Q_88_84)
      av_diff <- mean(diff)
      if(av_diff < 1e-17) {
        break
      } else {counter <- 1 + i}
      
    }
    
    # Outputs of Interest
    results <- list()
    results[["delta"]] <- mean(((new_88_86/8.375209) - 1)*1000)
    results[["delta_2se"]] <- 2*(sd(((new_88_86/8.375209) - 1)*1000))/sqrt(length(new_88_86))
    results[["sample_87_86"]] <- mean(new_87_84)/mean(new_86_84)
    results[["mix_84_86"]] <- 1/mean(calc_86_84)
    
    return(results)
    
  }
}


############# Batch.Reduce Function ##############
# This function applies get.stable to a batch of samples
# The code uses the given number of samples to index through a single table input with all the data

### Inputs
# n_samples > this is the number of samples your reducing - it's very important!!
# id_data > this is a tibble with all of your data in it. It should be have n_samples*3 columns, and the row lengths should be consistent with a sample but may vary between
#        depending on how many cycles you have for each measurement. For each sample, the columns should be in the order of 86/84, 88/84, 87/84
# ic_data > this is a tibble with mean unspiked 87Sr/86Sr ratios. It should have two columns and n_samples rows. The first column should be sample names as strings, 
#           and the second columns should contain the ratios. It is is critical that the first three columns in id_data are from the same sample as row one of ic_data, 
#           and so on.
# spike_comp > this is the spike composition that will be applied to all samples. It should be a tibble with three columns and one row, where the columns are 
#              sp_86_84, sp_87_84, and sp_88_84 respectively


batch.reduce <- function(n_samples, id_data, ic_data, spike_comp) {
  
  # Set up a table for results
  results <- tibble(sample_ID = character(length = n_samples), delta = numeric(length = n_samples), delta_2se = numeric(length = n_samples), sample_87_86 = numeric(length = n_samples), mix_84_86 = numeric(length = n_samples))
  
  
  # For-Loop
  for (i in 1:n_samples) {
    
    # i is the index for the current sample
    
    current_output <- get.stable(id_86_84 = id_data[,(((i-1)*3) + 1)], id_87_84 = id_data[,(((i-1)*3) + 3)], id_88_84 = id_data[,(((i-1)*3) + 2)], ic_87_86 = ic_data[i,2], spike_comp = spike_comp)
    
    results$sample_ID[i] <- ic_data[i,1]
    results$delta[i] <- current_output[["delta"]]
    results$delta_2se[i] <- current_output[["delta_2se"]]
    results$sample_87_86[i] <- current_output[["sample_87_86"]]
    results$mix_84_86[i] <- current_output[["mix_84_86"]]
    
  }
  
 return(results)
   
}



























