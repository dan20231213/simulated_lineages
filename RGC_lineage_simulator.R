#-------------------------------parameter setting-----------------------------#
set.seed(789)
sample_size <- 254
simulated_lineage <- list()
division_mode <- c("IPP", "IP", "N")

# Rule 1: The initial energy of the first generation RGC follows a mixing Poisson distribution.
# Energy ~ alpha1 *Poisson(lambda1) + alpha2 *Poisson(lambda2)
# Energy indicates the number of neurons that each progenitor cell can eventually produce.

# 0.35 * P(lambda = 4.09) + 0.65 * P(lambda = 7.62)
lambda1 = 4.09
lambda2 = 7.62
alpha1 = 0.35
alpha2 = 0.65

# Rule 2: Each asymmetric split produces the next generation of RGC and one IPP or IP or N.
# IPP: progenitor of intermediate progenitor cell
# IP: intermediate progenitor cell
# N: neuron
# Energy is reduced by IPP-4; IP-2; N-1 respectively.

eloss <- c(4, 2, 1)

# Rule 3: The probability ratio of the first generation RGC to generate IPP or IP or N is 1:1:1.
# From the second generation onwards, the probability of RGC_n's division mode depends only on RGC_{n-1}'s.

odd0 <- c(1, 1, 1)
odd1 <- data.frame(
  IPP = c(0.59, 0.21, 0.2),
  IP  = c(0.2, 0.58, 0.22),
  N   = c(0.13, 0.37, 0.5)
)

# Rule 4: The termination condition is that the generation energy is exhausted.
# When energy is less than 4, IPP cannot be produced.
# When energy is less than 2, IPP and IP cannot be produced.


#-------------------------------functions outline-----------------------------#

# Function 1: generates the initial potential distribution (energy)
# mixing_Poisson

mixing_Poisson <- function(size, a1, a2, l1, l2){
  p1 = dpois(0:20, lambda = l1)
  p1 = p1[3:21]
  p1 = p1/sum(p1)
  
  p2 = dpois(0:20, lambda = l2)
  p2 = p2[3:21]
  p2 = p2/sum(p2)
  
  p = a1*p1 + a2*p2
  energy = sample(2:20, size, replace = TRUE, prob = p)
  
  return(energy)
}

# Function 2: determines remaining energy
# E_residual
# input: e (energy)
# output: a signal passed to conditional_P

E_residual <- function(e){
  if(e>=4)
    signal = c(1, 1, 1)
  else if(e>=2)
    signal = c(0, 1, 3)
  else if(e>=1)
    signal = c(0, 0, 1)
  
  return(signal)
}

# Function 3: generates the conditional probability of division mode
# conditional_P 
# input: signal; RGC chain
# output: a set of probability values c(p1,p2,p3) corresponding to c("IPP","IP","N")

conditional_P <- function(signal, chain){
  prior <- tail(chain, 1)
  if(is.null(prior))
    odd = odd0
  else
    odd = odd1[1:3, prior]
  odd = odd * signal
  p = odd/sum(odd)
  
  return(p)
}

# Function 4: gets sample from division_mode according to conditional probability
# divide
# input: The probability value p from conditional_P; RGC chain
# output: a new RGC chain

divide <- function(p, chain){
  cell <- sample(x = division_mode, size = 1, replace = FALSE, prob = p)
  chain = append(chain, cell)
  
  return(chain)
}

# Function 5: deducts energy
# E_subtraction
# input: e (residual energy); RGC chain (an RGC lineage being generated)
# output: the remaining energy

E_subtraction <- function(e, chain){
  idx = which(division_mode == tail(chain, 1))
  e = e - eloss[idx]
  
  return(e)
}

#-------------------------------simulation process----------------------------#

energy = mixing_Poisson(sample_size, alpha1, alpha2, lambda1, lambda2)

for(e in energy){
  chain <- c()
  while(e>0){
    signal = E_residual(e)
    # calculate probability according to the tail of RGC chain
    p = conditional_P(signal, chain)
    # sample
    chain = divide(p, chain)
    # deduct energy
    e = E_subtraction(e, chain)
  }
  simulated_lineage <- append(simulated_lineage,list(chain))
}

#----------------------------descriptive statistics---------------------------#

# the frequency of generation in simulated_lineages

generation <- c()
for(i in 1:length(energy)){
  generation = append(generation, length(simulated_lineage[[i]]))
}
table(generation)

# proportion of IPP, IP, N division mode per generation

counter1 <- matrix(0, nrow = 3, ncol = max(generation))
rownames(counter1) <- c("IPP", "IP", "N")
for(i in 1:length(energy)){
  for(j in 1:length(simulated_lineage[[i]])){
    dt = simulated_lineage[[i]][j]
    counter1[dt,j] = counter1[dt,j] + 1
  }
}

# count division modes of RGC_n and RGC_{n-1}

counter2 <- matrix(0, nrow = 9, ncol = max(generation)-1)
rownames(counter2) <- c("IPP-IPP", "IPP-IP", "IPP-N",
                        "IP-IPP", "IP-IP", "IP-N",
                        "N-IPP", "N-IP", "N-N")
# depends on table(generation)
# needs being geared manually
colnames(counter2) <- c("1-2", "2-3", "3-4", "4-5", "5-6", "6-7",
                        "7-8", "8-9", "9-10")
for(i in 1:length(energy)){
  chain =  simulated_lineage[[i]]
  if(length(chain)>1){
    for(j in 2:length(chain)){
      dt = paste0(chain[j], "-", chain[j-1])
      counter2[dt,j-1] = counter2[dt,j-1] +1
    }
  }
}

#---------------------------repeat 100 times randomly-------------------------#
path = "C:/Users/23111/Desktop/simulated_data/"
for(i in 1:100){
  energy = mixing_Poisson(sample_size, alpha1, alpha2, lambda1, lambda2)
  simulated_lineage <- list()
  for(e in energy){
    chain <- c()
    while(e>0){
      signal = E_residual(e)
      p = conditional_P(signal, chain)
      chain = divide(p, chain)
      e = E_subtraction(e, chain)
    }
    simulated_lineage <- append(simulated_lineage,list(chain))
  }
  
  max_length <- max(sapply(simulated_lineage, length))
  simulated_lineage <- lapply(simulated_lineage, function(x) c(x, rep(NA, max_length - length(x))))
  df <- do.call(rbind, lapply(simulated_lineage, function(x) as.data.frame(t(x))))
  write.csv(df, paste0(path, i, ".csv"))
}

#------------------------------------summary----------------------------------#
# how many generations each lineage divide
s1 = as.data.frame(table(apply(df, 1, function(x) sum(!is.na(x)))))

# neuron number produced in each generation
s2 <- matrix(0, nrow = dim(df)[2], ncol = 1)
i = 1
for(col in df){
  PyN = 0
  for(char in col){
    if(!is.na(char))
      PyN = PyN + eloss[which(division_mode==char)]
  }
  s2[i,1] = PyN
  i = i + 1
}


# N/IP/IPP mode percentage in each generation

s3 <- matrix(0, nrow = 3, ncol = dim(df)[2])
rownames(s3) <- c("IPP", "IP", "N")
for(i in 1:dim(df)[1]){
  for(j in 1:dim(df)[2]){
    dt = df[i,j]
    if(!is.na(dt))
      s3[dt,j] = s3[dt,j] + 1
  }
}

# division modes and contribution

s4 <- as.data.frame(table(unlist(df)))
s4$con <- s4$Freq*c(2, 4, 1)

# all
combined_s1 <- data.frame()
combined_s2 <- data.frame()
combined_s3_N <- data.frame()
combined_s3_IP <- data.frame()
combined_s3_IPP <- data.frame()
combined_s4 <- data.frame()
combined_s5 <- data.frame()

for(f in 1:100){
  df <- read.csv(paste0(path, f, ".csv"))
  df = df[, -1]
  
  s1 = as.data.frame(table(apply(df, 1, function(x) sum(!is.na(x)))))
  colnames(s1) <- c("G", paste0("lineage", f))
  if(f == 1)
    combined_s1 <- s1
  else{
    combined_s1 <- merge(combined_s1, s1, by = "G", all = TRUE)
    combined_s1[is.na(combined_s1)] <- 0
  }
  
  s2 <- matrix(0, nrow = dim(df)[2], ncol = 1)
  i = 0
  for(col in df){
    PyN = 0
    for(char in col){
      if(!is.na(char))
        PyN = PyN + eloss[which(division_mode==char)]
    }
    i = i + 1
    s2[i,1] = PyN
  }
  s2 <- as.data.frame(s2)
  colnames(s2) <- c(paste0("lineage", f))
  s2$G <- 1:i
  if(f == 1)
    combined_s2 <- s2
  else{
    combined_s2 <- merge(combined_s2, s2, by = "G", all = TRUE)
    combined_s2[is.na(combined_s2)] <- 0
  }
  
  s3 <- matrix(0, nrow = 3, ncol = dim(df)[2])
  rownames(s3) <- c("IPP", "IP", "N")
  for(i in 1:dim(df)[1]){
    for(j in 1:dim(df)[2]){
      dt = df[i,j]
      if(!is.na(dt))
        s3[dt,j] = s3[dt,j] + 1
    }
  }
  s3 <- as.data.frame(s3)
  percentage_s3 <- as.data.frame(lapply(s3, function(x) x / sum(x) * 100))
  rownames(percentage_s3) <- rownames(s3)
  if(f==1){
    combined_s3_N = percentage_s3["N",1:6]
    combined_s3_IP = percentage_s3["IP",1:6]
    combined_s3_IPP = percentage_s3["IPP",1:6]
  }
  else{
    combined_s3_N = rbind(combined_s3_N, percentage_s3["N",1:6])
    combined_s3_IP = rbind(combined_s3_IP, percentage_s3["IP",1:6])
    combined_s3_IPP = rbind(combined_s3_IPP, percentage_s3["IPP",1:6])
  }
  
  s4 <- as.data.frame(table(unlist(df)))
  s5 <- data.frame(mode = c("IP", "IPP", "N"),
                   contribution = s4$Freq*c(2, 4, 1))
  colnames(s5) <- c("mode", paste0("contribution", f))
  colnames(s4) <- c("mode", paste0("freq", f))
  if(f == 1){
    combined_s4 <- s4
    combined_s5 <- s5
  }
  else{
    combined_s4 <- merge(combined_s4, s4, by = "mode", all = TRUE)
    combined_s5 <- merge(combined_s5, s5, by = "mode", all = TRUE)
  }
}

save_path = "C:/Users/23111/Desktop/results/"
write.csv(combined_s1, paste0(save_path, "g", ".csv"))
write.csv(combined_s2, paste0(save_path, "h", ".csv"))
write.csv(combined_s3_N, paste0(save_path, "iN", ".csv"))
write.csv(combined_s3_IP, paste0(save_path, "iIP", ".csv"))
write.csv(combined_s3_IPP, paste0(save_path, "iIPP", ".csv"))
write.csv(combined_s4, paste0(save_path, "j1", ".csv"))
write.csv(combined_s5, paste0(save_path, "j2", ".csv"))
