## Probably newer versions

Xcredible_cutoff <- function(target = 0.99,prior,object) {
  X0 <- object@X0
  pd <- postEst(prior,object)
  top <- quantile(pd,target)
  
  for (index in seq_along(pd)) {
    value <- pd[index]
    if (value >= top){
      if (X0[index] > 1) {
        co <- X0[index]
        return_list <- list(posterior_distribution = pd, 
          posterior_prob = top, cutoff = co)
        return(return_list)
      } else {
        return_list <- list(posterior_distribution = pd, 
          posterior_prob = top, cutoff = NA)
        return(return_list)
      }
    }
  }
}

Xcalculate_scores <- function(true,pred) {
  true_pos = 0
  false_pos = 0
  false_neg = 0
  true_neg = 0
  for (index in 1:length(true)) {
    real <- true[index]
    guess <- pred[index] 
    
    if (real == guess) {
      true_pos = true_pos + real
    }
    else if (guess > real) {
      true_pos = true_pos + real
      false_pos = false_pos + (guess - real)
    }
    else if (real > guess) {
      true_pos = true_pos + guess
      false_neg = false_neg + (real - guess)
    } else if (real == 0 & guess == 0) {
      true_neg = true_neg + 1
    }
  }
  specificity = true_neg / length(true == 0)
  precision = true_pos / (true_pos + false_pos)
  recall = true_pos / (true_pos + false_neg)
  f1 = (2 * precision * recall) / (precision + recall)
  return(list(precision, recall, f1, specificity))
}

Xbinary_calculate_scores <- function(true,pred) {
  true_pos = 0
  false_pos = 0
  false_neg = 0
  true_neg = 0
  for (index in 1:length(true)) {
    real <- true[index]
    guess <- pred[index]
    
    if (real == guess) {
      true_pos = true_pos + 1
    }
    else if (guess > real) {
      false_pos = false_pos + 1
    }
    else if (real > guess) {
      false_neg = false_neg + 1
    } else if (real == 0 & guess == 0) {
      true_neg = true_neg + 1
    }
  }
  specificity = true_neg / length(true == 0)
  precision = true_pos / (true_pos + false_pos)
  recall = true_pos / (true_pos + false_neg)
  f1 = (2 * precision * recall) / (precision + recall)
  return(list(precision, recall, f1, specificity))
}

# calculate_counts <- function(count, sig_count,pred_num,list) {
#   if(pred_num == sig_count) {
#     list[[1]] <- list[[1]] + sig_count # TP
#     list[[3]] <- list[[3]] + (count - sig_count) # TN
#     # list[[3]] <- if(pred_num == 0 & sig_count == 0) {
#     #   list[[3]] + 1
#     #  } else {
#     #   list[[3]] <- list[[3]]
#     # } # TN
#   }
#   else if (pred_num > sig_count) {
#     list[[1]] <- list[[1]] + (sig_count) # TP
#     list[[2]] <- list[[2]] + (pred_num - sig_count) # FP
#     list[[3]] <- list[[3]] + (count - pred_num) # TN
#   }
#   else if (sig_count > pred_num) {
#     list[[1]] <- list[[1]] + pred_num # TP
#     list[[3]] <- list[[3]] + (count - sig_count) # TN
#     list[[4]] <- list[[4]] + (sig_count - pred_num) # FN
    
#   }
#   return(list)
# }

Xcalculate_counts <- function(total, sig_loops, called_loops, ret_list) {
  if (called_loops == sig_loops) {
    TP <- sig_loops
    FN <- 0
    FP <- 0
    TN <- ifelse(called_loops == 0 & sig_loops == 0, 1, 0)
  } else if (called_loops > sig_loops) {
    TP <- sig_loops
    FN <- 0
    FP <- called_loops - sig_loops
    TN <- 0
  } else if (sig_loops > called_loops) {
    TP <- called_loops
    FN <- sig_loops - called_loops
    FP <- 0
    TN <- 0
  }
  ret_list[[1]] <- ret_list[[1]] + TP
  ret_list[[2]] <- ret_list[[2]] + FP
  ret_list[[3]] <- ret_list[[3]] + TN
  ret_list[[4]] <- ret_list[[4]] + FN
  return(ret_list)
}
