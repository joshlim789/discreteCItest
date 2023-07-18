# Setup -------------------------------------------------------------------
suppressMessages(library(bnlearn))
library(Rcpp)
sourceCpp("pDCorTest.cpp")

asia_mat <- as.matrix(asia)
asia_coln <- colnames(asia)

# Testing Same P-Values -------------------------------------------------------------------
test_that("Make sure you are obtaining the correct test results",{
  
  #Test all pairs of nodes
    for (i in 0:7){
      for (j in 0:7){
        
        #bnlearn ci.test won't accept when i == j
        if(i == j)
        {
          next
        }
        
        alpha <- 0.05
        set.seed(24601)
        
        # Comparing results with empty set
        k0 <- double()
        true_res <- ci.test(data = asia, asia_coln[i+1], asia_coln[j+1])
        est_res <- condInttestdis(asia_mat,i,j,k0,signif_level=alpha)
        # Compare p-values
        expect_equal(est_res$pval,true_res$p.value)
        
        #Comparing results with set of size 1
        # choose random conditioning set of size 1
        k1 <- sample(setdiff(0:7,c(i,j)),1,replace = FALSE)
        true_res <- ci.test(data = asia, asia_coln[i+1], asia_coln[j+1], asia_coln[k1+1])
        est_res <- condInttestdis(asia_mat,i,j,k1,signif_level=alpha)
        # Message if results do not match
        if (!all.equal(est_res$pval, true_res$p.value)){
          cat("\ni =",i,"| j =",j,"| k =",paste(k1,collapse = " "),
              "| pval (Est):",est_res$pval,"| pval (True): ",true_res$p.value)
        }
        # Compare results
        expect_equal(est_res$pval,true_res$p.value)
        
        #Comparing results with set of size 2
        # choose random conditioning set of size 2
        k2 <- sample(setdiff(0:7,c(i,j)),2,replace = FALSE)
        true_res <- ci.test(data = asia, asia_coln[i+1], asia_coln[j+1], asia_coln[k2+1])
        est_res <- condInttestdis(asia_mat,i,j,k2,signif_level=alpha)
        # Message if results do not match
        if (!all.equal(est_res$pval, true_res$p.value)){
          cat("\ni =",i,"| j =",j,"| k =",paste(k2,collapse = " "),
              "| pval (Est):",est_res$pval,"| pval (True): ",true_res$p.value)
        }
        # Compare results
        expect_equal(est_res$pval,true_res$p.value)
        
        #Comparing results with set of size 3
        # choose random conditioning set of size 3
        k3 <- sample(setdiff(0:7,c(i,j)),3,replace = FALSE)
        true_res <- ci.test(data = asia, asia_coln[i+1], asia_coln[j+1], asia_coln[k3+1])
        est_res <- condInttestdis(asia_mat,i,j,k3,signif_level=alpha)
        # Message if results do not match
        if (!all.equal(est_res$pval, true_res$p.value)){
          cat("\ni =",i,"| j =",j,"| k =",paste(k3,collapse = " "),
              "| pval (Est):",est_res$pval,"| pval (True): ",true_res$p.value)
        }
        # Compare results
        expect_equal(est_res$pval,true_res$p.value)
        
        #Comparing results with set of size 6
        k6 <- setdiff(0:7,c(i,j))
        true_res <- ci.test(data = asia, asia_coln[i+1], asia_coln[j+1], asia_coln[k6+1])
        est_res <- condInttestdis(asia_mat,i,j,k6,signif_level=alpha)
        # Message if results do not match
        if (!all.equal(est_res$pval, true_res$p.value)){
          cat("\ni =",i,"| j =",j,"| k =",paste(k6,collapse = " "),
              "| pval (Est):",est_res$pval,"| pval (True): ",true_res$p.value)
        }
        # Compare results
        expect_equal(est_res$pval,true_res$p.value)
      }
      
    }
})
