library(dplyr)
Final <- NULL

for (r in 1:Rep) {
  QFinal <- AFinal <- QFinal2 <- AFinal2 <- NULL
  for (nn in N) {
    for (d in dist) {
      cov <- paste0('X', 1:pp)
      set.seed(r)
      # X matrix by distribution
      if (d == 'bin') {
        Xmat <- matrix(rbinom(nn*pp, 1, prob=qq), nrow=nn, ncol=pp)
      } else if (d == 'norm') {
        Xmat <- matrix(rnorm(nn*pp, 0, 1), nrow=nn, ncol=pp)
      } else {
        Xmat <- cbind(matrix(rbinom(nn*2, 1, prob=qq), nrow=nn, ncol=2),
                      matrix(rnorm(nn*(pp-2), 0, 1), nrow=nn, ncol=pp-2))
      }
      
      colnames(Xmat) <- paste0('X', 1:pp)
      
      if (d == 'bin') {
        VXmat <- matrix(rbinom(Vnn*pp, 1, prob=qq), nrow=Vnn, ncol=pp)
      } else if (d == 'norm') {
        VXmat <- matrix(rnorm(Vnn*pp, 0, 1), nrow=Vnn, ncol=pp)
      } else {
        VXmat <- cbind(matrix(rbinom(Vnn*2, 1, prob=qq), nrow=Vnn, ncol=2),
                       matrix(rnorm(Vnn*(pp-2), 0, 1), nrow=Vnn, ncol=pp-2))
      }
      colnames(VXmat) <- paste0('X', 1:pp)
      
      ## for Linear pi(x), mu(x)
      Xmat2 <- cbind(rep(1,nn), Xmat)
      colnames(Xmat2) <- c("Intercept", paste0('X', 1:pp))
      
      VXmat2 <- cbind(rep(1,Vnn), VXmat)
      colnames(VXmat2) <- c("Intercept", paste0('X', 1:pp))
      
      ## for Linear tau(x)
      Tmat <- Tmat <- cbind(1, Xmat)
      colnames(Tmat) <- c("A", paste0('AX', 1:pp))
      
      VTmat <- cbind(1, VXmat)
      colnames(VTmat) <- c("A", paste0('AX', 1:pp))
      
      ## for Non-Linear pi(x)
      NL_Xmat2 <- cbind(2*(Xmat[,'X1']*Xmat[,'X5']),
                        ((Xmat[,'X1']>0.5) != (Xmat[,'X3']>0.5) ))
      
      VNL_Xmat2 <- cbind(2*(VXmat[,'X1']*VXmat[,'X5']),
                         ((VXmat[,'X1']>0.5) != (VXmat[,'X3']>0.5) ))
      
      ## for Non-Linear mu(x)
      NL_Xmat3 <- cbind((Xmat[,'X1']>0.5) != (Xmat[,'X2']>0.5)*1, 
                        sin(Xmat[,'X2']*Xmat[,'X5']+2*Xmat[,'X6'] ))
      
      VNL_Xmat3 <- cbind((VXmat[,'X1']>0.5) != (VXmat[,'X2']>0.5)*1, 
                         sin(VXmat[,'X2']*VXmat[,'X5']+2*VXmat[,'X6'] ))
      
      ## for Non-Linear tau(x)
      NL_Tmat <- cbind( 1, exp(Xmat[,'X1']), Xmat[,'X1']*Xmat[,'X2'], 
                        ((Xmat[,'X3']>0.5) != (Xmat[,'X4']>0.5))*1 )  
      
      VNL_Tmat <- cbind( 1, exp(VXmat[,'X1']), VXmat[,'X1']*VXmat[,'X2'], 
                         ((VXmat[,'X3']>0.5) != (VXmat[,'X4']>0.5))*1 )  
      
      for (s in setting) {
        print(sprintf("Setting %d - %s", s, d))
        # dataset 만들기
        # pi(x)
        # s %% 2 == 1, pi(x) is linear/ 0, pi(x) is non-linear
        if (s %% 2 == 1) {
          # pi(x) is Linear
          ps_beta <- c(1, -1,-1) # b0, b1, b3, b5
          ps_X <- c('Intercept', 'X1', 'X3')
          
          prob_A = c(exp(Xmat2[,ps_X] %*% ps_beta ) / (exp(Xmat2[,ps_X] %*% ps_beta ) + 1))
          A = rbinom(nn, 1, prob=prob_A)
          
          Vprob_A = exp(VXmat2[,ps_X] %*% ps_beta ) / (exp(VXmat2[,ps_X] %*% ps_beta ) + 1)
          VA = rbinom(Vnn, 1, prob=Vprob_A)
          VT_p_a <- c( VA*Vprob_A + (1-VA)*(1-Vprob_A) )
          
        } else if (s %% 2 == 0) {
          # pi(x) is Non-Linear
          ps_beta <- c(1, -1,1)
          
          prob_A = c( exp(NL_Xmat2 %*% ps_beta ) / (exp(NL_Xmat2 %*% ps_beta ) + 1) )
          A = rbinom(nn, 1, prob=prob_A)
          
          Vprob_A = exp(VNL_Xmat2 %*% ps_beta ) / (exp(VNL_Xmat2 %*% ps_beta ) + 1)
          VA = rbinom(Vnn, 1, prob=Vprob_A)
          VT_p_a <- c( VA*Vprob_A + (1-VA)*(1-Vprob_A) )
          
        }
        
        # tau(x)
        if (s %in% c(1, 2)) {
          # tau(x) is Linear, mu(x) is Linear
          Q0 <- 0.5 * c(1, 1, 1,-1, -1)  # b1, b2, b5, b6
          Q0_X <- c('Intercept', 'X1', 'X2', 'X5', 'X6')
          
          T0 <- c(1, -1, 1, -1) * alpha  # b1, b2, b3, b4
          T0_X <- c('AX1', 'AX2', 'AX3', 'AX4')
          
          data <- cbind(Xmat2[,1:11], Tmat[,1:11]*A)
          d1 <- cbind(Xmat, Tmat*A)
          
          mean_y = data[,c(Q0_X, T0_X)] %*% c(Q0, T0)
          Y = rnorm(nn, mean=mean_y, sd=1) # mean(mean_y)
          mu = data[,Q0_X] %*% Q0
          tau = Tmat[,T0_X] %*% T0
          
          Vdata <- cbind(VXmat2[,1:11], VTmat[,1:11]*VA)
          Vd1 <- cbind(VXmat, VTmat* VA)
          Vmean_y = Vdata[,c(Q0_X, T0_X)] %*% c(Q0, T0)
          VY = rnorm(Vnn, mean=Vmean_y, sd=1)
          T_tau <- c(VTmat[,T0_X] %*% T0)
          
        } else if (s %in% c(3, 4)){
          # tau(x) is Linear, mu(x) is Non Linear
          Q0 <- 0.5 * c( 2, -1)  # b1, b2, b5, b6
          T0 <- c(1, -1, 1, -1) * alpha  # b1, b2, b3, b4
          T0_X <- c('AX1', 'AX2', 'AX3', 'AX4')
          
          data <- cbind(NL_Xmat3, Tmat[,T0_X]*A)
          d1 <- cbind(Xmat, Tmat*A)
          
          mean_y = data %*% c(Q0, T0)
          Y = rnorm(nn, mean=mean_y, sd=1) # mean(mean_y)
          mu = NL_Xmat3 %*% Q0
          tau = Tmat[,T0_X] %*% T0
          
          Vdata <- cbind(VNL_Xmat3, VTmat[,T0_X]*VA)
          Vmean_y = Vdata %*% c(Q0, T0)
          VY = rnorm(Vnn, mean=Vmean_y, sd=1)
          T_tau <- c(VTmat[,T0_X] %*% T0)
          
        } else if (s %in% c(5, 6)){
          # tau(x) is Non Linear, mu(x) is Linear
          Q0 <- 0.5 * c(1, 1, 1,-1, -1)  # b1, b2, b5, b6
          Q0_X <- c('Intercept', 'X1', 'X2', 'X5', 'X6')
          T0 <- c(-2.07, 1, -1, 1) * alpha
          
          data <- cbind(Xmat2[,Q0_X], NL_Tmat*A)
          d1 <- cbind(Xmat, Tmat*A)
          
          mean_y = data %*% c(Q0, T0)
          Y = rnorm(nn, mean=mean_y, sd=1) # mean(mean_y)
          mu = data[,Q0_X] %*% Q0
          tau = NL_Tmat %*% T0
          
          Vdata <- cbind(VXmat2[,Q0_X], VNL_Tmat*VA)
          Vmean_y = Vdata %*% c(Q0, T0) 
          VY = rnorm(Vnn, mean=Vmean_y, sd=1)
          
          T_tau <- c(VNL_Tmat %*% T0)
          
        } else if (s %in% c(7, 8)){
          # tau(x) is Non Linear, mu(x) is Non Linear
          Q0 <- 0.5 * c( 2, -1)  # b1, b2, b5, b6
          T0 <- c(-2.07, 1, -1, 1) * alpha
          
          data <- cbind(NL_Xmat3, NL_Tmat*A)
          d1 <- cbind(Xmat, Tmat*A)
          
          mean_y = data %*% c(Q0, T0)
          Y = rnorm(nn, mean=mean_y, sd=1) # mean(mean_y) = 0.035
          mu = NL_Xmat3 %*% Q0
          tau = NL_Tmat %*% T0
          
          Vdata <- cbind(VNL_Xmat3, VNL_Tmat*VA)
          Vmean_y = Vdata %*% c(Q0, T0) 
          VY = rnorm(Vnn, mean=Vmean_y, sd=1)
          
          T_tau <- c(VNL_Tmat %*% T0)
        }
        
        muY = prob_A * (mu + tau) + (1-prob_A) * mu
        
        i_perm <-  sample.int(nn)
        
        T_d <- (sign(T_tau) + 1)/2
        T_d[T_d==0.5] <- 0
        T_ITR <-  (VA == T_d)
        T_MR <- sum( ((VY*T_ITR)/VT_p_a) ) / sum( (T_ITR / VT_p_a)  )
        
        QResults <- QResults2 <- AResults <- AResults2 <- NULL
        for (k in 1:numFolds) {
          ind_start <- floor((k-1)/numFolds*nn)+1
          ind_end <- floor(k/numFolds*nn)
          
          ind_val <- i_perm[ind_start:ind_end]
          ind_tr <- setdiff(1:nn, ind_val)
          
          X_valid <- d1[ind_val,]; X_train <- d1[ind_tr,]
          Y_valid <- Y[ind_val]; Y_train <- Y[ind_tr]
          A_valid <- X_valid[,'A']; A_train = X_train[,'A']
          
          mu_valid <- mu[ind_val]; mu_train <- mu[ind_tr]
          muY_valid <- muY[ind_val]; muY_train <- muY[ind_tr]
          pA_valid = prob_A[ind_val]; pA_train <- prob_A[ind_tr]
          
          ######################## Q-Learning ########################
          QLearn <- glmnet(x=X_train, y=Y_train, family='gaussian', lambda=lambda)
          
          prob_a = prob_A[ind_val]
          Qtau <- as.matrix(X_valid[,c('A', paste0('AX', 1:pp))] %*% coef(QLearn)[c('A', paste0('AX', 1:pp)), ])
          
          Qd <- (sign(Qtau)+1)/2
          Qd[Qd==0.5] <- 0
          Q_ITR <- (A_valid == Qd) # I_ITR must be boolean
          
          p_a = A_valid * pA_valid + (1-A_valid) * (1-pA_valid)
          
          QMR <- apply(Q_ITR, 2, function(t){
            return( sum( (Y_valid*t)/p_a ) /  sum( (t/p_a) ) )
          } )
          
          Qtmp <- data.frame(k_fold=k, lambda=QLearn$lambda, MeanReward=QMR)
          rownames(Qtmp) <- NULL
          QResults <- rbind(QResults, Qtmp)
          
          ######################## Q-Oracle ########################
          Qoracle <- glmnet(x=X_train, y=(Y_train-mu_train), family='gaussian', lambda=lambda)
          
          Qtau2 <- as.matrix(X_valid[,c('A', paste0('AX', 1:pp))] %*% coef(Qoracle)[c('A', paste0('AX', 1:pp)), ])
          
          Qd2 <- (sign(Qtau2)+1)/2
          Qd2[Qd2==0.5] <- 0
          Q_ITR2 <- (pA_valid == Qd2) # I_ITR must be boolean
          
          QMR2 <- apply(Q_ITR, 2, function(t){
            return( sum( (Y_valid*t)/pA_valid ) /  sum( (t/pA_valid) ) )
          } )
          
          Qtmp2 <- data.frame(k_fold=k, lambda=QLearn$lambda, MeanReward=QMR2)
          rownames(Qtmp2) <- NULL
          QResults2 <- rbind(QResults2, Qtmp2)
          
          ######################## A-Learning ########################
          # estimate My(Xi), Mt(Xi)
          tr_n = length(ind_tr); ind_halves <-  sample.int(tr_n)
          ind_1 <- ind_halves[1:(tr_n/2)]; ind_2 <- ind_halves[(tr_n/2+1):tr_n]
          
          X_train1 <- X_train[ind_1, cov]; X_train2 <- X_train[ind_2, cov]
          Y_train1 <- Y_train[ind_1]; Y_train2 <- Y_train[ind_2]
          A_train1 <- A_train[ind_1]; A_train2 <- A_train[ind_2]
          
          fit_Y1 <- randomForest(Y_train1 ~ ., data=X_train1)
          fit_A1 <- randomForest(factor(A_train1) ~ ., data=X_train1)
          
          fit_Y2 <- randomForest(Y_train2 ~ ., data=X_train2)
          fit_A2 <- randomForest(factor(A_train2) ~ ., data=X_train2)
          
          # predict all obs.
          pred_Y2 <- predict(fit_Y1, X_train2)
          pred_A2 <- predict(fit_A1, X_train2, type="prob")[,2]
          
          pred_Y1 <- predict(fit_Y2, X_train1)
          pred_A1 <- predict(fit_A2, X_train1, type="prob")[,2]
          
          pred_Y <- rep(0, tr_n); pred_A <- rep(0, tr_n)
          pred_Y[ind_1] <- pred_Y1; pred_Y[ind_2] <- pred_Y2
          pred_A[ind_1] <- pred_A1; pred_A[ind_2] <- pred_A2
          
          Awt <- (A_train-pred_A);
          Y2 <- Y_train-pred_Y;
          X_train2 <- cbind(1, X_train[, cov]) * Awt
          
          ALearn <- glmnet(x=X_train2, y=Y2, family = "gaussian", lambda=lambda, intercept = FALSE)
          
          X_valid2 <- cbind(1, X_valid[,cov])
          Atau <-  predict(ALearn, X_valid2)
          
          Ad <- (sign(Atau) + 1) /2 
          Ad[Ad==0.5] <- 0
          A_ITR <- (A_valid == Ad)
          
          AMR <- apply(A_ITR, 2, function(t){
            return( sum( (Y_valid*t)/pA_valid)  /  sum( (t/pA_valid) ) )
          } )
          
          Atmp <- data.frame(k_fold=k, lambda=ALearn$lambda, MeanReward=AMR)
          rownames(Atmp) <- NULL
          AResults <- rbind(AResults, Atmp)
          
          ######################## A-Oracle ########################
          # estimate My(Xi), Mt(Xi)
          Awt <- (A_train-pA_train);
          Y2 <- Y_train-muY_train;
          X_train2 <- cbind(1, X_train[, cov]) * Awt
          
          Aoracle <- glmnet(x=X_train2, y=Y2, family = "gaussian", lambda=lambda, intercept = FALSE)
          
          X_valid2 <- cbind(1, X_valid[,cov])
          Atau2 <-  predict(Aoracle, X_valid2)
          
          Ad2 <- (sign(Atau2) + 1) /2 
          Ad2[Ad2==0.5] <- 0
          A_ITR2 <- (A_valid == Ad2)
          
          AMR2 <- apply(A_ITR, 2, function(t){
            return( sum( (Y_valid*t)/pA_valid)  /  sum( (t/pA_valid) ) )
          } )
          
          Atmp2 <- data.frame(k_fold=k, lambda=Aoracle$lambda, MeanReward=AMR2)
          rownames(Atmp2) <- NULL
          AResults2 <- rbind(AResults2, Atmp2)
          
        }
        
        # Q Learning
        QResults_Final <- aggregate(MeanReward~lambda, QResults, "mean")
        
        Qbest_i <- order(-QResults_Final$MeanReward)[1]
        Qbest_lambda <- QResults_Final$lambda[order(-QResults_Final$MeanReward, -QResults_Final$lambda)[1]]
        Qbest_MR <- QResults_Final$MeanReward[order(-QResults_Final$MeanReward, -QResults_Final$lambda)[1]]
        
        Final_QLearn <-  glmnet(d1, Y,  family='gaussian', lambda=Qbest_lambda)
        Qbest_beta <- coef(Final_QLearn)[,'s0']
        
        Final_Qtau <- VTmat[,c('A', paste0('AX', 1:pp))] %*% Qbest_beta[c('A', paste0('AX', 1:pp))]
        Final_Qd <- ( sign( Final_Qtau )+1 )/2
        Final_Qd[Final_Qd==0.5] <- 0
        Final_QITR <- (VA == Final_Qd)
        Final_QMR <- sum(  ((VY*Final_QITR)/VT_p_a) ) /  sum( (Final_QITR/VT_p_a) ) 
        Final_QVR <- Final_QMR / T_MR 
        
        Qtmp2 <- data.frame('itr'=r, 'setting'=s, 'dist'=d, 'alpha'=alpha, 'pp'=pp, 'nn'=nn, 'method'='Q-Learning', 
                            'true_MR'=T_MR, 'true_tau'=mean(T_tau), 'best_MR'=Final_QMR, 'best_VR'=Final_QVR, 'lambda'=Qbest_lambda, 'Tau_MSE'=mean((T_tau-Final_Qtau)^2)) # Loss=MSE(R_valid, preds[,i]), 
        
        t <- (Qbest_beta[Qprint] != eff.mod)  
        print(Qbest_beta[Qprint])
        t2 <- c(t, 0, 0, 1, 1)
        eff.mod2 <- c(eff.mod, 0, 1, 0, 1)
        cm <- (table('actual'=1*eff.mod2, 'predicted'=1*t2)-1)
        TN <- cm[1,1]; FP <- cm[1,2]; FN <- cm[2,1]; TP <- cm[2,2]
        Sen <- TP/(TP+FN); Spe <- TN/(TN+FP);
        MCC = (TP*TN - FP*FN) / sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) ;
        Qtmp2 <- cbind(Qtmp2, t(Qbest_beta[Qprint]), 'Sensitivity'=Sen, 'Specificity'=Spe, 'MCC'=MCC )
        print(MCC)
        QFinal <- rbind(QFinal, Qtmp2)
        rm(list=c('Final_Qtau', 'Final_Qd', 'Final_QITR', 'Final_QMR', 'Final_QVR', 'Sen', 'Spe', 'MCC'))
        
        # A Learning
        AResults_Final <- aggregate(MeanReward~lambda, AResults, "mean")
        
        Abest_i <- order(-AResults_Final$MeanReward)[1]
        Abest_lambda <- AResults_Final$lambda[order(-AResults_Final$MeanReward, -AResults_Final$lambda)[1]]
        Abest_MR <- AResults_Final$MeanReward[order(-AResults_Final$MeanReward, -AResults_Final$lambda)[1]]
        
        tr_n = nrow(Xmat); ind_halves <-  sample.int(tr_n)
        ind_1 <- ind_halves[1:(tr_n/2)]; ind_2 <- ind_halves[(tr_n/2+1):tr_n]
        
        X1 <- Xmat[ind_1,]; X2 <- Xmat[ind_2,]
        Y1 <- Y[ind_1]; Y2 <- Y[ind_2]
        A1 <- A[ind_1]; A2 <- A[ind_2]
        
        fit_Y1 <- randomForest(Y1 ~ ., data=X1)
        fit_A1 <- randomForest(factor(A1) ~ ., data=X1)
        
        fit_Y2 <- randomForest(Y2 ~ ., data=X2)
        fit_A2 <- randomForest(factor(A2) ~ ., data=X2)
        
        # predict all obs.
        pred_Y2 <- predict(fit_Y1, X2)
        pred_A2 <- predict(fit_A1, X2, type="prob")[,2]
        
        pred_Y1 <- predict(fit_Y2, X1)
        pred_A1 <- predict(fit_A2, X1, type="prob")[,2]
        
        pred_Y <- rep(0, tr_n); pred_A <- rep(0, tr_n)
        pred_Y[ind_1] <- pred_Y1; pred_Y[ind_2] <- pred_Y2
        pred_A[ind_1] <- pred_A1; pred_A[ind_2] <- pred_A2
        
        ###
        Awt <- (A-pred_A);
        Y2 <- Y-pred_Y;
        X3 <- cbind(1, Xmat) * Awt
        
        Final_ALearn <- glmnet(X3, Y2, family='gaussian', lambda=Abest_lambda, intercept = FALSE)
        
        Final_Atau <-  predict(Final_ALearn, VTmat)
        Final_Ad <- (sign(Final_Atau)+1)/2
        Final_Ad[Final_Ad==0.5] <- 0
        Final_AITR <- c(VA == Final_Ad)
        Final_AMR <- sum( ( (VY*Final_AITR)/VT_p_a) ) /sum( ( Final_AITR/VT_p_a) )
        Final_AVR <- Final_AMR / T_MR
        
        Atmp2 <- data.frame('itr'=r, 'setting'=s, 'dist'=d, 'alpha'=alpha, 'pp'=pp, 'nn'=nn, 'method'='A-Learning', 
                            'true_MR'=T_MR, 'true_tau'=mean(T_tau), 'best_MR'=Final_AMR, 'best_VR'=Final_AVR, 'lambda'=Abest_lambda, 'Tau_MSE'=mean((T_tau-Final_Atau)^2)) # Loss=MSE(R_valid, preds[,i]), 
        
        Abest_beta <- coef(Final_ALearn)[,'s0']
        names(Abest_beta) <- c('zero', 'A', Qprint)
        
        t <- (Abest_beta[Qprint] != eff.mod)  
        t2 <- c(t, 0, 0, 1, 1)
        eff.mod2 <- c(eff.mod, 0, 1, 0, 1)
        cm <- (table('actual'=1*eff.mod2, 'predicted'=1*t2)-1)
        TN <- cm[1,1]; FP <- cm[1,2]; FN <- cm[2,1]; TP <- cm[2,2]
        Sen <- TP/(TP+FN); Spe <- TN/(TN+FP);
        MCC = (TP*TN - FP*FN) / sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) ;
        
        Atmp2 <- cbind(Atmp2, t( Abest_beta[Qprint] ) , 'Sensitivity'=Sen, 'Specificity'=Spe, 'MCC'=MCC )
        colnames(Atmp2) <- c(colnames(Atmp2))
        AFinal <- rbind(AFinal, Atmp2)
        
        rm(list=c('Final_Atau', 'Final_Ad', 'Final_AITR', 'Final_AMR', 'Final_AVR', 'Sen', 'Spe', 'MCC'))
        
        # Q Oracle
        QResults_Final2 <- aggregate(MeanReward~lambda, QResults2, "mean")
        
        Qbest_i2 <- order(-QResults_Final2$MeanReward)[1]
        Qbest_lambda2 <- QResults_Final2$lambda[order(-QResults_Final2$MeanReward, -QResults_Final2$lambda)[1]]
        Qbest_MR2 <- QResults_Final2$MeanReward[order(-QResults_Final2$MeanReward, -QResults_Final2$lambda)[1]]
        
        Final_QLearn2 <-  glmnet(d1, Y-mu,  family='gaussian', lambda=Qbest_lambda2)
        Qbest_beta2 <- coef(Final_QLearn2)[,'s0']
        
        Final_Qtau2 <- VTmat[,c('A', paste0('AX', 1:pp))] %*% Qbest_beta2[c('A', paste0('AX', 1:pp))]
        Final_Qd2 <- ( sign( Final_Qtau2 )+1 )/2
        Final_Qd2[Final_Qd2==0.5] <- 0
        Final_QITR2 <- (VA == Final_Qd2)
        Final_QMR2 <- sum(  ((VY*Final_QITR2)/VT_p_a) ) /  sum( (Final_QITR2/VT_p_a) ) 
        Final_QVR2 <- Final_QMR2 / T_MR 
        
        Qtmp2 <- data.frame('itr'=r, 'setting'=s, 'dist'=d, 'alpha'=alpha, 'pp'=pp, 'nn'=nn, 'method'='Q-Oracle', 
                            'true_MR'=T_MR, 'true_tau'=mean(T_tau), 'best_MR'=Final_QMR2, 'best_VR'=Final_QVR2, 'lambda'=Qbest_lambda2, 'Tau_MSE'=mean((T_tau-Final_Qtau2)^2))
        
        t <- (Qbest_beta2[Qprint] != eff.mod)  
        print(Qbest_beta2[Qprint])
        t2 <- c(t, 0, 0, 1, 1)
        eff.mod2 <- c(eff.mod, 0, 1, 0, 1)
        cm <- (table('actual'=1*eff.mod2, 'predicted'=1*t2)-1)
        TN <- cm[1,1]; FP <- cm[1,2]; FN <- cm[2,1]; TP <- cm[2,2]
        Sen <- TP/(TP+FN); Spe <- TN/(TN+FP);
        MCC = (TP*TN - FP*FN) / sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) ;
        Qtmp2 <- cbind(Qtmp2, t(Qbest_beta2[Qprint]), 'Sensitivity'=Sen, 'Specificity'=Spe, 'MCC'=MCC )
        print(MCC)
        QFinal2 <- rbind(QFinal2, Qtmp2)
        rm(list=c('Final_Qtau2', 'Final_Qd2', 'Final_QITR2', 'Final_QMR2', 'Final_QVR2', 'Sen', 'Spe', 'MCC'))
        
        # A Oracle
        AResults_Final2 <- aggregate(MeanReward~lambda, AResults2, "mean")
        
        Abest_i2 <- order(-AResults_Final2$MeanReward)[1]
        Abest_lambda2 <- AResults_Final2$lambda[order(-AResults_Final2$MeanReward, -AResults_Final2$lambda)[1]]
        Abest_MR2 <- AResults_Final2$MeanReward[order(-AResults_Final2$MeanReward, -AResults_Final2$lambda)[1]]
        
        Awt <- (A-prob_A);
        Y2 <- Y-muY;
        X3 <- cbind(1, Xmat) * Awt
        
        Final_ALearn2 <- glmnet(X3, Y2, family='gaussian', lambda=Abest_lambda2, intercept = FALSE)
        
        Final_Atau2 <-  predict(Final_ALearn2, VTmat)
        Final_Ad2 <- (sign(Final_Atau2)+1)/2
        Final_Ad2[Final_Ad2==0.5] <- 0
        Final_AITR2 <- c(VA == Final_Ad2)
        Final_AMR2 <- sum( ( (VY*Final_AITR2)/VT_p_a) )/sum( ( Final_AITR2/VT_p_a) )
        Final_AVR2 <- Final_AMR2 / T_MR
        
        Atmp2 <- data.frame('itr'=r, 'setting'=s, 'dist'=d, 'alpha'=alpha, 'pp'=pp, 'nn'=nn, 'method'='A-Oracle', 
                            'true_MR'=T_MR, 'true_tau'=mean(T_tau), 'best_MR'=Final_AMR2, 'best_VR'=Final_AVR2, 'lambda'=Abest_lambda2, 'Tau_MSE'=mean((T_tau-Final_Atau2)^2))
        
        Abest_beta2 <- coef(Final_ALearn2)[,'s0']
        names(Abest_beta2) <- c('zero', 'A', Qprint)
        
        t <- (Abest_beta2[Qprint] != eff.mod)  
        t2 <- c(t, 0, 0, 1, 1)
        eff.mod2 <- c(eff.mod, 0, 1, 0, 1)
        cm <- (table('actual'=1*eff.mod2, 'predicted'=1*t2)-1)
        TN <- cm[1,1]; FP <- cm[1,2]; FN <- cm[2,1]; TP <- cm[2,2]
        Sen <- TP/(TP+FN); Spe <- TN/(TN+FP);
        MCC = (TP*TN - FP*FN) / sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) ;
        
        Atmp2 <- cbind(Atmp2, t( Abest_beta2[Qprint] ) , 'Sensitivity'=Sen, 'Specificity'=Spe, 'MCC'=MCC )
        
        AFinal2 <- rbind(AFinal2, Atmp2)
        rm(list=c('Final_Atau2', 'Final_Ad2', 'Final_AITR2', 'Final_AMR2', 'Final_AVR2', 'Sen', 'Spe', 'MCC'))
      }
      
      
    }
  }
  
  
  Final <- rbind(Final, QFinal, QFinal2, AFinal, AFinal2)
  
  vars <- ls()
  vars2 <- vars[!vars %in% c('Rep', 
                             'alpha', 'lambda', 'numFolds', 'N',
                             'setting', 'dist', 'pp', 'nn', 'qq', 'alpha', 'Vnn',
                             'Qprint', 'Aprint', 'Final', 'eff.mod')]
  rm(list=vars2)
  write.csv(Final, './simulation/results2/12-31-10.csv', row.names=FALSE)
}

sum_results <- Final %>% group_by() %>% group_by(setting, dist, alpha, pp, nn, method) %>% summarise('true_MR'=mean(true_MR),
                                                                                      'true_tau'=mean(true_tau),
                                                                                      'best_MR'=mean(best_MR),
                                                                                      'best_VR'=mean(best_VR),
                                                                                      'tau_MSE'=mean(Tau_MSE),
                                                                                      'Sensitivity'=mean(Sensitivity),
                                                                                      'Specificity'=mean(Specificity),
                                                                                      'MCC'=mean(MCC)) %>% data.frame()

sum_results$best_VR[sum_results$best_VR > 1] <- 1
sum_results$MCC[is.na(sum_results$MCC)] <- 0

head(sum_results)

write.csv(sum_results, './simulation/results2/22-01-04-10.csv')

