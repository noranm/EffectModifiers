# Simulation Setting 1

# simulation
nn=100; pp=10;
qq=0.3
cov <- paste0('X', 1:pp)

lambda <-  2^seq(from=1, to=-10, length=100)
numFolds <- 5

id <- 1:nn

# True Coef #
ps_beta <- c(1,-1,0,1,0,-1,0,0,0,0,0)
#ps_beta <- c(0,0,0,0,0,0,0,0,0,0,0)
# mu
Q0 <- 0.3 * c(0, 1, 1, 0, 0, -1, -1, 0, 0, 0, 0) 
a = 1
# tau
T0 <- c(0, 1, -1, -1, 1, 0, 0, 0, 0, 0, 0) * a # a 약물은 처치 효과 Good / 하지만 x1, x2 의 하위군에는 효과 Bad

Rep=100
QFinal <- AFinal <- NULL
r=1
for (r in 1:Rep){
  if (r%%100==0) print(r)
  set.seed(r)
  # (1) data generation
  Xmat <- matrix(rbinom(nn*pp, 1, prob=qq), nrow=nn, ncol=pp)
  colnames(Xmat) <- paste0('X', 1:10)
  
  Xmat2 <- cbind(rep(1,nn), Xmat)
  colnames(Xmat2) <- c("Intercept", paste0('X', 1:10))
  cov <- paste0('X', 1:pp)
  
  prob_A = exp(Xmat2 %*% ps_beta ) / (exp(Xmat2 %*% ps_beta ) + 1)
  A = rbinom(nn, 1, prob=prob_A)
  
  Tmat <- cbind(A, Xmat * A)
  colnames(Tmat) <- c("A", paste0('AX', 1:10))
  data <- cbind(Xmat2, Tmat)
  d1 <- cbind(Xmat, Tmat)
  
  mean_y = data %*% c(Q0, T0) 
  Y = rnorm(nn, mean=mean_y, sd=1)
  i_perm <-  sample.int(nn)
  
  QResults <- QResults2 <- AResults <- AResults2 <- NULL
  
  # k=1
  # ps_model <- glm(A~Xmat2[,2:11], family='binomial')
  # ps <- ps_model$fitted.values
  for (k in 1:numFolds) {
    ind_start <- floor((k-1)/numFolds*nn)+1
    ind_end <- floor(k/numFolds*nn)
    
    ind_val <- i_perm[ind_start:ind_end]
    ind_tr <- setdiff(1:nn, ind_val)
    
    X_valid <- d1[ind_val,]; X_train <- d1[ind_tr,]
    Y_valid <- Y[ind_val]; Y_train <- Y[ind_tr]
    
    ######################## Q-Learning ########################
    QLearn <- glmnet(x=X_train, y=Y_train, family='gaussian', lambda=lambda)
    
    a = X_valid[,'A']
    # 나중에 validation set에서 추정 필요
    prob_a = prob_A[ind_val]
    Qtau <- as.matrix(X_valid[,11:21] %*% coef(QLearn)[c('A', 'AX1', 'AX2', 'AX3', 'AX4', 'AX5', 'AX6', 'AX7', 'AX8', 'AX9', 'AX10'), ])
    Qd <- (sign(Qtau)+1)/2
    Qd[Qd==0.5] <- 0
    Q_ITR <- (a == Qd) # I_ITR must be boolean
    
    # check point
    # dim(X_valid[,11:21]); dim(coef(QLearn)[c('A', 'AX1', 'AX2', 'AX3', 'AX4', 'AX5', 'AX6', 'AX7', 'AX8', 'AX9', 'AX10'), ])
    # dim(Qd); dim(Q_ITR);
    
    p_a = a * prob_a + (1-a) * (1-prob_a)
    # check
    # dim(Y_valid*Q_ITR); dim(p_a * Q_ITR) ; dim(p_a * Q_ITR)
    
    #  print(t)
    #  print(length(t))
    MR <- apply(Q_ITR, 2, function(t){
      return( sum( ( Y_valid[t] / p_a[t] ) /  ((1/p_a[t]) ) )/ (nn/5) )
    } )
    if (length(QLearn$lambda) != 100) {
      print("????")
      break
    }
    if (length(MR) != 100) {
      print("????")
      break
    }
    Qtmp <- data.frame(k_fold=k, lambda=QLearn$lambda, MeanReward=MR)
    rownames(Qtmp) <- NULL
    QResults <- rbind(QResults, Qtmp)
    
    ######################## A-Learning ########################
    # estimate My(Xi), Mt(Xi)
    tr_n = length(ind_tr); ind_halves <-  sample.int(tr_n)
    ind_1 <- ind_halves[1:(tr_n/2)]; ind_2 <- ind_halves[(tr_n/2+1):tr_n]
    
    X_train1 <- X_train[ind_1, cov]; X_train2 <- X_train[ind_2, cov]
    Y_train1 <- Y_train[ind_1]; Y_train2 <- Y_train[ind_2]
    
    A_train <- X_train[,'A']; A_valid <- X_valid[,'A']
    A_train1 <- A_train[ind_1]; A_train2 <- A_train[ind_2]
    
    fit_Y1 <- randomForest(factor(Y_train1) ~ ., data=X_train1)
    fit_A1 <- randomForest(factor(A_train1) ~ ., data=X_train1)
    
    fit_Y2 <- randomForest(factor(Y_train2) ~ ., data=X_train2)
    fit_A2 <- randomForest(factor(A_train2) ~ ., data=X_train2)
    
    # predict all obs.
    pred_Y2 <- predict(fit_Y1, X_train2, type="prob")[,2]
    pred_A2 <- predict(fit_A1, X_train2, type="prob")[,2]
    
    pred_Y1 <- predict(fit_Y2, X_train1, type="prob")[,2]
    pred_A1 <- predict(fit_A2, X_train1, type="prob")[,2]
    
    pred_Y <- rep(0, tr_n); pred_A <- rep(0, tr_n)
    pred_Y[ind_1] <- pred_Y1; pred_Y[ind_2] <- pred_Y2
    pred_A[ind_1] <- pred_A1; pred_A[ind_2] <- pred_A2
    
    Awt <- (A_train-pred_A);
    Y2 <- Y_train-pred_Y;
    X_train2 <- cbind(1, X_train) * Awt
    
    ALearn <- glmnet(x=X_train2, y=Y2, family = "gaussian", lambda=lambda, intercept = FALSE)
    
    X_valid2 <- cbind(1, X_valid)
    Atau <-  predict(ALearn, X_valid2)
    
    Ad <- (sign(Atau) + 1) /2 
    Ad[Ad==0.5] <- 0
    A_ITR <- (A_valid == Ad)
    MR <- apply(A_ITR, 2, function(t){
      return( sum( ( Y_valid[t] / p_a[t] ) /  ((1/p_a[t]) ) )/ (nn/5) )
    } )
    
    if (length(ALearn$lambda) != 100) {
      print("????")
      break
    }
    if (length(MR) != 100) {
      print("????")
      break
    }
    
    Atmp <- data.frame(k_fold=k, lambda=ALearn$lambda, MeanReward=MR)
    rownames(Atmp) <- NULL
    AResults <- rbind(AResults, Atmp)
  }
  
  # (왜?) dim 안맞음 100개가 아닌디...?
  QResults2 <- aggregate(MeanReward~lambda, QResults, "mean")
  AResults2 <- aggregate(MeanReward~lambda, AResults, "mean")
  
  Qbest_i <- order(QResults2$MeanReward,decreasing=TRUE)[1]
  Qbest_lambda <- QResults2$lambda[order(QResults2$MeanReward,decreasing=TRUE)[1]]
  Qbest_MR <- QResults2$MeanReward[order(QResults2$MeanReward,decreasing=TRUE)[1]]
  
  Abest_i <- order(AResults2$MeanReward,decreasing=TRUE)[1]
  Abest_lambda <- AResults2$lambda[order(AResults2$MeanReward,decreasing=TRUE)[1]]
  Abest_MR <- AResults2$MeanReward[order(AResults2$MeanReward,decreasing=TRUE)[1]]
  
  # Q 정리 
  Final_QLearn <-  glmnet(d1, Y,  family='gaussian', lambda=Qbest_lambda)
  
  Qtmp2 <- data.frame('itr'=r, Qbest_lambda, 'best_MR'=Qbest_MR) # Loss=MSE(R_valid, preds[,i]), 
  Qbest_beta <- coef(Final_QLearn)[,'s0']
  Qtmp2 <- cbind(Qtmp2, t(Qbest_beta) )
  QFinal <- rbind(QFinal, Qtmp2)
  
  # A 정리
  tr_n = nrow(Xmat); ind_halves <-  sample.int(tr_n)
  ind_1 <- ind_halves[1:(tr_n/2)]; ind_2 <- ind_halves[(tr_n/2+1):tr_n]
  
  X1 <- Xmat[ind_1,]; X2 <- Xmat[ind_2,]
  Y1 <- Y[ind_1]; Y2 <- Y[ind_2]
  A1 <- A[ind_1]; A2 <- A[ind_2]
  
  fit_Y1 <- randomForest(factor(Y1) ~ ., data=X1)
  fit_A1 <- randomForest(factor(A1) ~ ., data=X1)
  
  fit_Y2 <- randomForest(factor(Y2) ~ ., data=X2)
  fit_A2 <- randomForest(factor(A2) ~ ., data=X2)
  
  # predict all obs.
  pred_Y2 <- predict(fit_Y1, X2, type="prob")[,2]
  pred_A2 <- predict(fit_A1, X2, type="prob")[,2]
  
  pred_Y1 <- predict(fit_Y2, X1, type="prob")[,2]
  pred_A1 <- predict(fit_A2, X1, type="prob")[,2]
  
  pred_Y <- rep(0, tr_n); pred_A <- rep(0, tr_n)
  pred_Y[ind_1] <- pred_Y1; pred_Y[ind_2] <- pred_Y2
  pred_A[ind_1] <- pred_A1; pred_A[ind_2] <- pred_A2
  
  ###
  Awt <- (A-pred_A);
  Y2 <- Y-pred_Y;
  X3 <- cbind(1, Xmat) * Awt
  
  Final_ALearn <- glmnet(X3, Y2, family='gaussian', lambda=Abest_lambda, intercept = FALSE)
  Atmp2 <- data.frame('itr'=r, Abest_lambda, 'best_MR'=Abest_MR) # Loss=MSE(R_valid, preds[,i]), 
  Abest_beta <- coef(Final_ALearn)[,'s0']
  Atmp2 <- cbind(Atmp2, t( Abest_beta ) )
  
  AFinal <- rbind(AFinal, Atmp2)
}

print(mean(Final$best_MR))
print(median(Final$best_MR))
print(colSums(Final[,4:25])/1000)
