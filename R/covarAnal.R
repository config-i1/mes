mesVarAnal <- function(h, measurement, persistence, s2, Etype="M"){
    # This is currently developed only for pure multiplicative non-seasonal models
    QMatrix <- diag(as.vector(persistence),length(persistence),length(persistence));
    k <- length(persistence);
    varMat <- rep(1+s2, h);
    for(i in 2:h){
        varMat[i] <- varMat[i] * exp(measurement %*% log(matrixPowerWrap(diag(k)+matrixPowerWrap(QMatrix,2) * s2,i-1)
                                                         - diag(k)) %*% t(measurement));
    }
    varMat[1] <- varMat[1] - 1;
    varMat[-1] <- varMat[-1] + s2;

    return(varMat);
}

