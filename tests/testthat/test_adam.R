context("Tests for ADAM");

#### Basic ETS stuff ####
# Basic ADAM selection
testModel <- adam(Mcomp::M3[[1234]], "ZZZ");
test_that("Test ADAM(ZZZ) selection on N1234", {
    expect_match(modelType(testModel), "MAN");
})

# Basic ADAM selection on 2568
testModel <- adam(Mcomp::M3[[2568]], "ZZZ");
test_that("Test ADAM(ZZZ) selection on N2568", {
    expect_match(modelType(testModel), "MAM");
})

# Full ADAM selection
testModel <- adam(Mcomp::M3[[1234]], "FFF");
test_that("Test ADAM(FFF) selection on N1234", {
    expect_match(modelType(testModel), "MAN");
})

# ADAM with specified pool
testModel <- adam(Mcomp::M3[[1234]], c("AAA","ANN","MAN","MAM"));
test_that("Test ADAM selection with a pool on N1234", {
    expect_match(modelType(testModel), "MAN");
})

# Test ADAM forecasts with simulated interval
testForecast <- forecast(testModel,h=8,interval="sim",level=c(0.9,0.95));
test_that("Test ADAM forecast with simulated interval", {
    expect_equal(ncol(testForecast$lower), 2);
})

# ADAM combination
testModel <- adam(Mcomp::M3[[1234]], "CCC");
test_that("Test ADAM(CCC) on N1234", {
    expect_match(modelType(testModel), "CCC");
})

# Test ADAM forecasts with approximated interval
testForecast <- forecast(testModel,h=8,interval="app",level=c(0.9,0.95),side="upper");
test_that("Test ADAM forecast with simulated interval", {
    expect_equal(ncol(testForecast$lower), 2);
})


#### Advanced losses for ADAM ####
# ADAM with dalaplace
testModel <- adam(Mcomp::M3[[1234]], "MAN", distribution="dalaplace", lambda=0.05);
test_that("Test ADAM(MAN) with asymmetric Laplace on N1234", {
    expect_match(testModel$distribution, "dalaplace");
})

# ADAM with MSE
testModel <- adam(Mcomp::M3[[1234]], "MAN", loss="MSE");
test_that("Test ADAM(MAN) with MSE on N1234", {
    expect_match(testModel$loss, "MSE");
})

# ADAM with MSEh
testModel <- adam(Mcomp::M3[[1234]], "MAN", loss="MSEh");
test_that("Test ADAM(MAN) with MSE on N1234", {
    expect_match(testModel$loss, "MSEh");
})

# ADAM with GTMSE
testModel <- adam(Mcomp::M3[[1234]], "MAN", loss="GTMSE");
test_that("Test ADAM(MAN) with GTMSE on N1234", {
    expect_match(testModel$loss, "GTMSE");
})

# ADAM with GPL
testModel <- adam(Mcomp::M3[[1234]], "MAN", loss="GPL");
test_that("Test ADAM(MAN) with GPL on N1234", {
    expect_match(testModel$loss, "GPL");
})

# ADAM with LASSO
testModel <- adam(Mcomp::M3[[1234]], "MAN", loss="LASSO", lambda=0.5);
test_that("Test ADAM(MAN) with LASSO on N1234", {
    expect_match(testModel$loss, "LASSO");
})

# ADAM with custom loss function
loss <- function(actual, fitted, B){
    return(abs(actual-fitted)^3);
}
testModel <- adam(Mcomp::M3[[1234]], "AAN", loss=loss);
test_that("Test ADAM(AAN) with custom loss on N1234", {
    expect_match(testModel$loss, "custom");
})


#### ETS + occurrence model ####

#### ETS with several seasonalities ####

#### ETSX / Regression + formula ####

#### ETS + ARIMA / ARIMA ####

#### Provided initial / persistence / arma ####

#### auto.adam ####
