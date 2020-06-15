context("Tests for ADAM");

#### Basic ETS stuff ####
# Basic ADAM selection
testModel <- adam(Mcomp::M3[[1234]], "ZZZ");
test_that("ADAM ETS(ZZZ) selection on N1234", {
    expect_match(modelType(testModel), "MAN");
})

# Basic ADAM selection on 2568
testModel <- adam(Mcomp::M3[[2568]], "ZZZ");
test_that("ADAM ETS(ZZZ) selection on N2568", {
    expect_match(modelType(testModel), "MAM");
})

# Full ADAM selection
testModel <- adam(Mcomp::M3[[1234]], "FFF");
test_that("ADAM ETS(FFF) selection on N1234", {
    expect_match(modelType(testModel), "MAN");
})

# ADAM with specified pool
testModel <- adam(Mcomp::M3[[1234]], c("AAA","ANN","MAN","MAM"));
test_that("ADAM selection with a pool on N1234", {
    expect_match(modelType(testModel), "MAN");
})

# ADAM forecasts with simulated interval
testForecast <- forecast(testModel,h=8,interval="sim",level=c(0.9,0.95));
test_that("ADAM forecast with simulated interval", {
    expect_equal(ncol(testForecast$lower), 2);
})

# ADAM combination
testModel <- adam(Mcomp::M3[[1234]], "CCC");
test_that("ADAM ETS(CCC) on N1234", {
    expect_match(modelType(testModel), "CCC");
})

# ADAM forecasts with approximated interval
testForecast <- forecast(testModel,h=8,interval="app",level=c(0.9,0.95),side="upper");
test_that("ADAM forecast with simulated interval", {
    expect_equal(ncol(testForecast$lower), 2);
})


#### Advanced losses for ADAM ####
# ADAM with dalaplace
testModel <- adam(Mcomp::M3[[1234]], "MAN", distribution="dalaplace", lambda=0.05);
test_that("ADAM ETS(MAN) with asymmetric Laplace on N1234", {
    expect_match(testModel$distribution, "dalaplace");
})

# ADAM with MSE
testModel <- adam(Mcomp::M3[[1234]], "MAN", loss="MSE");
test_that("ADAM ETS(MAN) with MSE on N1234", {
    expect_match(testModel$loss, "MSE");
})

# ADAM with MSEh
testModel <- adam(Mcomp::M3[[1234]], "MAN", loss="MSEh");
test_that("ADAM ETS(MAN) with MSE on N1234", {
    expect_match(testModel$loss, "MSEh");
})

# ADAM with GTMSE
testModel <- adam(Mcomp::M3[[1234]], "MAN", loss="GTMSE");
test_that("ADAM ETS(MAN) with GTMSE on N1234", {
    expect_match(testModel$loss, "GTMSE");
})

# ADAM with GPL
testModel <- adam(Mcomp::M3[[1234]], "MAN", loss="GPL");
test_that("ADAM ETS(MAN) with GPL on N1234", {
    expect_match(testModel$loss, "GPL");
})

# ADAM with LASSO
testModel <- adam(Mcomp::M3[[1234]], "MAN", loss="LASSO", lambda=0.5);
test_that("ADAM ETS(MAN) with LASSO on N1234", {
    expect_match(testModel$loss, "LASSO");
})

# ADAM with custom loss function
loss <- function(actual, fitted, B){
    return(abs(actual-fitted)^3);
}
testModel <- adam(Mcomp::M3[[1234]], "AAN", loss=loss);
test_that("ADAM ETS(AAN) with custom loss on N1234", {
    expect_match(testModel$loss, "custom");
})


#### ETS + occurrence model ####
# Generate intermittent data
x <- sim.oes("MNN", 120, frequency=12, occurrence="general")
x <- sim.es("MNN", 120, frequency=12, probability=x$probability)

# iETS(M,N,N)_G
testModel <- adam(x$data, "MNN", occurrence="general")
test_that("ADAM iETS(MNN) with general occurrence", {
    expect_match(testModel$occurrence$occurrence, "general");
})

# iETS(M,M,M)_A
testModel <- adam(x$data, "MMM", occurrence="direct")
test_that("ADAM iETS(MMM) with direct occurrence", {
    expect_match(errorType(testModel), "M");
})

# iETS(M,M,N)_A
testModel <- adam(x$data, "MMN", occurrence="auto")
test_that("ADAM iETS(MMN) with auto occurrence", {
    expect_match(errorType(testModel), "M");
})

# iETS(Z,Z,N)_A
testModel <- adam(x$data, "ZZN", occurrence="auto")
test_that("ADAM iETS(MMN) with auto occurrence", {
    expect_true(is.occurrence(testModel$occurrence));
})

# Forecasts from the model
testForecast <- forecast(testModel, h=18, interval="semi")
test_that("Froecast from ADAM iETS(ZZZ)", {
    expect_true(is.adam(testForecast$model));
})


#### ETS with several seasonalities ####
# Double seasonality on N2568
testModel <- adam(Mcomp::M3[[2568]]$x, "YYY", lags=c(1,3,12), h=18);
test_that("ADAM ETS(YYY) with double seasonality on N2568", {
    expect_identical(testModel$lags, c(1,3,12));
})

# Double seasonality on N2568
testModel <- adam(Mcomp::M3[[2568]]$x, "FFF", lags=c(1,3,12), h=18, initial="backcasting");
test_that("ADAM ETS(FFF) + backcasting with double seasonality on N2568", {
    expect_identical(testModel$lags, c(1,3,12));
})

# Double seasonality on N2568
testModel <- adam(Mcomp::M3[[2568]]$x, "CCC", lags=c(1,3,12), h=18);
test_that("ADAM ETS(CCC) with double seasonality on N2568", {
    expect_identical(testModel$models[[1]]$lags, c(1,3,12));
})


#### ETSX / Regression + formula ####
# ETSX on N2568
xreg <- temporaldummy(Mcomp::M3[[2568]]$x);
testModel <- adam(Mcomp::M3[[2568]]$x, "MMN", h=18, holdout=TRUE, xreg=xreg);
test_that("ADAM ETSX(MMN) on N2568", {
    expect_true(is.matrix(testModel$xreg));
})

# ETSX selection on N2568
testModel <- adam(Mcomp::M3[[2568]]$x, "ZZZ", h=18, holdout=TRUE, xreg=xreg, xregDo="select");
test_that("ADAM ETSX(ZZZ) + xreg selection on N2568", {
    expect_match(testModel$xregDo, "select");
})

# ETSX adaption on N2568
testModel <- adam(Mcomp::M3[[2568]]$x, "MMN", h=18, holdout=TRUE, xreg=xreg, xregDo="adapt");
test_that("ADAM ETSX(MMN) + xreg adapt on N2568", {
    expect_match(testModel$xregDo, "adapt");
})

# Forecast from ETSX with formula
testForecast <- forecast(testModel, h=18, newxreg=tail(xreg, 18), interval="simulated");
test_that("Forecast for ADAM adaptive regression on N2568", {
    expect_equal(testForecast$level, 0.95);
})

# ETSX with formula
xreg <- data.frame(y=Mcomp::M3[[2568]]$x, x=factor(xreg %*% c(1:12)));
testModel <- adam(xreg, "MMN", h=18, holdout=TRUE, formula=y~x);
test_that("ADAM ETSX(MMN) + xreg formula on N2568", {
    expect_match(testModel$xregDo, "use");
})

# Forecast from ETSX with formula
testForecast <- forecast(testModel, h=18, newxreg=tail(xreg, 18), interval="nonp");
test_that("Forecast for ADAM ETSX(MMN) + xreg formula on N2568", {
    expect_equal(testForecast$level, 0.95);
})

# Adaptive regression
testModel <- adam(xreg, "NNN", h=18, holdout=TRUE, formula=y~x, xregDo="adapt");
test_that("ADAM adaptive regression on N2568", {
    expect_match(testModel$xregDo, "adapt");
})

# Forecast from ETSX with formula
testForecast <- forecast(testModel, h=18, newxreg=tail(xreg, 18), interval="nonp");
test_that("Forecast for ADAM adaptive regression on N2568", {
    expect_equal(testForecast$level, 0.95);
})

# Pure regression
testModel <- adam(xreg, "NNN", h=18, holdout=TRUE, formula=y~x);
test_that("ADAM regression (ALM) on N2568", {
    expect_true(is.alm(testModel));
})

#### ETS + ARIMA / ARIMA ####

#### ARIMAX ####

#### Provided initial / persistence / arma ####

#### auto.adam ####
