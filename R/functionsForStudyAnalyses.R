fitSPDE <- function(responseVec, covariateMatrix, coordinatesMatrix, timeVecNumeric, predCoordinatesMatrix, predCovariateMatrix, predTimeVecNumeric, numThreads = 1, dataCovarianceMatrix) {
  covariateFrame <- as.data.frame(covariateMatrix)
  predCovariateFrame <- as.data.frame(predCovariateMatrix)
  timeVecTraining <- timeVecNumeric - min(timeVecNumeric) + 1
  timeVecTest <- predTimeVecNumeric - min(predTimeVecNumeric) + 1
  knots <- seq(1, max(timeVecTraining), length = max(timeVecTraining))
  mesh1 <- INLA::inla.mesh.1d(loc = knots, degree = 2, boundary = "free")

  ## generate space mesh

  mesh2 <- INLA::inla.mesh.2d(loc = coordinatesMatrix[timeVecTraining == 1, ], cutoff = 0.01, offset = c(0.1, 0.2), max.n = 2000)

  # range0 and sigma0 control the prior means for the range and scale parameters.
  # See Lindgren INLA tutorial page 5.
  d <- 1
  alpha <- 2
  kappa <- 1 # = 1/(range parameter in my model)
  spatialSmoothness <- alpha - d/2 # cf p.3 INLA tutorial
  loghyperparaSDinMyModel <- log(10)
  # range0 and sigma0 seem to be the prior means...
  range0 <- sqrt(8 * spatialSmoothness)/kappa # sqrt(8 * spatial smoothness) / Kappa. In my model, I use 1 as prior mean for spatial range and fix smoothness at 1.5. This means Kappa = 1.
  sigma0 <- 1
  lkappa0 <- log(8 * spatialSmoothness)/2 - log(range0)
  ltau0 <- 0.5*log(gamma(spatialSmoothness)/(gamma(alpha)*(4*pi)^(d/2))) - log(sigma0) - spatialSmoothness * lkappa0

  ## build the spatial spde
  spde <- INLA::inla.spde2.matern(
    mesh2,
    B.tau = matrix(c(ltau0, -1, spatialSmoothness), 1, 3),
    B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
    theta.prior.mean = c(0,0), theta.prior.prec = c(1/loghyperparaSDinMyModel^2, 1/loghyperparaSDinMyModel^2))

  ## build the space time indices
  STindex <- INLA::inla.spde.make.index("space", n.spde = spde$n.spde, n.group = mesh1$m)

  ## Link data and process

  Atraining <- INLA::inla.spde.make.A(mesh2, loc = coordinatesMatrix, group = timeVecTraining, group.mesh = mesh1)
  Atest <- INLA::inla.spde.make.A(mesh2, loc = predCoordinatesMatrix, group = timeVecTest, group.mesh = mesh1)

  stackTraining <- INLA::inla.stack(
    data = list(y = responseVec),
    A = c(list(Atraining), lapply(rep(1, ncol(covariateMatrix)), identity)),
    effects = list(
      c(STindex, list(intercept = 1)),
      list(elevation = covariateFrame$elevation),
      list(May28 = covariateFrame$May28),
      list(May29 = covariateFrame$May29),
      list(EvergreenBroadleaf = covariateFrame$EvergreenBroadleaf),
      list(MixedForest = covariateFrame$MixedForest),
      list(ClosedShrublands = covariateFrame$ClosedShrublands),
      list(Savannas = covariateFrame$Savannas),
      list(Grasslands = covariateFrame$Grasslands),
      list(PermanentWetlands = covariateFrame$PermanentWetlands),
      list(Croplands = covariateFrame$Croplands),
      list(Urban = covariateFrame$Urban),
      list(CroplandNaturalMosaics = covariateFrame$CroplandNaturalMosaics),
      list(NonVegetated = covariateFrame$NonVegetated)
    ), tag = "est")

  stackTest <- INLA::inla.stack(
    data = list(y = NA),
    A = c(list(Atest), lapply(rep(1, ncol(predCovariateFrame)), identity)),
    effects = list(
      c(STindex, list(intercept = 1)),
      list(elevation = predCovariateFrame$elevation),
      list(May28 = predCovariateFrame$May28),
      list(May29 = predCovariateFrame$May29),
      list(EvergreenBroadleaf = predCovariateFrame$EvergreenBroadleaf),
      list(MixedForest = predCovariateFrame$MixedForest),
      list(ClosedShrublands = predCovariateFrame$ClosedShrublands),
      list(Savannas = predCovariateFrame$Savannas),
      list(Grasslands = predCovariateFrame$Grasslands),
      list(PermanentWetlands = predCovariateFrame$PermanentWetlands),
      list(Croplands = predCovariateFrame$Croplands),
      list(Urban = predCovariateFrame$Urban),
      list(CroplandNaturalMosaics = predCovariateFrame$CroplandNaturalMosaics),
      list(NonVegetated = predCovariateFrame$NonVegetated)
    ),
    tag = 'predictions')

  combinedStack <- INLA::inla.stack(stackTraining, stackTest)

  formulaForSPDE <- y ~ -1 + elevation + May28 + May29 + EvergreenBroadleaf + MixedForest + ClosedShrublands + Savannas + Grasslands + PermanentWetlands + Croplands + CroplandNaturalMosaics + NonVegetated + INLA::f(space, model = spde, group = space.group, control.group = list(model = "ar1"))

  SPDEresult <- INLA::inla(
    formulaForSPDE,
    data = INLA::inla.stack.data(combinedStack),
    control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(combinedStack)),
    num.threads = numThreads)
  list(fittedModel = SPDEresult, predictedValues = SPDEresult$summary.linear.predictor)
}

fitISMRA <- function(responseVec, coordinatesMatrix, predCoordinatesMatrix, covariateMatrix, predCovariateMatrix, timeVecNumeric, predTimeVecNumeric, numThreads = 1, dataCovarianceMatrix) {
  # FOR IS-MRA FIT (starting values):
  hyperStart <- list(
    space = c(rho = 0),
    time = c(rho = 0),
    scale = 0)

  errorSD <- 0.5 # Based on https://landval.gsfc.nasa.gov/Results.php?TitleID=mod11_valsup10
  fixedEffSD <- 10

  fixedHyperValues <- list(
    space = c(smoothness = log(1.5)),
    time = c(smoothness = log(0.5)),
    errorSD = log(errorSD),
    fixedEffSD = log(fixedEffSD)
  )

  logHyperpriorSD <- 2

  hyperNormalList <- list(
    space = list(
      smoothness = c(mu = fixedHyperValues$space[["smoothness"]], sigma = logHyperpriorSD),
      rho = c(mu = hyperStart$space[["rho"]], sigma = logHyperpriorSD)),
    time = list(
      smoothness = c(mu = fixedHyperValues$time[["smoothness"]], sigma = logHyperpriorSD),
      rho = c(mu = hyperStart$time[["rho"]], sigma = logHyperpriorSD)),
    scale = c(mu = hyperStart$scale, sigma = logHyperpriorSD * 2), # Prior should be more vague, see Lindgren INLA tutorial p. 12
    errorSD = c(mu = fixedHyperValues$errorSD , sigma = logHyperpriorSD),
    fixedEffSD = c(mu = fixedHyperValues$fixedEffSD, sigma = logHyperpriorSD)
  )
  ###

  ISMRAfit <- MRAinla::INLAMRA(
    responseVec = responseVec,
    covariateFrame = as.data.frame(covariateMatrix),
    spatialCoordMat = as.matrix(coordinatesMatrix),
    timePOSIXorNumericVec = timeVecNumeric,
    predCovariateFrame = as.data.frame(predCovariateMatrix),
    predSpatialCoordMat = predCoordinatesMatrix,
    predTimePOSIXorNumericVec = predTimeVecNumeric,
    spatialRangeList = list(start = hyperStart$space[["rho"]], hyperpars = hyperNormalList$space$rho),
    spatialSmoothnessList = list(start = fixedHyperValues$space[["smoothness"]]),
    timeRangeList = list(start = hyperStart$time[["rho"]], hyperpars = hyperNormalList$time$rho),
    timeSmoothnessList = list(start = fixedHyperValues$time[["smoothness"]]),
    scaleList = list(start = hyperStart$scale, hyperpars = hyperNormalList$scale),
    errorSDlist = list(start = fixedHyperValues$errorSD),
    fixedEffSDlist = list(start = fixedHyperValues$fixedEffSD),
    control = list(
      Mlon = 2,
      Mlat = 2,
      Mtime = 0,
      numKnotsRes0 = 20,
      numIterOptim = 20,
      tipKnotsThinningRate = 1,
      numOpenMPthreads = numThreads)
  )
  list(
    fittedModel = ISMRAfit,
    predictedValues = ISMRAfit$predictionMoments$predictMeans)
}

# This function is for GPVecchia. It must take two lists of (spatiotemporal) locations, and return a length-$k$ vector giving their covariances.

customCovFct <- function(locs1, locs2) {
  locsValues <- lapply(list(locs1, locs2), function(locsValue) {
    if (is.vector(locsValue)) {
      locsValue <- matrix(locsValue, 1)
    }
    locsValue
  })
  distValues <- geosphere::distHaversine(locsValues[[1]][ , 1:2], locsValues[[2]][ , 1:2])/1000
  timeDistValues <- abs(locsValues[[1]][ , 3] - locsValues[[2]][ , 3])
  MRAinla::maternCov(distValues, smoothness = 1.5, rho = 1, scale = 1) * MRAinla::maternCov(timeDistValues, smoothness = 0.5, rho = 1, scale = 1)
}

fitVecchia <- function(responseVec, covariateMatrix, coordinatesMatrix, predCovariateMatrix, predCoordinatesMatrix, timeVecNumeric, predTimeVecNumeric, dataCovarianceMatrix) {
  coordinatesMatrixIncremented <- cbind(coordinatesMatrix, time = timeVecNumeric)
  predCoordinatesMatrixIncremented <- cbind(predCoordinatesMatrix, time = predTimeVecNumeric)
  VecchiaModelFit <- GPvecchia::vecchia_estimate(data = responseVec, locs = coordinatesMatrixIncremented, X = covariateMatrix, theta.ini = c(1, .1, 0.5), covmodel = dataCovarianceMatrix, output.level = 0)
  vecchiaPredicted <- GPvecchia::vecchia_pred(VecchiaModelFit, locs.pred = predCoordinatesMatrixIncremented, X.pred = predCovariateMatrix)
  list(fittedModel = VecchiaModelFit, predictedValues = vecchiaPredicted)
}

fitModels <- function(responseVec, covariateMatrix, coordinatesMatrix, timeVecNumeric, obsIndicesForTraining, funToFitSPDE, funToFitVecchia, funToFitISMRA, dataCovarianceMatrix = dataCovarianceMatrix) {

  responseVecForTraining <- responseVec[obsIndicesForTraining]

  covariateMatrixForTraining <- covariateMatrix[obsIndicesForTraining, ]
  predCovariateMatrix <- covariateMatrix[!obsIndicesForTraining, ]

  coordinatesMatrixForTraining <- coordinatesMatrix[obsIndicesForTraining, ]
  predCoordinatesMatrix <- coordinatesMatrix[!obsIndicesForTraining, ]

  timeVecNumericForTraining <- timeVecNumeric[obsIndicesForTraining]
  predTimeVecNumeric <- timeVecNumeric[!obsIndicesForTraining]

  dataCovarianceMatrixSub <- dataCovarianceMatrix[obsIndicesForTraining, obsIndicesForTraining]

  # lapply(list(Vecchia = funToFitVecchia, SPDE = funToFitSPDE, ISMRA = funToFitISMRA),
  lapply(
    list(SPDE = funToFitSPDE, ISMRA = funToFitISMRA),
    FUN = function(fitFunctionName) fitFunctionName(responseVec = responseVecForTraining, covariateMatrix = covariateMatrixForTraining, coordinatesMatrix = coordinatesMatrixForTraining, predCovariateMatrix = predCovariateMatrix, predCoordinatesMatrix = predCoordinatesMatrix, timeVecNumeric = timeVecNumericForTraining, predTimeVecNumeric = predTimeVecNumeric, dataCovarianceMatrix = dataCovarianceMatrixSub))
}
