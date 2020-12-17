uniformiseLandCover <- function(landCoverPointsList) {
  landCovers <- do.call("c", lapply(landCoverPointsList, function(x) colnames(x@data)))
  uniqueLandCovers <- unique(landCovers)
  landCoverIndices <- as.numeric(substr(uniqueLandCovers, start = 10, stop = 100))
  uniqueLandCovers <- uniqueLandCovers[order(landCoverIndices)]
  lapply(landCoverPointsList, function(landCoverPoints) {
    if (length(missingCols <- setdiff(uniqueLandCovers, colnames(landCoverPoints@data))) > 0) {
      landCoverPoints@data[missingCols] <- 0
      landCoverPoints@data <- landCoverPoints@data[ , uniqueLandCovers]
    }
    landCoverPoints
  })
}

produceLandCover <- function(landCoverFiles, regionPolygon) {
  landCoverRasters <- lapply(landCoverFiles, function(filename) {
    landCoverSds <- MODIS::getSds(filename)
    landCover <- raster::raster(rgdal::readGDAL(landCoverSds$SDS4gdal[2], as.is = TRUE)) # Based on land type classification 2: https://lpdaac.usgs.gov/products/mcd12q1v006/
    landCover
  })
  landCoverRasters <- uniformiseLandCover(landCoverRasters)
  landCover <- do.call(raster::merge, landCoverRasters)
  smallerRaster <- raster::crop(x = landCover, y = regionPolygon)
  spObject <- raster::rasterToPoints(smallerRaster, spatial = TRUE)
  indiaValuesIndex <- sp::over(x = spObject, y = regionPolygon)
  pointsInIndia <- subset(spObject, subset = !is.na(indiaValuesIndex))
  raster::values(smallerRaster) <- rep(NA, raster::ncell(smallerRaster))
  output <- raster::rasterize(x = pointsInIndia, y = smallerRaster, field = "layer")
  output
}

prepareCovariateDataForISMRA <- function(elevationsRasterListWGS, landCoverRasterSinusoidal, collectionDatesPOSIX) {
  landCoverPoints <- raster::rasterToPoints(landCoverRasterSinusoidal, spatial = TRUE) # The last column indicates layer, which we don't need (output of rasterToPoints is a RasterLayer object).

  landCoverPointsWGS <- sp::spTransform(landCoverPoints, raster::crs(elevationsRasterListWGS[[1]]))
  elevationValues <- rep(NA, length(landCoverPointsWGS))
  elevationRasterIndex <- 0
  repeat {
    elevationRasterIndex <- elevationRasterIndex + 1
    indicesToExtract <- which(is.na(elevationValues))
    elevationValuesInListElement <- as.vector(raster::extract(elevationsRasterListWGS[[elevationRasterIndex]], landCoverPointsWGS[indicesToExtract, ])) # extract returns a matrix by default, each column for a different layer.
    elevationValues[indicesToExtract] <- elevationValuesInListElement
    if (!any(is.na(elevationValues)) | (elevationRasterIndex >= length(elevationsRasterListWGS))) break
  }

  combinedData <- cbind(landCover = landCoverPointsWGS@data[ , 1], elevation = elevationValues)

  coordinates <- landCoverPointsWGS@coords

  nonMissingLandCoverOrElevationIndices <- which(do.call("&" , lapply(1:ncol(combinedData), function(colIndex) !is.na(combinedData[ , colIndex]))))
  numNonMissing <- length(nonMissingLandCoverOrElevationIndices)
  timeValuesExtended <- rep(collectionDatesPOSIX, each = numNonMissing)
  coordinatesExtended <- coordinates[rep(nonMissingLandCoverOrElevationIndices, length(collectionDatesPOSIX)), ]
  rownames(coordinatesExtended) <- as.character(1:nrow(coordinatesExtended))
  combinedDataExtended <- combinedData[rep(nonMissingLandCoverOrElevationIndices, length(collectionDatesPOSIX)), ]

  spacetime::STIDF(sp = sp::SpatialPoints(coordinatesExtended, proj4string = raster::crs(elevationsRasterListWGS[[1]])), time = timeValuesExtended, data = as.data.frame(combinedDataExtended))
}

fitSPDE <- function(responseVec, covariateMatrix, coordinatesMatrix, timeVecNumeric, predCoordinatesMatrix, predCovariateMatrix, predTimeVecNumeric, numThreads = 1, dataCovarianceMatrix, control = list()) {
  control <- do.call("create.SPDE.control", control)
  covariateFrame <- as.data.frame(covariateMatrix)
  predCovariateFrame <- as.data.frame(predCovariateMatrix)
  timeVecTraining <- timeVecNumeric - min(timeVecNumeric) + 1
  timeVecTest <- predTimeVecNumeric - min(predTimeVecNumeric) + 1
  knots <- seq(1, max(timeVecTraining), length = max(timeVecTraining))
  mesh1 <- INLA::inla.mesh.1d(loc = knots, degree = 2, boundary = "free")

  ## generate space mesh

  mesh2 <- INLA::inla.mesh.2d(loc = coordinatesMatrix[timeVecTraining == 1, ], cutoff = control$mesh.2d.cutoff, offset = control$mesh.2d.offset, max.n = control$mesh.2d.max.n)

  # range0 and sigma0 control the prior means for the range and scale parameters.
  # See Lindgren INLA tutorial page 5.
  spatialSmoothness <- control$alpha - control$d/2 # cf p.3 INLA tutorial

  # range0 and sigma0 seem to be the prior means...
  range0 <- sqrt(8 * spatialSmoothness)/control$kappa # sqrt(8 * spatial smoothness) / Kappa. In my model, I use 1 as prior mean for spatial range and fix smoothness at 1.5. This means Kappa = 1.
  lkappa0 <- log(8 * spatialSmoothness)/2 - log(range0)
  ltau0 <- 0.5*log(gamma(spatialSmoothness)/(gamma(control$alpha)*(4*pi)^(control$d/2))) - log(control$sigma0) - spatialSmoothness * lkappa0

  ## build the spatial spde
  spde <- INLA::inla.spde2.matern(
    mesh2,
    B.tau = matrix(c(ltau0, -1, spatialSmoothness), 1, 3),
    B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
    theta.prior.mean = c(0,0), theta.prior.prec = c(1/control$loghyperparaSDinMyModel^2, 1/control$loghyperparaSDinMyModel^2))

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

  formulaForSPDE <- y ~ -1 + elevation + May28 + May29 + EvergreenBroadleaf + MixedForest + ClosedShrublands + Savannas + Grasslands + PermanentWetlands + Croplands + CroplandNaturalMosaics + NonVegetated + f(space, model = spde, group = space.group, control.group = list(model = "ar1"))

  SPDEresult <- INLA::inla(
    formulaForSPDE,
    data = INLA::inla.stack.data(combinedStack),
    control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(combinedStack)),
    num.threads = numThreads)
  list(fittedModel = SPDEresult, predictedValues = SPDEresult$summary.linear.predictor)
}

create.SPDE.control <- function(
  mesh.2d.cutoff = 0.01,
  mesh.2d.offset = c(0.1, 0.2),
  mesh.2d.max.n = 2000,
  d = 1,
  alpha = 2,
  kappa = 1,
  loghyperparaSDinMyModel = log(10),
  sigma0 = 1) {  # = 1/(range parameter in my model))
  list(mesh.2d.cutoff = mesh.2d.cutoff, mesh.2d.offset = mesh.2d.offset, mesh.2d.max.n = mesh.2d.max.n, d = d, alpha = alpha, kappa = kappa, loghyperparaSDinMyModel = loghyperparaSDinMyModel, sigma0 = sigma0)
}

fitISMRA <- function(responseVec, coordinatesMatrix, predCoordinatesMatrix, covariateMatrix, predCovariateMatrix, timeVecNumeric, predTimeVecNumeric, numThreads = 1, dataCovarianceMatrix, control) {
  control <- do.call("create.ISMRA.control", control)
  control$control$numOpenMPthreads <- numThreads
  hyperNormalList <- list(
    space = list(
      smoothness = c(mu = control$fixedHyperValues$space[["smoothness"]], sigma = control$logHyperpriorSD),
      rho = c(mu = control$hyperStart$space[["rho"]], sigma = control$logHyperpriorSD)),
    time = list(
      smoothness = c(mu = control$fixedHyperValues$time[["smoothness"]], sigma = control$logHyperpriorSD),
      rho = c(mu = control$hyperStart$time[["rho"]], sigma = control$logHyperpriorSD)),
    scale = c(mu = control$hyperStart$scale, sigma = control$logHyperpriorSD * 2), # Prior should be more vague, see Lindgren INLA tutorial p. 12
    errorSD = c(mu = control$fixedHyperValues$errorSD , sigma = control$logHyperpriorSD),
    fixedEffSD = c(mu = control$fixedHyperValues$fixedEffSD, sigma = control$logHyperpriorSD)
  )

  ISMRAfit <- MRAinla::INLAMRA(
    responseVec = responseVec,
    covariateFrame = as.data.frame(covariateMatrix),
    spatialCoordMat = as.matrix(coordinatesMatrix),
    timePOSIXorNumericVec = timeVecNumeric,
    predCovariateFrame = as.data.frame(predCovariateMatrix),
    predSpatialCoordMat = predCoordinatesMatrix,
    predTimePOSIXorNumericVec = predTimeVecNumeric,
    spatialRangeList = list(start = control$hyperStart$space[["rho"]], hyperpars = hyperNormalList$space$rho),
    spatialSmoothnessList = list(start = control$fixedHyperValues$space[["smoothness"]]),
    timeRangeList = list(start = control$hyperStart$time[["rho"]], hyperpars = hyperNormalList$time$rho),
    timeSmoothnessList = list(start = control$fixedHyperValues$time[["smoothness"]]),
    scaleList = list(start = control$hyperStart$scale, hyperpars = hyperNormalList$scale),
    errorSDlist = list(start = control$fixedHyperValues$errorSD),
    fixedEffSDlist = list(start = control$fixedHyperValues$fixedEffSD),
    control = control$control
  )
  list(
    fittedModel = ISMRAfit,
    predictedValues = ISMRAfit$predictionMoments$predictMeans)
}

create.ISMRA.control <- function(
  hyperStart = list(
    space = c(rho = 0),
    time = c(rho = 0),
    scale = 0),
  fixedHyperValues = list(
    space = c(smoothness = log(1.5)),
    time = c(smoothness = log(0.5)),
    errorSD = log(0.5),
    fixedEffSD = log(10)),
  logHyperpriorSD = 2,
  control = list(
    Mlon = 2,
    Mlat = 2,
    Mtime = 0,
    numKnotsRes0 = 20,
    numIterOptim = 20,
    tipKnotsThinningRate = 1,
    numValuesForIS = 100
    )
) {
  list(hyperStart = hyperStart, fixedHyperValues = fixedHyperValues, logHyperpriorSD = logHyperpriorSD, control = control)
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

fitModels <- function(responseVec, covariateMatrix, coordinatesMatrix, timeVecNumeric, obsIndicesForTraining, funToFitSPDE, funToFitVecchia, funToFitISMRA, dataCovarianceMatrix = dataCovarianceMatrix, controlForVecchia = list(), controlForISMRA = list(), controlForSPDE = list(), numThreads) {
  responseVecForTraining <- responseVec[obsIndicesForTraining]

  covariateMatrixForTraining <- covariateMatrix[obsIndicesForTraining, ]
  predCovariateMatrix <- covariateMatrix[!obsIndicesForTraining, ]

  coordinatesMatrixForTraining <- coordinatesMatrix[obsIndicesForTraining, ]
  predCoordinatesMatrix <- coordinatesMatrix[!obsIndicesForTraining, ]

  timeVecNumericForTraining <- timeVecNumeric[obsIndicesForTraining]
  predTimeVecNumeric <- timeVecNumeric[!obsIndicesForTraining]

  dataCovarianceMatrixSub <- dataCovarianceMatrix[obsIndicesForTraining, obsIndicesForTraining]

  controlAndFunToFitList <- list(
    # Vecchia = list(funToFit = fitVecchia, control = controlForVecchia),
    ISMRA = list(funToFit = fitISMRA, control = controlForISMRA),
    SPDE = list(funToFit = fitSPDE, control = controlForSPDE))
  lapply(
    X = controlAndFunToFitList,
    FUN = function(listElement) listElement$funToFit(responseVec = responseVecForTraining, covariateMatrix = covariateMatrixForTraining, coordinatesMatrix = coordinatesMatrixForTraining, predCovariateMatrix = predCovariateMatrix, predCoordinatesMatrix = predCoordinatesMatrix, timeVecNumeric = timeVecNumericForTraining, predTimeVecNumeric = predTimeVecNumeric, dataCovarianceMatrix = dataCovarianceMatrixSub, numThreads = numThreads, control = listElement$control))
}

# If saveDirectory is provided, simulationFun does not return anything.
# Might be preferable from a memory standpoint if outputs are large.

simulationFun <- function(datasetIndex, responseMatrix, covariateMatrix, coordinatesMatrix, timeVecNumeric, obsIndicesForTraining, funToFitSPDE, funToFitVecchia, funToFitISMRA, dataCovarianceMatrix, controlForVecchia = list(), controlForISMRA = list(), controlForSPDE = list(), saveDirectory = NULL, numThreads = 1) {
  cat("Processing simulated dataset", datasetIndex, "... ")
  fittedModel <- fitModels(responseVec = responseMatrix[ , datasetIndex], covariateMatrix = covariateMatrix, coordinatesMatrix = coordinatesMatrix, timeVecNumeric = timeVecNumeric, obsIndicesForTraining = obsIndicesForTraining, funToFitSPDE = funToFitSPDE, funToFitVecchia = funToFitVecchia, funToFitISMRA = funToFitISMRA, dataCovarianceMatrix = dataCovarianceMatrix, controlForVecchia = controlForVecchia, controlForISMRA = controlForISMRA, controlForSPDE = controlForSPDE, numThreads = numThreads)
  if (!is.null(saveDirectory)) {
    filename <- paste(saveDirectory, "/ISMRAsimulationResults_Dataset", datasetIndex, ".rds", sep = "")
    saveRDS(fittedModel, file = filename, compress = TRUE)
    return(NULL)
  } else {
    return(fittedModel)
  }
  cat("Done! \n")
  NULL
}

analysePartialResults <- function(folderForSimResults, pattern, datasetList, obsIndicesForTraining) {
  filesToImport <- list.files(folderForSimResults, pattern = pattern, full.names = TRUE)
  filesIndices <- as.numeric(stringr::str_extract(filesToImport, pattern = "(?<=Dataset).+(?=\\.rds)"))
  filesToImportInOrder <- filesToImport[order(filesIndices)]

  computePredStatsFromISMRA <- function(ISMRAoutput, realTestValues) {
      predictions <- ISMRAoutput$fittedModel$predMoments$Mean
      MSE <- mean((realTestValues - predictions)^2)
      MeanSD <- mean(ISMRAoutput$fittedModel$predMoments$SD)
      list(MSE = MSE, MeanSD = MeanSD)
  }

  computePredStatsFromSPDE <- function(SPDEoutput, realTestValues) {
    predictions <- ISMRAoutput$fittedModel$predMoments$Mean
    MSE <- mean((realTestValues - predictions)^2)
    MeanSD <- mean(ISMRAoutput$fittedModel$predMoments$SD)
    list(MSE = MSE, MeanSD = MeanSD)
  }

  computePredStats <- function(datasetIndex, obsIndicesForTraining, datasetList, filesToImportInOrder) {
    realTestValues <- datasetList$responses[!obsIndicesForTraining, datasetIndex]

  }

  predictionStatisticsByDataset <- lapply(sort(filesIndices), computePredStats, obsIndicesForTraining = obsIndicesForTraining, datasetList = datasetList, filesToImportInOrder = filesToImportInOrder)

}

getPredictionsFromINLAoutput <- function(INLAoutput, predictionLocations) {
  preds <- INLAoutput$summary.linear.predictor
  simdf <- testDataValidation@data
  simdf$pred <- tail(preds$mean, n = nrow(simdf))
  stackIndex <- inla::inla.stack.index(combinedStack, "predictions")$data
  preds$mean[stackIndex]
}

buildSpaceAndTimeMesh <- function(coordinatesMatrixTraining, timeVecNumericTraining, control) {
  timeVecTraining <- timeVecNumericTraining - min(timeVecNumericTraining) + 1
  knots <- seq(1, max(timeVecTraining), length = max(timeVecTraining))
  meshTime <- INLA::inla.mesh.1d(loc = knots, degree = 2, boundary = "free")

  ## generate space mesh

  meshSpace <- INLA::inla.mesh.2d(loc = coordinatesMatrixTraining[timeVecNumericTraining == 1, ], cutoff = control$mesh.2d.cutoff, offset = control$mesh.2d.offset, max.n = control$mesh.2d.max.n)
  list(meshTime = meshTime, meshSpace = meshSpace)
}

buildInlaStack <- function(coordinatesMatrixTraining, timeVecTraining, coordinatesMatrixTest, timeVecTest, meshForSpace, meshForTime, responseVec, covariateFrame, control) {
  # range0 and sigma0 control the prior means for the range and scale parameters.
  # See Lindgren INLA tutorial page 5.
  spatialSmoothness <- control$alpha - control$d/2 # cf p.3 INLA tutorial

  # range0 and sigma0 seem to be the prior means...
  range0 <- sqrt(8 * spatialSmoothness)/control$kappa # sqrt(8 * spatial smoothness) / Kappa. In my model, I use 1 as prior mean for spatial range and fix smoothness at 1.5. This means Kappa = 1.
  lkappa0 <- log(8 * spatialSmoothness)/2 - log(range0)
  ltau0 <- 0.5*log(gamma(spatialSmoothness)/(gamma(control$alpha)*(4*pi)^(control$d/2))) - log(control$sigma0) - spatialSmoothness * lkappa0

  ## build the spatial spde
  spde <- INLA::inla.spde2.matern(
    meshForSpace,
    B.tau = matrix(c(ltau0, -1, spatialSmoothness), 1, 3),
    B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
    theta.prior.mean = c(0,0), theta.prior.prec = c(1/control$loghyperparaSDinMyModel^2, 1/control$loghyperparaSDinMyModel^2))
  ## build the space time indices
  STindex <- INLA::inla.spde.make.index("space", n.spde = spde$n.spde, n.group = meshForTime$m)

  Atraining <- INLA::inla.spde.make.A(meshForSpace, loc = coordinatesMatrixTraining, group = timeVecTraining, group.mesh = meshForTime)
  Atest <- INLA::inla.spde.make.A(meshForSpace, loc = coordinatesMatrixTest, group = timeVecTest, group.mesh = meshForTime)

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

  INLA::inla.stack(stackTraining, stackTest)
}

