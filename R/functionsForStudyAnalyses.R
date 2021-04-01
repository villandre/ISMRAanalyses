listSDStoImport <- function(searchString, rawDataFilesLocation, dayOffset, dayRange, collectionDates) {
  temperatureFiles <- list.files(path = rawDataFilesLocation, pattern = searchString, full.names = TRUE)
  subFiles <- sapply(paste("A2012", dayOffset + dayRange, sep = ""), grep, x = temperatureFiles, value = TRUE)
  temperatures <- lapply(subFiles, MODIS::getSds)
  splitTemperatures <- split(temperatures, f = factor(substr(subFiles, start = 0, stop = gregexpr(pattern = ".h2", text = subFiles[[1]])[[1]] - 1)))
  names(splitTemperatures) <- collectionDates
  splitTemperatures
}

funToCreateRaster <- function(temperatureSdsList, polygonBound) {
  extractionFun <- function(x) {
    tempGrid <- rgdal::readGDAL(x$SDS4gdal[1], as.is = TRUE)
    hourGrid <- rgdal::readGDAL(x$SDS4gdal[3], as.is = TRUE)
    tempGrid$band1 <- tempGrid$band1 * 0.02 - 273.15 # See https://gis.stackexchange.com/questions/72524/how-do-i-convert-the-lst-values-on-the-modis-lst-image-to-degree-celsius
    # There's a 0.02 scaling factor applied to values in file to get the temperatures.
    # The -273.15 brings temperatures back in Celsius
    hourGrid@data[,1] <- hourGrid@data[,1] * 0.1

    list(temperatureRaster = raster::raster(tempGrid), hourRaster = raster::raster(hourGrid))
  }
  tempAndTimeRasters <- lapply(temperatureSdsList, extractionFun)

  createRaster <- function(rasterName) {
    rasterList <- lapply(tempAndTimeRasters, function(x) x[[rasterName]])
    mergedRasters <- do.call(raster::merge, rasterList)
    smallerRaster <- raster::crop(x = mergedRasters, y = polygonBound)
    spObject <- raster::rasterToPoints(smallerRaster, spatial = TRUE)
    polygonValuesIndex <- sp::over(x = spObject, y = polygonBound)
    pointsInPolygon <- subset(spObject, subset = !is.na(polygonValuesIndex))
    raster::values(smallerRaster) <- rep(NA, raster::ncell(smallerRaster))
    if (sum(!is.na(polygonValuesIndex)) == 0) {
      return(smallerRaster)
    }
    raster::rasterize(x = pointsInPolygon, y = smallerRaster, field = "layer")
  }
  rasterNames <- c("temperatureRaster", "hourRaster")
  tempAndTime <- lapply(rasterNames, FUN = createRaster)
  names(tempAndTime) <- rasterNames
  tempAndTime
}

funToGetDailyRastersAndSatelliteName <- function(var1, splitTemperaturesBySatellite, MaharashtraPolygonOtherCRS) {
  aquaRasters <- funToCreateRaster(splitTemperaturesBySatellite$Aqua[[var1]], polygonBound = MaharashtraPolygonOtherCRS)
  terraRasters <- funToCreateRaster(splitTemperaturesBySatellite$Terra[[var1]], polygonBound = MaharashtraPolygonOtherCRS)
  if (sum(!is.na(raster::values(aquaRasters$temperatureRaster))) >= sum(!is.na(raster::values(terraRasters$temperatureRaster)))) {
    cat("Returning Aqua!\n")
    c(aquaRasters, satellite = "Aqua")
  } else {
    cat("Returning Terra!\n")
    c(terraRasters, satellite = "Terra")
  }
}

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

# SPDEresult$misc$configs$config[[1]]$Q to access the Q matrix.
fitSPDE <- function(responseVec, covariateMatrix, coordinatesMatrix, timeVecNumeric, predCoordinatesMatrix, predCovariateMatrix, predTimeVecNumeric, numThreads = 1, control = list()) {
  coordinatesPoints <- sp::SpatialPoints(coords = coordinatesMatrix, proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
  newCoordinatesPoints <- sp::spTransform(coordinatesPoints, CRSobj = sp::CRS("+proj=utm +zone=43 +datum=WGS84 +units=km"))
  coordinatesMatrix <- newCoordinatesPoints@coords
  predCoordinatesMatrixPoints <- sp::SpatialPoints(coords = predCoordinatesMatrix, proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
  newPredCoordinatesPoints <- sp::spTransform(predCoordinatesMatrixPoints, CRSobj = sp::CRS("+proj=utm +zone=43 +datum=WGS84 +units=km"))
  predCoordinatesMatrix <- newPredCoordinatesPoints@coords
  timeVecTraining <- timeVecNumeric - min(timeVecNumeric) + 1
  spaceAndTimeMesh <- buildSpaceAndTimeMesh(coordinatesMatrixTraining = coordinatesMatrix, timeVecNumericTraining = timeVecTraining, control = control)

  inlaParameters <- produceINLAparameters(control)
  spde <- INLA::inla.spde2.matern(
    mesh = spaceAndTimeMesh$space,
    B.tau = matrix(c(inlaParameters$ltau0, -1, inlaParameters$spatialSmoothness), 1, 3),
    B.kappa = matrix(c(inlaParameters$lkappa0, 0, -1), 1, 3),
    theta.prior.mean = c(
      control$hyperStart$scale - log(control$sigma0),
      control$hyperStart$space[["rho"]] + log(2) - inlaParameters$lrange0),
    theta.prior.prec = c(1/control$logHyperpriorSDinISMRA^2, 1/control$logHyperpriorSDinISMRA^2))
  timeVecTest <- predTimeVecNumeric - min(timeVecNumeric) + 1
  combinedStack <- buildInlaStack(coordinatesMatrixTraining = coordinatesMatrix, timeVecTraining = timeVecTraining, coordinatesMatrixTest = predCoordinatesMatrix, timeVecTest = timeVecTest, meshForSpace = spaceAndTimeMesh$space, meshForTime = spaceAndTimeMesh$time, responseVecTraining = responseVec, covariateMatrixTraining = covariateMatrix, covariateMatrixTest = predCovariateMatrix, spdeObj = spde, control = control)
  error.prior.prec <- list(initial = 1/control$fixedHyperValues$errorSD^2, prior = "normal", fixed = TRUE)  # The precision in the Gaussian family is represented on the log-scale.
  control.family.value <- list(hyper = list(prec = error.prior.prec))
  randomValuesFromTimeRangePrior <- rnorm(10000, mean = control$hyperStart$time[["rho"]], sd = control$logHyperpriorSDinISMRA)
  transformedValues <- log(1 + exp(-1/exp(randomValuesFromTimeRangePrior))) - log(1 - exp(-1/exp(randomValuesFromTimeRangePrior)))
  meanForPrior <- mean(transformedValues)
  precForPrior <- 1/var(transformedValues)

  # formulaForSPDE <- y ~ 1 + elevation + May28 + May29 + EvergreenBroadleaf + MixedForest + ClosedShrublands + Savannas + Grasslands + PermanentWetlands + Croplands + Urban + CroplandNaturalMosaics + NonVegetated + f(space, model = spde, group = space.group, control.group = list(model = "ar1", hyper = list(theta = list(prior = "normal", param = c(mean = meanForPrior, precision = precForPrior), initial = meanForPrior, fixed = FALSE))))
  formulaForSPDE <- y ~ 1 + elevation + May28 + May29 + EvergreenBroadleaf + MixedForest + ClosedShrublands + Savannas + Grasslands + PermanentWetlands + Croplands + Urban + CroplandNaturalMosaics + NonVegetated + f(space, model = spde, group = space.group, control.group = list(model = "ar1", hyper = list("logit correlation" = list(prior = "normal", param = c(mean = meanForPrior, precision = precForPrior), initial = meanForPrior, fixed = FALSE))))

  SPDEresult <- tryCatch(
    expr = INLA::inla(
      formulaForSPDE,
      data = INLA::inla.stack.data(combinedStack),
      control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(combinedStack)),
      control.family = control.family.value, # Comment in to fix the precision for the error term.
      control.compute = list(config = TRUE, q = TRUE),
      control.fixed = list(mean = 0, prec = 1/exp(control$fixedHyperValues$fixedEffSD)^2),
      control.inla = list(h = 0.0075),
      num.threads = numThreads), error = function(e) e, finally = "Error in fitting SPDE! Return list will contain NAs.\n")
  returnResult <- predsAndSDs <- NULL
  if (!("simpleError" %in% class(SPDEresult))) {
    predsAndSDs <- getPredictionsAndSDsFromINLAoutputAlt(INLAoutput = SPDEresult, inlaStack = combinedStack, control = control)
    SPDEresult <- SPDEresult[grep(pattern = "summary", x = names(SPDEresult), value = TRUE)] # INLA objects can be huge. We only keep the elements we need.
    # We had to include in the vector of locations for which predictions are required two additional coordinates, one for day 1 (at the start) and the other for day 3 (at the end). The problem that else, the function to build A could not correctly identify the grouping. Some elements of SPDEresult will reflect that strategy. We will however subtract those unnecessary points from predictionMeans and predictionSDs.
    returnResult <- list(
      fittedModel = SPDEresult,
      predictionMeans = predsAndSDs$mean[2:(length(predsAndSDs$mean) - 1)],
      predictionSDs = predsAndSDs$sd[2:(length(predsAndSDs$sd) - 1)])
  } else {
    returnResult <- list(
      fittedModel = NA,
      predictionMeans = NA,
      predictionSDs = NA)
  }
  returnResult
}

create.SPDE.control <- function(
  mesh.2d.cutoff = 0.01,
  mesh.2d.offset = c(0.1, 0.2),
  mesh.2d.max.n = -1,
  mesh.2d.max.edge = 0.52,
  d = 1,
  alpha = 2,
  kappa0 = 1,
  sigma0 = 1,
  useFittedValues = FALSE) {  # = 1/(range parameter in my model))
  list(mesh.2d.cutoff = mesh.2d.cutoff, mesh.2d.max.edge = mesh.2d.max.edge, mesh.2d.offset = mesh.2d.offset, mesh.2d.max.n = mesh.2d.max.n, d = d, alpha = alpha, kappa0 = kappa0, sigma0 = sigma0, useFittedValues = useFittedValues)
}

fitISMRA <- function(responseVec, coordinatesMatrix, predCoordinatesMatrix, covariateMatrix, predCovariateMatrix, timeVecNumeric, predTimeVecNumeric, numThreads = 1, control) {
  control$control$numOpenMPthreads <- numThreads
  hyperNormalList <- list(
    space = list(
      smoothness = c(mu = control$fixedHyperValues$space[["smoothness"]], sigma = control$logHyperpriorSD),
      rho = c(mu = control$hyperStart$space[["rho"]], sigma = control$logHyperpriorSD)),
    time = list(
      smoothness = c(mu = control$fixedHyperValues$time[["smoothness"]], sigma = control$logHyperpriorSD),
      rho = c(mu = control$hyperStart$time[["rho"]], sigma = control$logHyperpriorSD)),
    scale = c(mu = control$hyperStart$scale, sigma = control$logHyperpriorSD),
    errorSD = c(mu = control$fixedHyperValues$errorSD , sigma = control$logHyperpriorSD),
    fixedEffSD = c(mu = control$fixedHyperValues$fixedEffSD, sigma = control$logHyperpriorSD)
  )

  ISMRAfit <- tryCatch(expr = MRAinla::INLAMRA(
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
  ), error = function(e) e)
  returnResult <- NULL
  if (!("simpleError" %in% class(ISMRAfit))) {
  returnResult <- list(
    fittedModel = ISMRAfit,
    predictionMeans = ISMRAfit$predMoments$Mean,
    predictionSDs = ISMRAfit$predMoments$SD)
  } else {
    returnResult <- list(
      fittedModel = NA,
      predictionMeans = NA,
      predictionSDs = NA)
  }
  returnResult
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

fitVecchia <- function(responseVec, covariateMatrix, coordinatesMatrix, predCovariateMatrix, predCoordinatesMatrix, timeVecNumeric, predTimeVecNumeric) {
  coordinatesMatrixIncremented <- cbind(coordinatesMatrix, time = timeVecNumeric)
  predCoordinatesMatrixIncremented <- cbind(predCoordinatesMatrix, time = predTimeVecNumeric)
  VecchiaModelFit <- GPvecchia::vecchia_estimate(data = responseVec, locs = coordinatesMatrixIncremented, X = covariateMatrix, theta.ini = c(1, .1, 0.5), covmodel = NULL, output.level = 0)
  vecchiaPredicted <- GPvecchia::vecchia_pred(VecchiaModelFit, locs.pred = predCoordinatesMatrixIncremented, X.pred = predCovariateMatrix)
  list(fittedModel = VecchiaModelFit, predictionMeans = vecchiaPredicted)
}

fitModels <- function(responseVec, covariateMatrix, coordinatesMatrix, timeVecNumeric, obsIndicesForTraining, funToFitSPDE, funToFitVecchia, funToFitISMRA, controlForVecchia = list(), controlForISMRA = list(), controlForSPDE = list(), numThreads) {
  responseVecForTraining <- responseVec[obsIndicesForTraining]
  controlForSPDE <- do.call("create.SPDE.control", controlForSPDE)
  controlForISMRA <- do.call("create.ISMRA.control", controlForISMRA)

  covariateMatrixForTraining <- covariateMatrix[obsIndicesForTraining, ]
  predCovariateMatrix <- covariateMatrix[!obsIndicesForTraining, ]

  coordinatesMatrixForTraining <- coordinatesMatrix[obsIndicesForTraining, ]
  predCoordinatesMatrix <- coordinatesMatrix[!obsIndicesForTraining, ]

  timeVecNumericForTraining <- timeVecNumeric[obsIndicesForTraining]
  predTimeVecNumeric <- timeVecNumeric[!obsIndicesForTraining]
  controlForSPDE$fixedHyperValues <- controlForISMRA$fixedHyperValues
  controlForSPDE$hyperStart <- controlForISMRA$hyperStart
  controlForSPDE$logHyperpriorSDinISMRA <- controlForISMRA$logHyperpriorSD
  controlAndFunToFitList <- list(
    # Vecchia = list(funToFit = fitVecchia, control = controlForVecchia),
    ISMRA = list(funToFit = fitISMRA, control = controlForISMRA),
    SPDE = list(funToFit = fitSPDE, control = controlForSPDE))

  lapply(
    X = controlAndFunToFitList,
    FUN = function(listElement) listElement$funToFit(responseVec = responseVecForTraining, covariateMatrix = covariateMatrixForTraining, coordinatesMatrix = coordinatesMatrixForTraining, predCovariateMatrix = predCovariateMatrix, predCoordinatesMatrix = predCoordinatesMatrix, timeVecNumeric = timeVecNumericForTraining, predTimeVecNumeric = predTimeVecNumeric, numThreads = numThreads, control = listElement$control))
}

# If saveDirectory is provided, simulationFun does not return anything.
# Might be preferable from a memory standpoint if outputs are large.

simulationFun <- function(datasetIndex, responseMatrix, covariateMatrix, coordinatesMatrix, timeVecNumeric, obsIndicesForTraining, funToFitSPDE = fitSPDE, funToFitVecchia = fitVecchia, funToFitISMRA = fitISMRA, controlForVecchia = list(), controlForISMRA = list(), controlForSPDE = list(), saveDirectory = NULL, numThreads = 1, resume = FALSE) {
  cat("Processing simulated dataset", datasetIndex, "... ")
  resultAvailable <- FALSE
  if (resume) {
    cat("Checking if dataset has already been handled... \n")
    filename <- paste(saveDirectory, "/ISMRAsimulationResults_Dataset", datasetIndex, ".rds", sep = "")
    resultAvailable <- file.exists(filename)
    if (resultAvailable) {
      cat("File ", filename, "already exists! Skipping...")
    }
  }
  fittedModel <- NULL
  if (!resultAvailable) {
    fittedModel <- fitModels(responseVec = responseMatrix[ , datasetIndex], covariateMatrix = covariateMatrix, coordinatesMatrix = coordinatesMatrix, timeVecNumeric = timeVecNumeric, obsIndicesForTraining = obsIndicesForTraining, funToFitSPDE = funToFitSPDE, funToFitVecchia = funToFitVecchia, funToFitISMRA = funToFitISMRA, controlForVecchia = controlForVecchia, controlForISMRA = controlForISMRA, controlForSPDE = controlForSPDE, numThreads = numThreads)
    if (!is.null(saveDirectory)) {
      filename <- paste(saveDirectory, "/ISMRAsimulationResults_Dataset", datasetIndex, ".rds", sep = "")
      saveRDS(fittedModel, file = filename, compress = TRUE)
      return(NULL)
    } else {
      return(fittedModel)
    }
    cat("Done! \n")
  }
  returnResult <- NULL
  if (is.null(saveDirectory)) {
    returnResult <- fittedModel
  }
  invisible(returnResult)
}

analysePredResults <- function(folderForSimResults, patternForFilename, simulatedDataList, obsIndicesForTraining, shiftISMRApostPredSDs = 0) {
  patternForExtractingNumber <- "[:digit:]+(?=\\.rds)"
  filesToImport <- list.files(folderForSimResults, pattern = patternForFilename, full.names = TRUE)
  filesIndices <- as.numeric(stringr::str_extract(filesToImport, pattern = patternForExtractingNumber))
  filesToImportInOrder <- filesToImport[order(filesIndices)]

  computePredStats <- function(filename, obsIndicesForTraining, simulatedDataList) {
    datasetIndex <- as.numeric(stringr::str_extract(filename, pattern = patternForExtractingNumber))
    simResults <- readRDS(filename)
    realValues <- simulatedDataList$responses[!obsIndicesForTraining, datasetIndex]
    modelNames <- names(simResults)
    names(modelNames) <- modelNames
    lapply(modelNames, function(modelName) {
      modelResult <- simResults[[modelName]]
      predictionSDs <- modelResult$predictionSDs
      if (modelName == "ISMRA") {
        predictionSDs <- modelResult$predictionSDs + shiftISMRApostPredSDs
      }
      coverageProb95 <- mean((realValues > modelResult$predictionMeans + predictionSDs * qnorm(0.025)) & (realValues < modelResult$predictionMeans + predictionSDs * qnorm(0.975)))
      c(MSPE = mean((modelResult$predictionMeans - realValues)^2), MedSPE = median((modelResult$predictionMeans - realValues)^2), MeanSD = mean(predictionSDs), MedSD = median(predictionSDs), CoverageProb_95 = coverageProb95)
    })
  }
  predStatsByDataset <- lapply(filesToImportInOrder, computePredStats, obsIndicesForTraining = obsIndicesForTraining, simulatedDataList = simulatedDataList)

  methodNames <- names(predStatsByDataset[[1]])
  names(methodNames) <- methodNames
  diffValues <- sapply(predStatsByDataset, function(x) x[["ISMRA"]] - x[["SPDE"]])
  frameForComparisonPlot <- data.frame(predStatName = rep(c("MSPE", "MedSPE"), each = ncol(diffValues)), Value = c(diffValues["MSPE", ], diffValues["MedSPE", ]))
  diffBoxPlot <- ggplot2::ggplot(data = frameForComparisonPlot, ggplot2::aes(predStatName, Value)) + ggplot2::geom_boxplot(outlier.colour = "red", outlier.shape = 1, colour = "blue") + ggplot2::theme_bw() + ggplot2::xlab("Prediction statistic") + ggplot2::theme(legend.position = "none", text = ggplot2::element_text(size = 16)) + ggplot2::scale_x_discrete(limits = c("MSPE", "MedSPE"))
  predStatsByMethod <- lapply(methodNames, function(methodName) {
    sapply(predStatsByDataset, "[[", methodName)
  })
  fancyMethodNames <- c("INLA-SPDE", "IS-MRA")
  names(fancyMethodNames) <- c("SPDE", "ISMRA")
  statsByMethodList <- lapply(methodNames, function(methodName) {
    data.frame(methodName = fancyMethodNames[[methodName]], predStatName = rep(c("MSPE", "MedSPE"), each = ncol(predStatsByMethod[[methodName]])), Value = c(predStatsByMethod[[methodName]]["MSPE", ], predStatsByMethod[[methodName]]["MedSPE", ]))
  })
  statsByMethodFrame <- do.call("rbind", statsByMethodList)
  boxPlotsByMethod <- ggplot2::ggplot(statsByMethodFrame, ggplot2::aes(x = predStatName, y = Value, colour = factor(methodName))) + ggplot2::geom_boxplot(outlier.colour = "red", outlier.shape = 1) + ggplot2::theme_bw() + ggplot2::xlab("Prediction statistic") + ggplot2::scale_colour_hue(name = "Method") + ggplot2::theme(legend.position = c(0.875, 0.875), text = ggplot2::element_text(size = 16)) + ggplot2::scale_x_discrete(limits = c("MSPE", "MedSPE"))
  getStatsByMethod <- function(resultMatrix) {
    t(apply(resultMatrix, 1, function(rowValues) c(Mean = mean(rowValues), SD = sd(rowValues), Median = median(rowValues), Minimum = min(rowValues), Maximum = max(rowValues))))
  }
  list(summaryStats = lapply(predStatsByMethod, getStatsByMethod), diffBoxPlots = diffBoxPlot, summaryBoxPlots = boxPlotsByMethod, summaryBoxPlotFrame = statsByMethodFrame, diffBoxPlotFrame = frameForComparisonPlot)
}

analyseParaEsts <- function(folderForSimResults, patternForFilename, simulatedDataList, obsIndicesForTraining, realFEs, realHyperparsLogScale, numSimsForGraphs = 50) {
  patternForExtractingNumber <- "[:digit:]+(?=\\.rds)"
  filesToImport <- list.files(folderForSimResults, pattern = patternForFilename, full.names = TRUE)
  filesIndices <- as.numeric(stringr::str_extract(filesToImport, pattern = patternForExtractingNumber))
  filesToImportInOrder <- filesToImport[order(filesIndices)]

  plotFrames <- .computePlotFrames(filesToImportInOrder = filesToImportInOrder, patternForExtractingNumber = patternForExtractingNumber, realFEs = realFEs)

  FEabsDiffByMethod <- lapply(unique(plotFrames$FE$Method), .getFEabsDiff, plotFrames = plotFrames, realFEs = realFEs)
  FEtableToPrint <- do.call("cbind", FEabsDiffByMethod)
  FEtableToPrint <- FEtableToPrint[ , c(rbind((1:ncol(FEabsDiffByMethod[[1]]) + ncol(FEabsDiffByMethod[[1]])), 1:ncol(FEabsDiffByMethod[[1]])))]

  FEcoverageByMethod <- sapply(unique(plotFrames$FE$Method), .getFEcoverageByMethod, plotFrames = plotFrames, realFEs = realFEs)
  colnames(FEcoverageByMethod) <- unique(plotFrames$FE$Method)
  rownames(FEcoverageByMethod) <- unique(plotFrames$FE$paraName)

  hyperparAbsDiffByMethod <- lapply(unique(plotFrames$hyperpar$Method), .funToGetHyperparAbsDiff, plotFrames = plotFrames, realHyperparsLogScale = realHyperparsLogScale)
  hyperparTableToPrint <- do.call("cbind", hyperparAbsDiffByMethod)
  hyperparTableToPrint <- hyperparTableToPrint[ , c(rbind((1:ncol(hyperparAbsDiffByMethod[[1]]) + ncol(hyperparAbsDiffByMethod[[1]])), 1:ncol(hyperparAbsDiffByMethod[[1]])))]

  hyperCoverageByMethod <- sapply(unique(plotFrames$hyperpar$Method), .getHyperCoverageByMethod, plotFrames = plotFrames, realHyperparsLogScale = realHyperparsLogScale)
  names(hyperCoverageByMethod) <- unique(plotFrames$hyperpar$Method)

  graphsForOutput <- .getAllGraphs(FEplotFrame = plotFrames$FE, hyperPlotFrame =  plotFrames$hyperpar, realHyperparsLogScale, realFEs = realFEs, numSimsForGraphs = numSimsForGraphs)
  list(hyperparGraphs = graphsForOutput$hyperparGraphs, FEgraphs = graphsForOutput$FEgraphs, coverageProbsMatrix = rbind(FEcoverageByMethod, hyperCoverageByMethod), parsAbsDiffSummary = rbind(FEtableToPrint, hyperparTableToPrint), parPlotFrames = plotFrames)
}

.getHyperCoverageByMethod <- function(methodName, plotFrames, realHyperparsLogScale) {
  getHyperCoverageByParaName <- function(hyperName) {
    subFrame <- subset(plotFrames$hyperpar, subset = (paraName == hyperName) & (Method == methodName))
    realValue <- realHyperparsLogScale[[hyperName]]
    testValues <- (realValue >= subFrame$CredInt2.5) & (realValue <= subFrame$CredInt97.5)
    mean(testValues)
  }
  hyperCoverageTest <- sapply(names(realHyperparsLogScale), getHyperCoverageByParaName)
  names(hyperCoverageTest) <- names(realHyperparsLogScale)
  hyperCoverageTest
}

.funToGetHyperparAbsDiff <- function(methodName, plotFrames, realHyperparsLogScale) {
  hyperparAbsDiffStatsByPara <- lapply(unique(plotFrames$hyperpar$paraName), function(parName) {
    output <- summary(abs(subset(plotFrames$hyperpar, subset = (paraName == parName) & (Method == methodName))$Mean - realHyperparsLogScale[[parName]]))
    # output <- output[c("Mean", "Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")]
    output <- output[c("Mean", "Min.", "Median", "Max.")]
    names(output) <- paste(names(output), methodName, sep = ".")
    output
  })
  hyperparAbsDiffStatsMat <- do.call("rbind", hyperparAbsDiffStatsByPara)
  rownames(hyperparAbsDiffStatsMat) <- unique(plotFrames$hyperpar$paraName)
  as.data.frame(hyperparAbsDiffStatsMat)
}

.getFEcoverageByMethod <- function(methodName, plotFrames, realFEs) {
  getFEcoverageByParaName <- function(FEname) {
    subFrame <- subset(plotFrames$FE, subset = (paraName == FEname) & (Method == methodName))
    realValue <- realFEs[[FEname]]
    testValues <- (realValue >= subFrame$CredInt2.5) & (realValue <= subFrame$CredInt97.5)
    mean(testValues)
  }
  FEcoverageTest <- sapply(unique(plotFrames$FE$paraName), getFEcoverageByParaName)
  names(FEcoverageTest) <- names(unique(plotFrames$FE$paraName))
  FEcoverageTest
}

.getFEabsDiff <- function(methodName, plotFrames, realFEs) {
  FEabsDiffStatsByPara <- lapply(unique(plotFrames$FE$paraName), function(parName) {
    output <- summary(abs(subset(plotFrames$FE, subset = (paraName == parName) & (Method == methodName))$Mean - realFEs[[parName]]))
    # output <- output[c("Mean", "Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")]
    output <- output[c("Mean", "Min.", "Median", "Max.")]
    names(output) <- paste(names(output), methodName, sep = ".")
    output
  })
  FEabsDiffStatsMat <- do.call("rbind", FEabsDiffStatsByPara)
  rownames(FEabsDiffStatsMat) <- unique(plotFrames$FE$paraName)
  as.data.frame(FEabsDiffStatsMat)
}

.getAllGraphs <- function(FEplotFrame, hyperPlotFrame, realHyperparsLogScale, realFEs, numSimsForGraphs) {
  hyperparGraphs <- lapply(unique(hyperPlotFrame$paraName), function(hyperparName) {
    ggplot2::ggplot(subset(hyperPlotFrame, subset = (dataIndex <= numSimsForGraphs) & (paraName == hyperparName)), ggplot2::aes(x = dataIndex, group = Method, colour = Method)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = CredInt2.5, ymax = CredInt97.5), width = .2, position = ggplot2::position_dodge(.9)) + ggplot2::theme_bw() + ggplot2::theme(legend.position = c(0.85, 0.9), text = ggplot2::element_text(size = 16)) + ggplot2::scale_colour_manual(values = c("goldenrod", "blue")) + ggplot2::geom_hline(yintercept = realHyperparsLogScale[[hyperparName]], linetype="dashed", color = "red") + ggplot2::xlab("Sim. dataset index") + ggplot2::ylab("Log-para. value")
  })
  names(hyperparGraphs) <- unique(hyperPlotFrame$paraName)
  FEgraphs <- lapply(unique(FEplotFrame$paraName), function(FEname) {
    ggplot2::ggplot(subset(FEplotFrame, subset = (dataIndex <= numSimsForGraphs) & (paraName == FEname)), ggplot2::aes(x = dataIndex, group = Method, colour = Method)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = CredInt2.5, ymax = CredInt97.5), width = .2, position = ggplot2::position_dodge(.9)) + ggplot2::theme_bw() + ggplot2::theme(legend.position = c(0.85, 0.9), text = ggplot2::element_text(size = 16)) + ggplot2::scale_colour_manual(values = c("goldenrod", "blue"), breaks = c("SPDE", "ISMRA"), labels = c("INLA-SPDE", "IS-MRA")) + ggplot2::geom_hline(yintercept = realFEs[[FEname]], linetype="dashed", color = "red") + ggplot2::xlab("Sim. dataset index") + ggplot2::ylab("Log-para. value")
  })
  names(FEgraphs) <- unique(FEplotFrame$paraName)
  list(hyperparGraphs = hyperparGraphs, FEgraphs = FEgraphs)
}

.computePlotFrames <- function(filesToImportInOrder, patternForExtractingNumber, realFEs) {
  computeFEplotFrames <- function(filename) {
    datasetIndex <- as.numeric(stringr::str_extract(filename, pattern = patternForExtractingNumber))
    simResults <- readRDS(filename)
    methodNames <- names(simResults)
    names(methodNames) <- methodNames
    getPlotFrame <- function(methodName) {
      FEandSDandCIs <- getFEmeansAndSDsAndCIs(simResults[[methodName]]$fittedModel)
      commonNames <- intersect(names(realFEs), rownames(FEandSDandCIs))
      output <- cbind(dataIndex = datasetIndex, Method = methodName, FEandSDandCIs[commonNames, ])
      rownames(output) <- NULL
      output
    }
    do.call("rbind", lapply(methodNames, getPlotFrame))
  }
  FEplotFrame <- do.call("rbind", lapply(filesToImportInOrder, FUN = computeFEplotFrames))
  hyperPlotFrame <- do.call("rbind", lapply(filesToImportInOrder, .funToGetHyperPlotFrame, patternForExtractingNumber = patternForExtractingNumber))
  list(FE = FEplotFrame, hyperpar = hyperPlotFrame)
}

.funToGetHyperPlotFrame <- function(filename, patternForExtractingNumber) {
  datasetIndex <- as.numeric(stringr::str_extract(filename, pattern = patternForExtractingNumber))
  simResults <- readRDS(filename)
  methodNames <- names(simResults)
  names(methodNames) <- methodNames
  ISMRAvalues <- simResults$ISMRA$fittedModel$hyperMarginalMoments[ , c("Mean", "CredInt_2.5%", "CredInt_97.5%")]
  SPDEparasToExtract <- c("Theta2 for space", "GroupRho for space", "Theta1 for space")
  SPDEvalues <- simResults$SPDE$fittedModel$summary.hyperpar[SPDEparasToExtract, c("mean", "0.025quant", "0.975quant")]
  # We adjust SPDE values to get them to match those in IS-MRA: parameterisations are different. We use medians because we'll need to transform GroupRho.
  SPDEvalues["Theta2 for space", ] <- SPDEvalues["Theta2 for space", ] - log(2) # Spatial range parameter in SPDE is twice that used in IS-MRA
  # SPDEvalues["GroupRho for space", ] <- log(-log((exp(SPDEvalues["GroupRho for space", ]) - 1)/(exp(SPDEvalues["GroupRho for space", ]) + 1)))
  simValuesForRhoTime <- rnorm(n = 5000, mean = simResults$SPDE$fittedModel$summary.hyperpar["GroupRho for space", "mean"], sd = simResults$SPDE$fittedModel$summary.hyperpar["GroupRho for space", "sd"]) ## In practice, skewness for GroupTheta is very small, hence the decision to simulate from a normal distribution.
  transformedSimValues <- -log(-log((exp(simValuesForRhoTime) - 1)/(exp(simValuesForRhoTime) + 1)))
  SPDEvalues["GroupRho for space", ] <- c(mean(transformedSimValues), quantile(x = transformedSimValues, probs = c(0.025, 0.975)))
  colnames(SPDEvalues) <- c("Mean", "CredInt_2.5%", "CredInt_97.5%")
  rownames(SPDEvalues) <- rownames(simResults$ISMRA$fittedModel$hyperMarginalMoments)
  SPDEresultFrame <- data.frame(dataIndex = datasetIndex, Mean = SPDEvalues[ , "Mean"], paraName = rownames(SPDEvalues), Method = "SPDE", CredInt2.5 = SPDEvalues[ , "CredInt_2.5%"], CredInt97.5 = SPDEvalues[ , "CredInt_97.5%"])
  ISMRAresultFrame <- data.frame(dataIndex = datasetIndex, Mean = ISMRAvalues[ , "Mean"], paraName = rownames(ISMRAvalues), Method = "ISMRA", CredInt2.5 = ISMRAvalues[ , "CredInt_2.5%"], CredInt97.5 = ISMRAvalues[ , "CredInt_97.5%"])
  rbind(SPDEresultFrame, ISMRAresultFrame)
}

getFEmeansAndSDsAndCIs <- function(outputObject) {
  output <- NULL
  if ("INLAMRA" %in% class(outputObject)) {
    output <- outputObject$FEmarginalMoments[ , c("Mean", "StdDev", "CredInt_2.5%", "CredInt_97.5%")]
  } else if (!is.null(outputObject$summary.fixed)) {
    output <- outputObject$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")]
  }
  colnames(output) <- c("Mean", "SD", "CredInt2.5", "CredInt97.5")
  cbind(paraName = rownames(output), output)
}

getPredictionsAndSDsFromINLAoutput <- function(INLAoutput, responseVecTraining, covariateMatrixTest, coordinatesMatrixTest, timeVecNumericTest, covariateMatrixTraining, coordinatesMatrixTraining, timeVecNumericTraining, control) {
  # Rebuilding components...
  control <- do.call("create.SPDE.control", control)
  timeVecNumericTraining <- timeVecNumericTraining - min(timeVecNumericTraining) + 1
  timeVecNumericTest <- timeVecNumericTest - min(timeVecNumericTest) + 1

  spaceAndTimeMesh <- buildSpaceAndTimeMesh(coordinatesMatrixTraining = coordinatesMatrixTraining, timeVecNumericTraining = timeVecNumericTraining, control = control)
  inlaParameters <- produceINLAparameters(control)
  ## build the spatial spde
  spde <- INLA::inla.spde2.matern(
    mesh = spaceAndTimeMesh$space,
    B.tau = matrix(c(inlaParameters$ltau0, -1, inlaParameters$spatialSmoothness), 1, 3),
    B.kappa = matrix(c(inlaParameters$lkappa0, 0, -1), 1, 3),
    theta.prior.mean = c(
      control$hyperStart$scale$mean - log(control$sigma0),
      control$hyperStart$space$range$mean + log(2) - inlaParameters$lrange0),
    theta.prior.prec = c(1/control$logHyperpriorSDinISMRA^2, 1/1/control$logHyperpriorSDinISMRA^2))

  combinedStack <- buildInlaStack(coordinatesMatrixTraining = coordinatesMatrixTraining, timeVecTraining = timeVecNumericTraining, coordinatesMatrixTest = coordinatesMatrixTest, timeVecTest = timeVecNumericTest, meshForSpace = spaceAndTimeMesh$space, meshForTime = spaceAndTimeMesh$time, responseVecTraining = responseVecTraining, covariateMatrixTraining = covariateMatrixTraining, covariateMatrixTest = covariateMatrixTest, spdeObj = spde, control = control)
  preds <- INLAoutput$summary.linear.predictor
  if (control$useFittedValues) {
    preds <- INLAoutput$summary.fitted.values
  }
  stackIndex <- INLA::inla.stack.index(combinedStack, "predictions")$data
  preds[stackIndex, c("mean", "sd")]
}

buildSpaceAndTimeMesh <- function(coordinatesMatrixTraining, timeVecNumericTraining, control) {

  # knots <- seq(1, max(timeVecNumericTraining), length = max(timeVecNumericTraining))
  knots <- range(timeVecNumericTraining)
  meshTime <- INLA::inla.mesh.1d(loc = knots, degree = 2, boundary = "free")

  ## generate space mesh

  meshSpace <- INLA::inla.mesh.2d(loc = coordinatesMatrixTraining[timeVecNumericTraining == 1, ], cutoff = control$mesh.2d.cutoff, offset = control$mesh.2d.offset, max.n = control$mesh.2d.max.n, max.edge = control$mesh.2d.max.edge)
  list(time = meshTime, space = meshSpace)
}

buildInlaStack <- function(coordinatesMatrixTraining, timeVecTraining, coordinatesMatrixTest, timeVecTest, meshForSpace, meshForTime, responseVecTraining, covariateMatrixTraining, covariateMatrixTest, spdeObj, control) {
  ## build the space time indices
  STindex <- INLA::inla.spde.make.index("space", n.spde = spdeObj$n.spde, n.group = length(unique(timeVecTraining)))
  # We add a dummy coordinate to coordinatesMatrixTest to account for a bug in inla.spde.make.A, where it cannot properly account for the prediction block being concentrated on one day in the middle.
  coordinatesMatrixTest <- rbind(coordinatesMatrixTraining[1, ], coordinatesMatrixTest, coordinatesMatrixTraining[nrow(coordinatesMatrixTraining), ])
  timeVecTest <- c(timeVecTraining[[1]], timeVecTest, tail(timeVecTraining, n = 1)[[1]])
  covariateMatrixTest <- rbind(covariateMatrixTraining[1, ], covariateMatrixTest, covariateMatrixTraining[nrow(covariateMatrixTraining), ])
  Atraining <- INLA::inla.spde.make.A(meshForSpace, loc = coordinatesMatrixTraining, group = timeVecTraining)
  Atest <- INLA::inla.spde.make.A(meshForSpace, loc = coordinatesMatrixTest, group = timeVecTest)
  covariateFrame <- as.data.frame(covariateMatrixTraining)
  predCovariateFrame <- as.data.frame(covariateMatrixTest)

  stackTraining <- INLA::inla.stack(
    data = list(y = responseVecTraining),
    A = c(list(Atraining), split(rep(1, ncol(covariateFrame)), 1:ncol(covariateFrame))),
    effects = list(
      STindex,
      elevation = covariateFrame$elevation,
      May28 = covariateFrame$May28,
      May29 = covariateFrame$May29,
      EvergreenBroadleaf = covariateFrame$EvergreenBroadleaf,
      MixedForest = covariateFrame$MixedForest,
      ClosedShrublands = covariateFrame$ClosedShrublands,
      Savannas = covariateFrame$Savannas,
      Grasslands = covariateFrame$Grasslands,
      PermanentWetlands = covariateFrame$PermanentWetlands,
      Croplands = covariateFrame$Croplands,
      Urban = covariateFrame$Urban,
      CroplandNaturalMosaics = covariateFrame$CroplandNaturalMosaics,
      NonVegetated = covariateFrame$NonVegetated), tag = "est")

  stackTest <- INLA::inla.stack(
    data = list(y = NA),
    A = c(list(Atest), split(rep(1, ncol(predCovariateFrame)), f = 1:ncol(predCovariateFrame))),
    effects = list(
      STindex,
      elevation = predCovariateFrame$elevation,
      May28 = predCovariateFrame$May28,
      May29 = predCovariateFrame$May29,
      EvergreenBroadleaf = predCovariateFrame$EvergreenBroadleaf,
      MixedForest = predCovariateFrame$MixedForest,
      ClosedShrublands = predCovariateFrame$ClosedShrublands,
      Savannas = predCovariateFrame$Savannas,
      Grasslands = predCovariateFrame$Grasslands,
      PermanentWetlands = predCovariateFrame$PermanentWetlands,
      Croplands = predCovariateFrame$Croplands,
      Urban = predCovariateFrame$Urban,
      CroplandNaturalMosaics = predCovariateFrame$CroplandNaturalMosaics,
     NonVegetated = predCovariateFrame$NonVegetated),
    tag = 'predictions')

  INLA::inla.stack(stackTraining, stackTest)
}

getPredictionsAndSDsFromINLAoutputAlt <- function(INLAoutput, inlaStack, control) {
  preds <- INLAoutput$summary.linear.predictor
  if (control$useFittedValues) {
    preds <- INLAoutput$summary.fitted.values
  }

  stackIndex <- INLA::inla.stack.index(inlaStack, "predictions")$data
  preds[stackIndex, c("mean", "sd")]
}

# This will update saved simulation results with predictedValues,
# which were missing due to a bug.
# The new version of the software should not require this function.

recomputePredictionsForSimOutputs <- function(folderForSimResults, patternForFilename = "Dataset.+\\.rds$", coordinatesMatrixTraining, coordinatesMatrixTest, timeVecNumericTraining, timeVecNumericTest, covariateMatrixTraining, covariateMatrixTest, responseMatrixTraining, saveResult = FALSE, controlForSPDE) {
  filesToImport <- list.files(folderForSimResults, pattern = patternForFilename, full.names = TRUE)

  fixFunction <- function(filename, responseMatrixTraining, covariateMatrixTraining, covariateMatrixTest, timeVecNumericTraining, timeVecNumericTest, coordinatesMatrixTraining, coordinatesMatrixTest) {
    dataIndex <- as.numeric(stringr::str_extract(filename, pattern = "[:digit:]+(?=\\.rds)"))
    simResults <- readRDS(filename)
    responseVec <- responseMatrixTraining[ , dataIndex]
    SPDEpredMeansAndSDs <- getPredictionsAndSDsFromINLAoutput(
      INLAoutput = simResults$SPDE$fittedModel,
      responseVecTraining = responseVec,
      covariateMatrixTest = covariateMatrixTest,
      coordinatesMatrixTest = coordinatesMatrixTest,
      timeVecNumericTest = timeVecNumericTest,
      covariateMatrixTraining = covariateMatrixTraining,
      coordinatesMatrixTraining =  coordinatesMatrixTraining,
      timeVecNumericTraining = timeVecNumericTraining,
      control = controlForSPDE)
    ISMRApredsAndSDs <- simResults$ISMRA$fittedModel$predMoments[ , c("Mean", "SD")]
    updatedSimResultsSPDE <- list(
      fittedModel = simResults$SPDE$fittedModel,
      predictionMeans = SPDEpredMeansAndSDs$mean,
      predictionSDs = SPDEpredMeansAndSDs$sd)
    updatedSimResultsISMRA <- list(
      fittedModel = simResults$ISMRA$fittedModel,
      predictionMeans = ISMRApredsAndSDs$Mean,
      predictionSDs = ISMRApredsAndSDs$SD)
    simResultsUpdate <- list(SPDE = updatedSimResultsSPDE, ISMRA = updatedSimResultsISMRA)
    if (saveResult) {
      saveRDS(simResultsUpdate, file = filename, compress = TRUE)
      cat("Updated predictions in ", filename, "\n")
      return(invisible(NULL))
    }
    simResultsUpdate
  }
  lapply(filesToImport, fixFunction , responseMatrixTraining = responseMatrixTraining, covariateMatrixTraining = covariateMatrixTraining,  timeVecNumericTraining = timeVecNumericTraining, covariateMatrixTest = covariateMatrixTest, timeVecNumericTest = timeVecNumericTest, coordinatesMatrixTraining = coordinatesMatrixTraining, coordinatesMatrixTest = coordinatesMatrixTest)
}

# SPDE objects are huge for nothing. We just need summaries. We therefore remove unneeded components
stripSPDEobjects <- function(folderForSimulationResults, patternForFilename = "Dataset") {
  filenames <- list.files(path = folderForSimulationResults, pattern = patternForFilename, full.names = TRUE)
  stripSPDEoutput <- function(filename) {
    bigSimResult <- readRDS(filename)
    bigSimResult$SPDE$fittedModel <- bigSimResult$SPDE$fittedModel[grep(pattern = "summary", x = names(bigSimResult$SPDE$fittedModel), value = TRUE)]
    saveRDS(bigSimResult, file = filename, compress = TRUE)
    cat("Stripped SPDE object in", filename, "\n")
    invisible(NULL)
  }
  lapply(filenames, stripSPDEoutput)
  invisible(NULL)
}

# This function will go through saved results and refit SPDE under a new set of control paramaters.
refitSPDE <- function(
  folderForSimulationResults,
  patternForFilename = "Dataset",
  responseMatrix,
  covariateMatrix,
  coordinatesMatrix,
  timeVecNumeric,
  funToFitSPDE = fitSPDE,
  controlForSPDE,
  controlForISMRA,
  numThreads,
  obsIndicesForTraining) {
  controlForSPDE <- do.call("create.SPDE.control", controlForSPDE)
  controlForISMRA <- do.call("create.ISMRA.control", controlForISMRA)
  controlForSPDE$fixedHyperValues <- controlForISMRA$fixedHyperValues
  controlForSPDE$hyperStart <- controlForISMRA$hyperStart
  controlForSPDE$logHyperpriorSDinISMRA <- controlForISMRA$logHyperpriorSD
  filenames <- list.files(path = folderForSimulationResults, pattern = patternForFilename, full.names = TRUE)
  refitSPDEinner <- function(filename) {
    dataIndex <- as.numeric(stringr::str_extract(filename, pattern = "[:digit:]+(?=\\.rds)"))
    responseVecForTraining <- responseMatrix[obsIndicesForTraining, dataIndex]
    covariateMatrixForTraining <- covariateMatrix[obsIndicesForTraining, ]
    coordinatesMatrixForTraining <- coordinatesMatrix[obsIndicesForTraining, ]
    timeVecNumericForTraining <- timeVecNumeric[obsIndicesForTraining]
    predCoordinatesMatrix <- coordinatesMatrix[!obsIndicesForTraining, ]
    predCovariateMatrix <- covariateMatrix[!obsIndicesForTraining, ]
    predTimeVecNumeric <- timeVecNumeric[!obsIndicesForTraining]

    listResult <- funToFitSPDE(responseVec = responseVecForTraining, covariateMatrix = covariateMatrixForTraining, coordinatesMatrix = coordinatesMatrixForTraining, predCovariateMatrix = predCovariateMatrix, predCoordinatesMatrix = predCoordinatesMatrix, timeVecNumeric = timeVecNumericForTraining, predTimeVecNumeric = predTimeVecNumeric, numThreads = numThreads, control = controlForSPDE)
    oldResult <- readRDS(filename)
    oldResult$SPDE <- listResult
    saveRDS(oldResult, filename)
    cat("Updated SPDE fit in", filename, "\n")
    invisible(NULL)
  }
  lapply(filenames, refitSPDEinner)
  invisible(NULL)
}

# This function will go through saved results and refit IS-MRA
refitISMRA <- function(
  folderForSimulationResults,
  patternForFilename = "Dataset",
  responseMatrix,
  covariateMatrix,
  coordinatesMatrix,
  timeVecNumeric,
  funToFitISMRA = fitISMRA,
  controlForISMRA,
  numThreads,
  obsIndicesForTraining) {
  controlForISMRA <- do.call("create.ISMRA.control", controlForISMRA)
  filenames <- list.files(path = folderForSimulationResults, pattern = patternForFilename, full.names = TRUE)
  refitISMRAinner <- function(filename) {
    dataIndex <- as.numeric(stringr::str_extract(filename, pattern = "[:digit:]+(?=\\.rds)"))
    responseVecForTraining <- responseMatrix[obsIndicesForTraining, dataIndex]
    covariateMatrixForTraining <- covariateMatrix[obsIndicesForTraining, ]
    coordinatesMatrixForTraining <- coordinatesMatrix[obsIndicesForTraining, ]
    timeVecNumericForTraining <- timeVecNumeric[obsIndicesForTraining]
    predCoordinatesMatrix <- coordinatesMatrix[!obsIndicesForTraining, ]
    predCovariateMatrix <- covariateMatrix[!obsIndicesForTraining, ]
    predTimeVecNumeric <- timeVecNumeric[!obsIndicesForTraining]

    listResult <- funToFitISMRA(responseVec = responseVecForTraining, covariateMatrix = covariateMatrixForTraining, coordinatesMatrix = coordinatesMatrixForTraining, predCovariateMatrix = predCovariateMatrix, predCoordinatesMatrix = predCoordinatesMatrix, timeVecNumeric = timeVecNumericForTraining, predTimeVecNumeric = predTimeVecNumeric, numThreads = numThreads, control = controlForISMRA)
    oldResult <- readRDS(filename)
    oldResult$ISMRA <- listResult
    saveRDS(oldResult, filename)
    cat("Updated ISMRA fit in", filename, "\n")
    invisible(NULL)
  }
  lapply(filenames, refitISMRAinner)
  invisible(NULL)
}

produceINLAparameters <- function(control) {
  # range0 and sigma0 control the prior means for the range and scale parameters.
  # See Lindgren INLA tutorial page 5.
  spatialSmoothness <- control$alpha - control$d/2 # cf p.3 INLA tutorial

  # range0 and sigma0 seem to be the prior means...
  range0 <- sqrt(8 * spatialSmoothness)/control$kappa0

  ltau0 <- 0.5 * log(gamma(spatialSmoothness)/(gamma(control$alpha)*(4*pi)^(control$d/2))) - log(control$sigma0) - spatialSmoothness * log(control$kappa0)
  list(spatialSmoothness = spatialSmoothness, lkappa0 = log(control$kappa0), ltau0 = ltau0, lrange0 = log(range0))
}

fitNewModel <- function(
  folderForSimulationResults,
  patternForFilename = "Dataset",
  responseMatrix,
  covariateMatrix,
  coordinatesMatrix,
  timeVecNumeric,
  funToFitNewModel = fitSPDE,
  newModelName = NULL,
  controlForNewModel,
  numThreads,
  obsIndicesForTraining) {
  if (is.null(newModelName)) {
    stop("Please input a name for the new model (argument newModelName)! \n")
  }
  filenames <- list.files(path = folderForSimulationResults, pattern = patternForFilename, full.names = TRUE)
  fitNewModel <- function(filename) {
    dataIndex <- as.numeric(stringr::str_extract(filename, pattern = "[:digit:]+(?=\\.rds)"))
    responseVecForTraining <- responseMatrix[obsIndicesForTraining, dataIndex]
    covariateMatrixForTraining <- covariateMatrix[obsIndicesForTraining, ]
    coordinatesMatrixForTraining <- coordinatesMatrix[obsIndicesForTraining, ]
    timeVecNumericForTraining <- timeVecNumeric[obsIndicesForTraining]
    predCoordinatesMatrix <- coordinatesMatrix[!obsIndicesForTraining, ]
    predCovariateMatrix <- covariateMatrix[!obsIndicesForTraining, ]
    predTimeVecNumeric <- timeVecNumeric[!obsIndicesForTraining]

    listResult <- funToFitNewModel(responseVec = responseVecForTraining, covariateMatrix = covariateMatrixForTraining, coordinatesMatrix = coordinatesMatrixForTraining, predCovariateMatrix = predCovariateMatrix, predCoordinatesMatrix = predCoordinatesMatrix, timeVecNumeric = timeVecNumericForTraining, predTimeVecNumeric = predTimeVecNumeric, numThreads = numThreads, control = controlForNewModel)
    oldResult <- readRDS(filename)
    oldResult[[newModelName]] <- listResult
    saveRDS(oldResult, filename)
    cat("Added", newModelName, "fit in", filename, "\n")
    invisible(NULL)
  }
  lapply(filenames, refitSPDE)
  invisible(NULL)
}

prepareDataForISMRA <- function(temperatures, elevations, landCover, satelliteNamesVec, collectionDatesPOSIX, completeDateVector = collectionDatesPOSIX) {
  if ("RasterLayer" %in% class(temperatures[[1]])) {
    temperaturePoints <- lapply(temperatures, FUN = raster::rasterToPoints, spatial = TRUE)
  } else {
    temperaturePoints <- temperatures
  }

  satelliteNamesList <- lapply(seq_along(satelliteNamesVec), function(dayIndex) {
    rep(satelliteNamesVec[[dayIndex]], nrow(temperaturePoints[[dayIndex]]@coords))
  })
  satellite <- do.call("c", satelliteNamesList)
  satellite <- as.numeric(factor(x = satellite, levels = c("Terra", "Aqua"))) - 1
  numTimePoints <- length(completeDateVector)

  timeValues <- do.call("c", lapply(seq_along(collectionDatesPOSIX), function(x) rep(collectionDatesPOSIX[[x]], length(temperaturePoints[[x]]))))

  timeLevels <- as.numeric(factor(as.character(timeValues), levels = as.character(completeDateVector))) - 1

  timeModelMatrix <- t(sapply(timeLevels, function(x) {
    unitVector <- rep(0, numTimePoints - 1)
    unitVector[x] <- 1 # When x takes value 0, the vector remains all 0s, which is what we want.
    unitVector
  }))
  colnames(timeModelMatrix) <- paste("time", 2:numTimePoints, sep = "")

  landCoverPoints <- lapply(temperaturePoints, function(tempPoints) {
    tempPointsReproj <- sp::spTransform(tempPoints, raster::crs(landCover))
    landCoverAtPoints <- raster::extract(landCover, tempPointsReproj)
    landCoverValues <- sort(unique(landCoverAtPoints))
    columnNames <- paste("landCover", landCoverValues, sep = "")
    landCoverMatrix <- t(sapply(landCoverAtPoints, function(x) {
      unitVec <- numeric(length(columnNames))
      unitVec[match(x, landCoverValues)] <- 1
      unitVec
    }))
    colnames(landCoverMatrix) <- columnNames
    sp::SpatialPointsDataFrame(coords = tempPoints@coords, data = as.data.frame(landCoverMatrix), proj4string = raster::crs(tempPoints))
  })
  landCoverPoints <- uniformiseLandCover(landCoverPoints)

  elevationPoints <- lapply(temperaturePoints, function(tempPoints) {
    tempPoints <- sp::spTransform(tempPoints, raster::crs(elevations[[1]]))
    elevationValues <- rep(0, length(tempPoints))
    lapply(elevations, function(elevationRaster) {
      extractedValues <- raster::extract(elevationRaster, tempPoints)
      elevationValues[!is.na(extractedValues)] <<- extractedValues[!is.na(extractedValues)]
      NULL
    })
    sp::SpatialPointsDataFrame(coords = tempPoints@coords, data = data.frame(elevation = elevationValues), proj4string = raster::crs(tempPoints))
  })

  latitudePoints <- lapply(temperaturePoints, function(x) {
    sp::SpatialPointsDataFrame(x@coords, data = data.frame(latitude = x@coords[, 2]), proj4string = raster::crs(x))
  })
  combinedData <- do.call("cbind", lapply(list(landCoverPoints, latitudePoints, elevationPoints), function(x) do.call("rbind", lapply(x, function(y) y@data))))
  combinedData <- cbind(combinedData, timeModelMatrix, Aqua = satellite)
  if ("RasterLayer" %in% class(temperatures[[1]])) {
    combinedData <- cbind(do.call("rbind", lapply(temperaturePoints, function(y) y@data)), combinedData)
    colnames(combinedData)[[1]] <- "y"
  }

  coordinates <- do.call("rbind", lapply(temperaturePoints, function(x) x@coords))
  rownames(coordinates) <- as.character(1:nrow(coordinates))
  missingLandCoverOrElevation <- (rowSums(combinedData[ , grep(colnames(combinedData), pattern = "landCover", value = TRUE)]) == 0) | is.na(combinedData[, "elevation"])

  spacetime::STIDF(sp = sp::SpatialPoints(coordinates[!missingLandCoverOrElevation, ], proj4string = raster::crs(temperaturePoints[[1]])), time = timeValues[!missingLandCoverOrElevation], data = as.data.frame(combinedData[!missingLandCoverOrElevation,]))
}

produceTestData <- function(indiaTemperatures, landCover, elevation, collectionDatesPOSIX, boundaryPolygon, satelliteNamesVec, dayIndex) {
  day28raster <- indiaTemperatures[[dayIndex]]
  raster::values(day28raster) <- replace(raster::values(day28raster), is.na(raster::values(day28raster)), -50)
  rasterCellsMidpoints <- raster::rasterToPoints(day28raster, spatial = TRUE)
  indiaValuesIndex <- sp::over(x = rasterCellsMidpoints, y = boundaryPolygon)
  pointsInIndia <- subset(rasterCellsMidpoints, subset = !is.na(indiaValuesIndex))
  missingIndianValueIndices <- pointsInIndia@data$layer == -50
  missingPoints <- subset(pointsInIndia, subset = missingIndianValueIndices)

  landCoverInMissingZones <- raster::extract(x = landCover, y = missingPoints)
  missingPointsNoWaterNoNA <- missingPoints[which(!is.na(landCoverInMissingZones) & !(landCoverInMissingZones == 0)), ]

  emptyRaster <- day28raster
  raster::values(emptyRaster) <- rep(NA, raster::ncell(emptyRaster))
  missingRaster <- raster::rasterize(x = missingPointsNoWaterNoNA, y = emptyRaster, field = "layer")

  testDataMay28 <- prepareDataForISMRA(landCover = landCover, elevations = elevation, temperatures = list(missingRaster), collectionDatesPOSIX = collectionDatesPOSIX[length(collectionDatesPOSIX) - 3], satelliteNamesVec = satelliteNamesVec, completeDateVector = collectionDatesPOSIX)
  testDataMay28
}

copySPDEresults <- function(originFolder, destinationFolder, patternForFilename) {
  filenames <- list.files(path = originFolder, pattern = patternForFilename, full.names = TRUE)
  filenamesWithoutFolder <- list.files(path = originFolder, pattern = patternForFilename, full.names = FALSE)
  dataIndices <- as.numeric(stringr::str_extract(filenames, pattern = "[:digit:]+(?=\\.rds)"))
  orderedFilenames <- filenames[order(dataIndices)]
  orderedFilenamesWithoutFolder <- filenamesWithoutFolder[order(dataIndices)]
  funToReplaceSPDEcomponent <- function(index) {
    oriSimResults <- readRDS(orderedFilenames[[index]])
    fileToLoad <- list.files(path = destinationFolder, pattern = orderedFilenamesWithoutFolder[[index]], full.names = TRUE)
    simResultsToUpdate <- readRDS(fileToLoad)
    simResultsToUpdate$SPDE <- oriSimResults$SPDE
    saveRDS(simResultsToUpdate, file = fileToLoad)
    cat("Updated", fileToLoad, "\n")
    NULL
  }
  lapply(sort(dataIndices), funToReplaceSPDEcomponent)
  cat("Done copying SPDE results.\n")
  invisible(NULL)
}
