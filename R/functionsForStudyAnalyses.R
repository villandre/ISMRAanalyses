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

fitSPDE <- function(responseVec, covariateMatrix, coordinatesMatrix, timeVecNumeric, predCoordinatesMatrix, predCovariateMatrix, predTimeVecNumeric, numThreads = 1, control = list()) {
  control <- do.call("create.SPDE.control", control)

  timeVecTraining <- timeVecNumeric - min(timeVecNumeric) + 1
  spaceAndTimeMesh <- buildSpaceAndTimeMesh(coordinatesMatrixTraining = coordinatesMatrix, timeVecNumericTraining = timeVecTraining, control = control)

  inlaParameters <- produceINLAparameters(control)
  spde <- INLA::inla.spde2.matern(
    mesh = spaceAndTimeMesh$space,
    B.tau = matrix(c(inlaParameters$ltau0, -1, inlaParameters$spatialSmoothness), 1, 3),
    B.kappa = matrix(c(inlaParameters$lkappa0, 0, -1), 1, 3),
    theta.prior.mean = c(0,0), theta.prior.prec = c(1/control$loghyperparaSDinMyModel^2, 1/control$loghyperparaSDinMyModel^2))

  timeVecTest <- predTimeVecNumeric - min(timeVecNumeric) + 1
  combinedStack <- buildInlaStack(coordinatesMatrixTraining = coordinatesMatrix, timeVecTraining = timeVecTraining, coordinatesMatrixTest = predCoordinatesMatrix, timeVecTest = timeVecTest, meshForSpace = spaceAndTimeMesh$space, meshForTime = spaceAndTimeMesh$time, responseVecTraining = responseVec, covariateMatrixTraining = covariateMatrix, covariateMatrixTest = predCovariateMatrix, spdeObj = spde, control = control)

  formulaForSPDE <- y ~ -1 + elevation + May28 + May29 + EvergreenBroadleaf + MixedForest + ClosedShrublands + Savannas + Grasslands + PermanentWetlands + Croplands + CroplandNaturalMosaics + NonVegetated + f(space, model = spde, group = space.group, control.group = list(model = "ar1"))

  SPDEresult <- tryCatch(
    expr = INLA::inla(
      formulaForSPDE,
      data = INLA::inla.stack.data(combinedStack),
      control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(combinedStack)),
      num.threads = numThreads), error = function(e) e, finally = "Error in fitting SPDE! Return list will contain NAs.\n")
  returnResult <- predsAndSDs <- NULL
  if (!("simpleError" %in% class(SPDEresult))) {
    predsAndSDs <- getPredictionsAndSDsFromINLAoutputAlt(INLAoutput = SPDEresult, inlaStack = combinedStack, control = control)
    SPDEresult <- SPDEresult[grep(pattern = "summary", x = names(SPDEresult), value = TRUE)] # INLA objects can be huge. We only keep the elements we need.
    returnResult <- list(
      fittedModel = SPDEresult,
      predictionMeans = predsAndSDs$mean,
      predictionSDs = predsAndSDs$sd)
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
  mesh.2d.max.n = 2000,
  d = 1,
  alpha = 2,
  kappa = 1,
  loghyperparaSDinMyModel = log(10),
  sigma0 = 1) {  # = 1/(range parameter in my model))
  list(mesh.2d.cutoff = mesh.2d.cutoff, mesh.2d.offset = mesh.2d.offset, mesh.2d.max.n = mesh.2d.max.n, d = d, alpha = alpha, kappa = kappa, loghyperparaSDinMyModel = loghyperparaSDinMyModel, sigma0 = sigma0)
}

fitISMRA <- function(responseVec, coordinatesMatrix, predCoordinatesMatrix, covariateMatrix, predCovariateMatrix, timeVecNumeric, predTimeVecNumeric, numThreads = 1, control) {
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

  ISMRAfit <- tryCatch(expr = MRAINLA::INLAMRA(
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
  ), error = function(e) e, finally = cat("ISMRA analysis failed! Return list will contain NAs.\n"))
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
  MRAINLA::maternCov(distValues, smoothness = 1.5, rho = 1, scale = 1) * MRAINLA::maternCov(timeDistValues, smoothness = 0.5, rho = 1, scale = 1)
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

  covariateMatrixForTraining <- covariateMatrix[obsIndicesForTraining, ]
  predCovariateMatrix <- covariateMatrix[!obsIndicesForTraining, ]

  coordinatesMatrixForTraining <- coordinatesMatrix[obsIndicesForTraining, ]
  predCoordinatesMatrix <- coordinatesMatrix[!obsIndicesForTraining, ]

  timeVecNumericForTraining <- timeVecNumeric[obsIndicesForTraining]
  predTimeVecNumeric <- timeVecNumeric[!obsIndicesForTraining]

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

analysePredResults <- function(folderForSimResults, patternForFilename, simulatedDataList, obsIndicesForTraining) {
  patternForExtractingNumber <- "[:digit:]+(?=\\.rds)"
  filesToImport <- list.files(folderForSimResults, pattern = patternForFilename, full.names = TRUE)
  filesIndices <- as.numeric(stringr::str_extract(filesToImport, pattern = patternForExtractingNumber))
  filesToImportInOrder <- filesToImport[order(filesIndices)]

  computePredStats <- function(filename, obsIndicesForTraining, simulatedDataList) {
    datasetIndex <- as.numeric(stringr::str_extract(filename, pattern = patternForExtractingNumber))
    simResults <- readRDS(filename)
    realValues <- simulatedDataList$responses[!obsIndicesForTraining, datasetIndex]
    lapply(simResults, function(modelResult) {
      c(MSPE = mean((modelResult$predictionMeans - realValues)^2), MedSPE = median((modelResult$predictionMeans - realValues)^2), MeanSD = mean(modelResult$predictionSDs), MedSD = median(modelResult$predictionSDs))
    })
  }
  predStatsByDataset <- lapply(filesToImportInOrder, computePredStats, obsIndicesForTraining = obsIndicesForTraining, simulatedDataList = simulatedDataList)
  methodNames <- names(predStatsByDataset[[1]])
  names(methodNames) <- methodNames
  predStatsByMethod <- lapply(methodNames, function(methodName) {
    sapply(predStatsByDataset, "[[", methodName)
  })
  getStatsByMethod <- function(resultMatrix) {
    t(apply(resultMatrix, 1, function(rowValues) c(Mean = mean(rowValues), SD = sd(rowValues), Median = median(rowValues), Minimum = min(rowValues), Maximum = max(rowValues))))
  }
  lapply(predStatsByMethod, getStatsByMethod)
}

analyseParaEsts <- function(folderForSimResults, patternForFilename, simulatedDataList, obsIndicesForTraining, realFEs) {
  patternForExtractingNumber <- "[:digit:]+(?=\\.rds)"
  filesToImport <- list.files(folderForSimResults, pattern = patternForFilename, full.names = TRUE)
  filesIndices <- as.numeric(stringr::str_extract(filesToImport, pattern = patternForExtractingNumber))
  filesToImportInOrder <- filesToImport[order(filesIndices)]
  computeFEerrorsAndSDs <- function(filename) {
    datasetIndex <- as.numeric(stringr::str_extract(filename, pattern = patternForExtractingNumber))
    simResults <- readRDS(filename)
    methodNames <- names(simResults)
    names(methodNames) <- methodNames
    lapply(methodNames, function(methodName) {
      FEandSD <- getFEmeansAndSDs(simResults[[methodName]]$fittedModel)
      commonNames <- intersect(names(realFEs), rownames(FEandSD))
      reorderedFEandSD <- FEandSD[commonNames, ]
      absError <- abs(realFEs[commonNames] - reorderedFEandSD$Mean)
      SDs <- reorderedFEandSD$SD
      data.frame(absError = absError, SD = SDs)
    })
  }
  FEerrorsAndSDsByDataset <- lapply(filesToImportInOrder, FUN = computeFEerrorsAndSDs)
  methodNames <- names(FEerrorsAndSDsByDataset[[1]])
  names(methodNames) <- methodNames
  meanFEerrorsByMethod <- lapply(methodNames, function(methodName) {
    FEerrorsAcrossDataset <- sapply(seq_along(FEerrorsAndSDsByDataset), function(datasetIndex) {
      FEerrorsAndSDsByDataset[[datasetIndex]][[methodName]]$absError
    })
    meanResult <- rowMeans(FEerrorsAcrossDataset)
    names(meanResult) <- rownames(FEerrorsAndSDsByDataset[[1]][[methodName]])
    meanResult
  })
  meanFEsdsByMethod <- lapply(methodNames, function(methodName) {
    FEsdsByDataset <- sapply(seq_along(FEerrorsAndSDsByDataset), function(datasetIndex) {
      FEerrorsAndSDsByDataset[[datasetIndex]][[methodName]]$SD
    })
    meanResult <- rowMeans(FEsdsByDataset)
    names(meanResult) <- rownames(FEerrorsAndSDsByDataset[[1]][[methodName]])
    meanResult
  })

  AbsErrorAndSDbyMethod <- lapply(methodNames, function(methodName) {
    data.frame(AbsError = meanFEerrorsByMethod[[methodName]], EstSD = meanFEsdsByMethod[[methodName]])
  })
  AbsErrorAndSDbyMethod
}

getFEmeansAndSDs <- function(outputObject) {
  output <- NULL
  if ("INLAMRA" %in% class(outputObject)) {
    output <- outputObject$FEmarginalMoments[ , c("Mean", "StdDev")]
  } else if (!is.null(outputObject$summary.fixed)) {
    output <- outputObject$summary.fixed[, c("mean", "sd")]
  }
  colnames(output) <- c("Mean", "SD")
  output
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
    B.kappa = matrix(c(inlaParameters$lkappa0, 0, -1), 1,3),
    theta.prior.mean = c(0,0), theta.prior.prec = c(1/control$loghyperparaSDinMyModel^2, 1/control$loghyperparaSDinMyModel^2))

  combinedStack <- buildInlaStack(coordinatesMatrixTraining = coordinatesMatrixTraining, timeVecTraining = timeVecNumericTraining, coordinatesMatrixTest = coordinatesMatrixTest, timeVecTest = timeVecNumericTest, meshForSpace = spaceAndTimeMesh$space, meshForTime = spaceAndTimeMesh$time, responseVecTraining = responseVecTraining, covariateMatrixTraining = covariateMatrixTraining, covariateMatrixTest = covariateMatrixTest, spdeObj = spde, control = control)

  preds <- INLAoutput$summary.linear.predictor
  stackIndex <- INLA::inla.stack.index(combinedStack, "predictions")$data
  preds[stackIndex, c("mean", "sd")]
}

buildSpaceAndTimeMesh <- function(coordinatesMatrixTraining, timeVecNumericTraining, control) {

  knots <- seq(1, max(timeVecNumericTraining), length = max(timeVecNumericTraining))
  meshTime <- INLA::inla.mesh.1d(loc = knots, degree = 2, boundary = "free")

  ## generate space mesh

  meshSpace <- INLA::inla.mesh.2d(loc = coordinatesMatrixTraining[timeVecNumericTraining == 1, ], cutoff = control$mesh.2d.cutoff, offset = control$mesh.2d.offset, max.n = control$mesh.2d.max.n)
  list(time = meshTime, space = meshSpace)
}

buildInlaStack <- function(coordinatesMatrixTraining, timeVecTraining, coordinatesMatrixTest, timeVecTest, meshForSpace, meshForTime, responseVecTraining, covariateMatrixTraining, covariateMatrixTest, spdeObj, control) {
  ## build the space time indices
  STindex <- INLA::inla.spde.make.index("space", n.spde = spdeObj$n.spde, n.group = meshForTime$m)

  Atraining <- INLA::inla.spde.make.A(meshForSpace, loc = coordinatesMatrixTraining, group = timeVecTraining, group.mesh = meshForTime)
  Atest <- INLA::inla.spde.make.A(meshForSpace, loc = coordinatesMatrixTest, group = timeVecTest, group.mesh = meshForTime)
  covariateFrame <- as.data.frame(covariateMatrixTraining)
  predCovariateFrame <- as.data.frame(covariateMatrixTest)

  stackTraining <- INLA::inla.stack(
    data = list(y = responseVecTraining),
    A = c(list(Atraining), lapply(rep(1, ncol(covariateFrame)), identity)),
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

getPredictionsAndSDsFromINLAoutputAlt <- function(INLAoutput, inlaStack, control) {
  preds <- INLAoutput$summary.linear.predictor
  stackIndex <- INLA::inla.stack.index(inlaStack, "predictions")$data
  preds[stackIndex, c("mean", "sd")]
}

# This will update saved simulation results with predictedValues,
# which were missing due to a bug.
# The new version of the software should not require this function.

recomputePredictionsForSimOutputs <- function(folderForSimResults, patternForFilename = "Dataset.+\\.rds$", coordinatesMatrixTraining, coordinatesMatrixTest, timeVecNumericTraining, timeVecNumericTest, covariateMatrixTraining, covariateMatrixTest, responseMatrixTraining, controlForSPDE) {
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
    saveRDS(simResultsUpdate, file = filename, compress = TRUE)
    cat("Updated predictions in ", filename, "\n")
    invisible(NULL)
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
  numThreads,
  obsIndicesForTraining) {
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

produceINLAparameters <- function(control) {
  # range0 and sigma0 control the prior means for the range and scale parameters.
  # See Lindgren INLA tutorial page 5.
  spatialSmoothness <- control$alpha - control$d/2 # cf p.3 INLA tutorial

  # range0 and sigma0 seem to be the prior means...
  range0 <- sqrt(8 * spatialSmoothness)/control$kappa # sqrt(8 * spatial smoothness) / Kappa. In my model, I use 1 as prior mean for spatial range and fix smoothness at 1.5. This means Kappa = 1.
  lkappa0 <- log(8 * spatialSmoothness)/2 - log(range0)
  ltau0 <- 0.5*log(gamma(spatialSmoothness)/(gamma(control$alpha)*(4*pi)^(control$d/2))) - log(control$sigma0) - spatialSmoothness * lkappa0
  list(spatialSmoothness = spatialSmoothness, lkappa0 = lkappa0, ltau0 = ltau0, range0 = range0)
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
