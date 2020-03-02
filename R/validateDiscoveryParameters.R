.validateDiscoveryParameters <- function(projectPath,
                                         debugFlag,
                                         threadNum,
                                         seedValue,
                                         archDescriptionList) {
  validateParameters::validateArchDescriptionListParameter(archDescriptionList)
  validateParameters::validateDebugFlagParameter(debugFlag)
  validateParameters::validateProjectPathParameter(projectPath)
  validateParameters::validateSeedValueParameter(seedValue)
  validateParameters::validateThreadNumberParameter(threadNum)

  return()
}
