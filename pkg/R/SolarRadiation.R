# Package dependencies ####

# R oxygen code for importing the proper classes used in this file
# Used for robust namespace management in R packages
#' @importFrom R6 R6Class
NULL

# R6 Class SolarRadiation ####

#' @export
#'
#' @title
#'   R6 class defining incoming solar radiation
#'
#' @description
#'   Provides the ability to make solar radiation calculations
#'   for a given location on Earth.
#'
SolarRadiation <- R6Class(
   classname = "SolarRadiation",
   public = list(
      
      ## Attributes ####
      
      #' @field latitudeAngle
      #'   Angle of the latitude
      latitudeAngle = NULL,
      
      #' @field longitudeAngle
      #'   Angle of the longitude
      longitudeAngle = NULL,
      
      #' @field solarConstantFlux
      #'   The solar constant for flux of radiation (for Earth)
      solarConstantFlux = NULL,
      
      #' @field earthAngularVelocity
      #'   The angular velocity of Earth's rotation on its axis
      earthAngularVelocity = NULL,
      
      ## Method: constructor ####

      #' @description
      #'   Construct a new instance of the class
      #'
      #' @param latitude
      #'   Latitude on earth in decimal degrees
      #' @param longitude
      #'   Longitude on earth in decimal degrees
      #' @param solarConstantFlux
      #'   Base flux of energy from the sun in Watts per square meter.
      #'   Default value is 1364.
      #' @param earthAngularVelocity
      #'   Angular velocity of earth's rotation in radians per hour.
      #'   Default value is 0.2618.
      #'
      initialize = function
      (
         latitude,
         longitude,
         solarConstantFlux = 1364,
         earthAngularVelocity = 0.2618
      )
      {
         self$latitudeAngle <- 2 * pi * (latitude / 360);
         self$longitudeAngle <- 2 * pi * (longitude / 360);
         self$solarConstantFlux <- solarConstantFlux;
         self$earthAngularVelocity <- earthAngularVelocity;
      },
      
      ## Method: getExtraterrestrialInsolation ####
      #
      #' @description
      #'   Calculates the incoming extraterrestrial solar radiation at
      #'   the times provided
      #'   (i.e. before any influence by atmosphere) for the location
      #'   on earth represented by the SolarRadiation object.
      #'
      #' @param timeCT
      #'   The time in POSIXct type.
      #' @param timeLT
      #'   The time in POSXlt type.
      #'   Default value is timeCT coerced into POSIXlt class.
      #' @param solarNoonCorrectionTime
      #'   The time basis for calculating solar noon corrections.
      #'   Defaults to a conversion to radians from a version of the
      #'   day angle.
      #' @param solarNoonCorrectionEcc
      #'   The correction for eccentricities in Earth's orbit.
      #'   Defaults to an empirical approximation from solarNoonCorrectionTime
      #' @param dayAngle
      #'   The angle representing the day of the year, if the year is represented
      #'   by 2 * pi radians
      #'   Default value is calculated based on the day of year form timeLT.
      #' @param eccCoefficient
      #'   The coefficient of correction for radiation based on eccentricities in
      #'   Earth's distance from the sun during its orbit.
      #'   Default value is an empirical calculation from the day angle.
      #' @param declinationAngle
      #'   The declination angle of the sun relative to Earth's equator.
      #'   Default value is an empirical calculation from the day angle.
      #'
      #' @return
      #'   The solar radiation flux at the top of Earth's atmosphere.
      #'   Units will be the same as those provided for the solar 
      #'   constant flux attribute.
      #'   Default solar constant flux in the constructor is in 
      #'   units of watts per square meter.
      #'
      getExtraterrestrialInsolation = function
      (
         timeCT,
         timeLT = as.POSIXlt(timeCT),
         solarNoonCorrectionTime =
            2 * pi * (timeLT$yday - 81) / 365,
         solarNoonCorrectionEcc =
            (
               9.87 * sin(2 * solarNoonCorrectionTime) -
                  7.53 * cos(solarNoonCorrectionTime) -
                  1.50 * sin(solarNoonCorrectionTime)
            ) /
            60,
         dayAngle = 2 * pi * (timeLT$yday / 365),
         eccCoefficient =
            1.000110 +
            0.034221 * cos(dayAngle) +
            0.001280 * sin(dayAngle) +
            0.000719 * cos(2 * dayAngle) +
            0.000077 * sin(2 * dayAngle),
         declinationAngle =
            0.006918 -
            0.399912 * cos(dayAngle) +
            0.070257 * sin(dayAngle) -
            0.006758 * cos(2 * dayAngle) +
            0.000907 * sin(2 * dayAngle) -
            0.002697 * cos(3 * dayAngle) +
            0.001480 * sin(3 * dayAngle)
      )
      {
         meridianAngle <-
            0.2617994 * (timeLT$gmtoff / 3600);
         
         solarNoonCorrection <-
            solarNoonCorrectionEcc +
            (self$longitudeAngle - meridianAngle) /
            self$earthAngularVelocity;
         
         solarTime <- timeCT + (solarNoonCorrection * 3600);
         
         solarNoon <- as.POSIXlt(solarTime);
         solarNoon$hour <- 12;
         solarNoon$min <- 0;
         solarNoon$sec <- 0;
         solarNoon <- as.POSIXct(solarNoon);
         
         timeAfterNoon <-
            (as.numeric(solarTime) - as.numeric(solarNoon)) / 3600;
         
         zenithAngle <-
            acos(
               sin(self$latitudeAngle) *
                  sin(declinationAngle) +
                  cos(self$latitudeAngle) *
                  cos(declinationAngle) *
                  cos(self$earthAngularVelocity * timeAfterNoon)
            );
         zenithAngle[zenithAngle > pi / 2] <- pi / 2;
         
         return(
            self$solarConstantFlux *
               eccCoefficient *
               cos(zenithAngle)
         );
      }
   )
)