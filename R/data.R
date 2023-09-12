#' @title ballField
#'
#' @description A wiffle ball is dropped from different heights (second column) and the
#'              time it took to hit the ground is measured (first column).
#'
#' @format      Matrix or data frame with 63 rows that represent observations and two
#'              columns that represent the response variable (time) and experimental
#'              variable (height) in that order
#' \describe{
#'   \item{time}{The time in seconds for the ball to hit the ground}
#'   \item{height}{Height of the ball in meters before dropping}
#' }
#' @source <https://bobby.gramacy.com/surrogates/ball.csv>
"ballField"

#' @title ballSim
#'
#' @description The time to hit the ground after dropping the wiffle ball is simulated
#'              using a mathematical model and is implemented as a code. The data
#'              represent code response (time) at different heights and gravity.
#'
#' @format      Matrix or data frame with 100 rows that represent simulation runs and
#'              three columns that represent response variable (time), experimental
#'              input (height) and calibration input (gravity).
#' \describe{
#'   \item{time}{The time in seconds for the ball to hit the ground}
#'   \item{height}{height of the ball in meters before dropping}
#'   \item{gravity}{Height of the ball (in meters) before dropping}
#' }
"ballSim"

#' @title analytic11F
#'
#' @description A simple analytic example with random simulated data, known functional
#'              components, and specified error terms. It has two columns. The first
#'              column represent the response, and second column represent a single
#'              experimental input.
#' @format      Matrix with 5 rows that represent simulation runs and two columns that
#'              represent field response and a single experimental input.
"analytic11F"


#' @title analytic11S
#'
#' @description A simple analytic example with random simulated data, known functional
#'              components, and specified error terms. It has two columns. The first
#'              column represents the the simulated code response, second column
#'              represents a single experimental input, and the third column represents a
#'              single calibration input.
#'
#' @format      Matrix with 15 rows that represent simulation runs and three columns that
#'              represent code response, experimental input, and calibration input.
"analytic11S"


#' @title Df2
#'
#' @description A simple analytic example with random simulated data, known functional
#'              components, and specified error terms. It has three columns. The first
#'              column represent the field response, and the next two columns represent
#'              two experimental inputs.
#' @format      Matrix with 50 rows that represent simulation runs and three columns that
#'              represent field response and two experimental inputs.
"Df2"


#' @title Ds2
#'
#' @description A simple analytic example with random simulated data, known functional
#'              components, and specified error terms. It has five columns. The first
#'              column represents the the simulated code response, the next two columns
#'              (second and third) represent two experimental inputs, and the last two
#'              columns (fourth and fifth) represent two calibration inputs.
#'
#' @format      Matrix with 100 rows that represent simulation runs and five columns that
#'              represent code response, 2 experimental inputs, and 2 calibration inputs.
"Ds2"


#' @title swField
#'
#' @description A simple analytic example with random simulated data, known functional
#'              components, and specified error terms. It has three columns. The first
#'              column represent the field response, and the next two columns represent
#'              two experimental inputs.
#' @format      Matrix with 50 rows that represent simulation runs and three columns that
#'              represent field response and two experimental inputs.
"swField"


#' @title swSim
#'
#' @description A simple analytic example with random simulated data, known functional
#'              components, and specified error terms. It has five columns. The first
#'              column represents the the simulated code response, the next two columns
#'              (second and third) represent two experimental inputs, and the last two
#'              columns (fourth and fifth) represent two calibration inputs.
#'
#' @format      Matrix with 100 rows that represent simulation runs and five columns that
#'              represent code response, 2 experimental inputs, and 2 calibration inputs.
"swSim"


