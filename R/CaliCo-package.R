#' Bayesian calibration for computational codes
#'
#' The CaliCo package provides three categories of important functions:
#' \code{\link{model}}, \code{\link{prior}}, \code{\link{calibrate}} and \code{\link{prediction}}.
#'
#' @useDynLib CaliCo
#' @importFrom R6 R6Class
#' @importFrom stats rnorm
#' @import ggplot2 DiceKriging DiceDesign FactoMineR coda parallel testthat MASS
#'
#' @details
#' Package: CaliCo
#'
#' Type:    Package
#'
#' Version: 0.1.0
#'
#' Date:    2017-11-07
#'
#' License: GPL-2 | GPL-3
#'
#' @docType package
#' @author Mathieu Carmassi
#' @author Maintainer: \email{mathieu.carmassi@gmail.com}
#' @references Bachoc, F., Blois, G., Garnie, J., and Martinez, J.-M. (2014). Calibration and improved prediction of computer models
#' by universal kriging. Computational Statistics and Data Analysis, pages 81–97
#' @references Bayarri, M., Berger, J., Sacks, P. R., Cafeo, J. A., Cavendish, J., Lin, C. H., and Tu, J. (2007 b). A framework for
#' validation of computer models. Technometrics.
#' @references Carmassi, M., Barbillon ,P., Chiodetti, M., Keller, M., Parent, E. (2018). Bayesian calibration of a numerical code for prediction,
#' arXiv preprint arXiv:1801.01810.
#' @references Cox, D., Park, J. S., and Singer, C. (2001). A statistical method for tuning a computer code to a data base. Computational
#' Statistics and Data Analysis.
#' @references Damblin, G. (2015). Contributions statistiques au calage et à la validation des codes de calculs. PhD thesis, University
#' Paris-Saclay
#' @references Hastings, W. K. (1970). Mont carlo sampling methods using marlov chains and their applications. Biometrika.
#' @references Higdon, D., Kennedy, M. C., Cavendish, J., Cafeo, J., and Ryne, R. (2004). Combining field data and computer
#' simulations for calibration and prediction. SIAM Journal on Scientific Computing.
#' @references Kennedy, M. C. and O’Hagan, A. (2001). Bayesian calibration of computer models. Journal of the Royal Statistical
#' Society, serie B, Methodological.
#' @references Kennedy, M. C. and O’Hagan, A. (2001b). Supplementary details on bayesian calibration of computer models. Journal
#' of the Royal Statistical Society, serie B, Methodological.
#' @references Liu, F., Bayarri, S., and Berger, J. (2009). Modularization in bayesian analysis, with emphasis on analysis of computer
#' models. Bayesian Analysis, pages 119–150.
#' @references Robert, C. (1996). Méthodes de monte carlo par chaines de markov. economica.
#' @references Roustant, O., Ginsbourger, D., and Devills, Y. (2012). Dicekriging, diceoptim : Two r packages for the analysis of
#' computer experiments by kriging-based metamodeling and optimization. Journal of Statistical Software.
#' @references Sacks, J., Welch, W. J., and Toby J. Mitchell, H. P. W. (1989). Design and analysis of computer experiments. Statistical
#' science, pages 409–423.
#' @references Santner, T., Williams, B., and Notz, W. (2003). The Design and Analysis of Computer Experiments. Springer-Verlab.
#' @examples
#' # Introduction to CaliCo
#' \dontrun{vignette("CaliCo-introduction")}
#' @name CaliCo
NULL
