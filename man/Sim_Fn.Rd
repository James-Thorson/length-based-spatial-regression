\name{Sim_Fn}
\alias{Sim_Fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for simulating data to test estimation performance
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Sim_Fn(n_bins = 10, n_stations = 20, SpatialScale = 0.1, SD_omega = 0.5, SD_nu = 0.2, SD_delta = 0.2, SD_extra = 0.1, rho = 0.8, logMeanDens = 1, Loc = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n_bins}{
%%     ~~Describe \code{n_bins} here~~
}
  \item{n_stations}{
%%     ~~Describe \code{n_stations} here~~
}
  \item{SpatialScale}{
%%     ~~Describe \code{SpatialScale} here~~
}
  \item{SD_omega}{
%%     ~~Describe \code{SD_omega} here~~
}
  \item{SD_nu}{
%%     ~~Describe \code{SD_nu} here~~
}
  \item{SD_delta}{
%%     ~~Describe \code{SD_delta} here~~
}
  \item{SD_extra}{
%%     ~~Describe \code{SD_extra} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{logMeanDens}{
%%     ~~Describe \code{logMeanDens} here~~
}
  \item{Loc}{
%%     ~~Describe \code{Loc} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n_bins = 10, n_stations = 20, SpatialScale = 0.1, SD_omega = 0.5, 
    SD_nu = 0.2, SD_delta = 0.2, SD_extra = 0.1, rho = 0.8, logMeanDens = 1, 
    Loc = NULL) 
{
    require(RandomFields)
    if (is.null(Loc)) 
        Loc = cbind(x = runif(n_stations, min = 0, max = 1), 
            y = runif(n_stations, min = 0, max = 1))
    model_omega <- RMgauss(var = SD_omega^2, scale = SpatialScale)
    model_delta <- RMgauss(var = SD_delta^2, scale = SpatialScale)
    Nu_b = logMeanDens - (1 - rho) * 1:n_bins
    Delta_s = RFsimulate(model = model_delta, x = Loc[, "x"], 
        y = Loc[, "y"])@data[, 1]
    Omega_sb = matrix(NA, nrow = n_stations, ncol = n_bins)
    Omega_sb[, 1] = RFsimulate(model = model_omega, x = Loc[, 
        "x"], y = Loc[, "y"])@data[, 1]
    for (b in 2:n_bins) {
        Omega_sb[, b] = rho * Omega_sb[, b - 1] + RFsimulate(model = model_omega, 
            x = Loc[, "x"], y = Loc[, "y"])@data[, 1]
    }
    log_chat_sb = outer(rep(1, n_stations), Nu_b) + outer(Delta_s, 
        rep(1, n_bins)) + Omega_sb
    c_sb = matrix(rbinom(n_stations * n_bins, size = 1, prob = 1 - 
        exp(-1 * exp(log_chat_sb))), nrow = n_stations, ncol = n_bins)
    c_sb = ifelse(c_sb == 1, rlnorm(n_stations * n_bins, meanlog = log_chat_sb, 
        sdlog = 1), 0)
    ReturnList = list(Loc = Loc, Nu_b = Nu_b, Delta_s = Delta_s, 
        Omega_sb = Omega_sb, c_sb = c_sb)
    return(ReturnList)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
