\name{MakeInput_Fn}
\alias{MakeInput_Fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to generate inputs for TMB based on data that is typically available, or from the simulator function Sim_Fn()
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MakeInput_Fn(Options_vec, c_i, s_i, b_i, X_ij)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Options_vec}{
%%     ~~Describe \code{Options_vec} here~~
}
  \item{c_i}{
%%     ~~Describe \code{c_i} here~~
}
  \item{s_i}{
%%     ~~Describe \code{s_i} here~~
}
  \item{b_i}{
%%     ~~Describe \code{b_i} here~~
}
  \item{X_ij}{
%%     ~~Describe \code{X_ij} here~~
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
function (Options_vec, c_i, s_i, b_i, X_ij) 
{
    Data = list(Options_vec = Options_vec, n_s = max(s_i) + 1, 
        n_b = max(b_i) + 1, c_i = c_i, x_ij = X_ij, s_i = s_i, 
        b_i = b_i)
    Data[["spde"]] = list(n_s = MeshList$spde$n.spde, n_tri = nrow(MeshList$mesh$graph$tv), 
        Tri_Area = MeshList$Tri_Area, E0 = MeshList$E0, E1 = MeshList$E1, 
        E2 = MeshList$E2, TV = MeshList$TV - 1, G0 = MeshList$spde$param.inla$M0, 
        G0_inv = inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)))
    Params = list(ln_H_input = rep(0, 2), gamma_j = rep(0, ncol(Data$x_ij)), 
        logeta_vec = rep(0, 3), rho_vec = rep(0, 2), logkappa_vec = rep(0, 
            2), theta_vec = rep(0, 3), Nuinput_b = rnorm(Data$n_b), 
        Omegainput_sb = matrix(rnorm(Data$n_b * Data$spde$n_s), 
            nrow = Data$spde$n_s, ncol = Data$n_b), Deltainput_s = rnorm(Data$spde$n_s))
    Random = c("Nuinput_b", "Omegainput_sb", "Deltainput_s")
    Map = NULL
    Map[["theta_vec"]] = factor(c(1, NA, 2))
    InputList = list(Data = Data, Params = Params, Random = Random, 
        Map = Map)
    return(InputList)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
