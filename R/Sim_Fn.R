Sim_Fn <-
function( n_bins=10, n_stations=20, SpatialScale=0.1, SD_omega=0.5, SD_nu=0.2, SD_delta=0.2, SD_extra=0.1, rho=0.8, logMeanDens=1, Loc=NULL ){

  # Spatial model
  require( RandomFields )
  if( is.null(Loc) ) Loc = cbind( "x"=runif(n_stations, min=0,max=1), "y"=runif(n_stations, min=0,max=1) )
  model_omega <- RMgauss(var=SD_omega^2, scale=SpatialScale)
  model_delta <- RMgauss(var=SD_delta^2, scale=SpatialScale)

  # Simulate Nu
  Nu_b = logMeanDens - (1-rho)*1:n_bins

  # Simulate Epsilon
  Delta_s = RFsimulate(model=model_delta, x=Loc[,'x'], y=Loc[,'y'])@data[,1]

  # Simulate Omega
  Omega_sb = matrix(NA, nrow=n_stations, ncol=n_bins)
  Omega_sb[,1] = RFsimulate(model=model_omega, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
  for(b in 2:n_bins){
    Omega_sb[,b] = rho*Omega_sb[,b-1] + RFsimulate(model=model_omega, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
  }

  # Simulate data
  log_chat_sb = outer(rep(1,n_stations),Nu_b) + outer(Delta_s,rep(1,n_bins)) + Omega_sb
  c_sb = matrix( rbinom(n_stations*n_bins,size=1,prob=1-exp(-1*exp(log_chat_sb))), nrow=n_stations, ncol=n_bins)
  c_sb = ifelse( c_sb==1, rlnorm(n_stations*n_bins,meanlog=log_chat_sb,sdlog=1), 0 )

  # return
  ReturnList = list("Loc"=Loc, "Nu_b"=Nu_b, "Delta_s"=Delta_s, "Omega_sb"=Omega_sb, "c_sb"=c_sb)
  return( ReturnList )
}
