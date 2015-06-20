MakeInput_Fn <-
function( Options_vec, c_i, s_i, b_i, X_ij ){

  # Data
  Data = list("Options_vec"=Options_vec, "n_s"=max(s_i)+1, "n_b"=max(b_i)+1, "c_i"=c_i, "x_ij"=X_ij, "s_i"=s_i, "b_i"=b_i)
  Data[['spde']] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )

  # Parameters
  Params = list("ln_H_input"=rep(0,2), "gamma_j"=rep(0,ncol(Data$x_ij)), "logeta_vec"=rep(0,3), "rho_vec"=rep(0,2), "logkappa_vec"=rep(0,2), "theta_vec"=rep(0,3), "Nuinput_b"=rnorm(Data$n_b), "Omegainput_sb"=matrix(rnorm(Data$n_b*Data$spde$n_s),nrow=Data$spde$n_s,ncol=Data$n_b), "Deltainput_s"=rnorm(Data$spde$n_s))

  # Random
  Random = c("Nuinput_b", "Omegainput_sb", "Deltainput_s")

  # Map
  Map = NULL
  Map[["theta_vec"]] = factor( c(1,NA,2) )

  # Return
  InputList = list( "Data"=Data, "Params"=Params, "Random"=Random, "Map"=Map)
  return(InputList)
}
