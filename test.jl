using MPI
if !MPI.Initialized()
  MPI.Init()
end

include("p4est.jl")
using p4est

let
  conn = p4est.Connectivity(5,7)
  pxest = p4est.PXEST(conn; min_lvl=2)
  p4est.balance(pxest)
  p4est.partition(pxest)
end

if !isinteractive()
  # Run gc to make sure cleanup happens before MPI is finalized
  gc()
  MPI.Finalize()
end
