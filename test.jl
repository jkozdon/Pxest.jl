using MPI
if !MPI.Initialized()
  MPI.Init()
  const mpicomm = MPI.COMM_WORLD
  const mpirank = MPI.Comm_rank(mpicomm)
end


include("p4est.jl")
using p4est
let
  conn = p4est.Connectivity(5,7)
  pxest = p4est.PXEST(conn; min_lvl=2)
  p4est.balance(pxest)
  p4est.partition(pxest)

  vtk_dir = "vtk_files"
  vtk_base = "mesh_p4est"
  mpirank == 0 ? mkpath(vtk_dir) : nothing
  MPI.Barrier(mpicomm)
  p4est.vtk_write_file(pxest, string(vtk_dir, "/", vtk_base))
  if mpirank == 0
    mv(string(vtk_dir, "/", vtk_base, ".pvtu"), string(vtk_base, ".pvtu"),
      remove_destination=true)
    mv(string(vtk_dir, "/", vtk_base, ".visit"), string(vtk_base, ".visit"),
      remove_destination=true)
  end
end

include("p8est.jl")
using p8est
let
  conn = p8est.Connectivity(5,7,2)
  pxest = p8est.PXEST(conn; min_lvl=1)
  p8est.balance(pxest)
  p8est.partition(pxest)

  vtk_dir = "vtk_files"
  vtk_base = "mesh_p8est"
  mpirank == 0 ? mkpath(vtk_dir) : nothing
  MPI.Barrier(mpicomm)
  p8est.vtk_write_file(pxest, string(vtk_dir, "/", vtk_base))
  if mpirank == 0
    mv(string(vtk_dir, "/", vtk_base, ".pvtu"), string(vtk_base, ".pvtu"),
      remove_destination=true)
    mv(string(vtk_dir, "/", vtk_base, ".visit"), string(vtk_base, ".visit"),
      remove_destination=true)
  end
end

if !isinteractive()
  # Run gc to make sure cleanup happens before MPI is finalized
  gc()
  MPI.Finalize()
end
