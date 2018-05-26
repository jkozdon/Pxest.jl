using MPI
using Compat.Random
using Compat.GC

if !MPI.Initialized()
  MPI.Init()
  const mpicomm = MPI.COMM_WORLD
  const mpirank = MPI.Comm_rank(mpicomm)
  srand(mpirank)
end


function random_refinement(pxest, which_tree, quadrant)::Cint
  if rand() > 0.9
    return Cint(1)
  else
    return Cint(0)
  end
end


using Pxest.p4est
let
  conn = p4est.Connectivity(5,7)
  pxest = p4est.PXEST(conn; min_lvl=0)
  p4est.refine!(pxest, random_refinement, 3)
  p4est.balance!(pxest)
  p4est.partition!(pxest)
  p4est.ghost!(pxest)
  mesh = p4est.Mesh(pxest)

  # dump VTK
  vtk_dir = "vtk_files"
  vtk_base = "mesh_p4est"
  mpirank == 0 ? mkpath(vtk_dir) : nothing
  MPI.Barrier(mpicomm)
  p4est.vtk_write_file(pxest, string(vtk_dir, "/", vtk_base))
  if mpirank == 0
    mv(string(vtk_dir, "/", vtk_base, ".pvtu"), string(vtk_base, ".pvtu"),
      force=true)
    mv(string(vtk_dir, "/", vtk_base, ".visit"), string(vtk_base, ".visit"),
      force=true)
  end
end

using Pxest.p8est
let
  conn = p8est.Connectivity(5,7,2)
  pxest = p8est.PXEST(conn; min_lvl=0)
  p8est.refine!(pxest, random_refinement, 3)
  p8est.balance!(pxest)
  p8est.partition!(pxest)
  p8est.ghost!(pxest)
  mesh = p8est.Mesh(pxest)

  # dump VTK
  vtk_dir = "vtk_files"
  vtk_base = "mesh_p8est"
  mpirank == 0 ? mkpath(vtk_dir) : nothing
  MPI.Barrier(mpicomm)
  p8est.vtk_write_file(pxest, string(vtk_dir, "/", vtk_base))
  if mpirank == 0
    mv(string(vtk_dir, "/", vtk_base, ".pvtu"), string(vtk_base, ".pvtu"),
      force=true)
    mv(string(vtk_dir, "/", vtk_base, ".visit"), string(vtk_base, ".visit"),
      force=true)
  end
end

if !isinteractive()
  # Run gc to make sure cleanup happens before MPI is finalized
  GC.gc()
  MPI.Finalize()
end
