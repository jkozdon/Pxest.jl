using MPI
using Random

if !MPI.Initialized()
  MPI.Init()
  const mpicomm = MPI.COMM_WORLD
  const mpirank = MPI.Comm_rank(mpicomm)
  srand(mpirank)
end




using Pxest.p4est
let
  conn = p4est.Connectivity(5,7)
  pxest = p4est.PXEST(conn; min_lvl=0)
  p4est.refine!(pxest; maxlevel=3) do which_tree, quadrant
    if rand() > 0.9
      return Cint(1)
    else
      return Cint(0)
    end
  end
  p4est.balance!(pxest)
  p4est.partition!(pxest)
  p4est.ghost!(pxest)
  p4est.lnodes!(pxest, 2)
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

  #=
  # Count quadrants, faces, and corners separately
  q = 0
  p4est.quadrants(pxest) do quad
    q = q+1
  end
  f = 0
  p4est.faces(pxest) do face
    f = f+1
  end
  e = 0
  c = 0
  p4est.corners(pxest) do corn
    c = c+1
  end
  println((q, f, e, c))

  # Count quadrants, faces, and corners together
  q = 0
  f = 0
  c = 0
  p4est.iterator(pxest,
                 (quad,)->begin
                   q+=1
                 end,
                 (face,)->begin
                   f+=1
                 end,
                 (corn,)->begin
                   c+=1
                 end)
  println((q, f, c))
  =#
end

using Pxest.p8est
let
  conn = p8est.Connectivity(5,7,2)
  pxest = p8est.PXEST(conn; min_lvl=0)
  p8est.refine!(pxest; maxlevel=3) do which_tree, quadrant
    if rand() > 0.9
      return Cint(1)
    else
      return Cint(0)
    end
  end
  p8est.balance!(pxest)
  p8est.partition!(pxest)
  p8est.ghost!(pxest)
  p8est.lnodes!(pxest, 2)
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

  #=
  # Count quadrants, faces, edges, and corners separately
  q = 0
  p8est.quadrants(pxest) do quad
    q = q+1
  end
  f = 0
  p8est.faces(pxest) do face
    f = f+1
  end
  e = 0
  p8est.edges(pxest) do edge
    e = e+1
  end
  c = 0
  p8est.corners(pxest) do corn
    c = c+1
  end
  println((q, f, e, c))

  # Count quadrants, faces, edges, and corners together
  q = 0
  f = 0
  e = 0
  c = 0
  p8est.iterator(pxest,
                 (quad,)->begin
                   q+=1
                 end,
                 (face,)->begin
                   f+=1
                 end,
                 (edge,)->begin
                   e+=1
                 end,
                 (corn,)->begin
                   c+=1
                 end)
  println(length(pxest))
  =#
end

if !isinteractive()
  # Run gc to make sure cleanup happens before MPI is finalized
  GC.gc()
  MPI.Finalize()
end
