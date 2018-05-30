mutable struct Mesh
  K::Int32
  Kintra::Int32
  Kuniqmirror::Int32
  Kmirror::Int32
  EToL ::Array{Int8 , 1} # element to p4est level
  EToT ::Array{Int32, 1} # element to p4est treeid
  EToX ::Array{Int32, 1} # element to p4est x-qcoord
  EToY ::Array{Int32, 1} # element to p4est y-qcoord
  EToZ ::Array{Int32, 1} # element to p4est z-qcoord
  IToE ::Array{Int32, 1} # interior to element number
  UMToE::Array{Int32, 1} # unique mirror elements
  MToE ::Array{Int32, 1} # mirror elements arranged for MPI sends
  PToM ::Array{Int32, 1} # processors offsets for the mirrors (size: mpisize+1)

  EToB ::Array{Cint , 2} # element to boundary condition
  EToE ::Array{Cint , 2} # element to neighbor element
  EToF ::Array{Cint , 2} # element to neighbor face
  EToO ::Array{Cint , 2} # element to neighbor orientation
  EToP ::Array{Cint , 2} # element to periodicity mask (filled only for brick)

  function Mesh(pxest)
    mesh = new()

    ghost = unsafe_load(pxest.ghost)
    mirrors = unsafe_wrap(Array{pxest_quadrant_t, 1}, ghost.mirrors.array,
                          (Int(ghost.mirrors.elem_count),))
    mirror_lid = Dict{pxest_locidx_t, Bool}()
    for m in mirrors
      # +1 needed since julia counts from 1 not 0
      local_num = ccall(PXEST_JULIA_QUADRANT_P_PIGGY3_LOCAL_NUM,
                         pxest_locidx_t, (Ref{pxest_quadrant_t},), m) + 1
      mirror_lid[local_num] = true
    end

    mirror_proc_offsets = unsafe_wrap(Array{pxest_locidx_t, 1},
                                      ghost.mirror_proc_offsets,
                                      (ghost.mpisize+1,))
    mirror_proc_mirrors = unsafe_wrap(Array{pxest_locidx_t, 1},
                                      ghost.mirror_proc_mirrors,
                                      (mirror_proc_offsets[end],))
    PToM = mesh.PToM = zeros(Int32, ghost.mpisize)
    for p = 1:ghost.mpisize
      PToM[p] = mirror_proc_offsets[p] + 1
    end

    K           = mesh.K           = Int32(length(pxest))
    Kuniqmirror = mesh.Kuniqmirror = Int32(ghost.mirrors.elem_count)
    Kintra      = mesh.Kintra      = Int32(K - Kuniqmirror)
    Kmirror     = mesh.Kmirror     = Int32(PToM[end] - 1)
    @assert length(mirror_lid) == Kuniqmirror

    EToL = mesh.EToL = zeros(Int8, K)
    EToT = mesh.EToT = zeros(Int8, K)
    EToX = mesh.EToX = zeros(Int32, K)
    EToY = mesh.EToY = zeros(Int32, K)
    EToZ = mesh.EToZ = zeros(Int32, K)

    IToE  = mesh.IToE  = zeros(Int32, Kintra)
    UMToE = mesh.UMToE = zeros(Int32, Kuniqmirror)
    MToE  = mesh.MToE  = zeros(Int32, Kmirror)

    Ie  = 0
    UMe = 0
    quadrants(pxest) do quad
      qid = quadLID(quad)
      EToL[qid] = quadlevel(quad)
      EToX[qid] = quadX(quad)
      EToY[qid] = quadY(quad)
      @p8est EToZ[qid] = quadZ(quad)
      @p4est EToZ[qid] = 0
      if haskey(mirror_lid, qid)
        UMe += 1
        @assert UMe <= Kuniqmirror
        UMToE[UMe] = qid
      else
        Ie += 1
        @assert Ie <= Kintra
        IToE[Ie] = qid
      end
    end
    @assert Ie == Kintra
    @assert UMe == Kuniqmirror

    Me = 0
    for Me = 1:Kmirror
      MToE[Me] = UMToE[mirror_proc_mirrors[Me]+1]
    end

    faces(pxest) do face
      if length(face) == 1
        # Boundary
        side1 = face[1]
        println((sidetreeid(side1), sideface(side1), sideishanging(side1)))
      else
        side1 = face[1]
        side2 = face[2]
      end
    end

    mesh
  end
end
