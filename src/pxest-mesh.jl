mutable struct Mesh
  Klocal::Int32
  Kintra::Int32
  Kuniqmirror::Int32
  Kmirror::Int32
  Kghost::Int32
  Ktotal::Int32 # Klocal + Kghost
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
  #FIXME EToP ::Array{Bool , 2} # element to periodicity mask (filled only for brick)

  EToFC :: Array{Cint, 1} # element to p4est face code
  CToD_starts::Array{Cint, 1}  # start indices for each continuous node
  CToD_indices::Array{Cint, 1} # indices for continuous to discontinuous map
  DToC::Array{Int64, 2}  # start indices for each continuous node

  function Mesh(pxest; default_bc=1)
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

    Klocal      = mesh.Klocal      = Int32(length(pxest))
    Kuniqmirror = mesh.Kuniqmirror = Int32(ghost.mirrors.elem_count)
    Kintra      = mesh.Kintra      = Int32(Klocal - Kuniqmirror)
    Kmirror     = mesh.Kmirror     = Int32(PToM[end] - 1)
    Kghost      = mesh.Kghost      = Int32(ghost.ghosts.elem_count)
    @assert length(mirror_lid) == Kuniqmirror
    Ktotal      = mesh.Ktotal      = Int32(Klocal + Kghost)

    # elememt info
    EToL = mesh.EToL = zeros(Int8,  Klocal)
    EToT = mesh.EToT = zeros(Int8,  Klocal)
    EToX = mesh.EToX = zeros(Int32, Klocal)
    EToY = mesh.EToY = zeros(Int32, Klocal)
    EToZ = mesh.EToZ = zeros(Int32, Klocal)

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

    # Element to Element connectivity info
    EToB  = mesh.EToB = zeros(Cint, PXEST_FACES, Ktotal)
    EToE  = mesh.EToE = zeros(Cint, PXEST_FACES, Ktotal)
    EToF  = mesh.EToF = zeros(Cint, PXEST_FACES, Ktotal)
    EToO  = mesh.EToO = zeros(Cint, PXEST_FACES, Ktotal)
    #FIXME EToP  = mesh.EToP = zeros(Bool, PXEST_FACES, Ktotal or Klocal)
    faces(pxest) do face
      if length(face) == 1
        # Boundary
        side1 = face[1]
        f = sideface(side1)
        tid = sidetreeid(side1)
        (qid,) = sidequadid(side1)
        elm = pxest[qid, tid]
        EToB[f, elm] = default_bc
        EToE[f, elm] = elm
        EToF[f, elm] = f
        EToO[f, elm] = 0
        #FIXME EToP[f, elm] = false
      else
        side1 = face[1]
        side2 = face[2]
        o1 = o2 = faceorientation(face) # FIXME
        if !sideishanging(side1) && !sideishanging(side2)
          # Full Face
          (f1,) = sideface(side1)
          (g1,) = sideisghost(side1)
          (q1,) = sidequadid(side1)
          if !g1
            (t1,) = sidetreeid(side1)
            e1 = pxest[q1, t1]
          else
            e1 = q1 + Klocal
          end

          # Full Face
          (f2,) = sideface(side2)
          (g2,) = sideisghost(side2)
          (q2,) = sidequadid(side2)
          if !g2
            (t2,) = sidetreeid(side2)
            e2 = pxest[q2, t2]
          else
            e2 = q2 + Klocal
          end

          EToB[f1, e1] = -1
          EToE[f1, e1] = e2
          EToF[f1, e1] = f2
          EToO[f1, e1] = o2 # FIXME
          #FIXME EToP[f1, e1] =

          EToB[f2, e2] = -1
          EToE[f2, e2] = e1
          EToF[f2, e2] = f1
          EToO[f2, e2] = o1 # FIXME
          #FIXME EToP[f2, e2] =

        else # hanging face
          # make sure side2 is the hanging face
          sideishanging(side1) && ((side1, side2) = (side2, side1))
          @assert !sideishanging(side1) && sideishanging(side2)
        end
      end
    end

    # Continuous to discontinuous info
    if pxest.lnodes != C_NULL
      lnodes = unsafe_load(pxest.lnodes)
      EToFC = mesh.EToFC = zeros(Cint, Ktotal)

      #{{{ Get EToFc (based on p4est_plex.c)
      let
        # Get the face_code for each element
        for e = 1:Klocal
          EToFC[e] = unsafe_load(lnodes.face_code, e)
        end

        # Copy the face codes to a mirror array to send out to ghosts (note we
        # copy the address of the data to send not the data itself!)
        mirror_EToFC = Array{Ptr{Cint}}(undef, Kuniqmirror)
        for m = 1:Kuniqmirror
          mirror_EToFC[m] = pointer(EToFC, UMToE[m])
        end

        # Send the data out from the mirror array to the ghosts
        ccall(PXEST_GHOST_EXCHANGE_CUSTOM, Cvoid,
              (Ref{pxest_t}, Ptr{pxest_ghost_t}, Csize_t, Ref{Ptr{Cint}},
               Ptr{Cint}),
              pxest.pxest, pxest.ghost, sizeof(Cint), mirror_EToFC,
              pointer(EToFC, Klocal))
      end
      #}}}

      Np = lnodes.vnodes
      DToC = mesh.DToC = zeros(Int64, Np, Ktotal)
      owned_count = lnodes.owned_count
      global_offset = lnodes.global_offset
      #{{{ Get DtoC (based on p4est_plex.c)
      let
        # For all the local elements, get their nodes global continuous node
        # number
        for n = 1:(Np * Klocal)
          lidx = unsafe_load(lnodes.element_nodes, n) + 1
          DToC[n] = (lidx <= owned_count) ? (lidx + global_offset) :
          unsafe_load(lnodes.nonlocal_nodes, lidx - owned_count)
        end

        # Copy the discontinuous to continuous for my local elements to the
        # mirror out to ghosts. We will copy all of an elements numbers, so we
        # only copy the starting address for each element (since we copy the
        # address of the data to send not the data itself!)
        mirror_DToC = Array{Ptr{Int64}}(undef, Kuniqmirror)
        for m = 1:Kuniqmirror
          mirror_DToC[m] = pointer(DToC, (UMToE[m]-1) * Np + 1)
        end

        # Send the data out from the mirror array to the ghosts (i.e., figure
        # out my neighbors DToC map)
        ccall(PXEST_GHOST_EXCHANGE_CUSTOM, Cvoid,
              (Ref{pxest_t}, Ptr{pxest_ghost_t}, Csize_t, Ref{Ptr{Int64}},
               Ptr{Int64}),
              pxest.pxest, pxest.ghost, Np * sizeof(Int64), mirror_DToC,
              pointer(DToC, Klocal * Np + 1))
      end
      #}}}

      CToD = Dict{Tuple{Int64, Int8, Cint, Int32, Int32}, Int32}()

      # storage to mark that we have a continuous node locally
      hasC = Dict{Int64, Bool}()

      #{{{ Fill CToD Dictionary for local elements
      hanging = zeros(Int8, Np)
      for e = 1:Klocal
        # TODO: gethanging!(hanging, EToC[e])
        for n = 1:Np
          CToD[(DToC[n, e], hanging[n], pxest.pxest.mpirank, e, n)] = Np * (e-1) + n
          hasC[DToC[n, e]] = true
        end
      end
      #}}}

      #{{{ Fill CToD Dictionary for ghost
      let
        st = en = unsafe_load(ghost.proc_offsets, 1)+1
        # Go through all the ranks and grab the ghost elements they have
        for grank = 1:pxest.pxest.mpisize
          # my start is my neighbors en
          (st, en) = (en, unsafe_load(ghost.proc_offsets, grank+1)+1)
          for g = st:en-1
            # TODO: gethanging!(hanging, EToC[Klocal + g])
            for n = 1:Np
              if haskey(hasC, DToC[n, g])
                CToD[(DToC[n, g], hanging[n], grank, g, n)] = Np * (g-1) + n
              end
            end
          end
        end
      end
      #}}}

      ## TODO: FINISH THIS!
    end

    mesh
  end
end
