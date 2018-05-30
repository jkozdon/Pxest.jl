run(`make`)

using Libdl

const libsc = "p4est/local/lib/libsc.dylib"
const libpxest = "p4est/local/lib/libp4est.dylib"

const libsc_handle =
Libdl.dlopen(libsc, Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)

const libpxest_handle =
Libdl.dlopen(libpxest, Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)

const pxest_qcoord_t = Int32
struct p4est_quadrant_t
  # coordinates
  x::pxest_qcoord_t
  y::pxest_qcoord_t

  level::Int8    # level of refinement
  pad8::Int8     # padding
  pad16::Int16    # padding
end

struct p8est_quadrant_t
  # coordinates
  x::pxest_qcoord_t
  y::pxest_qcoord_t
  z::pxest_qcoord_t

  level::Int8    # level of refinement
  pad8::Int8     # padding
  pad16::Int16    # padding
end

const pxest_topidx_t = Int32
struct pxest_iter_face_side_t
  treeid::pxest_topidx_t # the tree on this side
  face::Int8             # which quadrant side the face touches
  is_hanging::Int8       # boolean: one full quad (0) or two smaller quads (1)
end

struct pxest_iter_edge_side_t
  treeid::pxest_topidx_t # the tree on this side
  edge::Int8             # which quadrant side the edge touches
  orientation::Int8      # the orientation of each quadrant relative to this
                         # edge, as in the definition of p8est_connectivity_t
  is_hanging::Int8       # boolean: one full quad (0) or two smaller quads (1)
  faces::NTuple{2, Int8} # FIXME: Check that this is correct
end

let
  p4est_quadrant_t_size = ccall(("p4est_julia_sizeof_quadrant_t", libpxest),
                                Csize_t, ())
  p4est_iter_face_side_size = ccall(("p4est_julia_sizeof_iter_face_side_t",
                                     libpxest), Csize_t, ())
  p8est_quadrant_t_size = ccall(("p8est_julia_sizeof_quadrant_t", libpxest),
                                Csize_t, ())
  p8est_iter_face_side_size = ccall(("p8est_julia_sizeof_iter_face_side_t",
                                     libpxest), Csize_t, ())
  p8est_iter_edge_side_size = ccall(("p8est_julia_sizeof_iter_edge_side_t",
                                     libpxest), Csize_t, ())

  p4est_quadrant_piggy_size = p4est_quadrant_t_size - sizeof(p4est_quadrant_t)
  p4est_iter_face_side_is_size = p4est_iter_face_side_size -
                                 sizeof(pxest_iter_face_side_t)

  p8est_quadrant_piggy_size = p8est_quadrant_t_size - sizeof(p8est_quadrant_t)
  p8est_iter_face_side_is_size = p8est_iter_face_side_size -
                                 sizeof(pxest_iter_face_side_t)
  p8est_iter_edge_side_is_size = p8est_iter_edge_side_size -
                                 sizeof(pxest_iter_edge_side_t)

  open("compile-time.jl", "w") do io
    write(io, "@p4est const PXEST_QUADRANT_PIGGY_SIZE = ",
          string(p4est_quadrant_piggy_size), "\n")
    write(io, "@p4est const PXEST_ITER_FACE_SIDE_IS_SIZE = ",
          string(p4est_iter_face_side_is_size), "\n")
    write(io, "@p8est const PXEST_QUADRANT_PIGGY_SIZE = ",
          string(p8est_quadrant_piggy_size), "\n")
    write(io, "@p8est const PXEST_ITER_FACE_SIDE_IS_SIZE = ",
          string(p8est_iter_face_side_is_size), "\n")
    write(io, "@p8est const PXEST_ITER_EDGE_SIDE_IS_SIZE = ",
          string(p8est_iter_edge_side_is_size), "\n")
  end
end
