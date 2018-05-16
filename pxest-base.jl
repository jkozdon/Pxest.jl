using MPI

# Path to the library
const libsc =
  "/Users/jekozdon/codes/p4est_julia/p4est/local/lib/libsc.dylib"
const libpxest =
  "/Users/jekozdon/codes/p4est_julia/p4est/local/lib/libp4est.dylib"

# Functions to load
@p4est const _pxest_functions = Dict{Symbol, String}(
    :PXEST_CONNECTIVITY_NEW_BRICK  => "p4est_connectivity_new_brick",
    :PXEST_CONNECTIVITY_READ_INP   => "p4est_connectivity_read_inp",
    :PXEST_CONNECTIVITY_NEW_BYNAME => "p4est_connectivity_new_byname",
    :PXEST_CONNECTIVITY_DESTROY    => "p4est_connectivity_destroy",
    :PXEST_NEW_EXT                 => "p4est_new_ext",
    :PXEST_DESTROY                 => "p4est_destroy",
   )

@p8est const _pxest_functions = Dict{Symbol, String}(
    :PXEST_CONNECTIVITY_NEW_BRICK  => "p8est_connectivity_new_brick",
    :PXEST_CONNECTIVITY_READ_INP   => "p8est_connectivity_read_inp",
    :PXEST_CONNECTIVITY_NEW_BYNAME => "p8est_connectivity_new_byname",
    :PXEST_CONNECTIVITY_DESTROY    => "p8est_connectivity_destroy",
    :PXEST_NEW_EXT                 => "p8est_new_ext",
    :PXEST_DESTROY                 => "p8est_destroy",
   )

# Build symbols
function __init__()

  @eval const libsc_handle =
    Libdl.dlopen(libsc, Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)

  @eval const libpxest_handle =
    Libdl.dlopen(libpxest, Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)

    # look up all symbols ahead of time
    for (jname, fname) in _pxest_functions
      eval(:(const $jname = Libdl.dlsym(libpxest_handle, $fname)))
    end

end

#
@p4est const PXEST_CHILDREN = 4
@p8est const PXEST_CHILDREN = 8

# Connectivity data structures
const pxest_topidx_t = Int32
const pxest_locidx_t = Int32
struct pxest_connectivity_t
  # the number of vertices that define the \a embedding of the forest (not the
  # topology)
  num_vertices::pxest_topidx_t

  # the number of trees
  num_trees::pxest_topidx_t

  # the number of edges that help define the topology
  @p8est num_edges::pxest_topidx_t

  # the number of corners that help define topology
  num_corners::pxest_topidx_t

  # an array of size (3 * \a num_vertices)
  vertices::Ptr{Cdouble}

  # embed each tree into \f$R^3\f$ for e.g. visualization (see p4est_vtk.h)
  tree_to_vertex::Ptr{pxest_topidx_t}

  # bytes per tree in tree_to_attr
  tree_attr_bytes::Csize_t

  # not touched by p4est
  tree_to_attr::Ptr{Cchar}

  # (4 * \a num_trees) neighbors across faces
  tree_to_tree::Ptr{pxest_topidx_t}

  # (4 * \a num_trees) face to face+orientation (see description)
  tree_to_face::Ptr{Int8}

  #(12 * \a num_trees) or NULL (see description) */
  @p8est tree_to_edge::Ptr{pxest_topidx_t}

  #edge to offset in \a edge_to_tree and \a edge_to_edge */
  @p8est ett_offset::Ptr{pxest_topidx_t}

  #list of trees that meet at an edge */
  @p8est edge_to_tree::Ptr{pxest_topidx_t}

  #list of tree-edges+orientations that meet at an edge (see description) */
  @p8est edge_to_edge::Ptr{Int8}

  #(4 * \a num_trees) or NULL (see description)
  tree_to_corner::Ptr{pxest_topidx_t}

  #corner to offset in \a corner_to_tree and \a corner_to_corner
  ctt_offset::Ptr{pxest_topidx_t}

  #list of trees that meet at a corner
  corner_to_tree::Ptr{pxest_topidx_t}

  #list of tree-corners that meet at a corner
  corner_to_corner::Ptr{Int8}
end

function pxest_connectivity_read_inp(filename)
  ccall(PXEST_CONNECTIVITY_READ_INP, Ptr{pxest_connectivity_t},
        (Ptr{Cchar},), filename)
end

function pxest_connectivity_new_byname(name)
  ccall(PXEST_CONNECTIVITY_NEW_BYNAME, Ptr{pxest_connectivity_t},
        (Ptr{Cchar},), name)
end

@p4est function pxest_connectivity_new_brick(nx, ny; px=false, py=false)
  ccall(PXEST_CONNECTIVITY_NEW_BRICK, Ptr{pxest_connectivity_t},
        (Cint, Cint, Cint, Cint), nx, ny, px, py)
end

@p8est function pxest_connectivity_new_brick(nx, ny, nz;
                                             px=false, py=false, pz=false)
  ccall(PXEST_CONNECTIVITY_NEW_BRICK, Ptr{pxest_connectivity_t},
        (Cint, Cint, Cint, Cint, Cint, Cint), nx, ny, nz, px, py, pz)
end

function pxest_connectivity_destroy(conn)
  ccall(PXEST_CONNECTIVITY_DESTROY, Void, (Ptr{Void},), conn)
end

mutable struct Connectivity
  # the number of vertices that define the \a embedding of the forest (not the
  # topology)
  num_vertices::pxest_topidx_t

  # the number of trees
  num_trees::pxest_topidx_t

  # an array of size (3 * \a num_vertices)
  vertices::Array{Cdouble, 2}

  # embed each tree into \f$R^3\f$ for e.g. visualization (see p4est_vtk.h)
  tree_to_vertex::Array{pxest_topidx_t, 2}

  pxest_conn_ptr::Ptr{pxest_connectivity_t}

  @p4est function Connectivity(nx, ny; px = false, py = false)
    Connectivity(pxest_connectivity_new_brick(nx, ny; px=px, py=py))
  end
  @p8est function Connectivity(nx, ny, nz; px = false, py = false, pz = false)
    Connectivity(pxest_connectivity_new_brick(nx, ny, nz; px=px, py=py, pz=pz))
  end
  function Connectivity(filename)
    Connectivity(isfile(filename) ?
                 pxest_connectivity_read_inp(filename) :
                 pxest_connectivity_new_byname(filename))
  end

  function Connectivity(pxest_conn_ptr::Ptr{pxest_connectivity_t})
    C = unsafe_load(Ptr{pxest_connectivity_t}(pxest_conn_ptr))
    num_vertices = C.num_vertices
    num_trees = C.num_trees
    vertices = unsafe_wrap(Array{Cdouble, 2}, C.vertices,
                           (3, Int(num_vertices)), false)
    tree_to_vertex = unsafe_wrap(Array{pxest_topidx_t, 2}, C.tree_to_vertex,
                                 (PXEST_CHILDREN, Int(num_trees)), false)
    this = new(num_vertices, num_trees, vertices, tree_to_vertex,
               pxest_conn_ptr)
    finalizer(this, connectivity_destroy)
    return this
  end
end

function connectivity_destroy(conn)
  pxest_connectivity_destroy(conn.pxest_conn_ptr)
  conn.pxest_conn_ptr = C_NULL
  conn.num_vertices = 0
  conn.num_trees = 0
  conn.vertices = zeros(Cdouble, 0,0)
  conn.tree_to_vertex = zeros(pxest_topidx_t, 0,0)
end

struct pxest_t
end

mutable struct PXEST
  pxest_ptr::Ptr{Void}

  function PXEST(conn ;mpicomm=MPI.COMM_WORLD, min_lvl = 0)
    pxest = ccall(PXEST_NEW_EXT, Ptr{Void},
                  (MPI.CComm, Ptr{pxest_connectivity_t}, pxest_locidx_t,
                   Cint, Cint, Csize_t, Ptr{Void}, Ptr{Void}),
                  MPI.CComm(mpicomm), conn.pxest_conn_ptr, 0, min_lvl, 1, 0,
                  C_NULL, C_NULL)

    this = new(pxest)
    finalizer(this, pxest_destroy)
    return this
  end
end

function pxest_destroy(pxest)
  ccall(PXEST_DESTROY, Void, (Ptr{Void},), pxest.pxest_ptr)
  pxest.pxest_ptr = C_NULL
end
