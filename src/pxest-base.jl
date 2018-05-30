using MPI
using Libdl

# Path to the library
const libsc = joinpath(dirname(@__FILE__), "../deps/p4est/local/lib/libsc.dylib")
const libpxest = joinpath(dirname(@__FILE__), "../deps/p4est/local/lib/libp4est.dylib")

# Functions to load
@p4est const _pxest_functions = Dict{Symbol, String}(
    :PXEST_CONNECTIVITY_NEW_BRICK     => "p4est_connectivity_new_brick",
    :PXEST_CONNECTIVITY_READ_INP      => "p4est_connectivity_read_inp",
    :PXEST_CONNECTIVITY_NEW_BYNAME    => "p4est_connectivity_new_byname",
    :PXEST_CONNECTIVITY_DESTROY       => "p4est_connectivity_destroy",
    :PXEST_INIT_EXT                   => "p4est_init_ext",
    :PXEST_DESTROY_EXT                => "p4est_destroy_ext",
    :PXEST_BALANCE_EXT                => "p4est_balance_ext",
    :PXEST_PARTITION                  => "p4est_partition",
    :PXEST_GHOST_NEW                  => "p4est_ghost_new",
    :PXEST_GHOST_DESTROY              => "p4est_ghost_destroy",
    :PXEST_REFINE_EXT                 => "p4est_refine_ext",
    :PXEST_VTK_CONTEXT_NEW            => "p4est_vtk_context_new",
    :PXEST_VTK_CONTEXT_SET_SCALE      => "p4est_vtk_context_set_scale",
    :PXEST_VTK_CONTEXT_SET_CONTINUOUS => "p4est_vtk_context_set_continuous",
    :PXEST_VTK_WRITE_HEADER           => "p4est_vtk_write_header",
    :PXEST_VTK_WRITE_CELL_DATAF       => "p4est_vtk_write_cell_dataf",
    :PXEST_VTK_WRITE_FOOTER           => "p4est_vtk_write_footer",
    :PXEST_ITERATE                    => "p4est_iterate",
    :PXEST_JULIA_QUADRANT_P_WHICH_TREE =>
      "p4est_julia_quadrant_p_which_tree",
    :PXEST_JULIA_QUADRANT_P_PIGGY1_WHICH_TREE =>
      "p4est_julia_quadrant_p_piggy1_which_tree",
    :PXEST_JULIA_QUADRANT_P_PIGGY1_OWNER_RANK =>
      "p4est_julia_quadrant_p_piggy1_owner_rank",
    :PXEST_JULIA_QUADRANT_P_PIGGY2_WHICH_TREE =>
      "p4est_julia_quadrant_p_piggy2_which_tree",
    :PXEST_JULIA_QUADRANT_P_PIGGY2_FROM_TREE  =>
      "p4est_julia_quadrant_p_piggy2_from_tree",
    :PXEST_JULIA_QUADRANT_P_PIGGY3_WHICH_TREE =>
      "p4est_julia_quadrant_p_piggy3_which_tree",
    :PXEST_JULIA_QUADRANT_P_PIGGY3_LOCAL_NUM  =>
      "p4est_julia_quadrant_p_piggy3_local_num",
    #TODO: Use this to set piggy size
    :PXEST_JULIA_SIZEOF_QUADRANT_T    => "p4est_julia_sizeof_quadrant_t",
   )

@p8est const _pxest_functions = Dict{Symbol, String}(
    :PXEST_CONNECTIVITY_NEW_BRICK     => "p8est_connectivity_new_brick",
    :PXEST_CONNECTIVITY_READ_INP      => "p8est_connectivity_read_inp",
    :PXEST_CONNECTIVITY_NEW_BYNAME    => "p8est_connectivity_new_byname",
    :PXEST_CONNECTIVITY_DESTROY       => "p8est_connectivity_destroy",
    :PXEST_INIT_EXT                   => "p8est_init_ext",
    :PXEST_DESTROY_EXT                => "p8est_destroy_ext",
    :PXEST_BALANCE_EXT                => "p8est_balance_ext",
    :PXEST_PARTITION                  => "p8est_partition",
    :PXEST_GHOST_NEW                  => "p8est_ghost_new",
    :PXEST_GHOST_DESTROY              => "p8est_ghost_destroy",
    :PXEST_REFINE_EXT                 => "p8est_refine_ext",
    :PXEST_VTK_CONTEXT_NEW            => "p8est_vtk_context_new",
    :PXEST_VTK_CONTEXT_SET_SCALE      => "p8est_vtk_context_set_scale",
    :PXEST_VTK_CONTEXT_SET_CONTINUOUS => "p8est_vtk_context_set_continuous",
    :PXEST_VTK_WRITE_HEADER           => "p8est_vtk_write_header",
    :PXEST_VTK_WRITE_CELL_DATAF       => "p8est_vtk_write_cell_dataf",
    :PXEST_VTK_WRITE_FOOTER           => "p8est_vtk_write_footer",
    :PXEST_ITERATE                    => "p8est_iterate",
    :PXEST_JULIA_QUADRANT_P_WHICH_TREE =>
      "p8est_julia_quadrant_p_which_tree",
    :PXEST_JULIA_QUADRANT_P_PIGGY1_WHICH_TREE =>
      "p8est_julia_quadrant_p_piggy1_which_tree",
    :PXEST_JULIA_QUADRANT_P_PIGGY1_OWNER_RANK =>
      "p8est_julia_quadrant_p_piggy1_owner_rank",
    :PXEST_JULIA_QUADRANT_P_PIGGY2_WHICH_TREE =>
      "p8est_julia_quadrant_p_piggy2_which_tree",
    :PXEST_JULIA_QUADRANT_P_PIGGY2_FROM_TREE  =>
      "p8est_julia_quadrant_p_piggy2_from_tree",
    :PXEST_JULIA_QUADRANT_P_PIGGY3_WHICH_TREE =>
      "p8est_julia_quadrant_p_piggy3_which_tree",
    :PXEST_JULIA_QUADRANT_P_PIGGY3_LOCAL_NUM  =>
      "p8est_julia_quadrant_p_piggy3_local_num",
    #TODO: Use this to set piggy size
    :PXEST_JULIA_SIZEOF_QUADRANT_T    => "p8est_julia_sizeof_quadrant_t",
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

#{{{ p4est constants and types
@p4est const PXEST_CHILDREN = 4
@p8est const PXEST_CHILDREN = 8
@p4est const PXEST_MAXLEVEL = 30
@p8est const PXEST_MAXLEVEL = 19
const PXEST_ROOT_LEN = 1 << PXEST_MAXLEVEL

struct sc_array_t{T}
  # interface variables
  elem_size::Csize_t  # size of a single element
  elem_count::Csize_t # number of valid elements

  # implementation variables
  byte_alloc::Cssize_t # number of allocated bytes or
                      # -(number of viewed bytes + 1)
                      # if this is a view: the "+ 1"
                      # distinguishes an array of size 0
                      # from a view of size 0

  array::Ptr{T}       # linear array to store elements
end

# Connectivity data structures
const pxest_topidx_t = Int32
const pxest_locidx_t = Int32
const pxest_gloidx_t = Int64
const pxest_qcoord_t = Int32
const sc_mempool_t = Cvoid
const pxest_inspect_t = Cvoid
const pxest_connect_type_t = Cuint
@p4est const PXEST_CONNECT_FACE = pxest_connect_type_t(21)
@p4est const PXEST_CONNECT_CORNER = pxest_connect_type_t(22)

@p8est const PXEST_CONNECT_FACE = pxest_connect_type_t(31)
@p8est const PXEST_CONNECT_EDGE = pxest_connect_type_t(32)
@p8est const PXEST_CONNECT_CORNER = pxest_connect_type_t(33)

const PXEST_CONNECT_FULL = PXEST_CONNECT_CORNER
#}}}

#{{{ connectivity
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
  ccall(PXEST_CONNECTIVITY_DESTROY, Cvoid, (Ptr{Cvoid},), conn)
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
                           (3, Int(num_vertices)))
    tree_to_vertex = unsafe_wrap(Array{pxest_topidx_t, 2}, C.tree_to_vertex,
                                 (PXEST_CHILDREN, Int(num_trees)))
    this = new(num_vertices, num_trees, vertices, tree_to_vertex,
               pxest_conn_ptr)
    finalizer(connectivity_destroy, this)
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
#}}}

#{{{ pxest data structure
const piggy_size = max(sizeof(Ptr{Cvoid}),
                       sizeof(Clong),
                       sizeof(Cint),
                       sizeof(pxest_topidx_t),
                       sizeof(pxest_topidx_t) + sizeof(Cint),
                       sizeof(pxest_topidx_t) + sizeof(pxest_topidx_t),
                       sizeof(pxest_topidx_t) + sizeof(pxest_locidx_t))

struct pxest_quadrant_t
  # coordinates
  x::pxest_qcoord_t
  y::pxest_qcoord_t
  @p8est z::pxest_qcoord_t

  level::Int8    # level of refinement
  pad8::Int8     # padding
  pad16::Int16    # padding

  # TODO: Figure out how to handle this better (How to cast?)
  piggy_data::NTuple{piggy_size, Cchar}
  @p4est dummy::pxest_qcoord_t # FIXME: Is this correct???
end

struct pxest_ghost_t
  mpisize::Cint
  num_trees::pxest_topidx_t
  btype::pxest_connect_type_t # which neighbors are in the ghost layer

  #=
  # An array of quadrants which make up the ghost layer around \a
  # p4est.  Their piggy3 data member is filled with their owner's tree
  # and local number (cumulative over trees).  Quadrants are ordered in \c
  # p4est_quadrant_compare_piggy order.  These are quadrants inside the
  # neighboring tree, i.e., \c p4est_quadrant_is_inside() is true for the
  # quadrant and the neighboring tree.
  =#

  ghosts::sc_array_t{pxest_quadrant_t} # array of p4est_quadrant_t type
  tree_offsets::Ptr{pxest_locidx_t} # num_trees + 1 ghost indices
  proc_offsets::Ptr{pxest_locidx_t} # mpisize + 1 ghost indices

  #=
  # An array of local quadrants that touch the parallel boundary from the
  # inside, i.e., that are ghosts in the perspective of at least one other
  # processor.  The storage convention is the same as for \c ghosts above.
  =#
  mirrors::sc_array_t{pxest_quadrant_t} # array of p4est_quadrant_t type
  mirror_tree_offsets::Ptr{pxest_locidx_t} # num_trees + 1 mirror indices
  #=
  # indices into mirrors grouped by outside processor rank and ascending within
  # each rank
  =#
  mirror_proc_mirrors::Ptr{pxest_locidx_t}

  # mpisize + 1 indices into mirror_proc_mirrors
  mirror_proc_offsets::Ptr{pxest_locidx_t}

  #=
  # like mirror_proc_mirrors, but limited to the outermost octants.  This is
  # NULL until p4est_ghost_expand is called
  =#
  mirror_proc_fronts::Ptr{pxest_locidx_t}

  # NULL until p4est_ghost_expand is called
  mirror_proc_front_offsets::Ptr{pxest_locidx_t}
end


function quadX(quad::pxest_quadrant_t)
  quad.x
end
function quadY(quad::pxest_quadrant_t)
  quad.y
end
@p8est function quadZ(quad::pxest_quadrant_t)
  quad.z
end
function quadlevel(quad::pxest_quadrant_t)
  quad.level
end

struct pxest_tree_t
  quadrants::sc_array_t{pxest_quadrant_t} # locally stored quadrants
  first_desc::pxest_quadrant_t            # first local descendant */
  last_desc::pxest_quadrant_t             # last local descendant */

  # cumulative sum over earlier trees on this processor (locals only)
  quadrants_offset::pxest_locidx_t

  quadrants_per_level::NTuple{PXEST_MAXLEVEL+1, pxest_locidx_t} # locals only
  maxlevel::Int8                   # highest local quadrant level
end

mutable struct pxest_t
  # MPI communicator
  mpicomm::MPI.CComm

  # number of MPI processes
  mpisize::Cint

  # this process's MPI rank
  mpirank::Cint

  # flag if communicator is owned
  mpicomm_owned::Cint

  # size of per-quadrant p.user_data (see
  # p4est_quadrant_t::p4est_quadrant_data::user_data)
  data_size::Csize_t

  # convenience pointer for users, never touched by p4est
  user_pointer::Ptr{Cvoid}

  # Gets bumped on mesh change
  revision::Clong

  # 0-based index of first local tree, must be -1 for an empty processor
  first_local_tree::pxest_topidx_t

  # 0-based index of last local tree, must be -2 for an empty processor
  last_local_tree::pxest_topidx_t

  # number of quadrants on all trees on this processor
  local_num_quadrants::pxest_locidx_t

  # number of quadrants on all trees on all processors
  global_num_quadrants::pxest_gloidx_t

  # first global quadrant index for each process and 1 beyond
  global_first_quadrant::Ptr{pxest_gloidx_t}

  # first smallest possible quad for each process and 1 beyond
  global_first_position::Ptr{pxest_quadrant_t}

  # connectivity structure, not owned
  connectivity::Ptr{pxest_connectivity_t}

  # array of all trees
  trees::Ptr{sc_array_t{pxest_tree_t}}

  # memory allocator for user data
  #   WARNING: This is NULL if data size equals zero.
  user_data_pool::Ptr{sc_mempool_t}

  # memory allocator for temporary quadrants
  quadrant_pool::Ptr{sc_mempool_t}

  # algorithmic switches
  inspect::Ptr{pxest_inspect_t}

  function pxest_t()
    new()
  end
end

mutable struct PXEST
  pxest::pxest_t
  conn::Connectivity
  ghost::Ptr{pxest_ghost_t}

  function PXEST(conn ;mpicomm=MPI.COMM_WORLD, min_lvl = 0)
    pxest = pxest_t()
    ccall(PXEST_INIT_EXT, Cvoid,
          (Ref{pxest_t}, MPI.CComm, Ptr{pxest_connectivity_t}, pxest_locidx_t,
           Cint, Cint, Csize_t, Ptr{Cvoid}, Ptr{Cvoid}),
          pxest, MPI.CComm(mpicomm), conn.pxest_conn_ptr, 0, min_lvl, 1, 0,
          C_NULL, C_NULL)

    this = new(pxest, conn, C_NULL)
    finalizer(pxest_destroy, this)
    return this
  end
end

function Base.length(pxest::PXEST)
  return pxest.pxest.local_num_quadrants
end

function pxest_destroy(pxest)
  if pxest.ghost != C_NULL
    ghost_destroy!(pxest)
  end
  ccall(PXEST_DESTROY_EXT, Cvoid, (Ref{pxest_t}, Cint), pxest.pxest, 0)
end

function balance!(pxest; connect = PXEST_CONNECT_FULL)
  ccall(PXEST_BALANCE_EXT, Cvoid,
        (Ref{pxest_t}, pxest_connect_type_t, Ptr{Cvoid}, Ptr{Cvoid}),
        pxest.pxest, connect, C_NULL, C_NULL)
end

function partition!(pxest; allow_for_coarsening=true)
  ccall(PXEST_PARTITION, Cvoid, (Ref{pxest_t}, Cint , Ptr{Cvoid}),
        pxest.pxest, allow_for_coarsening, C_NULL)
end

function refine!(refine_fn::Function, pxest::PXEST; maxlevel=-1, refine_recursive=1)
  refine_call!(pxest, (x,y,z)->refine_fn(y,z), maxlevel, refine_recursive)
end

function refine!(pxest::PXEST, refine_fn::Function; maxlevel=-1, refine_recursive=1)
  refine_call!(pxest, (x,y,z)->refine_fn(y,z), maxlevel, refine_recursive)
end

function refine_call!(pxest::PXEST, refine_fn, maxlevel, refine_recursive)
  refine_fn_c = @cfunction($refine_fn, Cint, (Ref{pxest_t}, pxest_topidx_t,
                                                  Ref{pxest_quadrant_t}))
  ccall(PXEST_REFINE_EXT, Cvoid,
        (Ref{pxest_t}, Cint, Cint, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        pxest.pxest, refine_recursive, maxlevel, refine_fn_c, C_NULL,
        C_NULL)
end

function ghost!(pxest; connect = PXEST_CONNECT_FULL)
  pxest.ghost == C_NULL || error("ghost is not C_NULL")
  pxest.ghost = ccall(PXEST_GHOST_NEW, Ptr{pxest_ghost_t},
                      (Ref{pxest_t}, pxest_connect_type_t),
                      pxest.pxest, connect)
end

function ghost_destroy!(pxest; connect = PXEST_CONNECT_FULL)
  pxest.ghost != C_NULL || error("ghost is C_NULL")
  ccall(PXEST_GHOST_DESTROY, Cvoid, (Ptr{pxest_ghost_t},), pxest.ghost)
  pxest.ghost = C_NULL
end

struct pxest_iter_volume_info_t
  pxest::Ptr{pxest_t}
  ghost_layer::Ptr{pxest_ghost_t}
  quad::Ptr{pxest_quadrant_t} # the quadrant of the callback
  quadid::pxest_locidx_t      # id in \a quad's tree array (see p4est_tree_t)
  treeid::pxest_topidx_t      # the tree containing \a quad
end

function quadX(volinfo::pxest_iter_volume_info_t)
  quadX(unsafe_load(volinfo.quad))
end
function quadY(volinfo::pxest_iter_volume_info_t)
  quadY(unsafe_load(volinfo.quad))
end
@p8est function quadZ(volinfo::pxest_iter_volume_info_t)
  quadZ(unsafe_load(volinfo.quad))
end
function quadlevel(volinfo::pxest_iter_volume_info_t)
  quadlevel(unsafe_load(volinfo.quad))
end
function quadlevel(volinfo::pxest_iter_volume_info_t)
  volinfo.treeid
end
function quadLID(volinfo::pxest_iter_volume_info_t)
  trees = unsafe_load(unsafe_load(volinfo.pxest).trees)
  # FIXME: If this fails, then likely something wrong with the union in
  # pxest_quadrant_t
  @assert trees.elem_size == sizeof(pxest_tree_t)
  tree  = unsafe_load(trees.array, volinfo.treeid+1)
  tree.quadrants_offset + volinfo.quadid + 1
end

const pxest_iter_face_info_t = Cvoid
@p8est const pxest_iter_edge_info_t = Cvoid
const pxest_iter_corner_info_t = Cvoid

function quad_fn_wrapper(pxest_ptr, (pxest, quad_fn))::Cvoid
  quad_fn()
  nothing
end

function quadrants(quad_fn::Function, pxest)
  quadrants(pxest, quad_fn)
end

function quadrants(pxest, quad_fn::Function)
  @p8est iterator(pxest, quad_fn, C_NULL, C_NULL, C_NULL)
  @p4est iterator(pxest, quad_fn, C_NULL, C_NULL)
end

function faces(face_fn::Function, pxest)
  faces(pxest, face_fn)
end

function faces(pxest, face_fn::Function)
  @p8est iterator(pxest, C_NULL, face_fn, C_NULL, C_NULL)
  @p4est iterator(pxest, C_NULL, face_fn, C_NULL)
end

@p8est function edges(edge_fn::Function, pxest)
  edges(pxest, edge_fn)
end

@p8est function edges(pxest, edge_fn::Function)
  iterator(pxest, C_NULL, C_NULL, edge_fn, C_NULL)
end

function corners(corn_fn::Function, pxest)
  corners(pxest, corn_fn)
end

function corners(pxest, corn_fn::Function)
  @p8est iterator(pxest, C_NULL, C_NULL, C_NULL, corn_fn)
  @p4est iterator(pxest, C_NULL, C_NULL, corn_fn)
end

@p4est function iterator(pxest, quad_fn, face_fn, corn_fn)
  println(typeof(pxest))
  iterator_call(pxest,
                quad_fn == C_NULL ? C_NULL : (x, y)->begin
                  quad_fn(x)
                  nothing
                end,
                face_fn == C_NULL ? C_NULL : (x, y)->begin
                  face_fn(x)
                  nothing
                end,
                C_NULL,
                corn_fn == C_NULL ? C_NULL : (x, y)->begin
                  corn_fn(x)
                  nothing
                end)
end

@p8est function iterator(pxest, quad_fn, face_fn, edge_fn, corn_fn)
  iterator_call(pxest,
                quad_fn == C_NULL ? C_NULL : (x, y)->begin
                  quad_fn(x)
                  nothing
                end,
                face_fn == C_NULL ? C_NULL : (x, y)->begin
                  face_fn(x)
                  nothing
                end,
                edge_fn == C_NULL ? C_NULL : (x, y)->begin
                  edge_fn(x)
                  nothing
                end,
                corn_fn == C_NULL ? C_NULL : (x, y)->begin
                  corn_fn(x)
                  nothing
                end)
end

function iterator_call(pxest, quad_fn, face_fn, edge_fn, corn_fn)
  if quad_fn != C_NULL
    quad_fn_ptr = @cfunction($quad_fn, Cvoid, (Ref{pxest_iter_volume_info_t},
                                               Ptr{Cvoid}))
  else
    quad_fn_ptr = C_NULL
  end
  if face_fn != C_NULL
    face_fn_ptr = @cfunction($face_fn, Cvoid, (Ref{pxest_iter_face_info_t},
                                               Ptr{Cvoid}))
  else
    face_fn_ptr = C_NULL
  end
  @p8est if edge_fn != C_NULL
    edge_fn_ptr = @cfunction($edge_fn, Cvoid, (Ref{pxest_iter_edge_info_t},
                                               Ptr{Cvoid}))
  else
    edge_fn_ptr = C_NULL
  end
  if corn_fn != C_NULL
    corn_fn_ptr = @cfunction($corn_fn, Cvoid, (Ref{pxest_iter_corner_info_t},
                                               Ptr{Cvoid}))
  else
    corn_fn_ptr = C_NULL
  end

  @p4est ccall(PXEST_ITERATE, Cvoid, (Ref{pxest_t}, Ptr{pxest_ghost_t},
                                      Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid},
                                      Ptr{Cvoid}),
               pxest.pxest, C_NULL, C_NULL, quad_fn_ptr, face_fn_ptr,
               corn_fn_ptr)

  @p8est ccall(PXEST_ITERATE, Cvoid, (Ref{pxest_t}, Ptr{pxest_ghost_t},
                                      Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid},
                                      Ptr{Cvoid}, Ptr{Cvoid}),
               pxest.pxest, C_NULL, C_NULL, quad_fn_ptr, face_fn_ptr,
               edge_fn_ptr, corn_fn_ptr)
end



#}}}

#{{{vtk
function vtk_write_file(pxest, fn; scale = 1)
  cont = ccall(PXEST_VTK_CONTEXT_NEW, Ptr{Cvoid}, (Ref{pxest_t}, Ptr{Cchar},),
               pxest.pxest, fn)

  ccall(PXEST_VTK_CONTEXT_SET_SCALE, Cvoid, (Ptr{Cvoid}, Cdouble), cont, scale)

  ccall(PXEST_VTK_CONTEXT_SET_CONTINUOUS, Cvoid, (Ptr{Cvoid}, Cint), cont, 1)

  cont = ccall(PXEST_VTK_WRITE_HEADER, Ptr{Cvoid}, (Ptr{Cvoid},), cont)
  cont != C_NULL || error(string(fn, "_vtk: Error writing header"))

  cont = ccall(PXEST_VTK_WRITE_CELL_DATAF, Ptr{Cvoid},
               (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Cint, Cint, Ptr{Cvoid}...),
               cont, 1, 1, 1, 0, 0, 0, cont)
  cont != C_NULL || error(string(fn, "_vtk: Error writing cell data"))

  retval = ccall(PXEST_VTK_WRITE_FOOTER, Cint, (Ptr{Cvoid},), cont)
  retval == 0 || error(string(fn, "_vtk: Error writing footer"))
end
#}}}
