using MPI
using Libdl

include(joinpath(dirname(@__FILE__), "../deps/compile-time.jl"))

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
    :PXEST_GHOST_EXCHANGE_CUSTOM      => "p4est_ghost_exchange_custom",
    :PXEST_LNODES_NEW                 => "p4est_lnodes_new",
    :PXEST_LNODES_DESTROY             => "p4est_lnodes_destroy",
    :PXEST_GHOST_SUPPORT_LNODES       => "p4est_ghost_support_lnodes",
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
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_FULL_IS_GHOST =>
      "p4est_julia_iter_face_side_t_is_full_is_ghost",
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_FULL_QUAD =>
      "p4est_julia_iter_face_side_t_is_full_quad",
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_FULL_QUADID =>
      "p4est_julia_iter_face_side_t_is_full_quadid",
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_IS_GHOST =>
      "p4est_julia_iter_face_side_t_is_hanging_is_ghost",
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUAD =>
      "p4est_julia_iter_face_side_t_is_hanging_quad",
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUADID =>
      "p4est_julia_iter_face_side_t_is_hanging_quadid",
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
    :PXEST_GHOST_EXCHANGE_CUSTOM      => "p8est_ghost_exchange_custom",
    :PXEST_LNODES_NEW                 => "p8est_lnodes_new",
    :PXEST_LNODES_DESTROY             => "p8est_lnodes_destroy",
    :PXEST_GHOST_SUPPORT_LNODES       => "p8est_ghost_support_lnodes",
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
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_FULL_IS_GHOST =>
      "p8est_julia_iter_face_side_t_is_full_is_ghost",
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_FULL_QUAD =>
      "p8est_julia_iter_face_side_t_is_full_quad",
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_FULL_QUADID =>
      "p8est_julia_iter_face_side_t_is_full_quadid",
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_IS_GHOST =>
      "p8est_julia_iter_face_side_t_is_hanging_is_ghost",
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUAD =>
      "p8est_julia_iter_face_side_t_is_hanging_quad",
    :PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUADID =>
      "p8est_julia_iter_face_side_t_is_hanging_quadid",
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
const PXEST_HALF = div(PXEST_CHILDREN, 2)
@p4est const PXEST_MAXLEVEL = 30
@p8est const PXEST_MAXLEVEL = 19
const PXEST_ROOT_LEN = 1 << PXEST_MAXLEVEL
@p4est const PXEST_FACES = 4
@p8est const PXEST_FACES = 6

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

# Since no guarantee on GC order, we keep track of references to the
# connectivity pointer and only free on the last reference
const pxest_conn_ptr_refs = Dict{Ptr{pxest_connectivity_t}, Int}()

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

function pxest_connectivity_destroy(conn_ptr)
  if pxest_conn_ptr_refs[conn_ptr] == 1
    ccall(PXEST_CONNECTIVITY_DESTROY, Cvoid, (Ptr{Cvoid},), conn_ptr)
    delete!(pxest_conn_ptr_refs, conn_ptr)
  else
    pxest_conn_ptr_refs[conn_ptr] -= 1
  end
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
    @assert !haskey(pxest_conn_ptr_refs, pxest_conn_ptr)
    pxest_conn_ptr_refs[pxest_conn_ptr] = 1
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

  piggy_data::NTuple{PXEST_QUADRANT_PIGGY_SIZE, Cchar}
end

@p4est const pxest_lnodes_code_t = Int8
@p8est const pxest_lnodes_code_t = Int16
struct pxest_lnodes_rank_t
  rank::Cint
  shared_nodes::sc_array_t{pxest_locidx_t}
  shared_mine_offset::pxest_locidx_t
  shared_mine_count::pxest_locidx_t
  owned_offset::pxest_locidx_t
  owned_count::pxest_locidx_t
end

struct pxest_lnodes_t
  mpicomm::MPI.CComm
  num_local_nodes::pxest_locidx_t
  owned_count::pxest_locidx_t
  global_offset::pxest_gloidx_t
  nonlocal_nodes::Ptr{pxest_gloidx_t}
  sharers::Ptr{sc_array_t{pxest_lnodes_rank_t}}
  global_owned_count::Ptr{pxest_locidx_t}

  degree::Cint
  vnodes::Cint
  num_local_elements::pxest_locidx_t
  face_code::Ptr{pxest_lnodes_code_t}
  element_nodes::Ptr{pxest_locidx_t}
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

function Base.getindex(pxest::pxest_t, quadid, treeid)
  trees = unsafe_load(pxest.trees)
  # FIXME: If this fails, then likely something wrong with the union in
  # pxest_quadrant_t
  @assert trees.elem_size == sizeof(pxest_tree_t)
  tree  = unsafe_load(trees.array, treeid)
  tree.quadrants_offset + quadid
end
function Base.getindex(pxest_ptr::Ptr{pxest_t}, quadid, treeid)
  pxest = unsafe_load(pxest_ptr)
  pxest[quadid, treeid]
end

mutable struct PXEST
  pxest::pxest_t
  conn::Connectivity
  ghost::Ptr{pxest_ghost_t}
  lnodes::Ptr{pxest_lnodes_t}

  function PXEST(conn ;mpicomm=MPI.COMM_WORLD, min_lvl = 0)
    pxest = pxest_t()
    ccall(PXEST_INIT_EXT, Cvoid,
          (Ref{pxest_t}, MPI.CComm, Ptr{pxest_connectivity_t}, pxest_locidx_t,
           Cint, Cint, Csize_t, Ptr{Cvoid}, Ptr{Cvoid}),
          pxest, MPI.CComm(mpicomm), conn.pxest_conn_ptr, 0, min_lvl, 1, 0,
          C_NULL, C_NULL)

    @assert haskey(pxest_conn_ptr_refs, conn.pxest_conn_ptr)
    pxest_conn_ptr_refs[conn.pxest_conn_ptr] +=1
    this = new(pxest, conn, C_NULL, C_NULL)
    finalizer(pxest_destroy, this)
    return this
  end
end

function Base.length(pxest::PXEST)
  return pxest.pxest.local_num_quadrants
end
function Base.getindex(pxest::PXEST, quadid, treeid)
  pxest.pxest[quadid, treeid]
end

function pxest_destroy(pxest)
  if pxest.ghost != C_NULL
    ghost_destroy!(pxest)
  end
  if pxest.lnodes != C_NULL
    lnodes_destroy!(pxest)
  end
  conn_ptr = pxest.pxest.connectivity
  ccall(PXEST_DESTROY_EXT, Cvoid, (Ref{pxest_t}, Cint), pxest.pxest, 0)
  pxest_connectivity_destroy(conn_ptr)
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

function lnodes!(pxest, degree)
  pxest.lnodes == C_NULL || error("lnodes is not C_NULL")
  pxest.ghost  != C_NULL ||
   error("ghost is C_NULL (ghost must be set before lnodes! is called")

   # set up lnodes
  pxest.lnodes = ccall(PXEST_LNODES_NEW, Ptr{pxest_lnodes_t},
                       (Ref{pxest_t}, Ptr{pxest_ghost_t}, Cint),
                      pxest.pxest, pxest.ghost, degree)

  # expand ghost to include lnodes
  ccall(PXEST_GHOST_SUPPORT_LNODES, Cvoid,
        (Ref{pxest_t}, Ptr{pxest_lnodes_t}, Ptr{pxest_ghost_t}),
        pxest.pxest, pxest.lnodes, pxest.ghost)
end

function ghost_destroy!(pxest)
  pxest.ghost != C_NULL || error("ghost is C_NULL")
  ccall(PXEST_GHOST_DESTROY, Cvoid, (Ptr{pxest_ghost_t},), pxest.ghost)
  pxest.ghost = C_NULL
end

function lnodes_destroy!(pxest)
  pxest.lnodes != C_NULL || error("lnodes is C_NULL")
  ccall(PXEST_LNODES_DESTROY, Cvoid, (Ptr{pxest_lnodes_t},), pxest.lnodes)
  pxest.lnodes = C_NULL
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
  volinfo.pxest[volinfo.quadid+1, volinfo.treeid+1]
end

struct pxest_iter_face_side_t
  treeid::pxest_topidx_t # the tree on this side
  face::Int8             # which quadrant side the face touches
  is_hanging::Int8       # boolean: one full quad (0) or two smaller quads (1)
  pxest_iter_face_side_data::NTuple{PXEST_ITER_FACE_SIDE_IS_SIZE, Cchar}
end

struct pxest_iter_face_info_t
  pxest::Ptr{pxest_t}
  ghost_layer::Ptr{pxest_ghost_t}
  orientation::Int8   # the orientation of the sides to each other, as in the
                      # definition of p8est_connectivity_t
  tree_boundary::Int8 # boolean: interior face (0), boundary face (1)
  # array of p4est_iter_face_side_t type
  sides::sc_array_t{pxest_iter_face_side_t}
end
function Base.length(face::pxest_iter_face_info_t)
  face.sides.elem_count
end
function Base.getindex(face::pxest_iter_face_info_t, i)
  @assert i == 1 || i == 2
  unsafe_load(face.sides.array, i)
end
function sidetreeid(side::pxest_iter_face_side_t)
  side.treeid+1
end
function sideface(side::pxest_iter_face_side_t)
  side.face+1
end
function sideishanging(side::pxest_iter_face_side_t)
  side.is_hanging!=0
end
function faceorientation(face::pxest_iter_face_info_t)
  face.orientation
end
function sideisghost(side::pxest_iter_face_side_t) :: NTuple{PXEST_HALF, Bool}
  if side.is_hanging != 0
    g0 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_IS_GHOST, Int8,
               (Ref{pxest_iter_face_side_t}, Int), side, 0)
    g1 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_IS_GHOST, Int8,
               (Ref{pxest_iter_face_side_t}, Int), side, 1)
    @p4est return (g0==1, g1==1)
    @p8est g2 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_IS_GHOST, Int8,
                      (Ref{pxest_iter_face_side_t}, Int), side, 2)
    @p8est g3 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_IS_GHOST, Int8,
                      (Ref{pxest_iter_face_side_t}, Int), side, 3)
    @p8est return (g0==1, g1==1, g2==1, g3==1)
  else
    g = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_FULL_IS_GHOST, Int8,
              (Ref{pxest_iter_face_side_t},), side)
    @p4est return (g==1, false)
    @p8est return (g==1, false, false, false)
  end
end

function sidequadid(side::pxest_iter_face_side_t) :: NTuple{PXEST_HALF,
                                                            pxest_locidx_t}
  if side.is_hanging != 0
    q0 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUADID, pxest_locidx_t,
               (Ref{pxest_iter_face_side_t}, Int), side, 0)
    q1 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUADID, pxest_locidx_t,
               (Ref{pxest_iter_face_side_t}, Int), side, 1)
    @p4est return (q0+1, q1+1)
    @p8est q2 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUADID,
                      pxest_locidx_t,
                      (Ref{pxest_iter_face_side_t}, Int), side, 2)
    @p8est q3 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUADID,
                      pxest_locidx_t,
                      (Ref{pxest_iter_face_side_t}, Int), side, 3)
    @p8est return (q0+1, q1+1, q2+1, q3+1)
  else
    q = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_FULL_QUADID, pxest_locidx_t,
              (Ref{pxest_iter_face_side_t},), side)
    @p4est return (q+1, 0)
    @p8est return (q+1, 0, 0, 0)
  end
end
function sidequad(side::pxest_iter_face_side_t) :: NTuple{PXEST_HALF,
                                                          pxest_quadrant_t}
  if side.is_hanging != 0
    q0 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUAD,
               Ptr{pxest_quadrant_t},
               (Ref{pxest_iter_face_side_t}, Int), side, 0)
    q1 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUAD,
               Ptr{pxest_quadrant_t},
               (Ref{pxest_iter_face_side_t}, Int), side, 1)
    @p4est return (q0, q1)
    @p8est q2 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUAD,
                      Ptr{pxest_quadrant_t},
                      (Ref{pxest_iter_face_side_t}, Int), side, 2)
    @p8est q3 = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_HANGING_QUAD,
                      Ptr{pxest_quadrant_t},
                      (Ref{pxest_iter_face_side_t}, Int), side, 3)
    @p8est return (q0, q1, q2, q3)
  else
    q = ccall(PXEST_JULIA_ITER_FACE_SIDE_T_IS_FULL_QUAD, Ptr{pxest_quadrant_t},
              (Ref{pxest_iter_face_side_t},), side)
    @p4est return (q, 0)
    @p8est return (q, 0, 0, 0)
  end
end

@p8est struct pxest_iter_edge_side_t
  treeid::pxest_topidx_t # the tree on this side
  edge::Int8             # which quadrant side the edge touches
  orientation::Int8      # the orientation of each quadrant relative to this
                         # edge, as in the definition of p8est_connectivity_t
  is_hanging::Int8       # boolean: one full quad (0) or two smaller quads (1)
  pxest_iter_edge_side_data::NTuple{PXEST_ITER_EDGE_SIDE_IS_SIZE, Cchar}
  faces::NTuple{2, Int8} # FIXME: Check that this is correct
end

@p8est struct pxest_iter_edge_info_t
  pxest::Ptr{pxest_t}
  ghost_layer::Ptr{pxest_ghost_t}
  tree_boundary::Int8 # boolean: interior face (0), boundary face (1)
  # array of p8est_iter_edge_side_t type
  sides::sc_array_t{pxest_iter_edge_side_t}
end

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
  faces(pxest, face_fn, pxest.ghost)
end

function faces(pxest, face_fn::Function)
  faces(pxest, face_fn, pxest.ghost)
end


function faces(pxest, face_fn::Function, ghost)
  @p8est iterator(pxest, C_NULL, face_fn, C_NULL, C_NULL, ghost)
  @p4est iterator(pxest, C_NULL, face_fn, C_NULL, ghost)
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

@p4est function iterator(pxest, quad_fn, face_fn, corn_fn, ghost=C_NULL)
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
                end, ghost)
end

@p8est function iterator(pxest, quad_fn, face_fn, edge_fn, corn_fn,
                         ghost=C_NULL)
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
                end, ghost)
end

function iterator_call(pxest, quad_fn, face_fn, edge_fn, corn_fn, ghost=C_NULL)
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
               pxest.pxest, ghost, C_NULL, quad_fn_ptr, face_fn_ptr,
               corn_fn_ptr)

  @p8est ccall(PXEST_ITERATE, Cvoid, (Ref{pxest_t}, Ptr{pxest_ghost_t},
                                      Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid},
                                      Ptr{Cvoid}, Ptr{Cvoid}),
               pxest.pxest, ghost, C_NULL, quad_fn_ptr, face_fn_ptr,
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
