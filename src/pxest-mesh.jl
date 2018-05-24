mutable struct mesh_connectivity
  EToL::Array{Cint, 1} # element to p4est level
  EToT::Array{Cint, 1} # element to p4est treeid
  EToX::Array{Cint, 1} # element to p4est x-qcoord
  EToY::Array{Cint, 1} # element to p4est y-qcoord
  EToZ::Array{Cint, 1} # element to p4est z-qcoord
  EToB::Array{Cint, 2} # element to boundary condition
  EToE::Array{Cint, 2} # element to neighbor element
  EToF::Array{Cint, 2} # element to neighbor face
  EToO::Array{Cint, 2} # element to neighbor orientation
  EToP::Array{Cint, 2} # element to periodicity mask (filled only for brick)

  function mesh_connectivity(pxest)
    mesh = new()
  end
end
