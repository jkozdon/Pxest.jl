module Pxest

module p4est
  macro p4est(ex)
    return :($(esc(ex)))
  end

  macro p8est(ex)
    return :()
  end

  include("pxest-base.jl")

end

module p8est

  macro p4est(ex)
    return :()
  end

  macro p8est(ex)
    return :($(esc(ex)))
  end

  include("pxest-base.jl")

end

end
