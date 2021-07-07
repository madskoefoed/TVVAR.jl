#module TVVAR

# Code based on Forecasting with time-varying vector autoregressive models (2008), K. Triantafyllopoulos

# Import
#using Base: AbstractFloat
using LinearAlgebra: diag, kron, I, cholesky
using Distributions: Normal, MvNormal, MvTDist, InverseWishart, Wishart, MatrixNormal
using Parameters: @unpack

# Include scripts
include("./src/types.jl")
include("./src/utils.jl")
include("./src/estimation.jl")

include("./example/univariate.jl")
include("./example/multivariate.jl")

# Exported types
export TVVAR, KF

# Exported functions
export estimate

#end