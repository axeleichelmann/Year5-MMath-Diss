mutable struct D__type
    ms::ms_type                              # data of master problem sets
    mp::mp_type                              # data of master problem parameters
    ps::ps_type                              # data of subproblem sets
    pp::pp_type                              # data of subproblem parameters
   unc::u_type                               # data of uncertain parameters
  case::String                               # case to solve (0, 1, 2, or 3)
  algm::Int64                                # algorithm used (0, 1, 2, or 3)
     q::Float64                              # tolerance damping parameter                            -> valid only for algorithm 3 (Zaker_Bend)
     w::Int64                                # number of subproblems solved each iteration            -> valid only for algorithm 2 (Adapt_Bend)
  kidx::Array{BitArray{1},1}                 # array of binary association with subset ℰₖ             -> valid only for algorithm 2 (Adapt_Bend)
  iidx::Array{Array{Int64,1},1}              # array of subsets {ℰₖ, 𝑘=1,..,w}                        -> valid only for algorithm 2 (Adapt_Bend)
     J::Int64                                # maximum number of iteration
     δ::Float64                              # convergence tolerance
   emb::Int64                                # embedded LB oracle?
   adj::Float64                              # time adjust (for reporting), needs to be moved to proper location
  stab::Int64                                # stabilisation on or off
    γs::Float64                              # stabilisation factor
    ad::Int64                                # adaptive number of subproblems to solve
    ac::Int64                                # set if perform avoid adding cuts
    pg::Int64                                # new point update gap or old point
   i2n::Dict{Any,Any}                        # index2node
end

mutable struct R__type
     m::JuMP.Model                           # relaxed master problem JuMP model
  vars::Array{JuMP.VariableRef,2}            # reference to variablems {(βᵢ,xᵢ), ∀i}
  cons::Array{JuMP.ConstraintRef,2}          # reference to Benders cuts
  objf::JuMP.AffExpr                         # rmp original objective
end

mutable struct L__type
     m::JuMP.Model                           # level method problem JuMP model
end

mutable struct T1_type
     x::Array{Float64,2}                     # decisions {xᵢ, ∀i} set by the master problem @ iter k
     c::Array{Float64,2}                     # parameters {cᵢ, ∀i}
     β::Array{Float64,1}                     # variables value {βᵢ, ∀i} of the master problem @ iter k
     θ::Array{Float64,1}                     # optimal objectives {θᵢ, ∀i} of subproblems @ iter k
     λ::Array{Float64,2}                     # subgradient {λᵢ, ∀i} wrt xᵢ of subproblems @ iter k
     Δ::Float64                              # optimality gap @ iter k
     xr::Array{Float64,2}                     # decisions of the RMP
     fr::Float64                              # RMP investment only cost
     lf::Float64                              # stabilisation factor
   aux::Dict{Any,Any}                       # auxiliary information, to be held (depends on the instance, in here for LMP)
end

mutable struct T2_type
     x::Array{Float64,2}                     # decisions {xᵢ, ∀i} set by the master problem @ iter k
     c::Array{Float64,2}                     # parameters {cᵢ, ∀i}
     β::Array{Float64,1}                     # variables value {βᵢ, ∀i} of the master problem @ iter k
    θl::Array{Float64,1}                     # valid lower bounds on objectives {θᵢ, ∀i} of subproblems @ iter k
    λl::Array{Float64,2}                     # valid subgradients {λᵢ, ∀i} wrt xᵢ of subproblems @ iter k
    θu::Array{Float64,1}                     # valid upper bounds on objectives {θᵢ, ∀i} of subproblems @ iter k
    ϕu::Array{Float64,2}                     # valid subgradients {ϕᵢ, ∀i} wrt cᵢ of subproblems @ iter k
    Ie::Array{Int64,1}                       # index of w subproblems solved exactly @ iter k
     Δ::Float64                              # optimality gap @ iter k
    xr::Array{Float64,2}                     # decisions of the RMP
    fr::Float64                              # RMP investment only cost
    lf::Float64                              # stabilisation factor
    ne::Int64                                # number of subproblems solved
    lb::Float64                              # lb for the whole problem
    ub::Float64                              # ub for the whole problem
    δl::Float64                              # relative gap between bounds lb and ub
    bub::Float64                             # best upper bound
    aux::Dict{Any,Any}                       # auxiliary information, to be held (depends on the instance, in here for LMP)

end

mutable struct T3_type
     x::Array{Float64,2}                     # decisions {xᵢ, ∀i} set by the master problem @ iter k
     c::Array{Float64,2}                     # parameters {cᵢ, ∀i}
     β::Array{Float64,1}                     # variables value {βᵢ, ∀i} of the master problem @ iter k
    θp::Array{Float64,1}                     # valid upper bounds on objectives {θᵢ, ∀i} of subproblems @ iter k
    θd::Array{Float64,1}                     # valid lower bounds on objectives {θᵢ, ∀i} of subproblems @ iter k
     λ::Array{Float64,2}                     # valid subgradients {λᵢ, ∀i} wrt xᵢ of subproblems @ iter k
     Δ::Float64                              # optimality gap @ iter k
end

mutable struct H__type
     L::Array{Float64,1}                     # lower bound values (up to iter k)
     U::Array{Float64,1}                     # upper bound values (up to iter k)
     t::Array{Float64,1}                     # time spend to perform each iteration (up to iter k)
     T::Array{Float64,2}                     # time spend to perform each iteration, split by steps (up to iter k)
     k::Int64                                # current iteration number
end

mutable struct M__type
     n::Int64                                # number of exact solution stored
     θ::Array{Float64,1}                     # optimal objective of exact solutions stored
     λ::Array{Float64,2}                     # subgradient wrt x of exact solutions stored
     ϕ::Array{Float64,2}                     # subgradient wrt c of exact solutions stored
     x::Array{Float64,2}                     # values of x of exact solutions stored
     c::Array{Float64,2}                     # values of c of exact solutions stored
end

mutable struct B1_type
   rmp::R__type                              # relaxed master problem
  temp::T1_type                              # temporary data (updated each iteration)
  hist::H__type                              # data stored during the solution
  data::D__type                              # problem data
   lmp::L__type                              # level set method problem
end

mutable struct B2_type
   rmp::R__type                              # relaxed master problem
  temp::T2_type                              # temporary data (updated each iteration)
  hist::H__type                              # data stored during the solution
    m0::M__type                              # exact solutions stored
  data::D__type                              # problem data
   lmp::L__type                              # level set method problem
 end

mutable struct B3_type
   rmp::R__type                              # relaxed master problem
  temp::T3_type                              # temporary data (updated each iteration)
  hist::H__type                              # data stored during the solution
  data::D__type                              # problem data
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- structures of data : subproblem --------------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

export E__type                               # structure of exact subproblem
export O__type                               # structure of oracles
export OP_type                               # structure of oracles (primal version)
export Q__type                               # structure of extended exact solutions stored during the solution
export d__type                               # structure of subproblem data
export t1_type                               # structure of temporary data (updated each iteration)    -> valid only for algorithm 1 (Stand_Bend)
export t2_type                               # structure of temporary data (updated each iteration)    -> valid only for algorithm 2 (Adapt_Bend)
export t3_type                               # structure of temporary data (updated each iteration)    -> valid only for algorithm 3 (Zaker_Bend)
export S1_type                               # structure of Benders subproblem                         -> valid only for algorithm 1 (Stand_Bend)
export S2_type                               # structure of Benders subproblem                         -> valid only for algorithm 2 (Adapt_Bend)
export S3_type                               # structure of Benders subproblem                         -> valid only for algorithm 3 (Zaker_Bend)

mutable struct E__type
      m::JuMP.Model                          # exact subproblem JuMP model
   vars::Array{JuMP.VariableRef,1}           # reference to variables (c₀,ϕ)
   objf::GenericAffExpr{Float64,VariableRef} # reference to objective function
   objfor::GenericAffExpr{Float64,VariableRef} # reference to objective function (original)
end

mutable struct O__type
      m::JuMP.Model                          # oracle JuMP model
   cons::Array{JuMP.ConstraintRef,1}         # reference to oracle new constraints (new exact solutions)
   vars::Array{JuMP.VariableRef,1}           # reference to variables (ϕ,γ)
   objf::GenericAffExpr{Float64,VariableRef} # reference to objective function
   cexp::GenericAffExpr{Float64,VariableRef} # reference to objective function
   help::Array{Float64,1}                    # rhs values of old constraints (old exact solutions)
      n::Int64                               # number of exact solutions added to the oracle
end

mutable struct OP_type
      m::JuMP.Model                          # oracle JuMP model
   objf::GenericAffExpr{Float64,VariableRef} # reference to objective function
      n::Int64                               # number of exact solutions added to the oracle
end

mutable struct Q__type
     ϕ::Array{Float64,1}                     # operational cost dependent on uncertain cost cᵢ
    c0::Float64                              # operational cost independent of uncertain costs c
    yG::Array{Float64,3}                     # generation level of conventional generation
    yI::Array{Float64,3}                     # charging power of storage generation
    yO::Array{Float64,3}                     # discharging power of storage generation
    yL::Array{Float64,3}                     # energy level of storage generation
    yS::Array{Float64,2}                     # load shedding
    x0::Array{Float64,1}                     # values of rhs parameters in the subproblem (eg, capacity)
 end

mutable struct d__type
     ps::ps_type                             # data of subproblem sets
     pp::pp_type                             # data of subproblem parameters
    unc::u_type                              # data of uncertain parameters
   algm::Int64                               # algorithm used (0, 1, 2, or 3)
      J::Int64                               # maximum number of exact solution to store
end

mutable struct t1_type
      x::Array{Float64,2}                    # decisions {xᵢ, ∀i} set by the master problem @ iter k
      c::Array{Float64,2}                    # parameters {cᵢ, ∀i}
      θ::Array{Float64,1}                    # optimal objectives {θᵢ, ∀i} of subproblems @ iter k
      λ::Array{Float64,2}                    # subgradient {λᵢ, ∀i} wrt xᵢ of subproblems @ iter k
      i::Int64                               # index of subproblem to solve
end

mutable struct t2_type
     x::Array{Float64,2}                     # decisions {xᵢ, ∀i} set by the master problem @ iter k
     c::Array{Float64,2}                     # parameters {cᵢ, ∀i}
    θl::Array{Float64,1}                     # valid lower bounds on objectives {θᵢ, ∀i} of subproblems @ iter k
    θu::Array{Float64,1}                     # valid upper bounds on objectives {θᵢ, ∀i} of subproblems @ iter k
    λl::Array{Float64,2}                     # valid subgradients {λᵢ, ∀i} wrt xᵢ of subproblems @ iter k
    ϕu::Array{Float64,2}                     # valid subgradients {ϕᵢ, ∀i} wrt cᵢ of subproblems @ iter k
     i::Int64                                # index of subproblem to solve
end

mutable struct t3_type
     x::Array{Float64,2}                     # decisions {xᵢ, ∀i} set by the master problem @ iter k
     c::Array{Float64,2}                     # parameters {cᵢ, ∀i}
    θp::Array{Float64,1}                     # valid upper bounds on objectives {θᵢ, ∀i} of subproblems @ iter k
    θd::Array{Float64,1}                     # valid lower bounds on objectives {θᵢ, ∀i} of subproblems @ iter k
     λ::Array{Float64,2}                     # valid subgradients {λᵢ, ∀i} wrt xᵢ of subproblems @ iter k
     i::Int64                                # optimality gap @ iter k
end

mutable struct S1_type
     ex::E__type                             # exact subproblem
   data::d__type                             # subproblem data
   temp::t1_type                             # temporary data (updated each iteration)
end

mutable struct S2_type
     ex::E__type                             # exact subproblem
    olb::O__type                             # oracle for lower bound
    oub::O__type                             # oracle for upper bound
      m::M__type                             # exact solutions stored
      q::Array{Q__type,1}                    # extended exact solutions stored
   data::d__type                             # subproblem data
   temp::t2_type                             # temporary data (updated each iteration)
   olbP::OP_type                             # oracle for lower bound (primal version)
   oubP::OP_type                             # oracle for upper bound (primal version)
end

mutable struct S3_type
     ex::E__type                             # exact subproblem
   data::d__type                             # subproblem data
   temp::t3_type                             # temporary data (updated each iteration)
end


mutable struct N1_type
     UB::Array{Float64,1}                      # upper bound
     LB::Array{Float64,1}                    # lower bound
     λc::Array{Array{Float64,2},1}             # subgradient of the cuts (used for cut signature)
     hc::Array{Array{Float64,1},1}             # height of the cuts (used for cut signature)
     sc::Array{Array{Float64,1},1}             # signature of the cuts (used for cut signature)
     rv::Array{Float64,1}                      # random vector (used for cut signature)

end

mutable struct N2_type
   xrmp::Array{Array{Float64,2},1}           # accumulated capacity
  x0rmp::Array{Array{Float64,2},1}           # newly invested capacity
     LB::Array{Float64,1}                    # lower bound
 LMP_target::Array{Float64,1}                # LMP target
 LMP_distance::Array{Float64,1}              # distance moved from one iter to the next
  xlmp::Array{Array{Float64,2},1}            # LMP accumulated capacity
  x0lmp::Array{Array{Float64,2},1}           # LMP newly invested capacity
  LBOα::Array{Float64,1}                     # lower bound oracle before solving SP
  UBOα::Array{Float64,1}                     # upper bound oracle before solving SP
   LBO::Array{Float64,1}                     # lower bound oracle as solving SP
   UBO::Array{Float64,1}                     # upper bound oracle as solving SP
 LBOiα::Array{Array{Float64,1},1}                      # lower bound oracle before solving SP
 UBOiα::Array{Array{Float64,1},1}                     # upper bound oracle before solving SP
  LBOi::Array{Array{Float64,1},1}                     # lower bound oracle as solving SP
  UBOi::Array{Array{Float64,1},1}                     # upper bound oracle as solving SP
   nSP::Array{Int64,1}                       # number of SP solved at an iteration
   UB::Array{Float64,1}                      # upper bound
   tRMP::Array{Float64,1}                    # time spent solving RMP
   tLMP::Array{Float64,1}                    # time spent solving LMP
   tSP::Array{Float64,1}                     # time spent solving SP
   ttotal::Array{Float64,1}                  # total time spent
   ubo_matrix::Array{Array{Float64,2},1}     # ubo matrix
   θc::Array{Array{Float64,1},1}             # height of the cuts
   λc::Array{Array{Float64,2},1}             # subgradient of the cuts
   xc::Array{Array{Float64,2},1}             # point in cuts
   hc::Array{Array{Float64,1},1}             # height of the cuts
   sc::Array{Array{Float64,1},1}             # signature of the cuts
   rv::Array{Float64,1}                      # random vector
 gapi::Array{Array{Float64,1},1}             # gap of each node i
gapis::Array{Array{Float64,1},1}             # gap of each node i after solving SP
    δ::Array{Float64,1}                      # convergence
  sSP::Array{Int64,1}                        # store the solves SP
 UBOm::Array{DataFrame,1}                    # array of UBO
   vSP::Array{Array{Int64,1},1}              # list of SP solved in the iteration
  Rval::Array{Float64,1}                     # real function value if all SP are solved exactly
 ubo_sol::Array{Array{Float64,1},1}          # variable values in ubo
 ubo_dual::Array{Array{Float64,1},1}         # dual of constraints in ubo
end
