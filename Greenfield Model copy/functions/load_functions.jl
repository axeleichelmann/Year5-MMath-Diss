print("Load decomp. functions...    "); a=time();


# */ --------------------------------------------------------------------------------------------- /* #
# */ --- optimization models --------------------------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function LBoracle!(m::JuMP.Model,unc::u_type,J::Int64)::JuMP.Model

    # problem structure
    # min   ϕ + γᵀcᵢ
    # s.t.  ϕ + γᵀcₛ ≧ θₛ + λᵀₛ(x-xₛ),   ∀s
    #       x = xᵢ

    # */ -- variables ---------------------------------------------- /* #
    @variable(m, phi, container=Array)
    @variable(m, gamma[1:unc.nc] >= .0, container=Array)
    @variable(m, x[1:unc.nx], container=Array)

    # @constraint(m,csr_nonneg,phi+gamma'*ones((unc.nc))>=0.0)
    # */ ----------------------------------------------------------- /* #

    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Min, 0. )
    # */ ----------------------------------------------------------- /* #

    # */ -- constraints -------------------------------------------- /* #
    # */ ----------------------------------------------------------- /* #
    #optimize!(m)

    return m
end

function UBoracle!(m::JuMP.Model,unc::u_type,J::Int64)::JuMP.Model

    # problem structure
    # max   ϕ + γᵀxᵢ
    # s.t.  ϕ + γᵀxₛ ≦ θₛ + ϕᵀₛ(c-cₛ),   ∀s
    #       c = cᵢ

    # */ -- variables ---------------------------------------------- /* #
    @variable(m, phi, container=Array)
    @variable(m, gamma[1:unc.nx] <= .0, container=Array)
    @variable(m, c[1:unc.nc], container=Array)
    # */ ----------------------------------------------------------- /* #

    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Max, 0. )
    # */ ----------------------------------------------------------- /* #

    # */ -- constraints -------------------------------------------- /* #
    # */ ----------------------------------------------------------- /* #
    #optimize!(m)

    return m
end




function UBoraclePrimal!(m::JuMP.Model,unc::u_type,J::Int64)::JuMP.Model


    # */ -- variables ---------------------------------------------- /* #
    @variable(m,ν[1:J]==0.0)

    # */ -- obj function ------------------------------------------- /* #
    @objective(m,Min,0. )
    # */ ----------------------------------------------------------- /* #

    # */ -- constraints -------------------------------------------- /* #
    @constraint(m,csr_lim,zeros(unc.nx,J)*ν.<=zeros(unc.nx));
    @constraint(m,csr_conv,sum(ν)==1.0);
    # */ ----------------------------------------------------------- /* #

    #optimize!(m)

    return m
end


function LBoraclePrimal!(m::JuMP.Model,unc::u_type,J::Int64)::JuMP.Model

    # */ -- variables ---------------------------------------------- /* #
    @variable(m,μ[1:J]==0.0)

    # */ -- obj function ------------------------------------------- /* #
    @objective(m,Max,0. )
    # */ ----------------------------------------------------------- /* #

    # */ -- constraints -------------------------------------------- /* #
    @constraint(m,csr_lim,zeros(unc.nc,J)*μ.<=zeros(unc.nc));
    @constraint(m,csr_conv,sum(μ)==1.0);
    # */ ----------------------------------------------------------- /* #

    #optimize!(m)

    return m
end


# */ --------------------------------------------------------------------------------------------- /* #
# */ --- functions to generate structs : master -------------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function gen_kiidx(u::u_type,w::Int64)::Tuple{Array{BitArray{1},1},Array{Array{Int64,1},1}}

    # function to generate (kidx,iidx) for given number of subproblems w solved at each iteration
    # iidx -> array of subsets {ℰₖ, 𝑘=1,..,w}
    # kidx -> array of 0/1 association with subset ℰₖ

    if w == 1 # (trivial case with one subproblem) ℰₖ := ℰ
        kidx = Array{BitArray{1},1}(undef,0)
        push!(kidx,zeros(Int64,u.ni).==0)
        iidx = Array{Array{Int64,1},1}(undef,0)
        push!(iidx,[i for i in 1:u.ni])
    else
        x = zeros(u.nh+u.nc,u.ni) # generate empty matrix x for (H,C) values of each subproblem
        for j in 1:u.nh
            uh = u.h[:,j] # collect values of H (rhs uncertain params)
            hmin,hmax = minimum(uh),maximum(uh) # evaluate min & max values
            (hmax > hmin) ? x[j,:] .= (uh.-hmin)./(hmax-hmin) :  x[j,:] .= zeros(u.ni) # set normalized H values into x
        end
        for j in 1:u.nc
            uc = u.c[:,j] # collect values of C (uncertain cost coeff)
            cmin,cmax = minimum(uc),maximum(uc) # evaluate min & max values
            (cmax > cmin) ? x[j+u.nh,:] .= (uc.-cmin)./(cmax-cmin) : x[j+u.nh,:] .= zeros(u.ni) # set normalized C values into x
        end
        kass = kmeans(x,w).assignments # run kmeans clustering on x to obtain w subsets
        kidx = Array{BitArray{1},1}(undef,0) # generate empty kidx
        for j in 1:w
            push!(kidx,kass.==j) # for each subset j = 1,..,w generate a bitarry with 1 if subproblem i belongs to subset j (and 0 otherwise)
        end
        idx = [i for i in 1:u.ni] # generate array of indices of subproblems
        iidx = Array{Array{Int64,1},1}(undef,0) # generate empty iidx
        for j in 1:w
            push!(iidx,idx[kidx[j]]) # for each subset j = 1,..,w generate an array with all the subproblems that belongs to subset j
        end
    end

    return kidx,iidx
end

function gen_D__type(ms::ms_type,mp::mp_type,ps::ps_type,pp::pp_type,u::u_type,c::String,a::Int64,q::Float64,w::Int64,j::Int64,e::Float64,emb::Int64,stab::Int64,γs::Float64,ad::Int64,ac::Int64,pg::Int64,i2n)::D__type

    # function to generate D__type (structure of problem data) for given ms (master problem sets), mp (master problem params),
    # ps (subproblem sets), pp (subproblem params), u (uncertainty params), c (case), a (algorithm), q (tolerance damping),
    # w (number of exact subproblems), j (maximum number of iteration), and e (convergence tolerance)

    if ad==1 w=1 end

    if a==2
        k,i = gen_kiidx(u,w) # generate (kidx,iidx) for given uncertain data u and given number of subproblems w solved at each iteration
    else
        k,i=Array{BitArray{1},1}(),Array{Array{Int64,1},1}()
    end

    return D__type(ms,mp,ps,pp,u,c,a,q,w,k,i,j,e,emb,0.0,stab,γs,ad,ac,pg,i2n)
end



function gen_R__type(d::D__type,m::Model)::R__type

    # function to generate R__type (structure of relax master problem) given data D__type (problem data)

    o = objective_function(m);
    v = Array{JuMP.VariableRef,2}(undef,d.unc.nx+1,d.unc.ni) # generate empty array of vars
    for i in 1:d.unc.ni v[:,i] .= vcat([m[:beta][d.i2n[i]]],m[:x][1:d.unc.nx,i]) end # set reference to variablems {(βᵢ,xᵢ), ∀i}
    c = Array{JuMP.ConstraintRef,2}(undef,d.J,d.unc.ni) # generate empty array of cons

    return R__type(m,v,c,o)
end

function gen_L__type(d::D__type,m::Model)::L__type

    # function to generate L__type (structure of level method problem) given data D__type (problem data)
    return L__type(m)
end

function gen_T1_type(d::D__type)::T1_type

    # function to generate T1_type (structure of temporary data, algorithm 1 ) given data D__type (problem data)

    x = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions {xᵢ, ∀i}
    c = convert(Array{Float64,2},d.unc.c') # set values of parameters {cᵢ, ∀i}
    β = zeros(d.unc.ni) # generate empty array for variables value {βᵢ, ∀i}
    θ = zeros(d.unc.ni) # generate empty array for optimal objectives {θᵢ, ∀i}
    λ = zeros(d.unc.nx,d.unc.ni) # generate empty array for subgradient {λᵢ, ∀i} wrt xᵢ
    xr  = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions RMP {xᵢ, ∀i}


    return T1_type(x,c,β,θ,λ,0.,xr,0.0,0.0,Dict{Any,Any}())
end

function gen_T2_type(d::D__type)::T2_type

    # function to generate T1_type (structure of temporary data, algorithm 2) given data D__type (problem data)

    x  = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions current {xᵢ, ∀i}
    c  = convert(Array{Float64,2},d.unc.c') # set values of parameters {cᵢ, ∀i}
    β  = zeros(d.unc.ni) # generate empty array for variables value {βᵢ, ∀i}
    θl = zeros(d.unc.ni) # generate empty array for valid lower bounds on objectives {θᵢ, ∀i}
    λl = zeros(d.unc.nx,d.unc.ni) # generate empty array for valid subgradient {λᵢ, ∀i} wrt xᵢ
    θu = zeros(d.unc.ni) # generate empty array for valid upper bounds on objectives {θᵢ, ∀i}
    ϕu = zeros(d.unc.nc,d.unc.ni) # generate empty array for valid subgradient {ϕᵢ, ∀i} wrt cᵢ
    Ie = zeros(Int64,d.w) # generate empty array for indices of w subproblems
    xr  = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions RMP {xᵢ, ∀i}

    return T2_type(x,c,β,θl,λl,θu,ϕu,Ie,0.,xr,Inf,0.0,0,0.0,Inf,Inf,Inf,Dict{Any,Any}())
end

function gen_H__type(d::D__type)::H__type

    # function to generate H_type (structure of stored data) given data D__type (problem data)

    L = zeros(d.J) # generate empty array for lower bound values
    U = zeros(d.J) # generate empty array for upper bound values
    t = zeros(d.J) # generate empty array for time spend to perform each iteration
    T = zeros(d.J,6) # generate empty array for time spend to perform each iteration, split by steps

    return H__type(L,U,t,T,0)
end

function gen_M__type(d::D__type)::M__type

    # function to generate M_type (structure of exact solutions stored) given data D__type (problem data)
    θ = zeros(1) # generate empty array of optimal objective
    λ = zeros(1,d.unc.nx) # generate empty array of subgradient wrt x
    ϕ = zeros(1,d.unc.nc) # generate empty array of subgradient wrt c
    x = zeros(1,d.unc.nx) # generate empty array of values of x
    c = zeros(1,d.unc.nc) # generate empty array of values of c

    return M__type(0,θ,λ,ϕ,x,c)
end


#function gen_B1_type(c::String,a::Int64,j::Int64,e::Float64,emb::Int64,stab::Int64,γs::Float64,ms,mp,ps,pp,unc,m,i2n,mLP)::B1_type
function gen_B1_type(c::String,a::Int64,w::Int64,j::Int64,e::Float64,emb::Int64,stab::Int64,γs::Float64,ad::Int64,ac::Int64,pg::Int64,ms,mp,ps,pp,unc,m,i2n,m3)::B1_type

    # function to generate B1_type (Benders master structure of algorithm 1 ) given c (case), a (algorithm),
    # j (maximum number of iteration), and e (convergence tolerance).

    d = gen_D__type(ms,mp,ps,pp,unc,c,a,.0,0,j,e,emb,stab,γs,0,ac,0,i2n)
    r = gen_R__type(d,m) # generate R__type (structure of relax master problem)
    t = gen_T1_type(d) # generate T1_type (structure of temporary data)
    h = gen_H__type(d) # generate H__type (structure of stored data)

    l = gen_L__type(d,m3) # generate R__type (structure of relax master problem)

    return B1_type(r,t,h,d,l)
end

function gen_B2_type(c::String,a::Int64,w::Int64,j::Int64,e::Float64,emb::Int64,stab::Int64,γs::Float64,ad::Int64,ac::Int64,pg::Int64,ms,mp,ps,pp,unc,m,i2n,m3)::B2_type

    # function to generate B2_type (Benders master structure of algorithm 2) given c (case), a (algorithm), w (number of exact subproblems),
    # j (maximum number of iteration), and e (convergence tolerance).

    #ms,mp,ps,pp,unc = load_data(c) # generate data structures ms,mp,ps,pp,unc for given value of c (case)
    d = gen_D__type(ms,mp,ps,pp,unc,c,a,.0,w,j,e,emb,stab,γs,ad,ac,pg,i2n)  # generate D__type (structure of problem data)
    r = gen_R__type(d,m) # generate R__type (structure of relax master problem)
    t = gen_T2_type(d) # generate T2_type (structure of temporary data)
    h = gen_H__type(d) # generate H__type (structure of stored data)
    m = gen_M__type(d) # generate M__type (structure of exact solutions stored)

    l = gen_L__type(d,m3) # generate R__type (structure of relax master problem)

    return B2_type(r,t,h,m,d,l)
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- functions to generate structs : subproblem ---------------------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function gen_d__type(d::D__type)::d__type

    # function to generate d__type (structure of subproblem data) for given D__type (structure of problem data)

    s = d.ps # set subproblem sets
    p = d.pp # set subproblem parameters
    u = d.unc # set uncertainty parameters
    a = d.algm # set algorithm
    j = Int64(round(d.J*d.w;digits=0)) # set maximum number of exact solutions to store


    return d__type(s,p,u,a,j)
end

function gen_Lb_type(d::d__type)::O__type

    # function to generate O__type (structure of oracle for lower bound) for given d__type (structure of subproblem data)

    #m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env),"OutputFlag"=>0)) # generate empty JuMP model
    #JuMP.set_optimizer_attribute(m,"Method",1)  # set solution method (barrier)
    #JuMP.set_optimizer_attribute(m,"Presolve",1)  # turn off crossover
    m = Model(HiGHS.Optimizer) # generate empty JuMP model
    MOI.set(m, MOI.Silent(), true)

    LBoracle!(m,d.unc,d.J+1) # build lower bound oracle
    v = vcat(m[:phi],m[:gamma]) # set reference to variables (ϕ,γ)
    c = Array{JuMP.ConstraintRef,1}(undef,d.J+1) # generate empty reference to new constraints
    o = v'*ones(length(v)) # generate reference for objective function
    e = v'*ones(length(v)) # generate reference for objective function
    h = zeros(d.J+1) # generate empty array for rhs values of old constraints (old exact solutions)

    return O__type(m,c,v,o,e,h,0)
end

function gen_Ub_type(d::d__type)::O__type

    # function to generate O__type (structure of oracle for upper bound) for given d__type (structure of subproblem data)

    #m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env),"OutputFlag"=>0))
    #JuMP.set_optimizer_attribute(m,"Method",1)  # set solution method (barrier)
    #JuMP.set_optimizer_attribute(m,"Presolve",1)  # turn off crossover
    m = Model(HiGHS.Optimizer) # generate empty JuMP model
    MOI.set(m, MOI.Silent(), true)
        
    UBoracle!(m,d.unc,d.J+1) # build upper bound oracle
    v = vcat(m[:phi],m[:gamma]) # set reference to variables (ϕ,γ)
    c = Array{JuMP.ConstraintRef,1}(undef,d.J+1) # generate empty reference to new constraints
    o = v'*ones(length(v)) # generate reference for objective function
    e = v'*ones(length(v)) # generate reference for objective function
    h = zeros(d.J+1) # generate empty array for rhs values of old constraints (old exact solutions)

    return O__type(m,c,v,o,e,h,0)
end


function gen_LbP_type(d::d__type)::OP_type

    # function to generate OP_type (structure of oracle for upper bound) for given d__type (structure of subproblem data)

    #m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env),"OutputFlag"=>0))

    m = Model(HiGHS.Optimizer) # generate empty JuMP model
    MOI.set(m, MOI.Silent(), true)


    LBoraclePrimal!(m,d.unc,d.J+1) # build upper bound oracle

    return OP_type(m,AffExpr(0.0),0)
end


function gen_UbP_type(d::d__type)::OP_type

    # function to generate OP_type (structure of oracle for upper bound) for given d__type (structure of subproblem data)

    #m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env),"OutputFlag"=>0))
    #UBoraclePrimal!(m,d.unc,d.J+1) # build upper bound oracle

    m = Model(HiGHS.Optimizer) # generate empty JuMP model
    MOI.set(m, MOI.Silent(), true)


    return OP_type(m,AffExpr(0.0),0)
end


function gen_E__type(d::d__type,m)::E__type

    # function to generate E__type (structure of exact subproblem) for given d__type (structure of subproblem data)

    d.algm == 3 ? JuMP.set_optimizer_attribute(m,"Crossover",0) : nothing # deactivate crossover (if algorithm 3)
    v = vcat(m[:phi][:]) # set reference to variables (c₀,ϕ)
    o = v'*ones(length(v)) # generate reference for objective function
    of= objective_function(m)

    return E__type(m,v,o,of)
end

function gen_M__type(d::d__type)::M__type

    # function to generate M__type (structure of exact solutions stored) for given d__type (structure of subproblem data)

    θ = zeros(d.J+1) # generate empty array for optimal objective of exact solutions
    λ = zeros(d.J+1,d.unc.nx) # generate empty array for subgradient wrt x of exact solutions
    ϕ = zeros(d.J+1,d.unc.nc) # generate empty array for subgradient wrt c of exact solutions
    x = zeros(d.J+1,d.unc.nx) # generate empty array for values of x of exact solutions
    c = zeros(d.J+1,d.unc.nc) # generate empty array for values of c of exact solutions

    return M__type(0,θ,λ,ϕ,x,c)
end

function gen_Q__type(d::d__type)::Array{Q__type,1}

    # function to generate Q__type (structure of extended exact solutions stored) for given d__type (structure of subproblem data)

    return Array{Q__type,1}(undef,0)
end

function gen_t1_type(d::d__type)::t1_type

    # function to generate t1_type (structure of temporary data, algorithm 1 ) for given d__type (structure of subproblem data)

    x = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions {xᵢ, ∀i}
    c = convert(Array{Float64,2},d.unc.c') # set values of parameters {cᵢ, ∀i}
    θ = zeros(d.unc.ni) # generate empty array for optimal objectives {θᵢ, ∀i}
    λ = zeros(d.unc.nx,d.unc.ni) # generate empty array for subgradient {λᵢ, ∀i} wrt xᵢ

   return t1_type(x,c,θ,λ,0)
end


function gen_t2_type(d::d__type)::t2_type

    # function to generate t2_type (structure of temporary data, algorithm 2) for given d__type (structure of subproblem data)

    x  = zeros(d.unc.nx,d.unc.ni) # generate empty array for decisions {xᵢ, ∀i}
    c  = convert(Array{Float64,2},d.unc.c') # set values of parameters {cᵢ, ∀i}
    θl = zeros(d.unc.ni) # generate empty array for valid power bounds on objectives {θᵢ, ∀i}
    θu = zeros(d.unc.ni) # generate empty array for valid upper bounds on objectives {θᵢ, ∀i}
    λl = zeros(d.unc.nx,d.unc.ni) # generate empty array for valid subgradients {λᵢ, ∀i} wrt xᵢ
    ϕu = zeros(d.unc.nc,d.unc.ni) # generate empty array for valid subgradients {ϕᵢ, ∀i} wrt cᵢ

   return t2_type(x,c,θl,θu,λl,ϕu,0,)
end

function gen_S1_type(D::D__type,m)::S1_type

    # function to generate S1_type (structure of Benders subproblem, algorithm 1 ) for given d__type (structure of subproblem data)

    d = gen_d__type(D) # generate d__type (structure of subproblem data)
    e = gen_E__type(d,m) # generate E__type (structure of exact subproblem)
    t = gen_t1_type(d) # generate t1_type (structure of temporary data)

    return S1_type(e,d,t)
end

function gen_S3_type(D::D__type,m)::S3_type

    # function to generate S3_type (structure of Benders subproblem, algorithm 3) for given d__type (structure of subproblem data)

    d = gen_d__type(D) # generate d__type (structure of subproblem data)
    e = gen_E__type(d,m) # generate E__type (structure of exact subproblem)
    t = gen_t3_type(d) # generate t3_type (structure of temporary data)

    return S3_type(e,d,t)
end

function gen_S2_type(D::D__type,m)::S2_type

    # function to generate S2_type (structure of Benders subproblem, algorithm 2) for given d__type (structure of subproblem data)

    d = gen_d__type(D) # generate d__type (structure of subproblem data)
    e = gen_E__type(d,m) # generate E__type (structure of exact subproblem)
    l = gen_Lb_type(d) # generate O__type (structure of oracle for lower bound)
    u = gen_Ub_type(d) # generate O__type (structure of oracle for upper bound)
    m = gen_M__type(d) # generate M__type (structure of exact solutions stored)
    q = gen_Q__type(d) # generate Q__type (structure of extended exact solutions stored)
    t = gen_t2_type(d) # generate t2_type (structure of temporary data)
    lp = gen_LbP_type(d) # generate OP_type (structure of oracle for lower bound - primal version)
    up = gen_UbP_type(d) # generate OP_type (structure of oracle for upper bound - primal version)


    return S2_type(e,l,u,m,q,d,t,lp,up)
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- functions to generate structs : master and subproblem ----------------------------------- /* #
# */ --------------------------------------------------------------------------------------------- /* #

function generate_structures(c::String,a::Int64,q::Float64,w::Int64,j::Int64,e::Float64,emb::Int64,stab::Int64,γs::Float64,ad::Int64,ac::Int64,pg::Int64,ms,mp,ps,pp,unc,mMP,mSP,i2n;mLP=Model(),mU=Model())

    # function to generate B__struct and S__struct for given values of c (case), a (algorithm), q (tolerance damping),
    # w (number of exact subproblems), j (maximum number of iteration), and e (convergence tolerance)

    if a == 0
        B = gen_B0_type(c,ms,mp,ps,pp,unc,mU,i2n,mU) # generate B__struct (Benders master structure)
        S = nothing # empty subproblem structure
    end
    if a == 1 # (if Stand_Bend)
        #B = gen_B1_type(c,a,j,e,emb,stab,γs,ms,mp,ps,pp,unc,mMP,i2n,mLP)  # generate B__struct (Benders master structure) - what is mLP?
        B = gen_B1_type(c,a,w,j,e,emb,stab,γs,ad,ac,pg,ms,mp,ps,pp,unc,mMP,i2n,mLP)
        S = gen_S1_type(B.data,mSP) # generate S__struct (Benders subproblem structure)
    end
    if a == 2 # (if Adapt_Bend)
        B = gen_B2_type(c,a,w,j,e,emb,stab,γs,ad,ac,pg,ms,mp,ps,pp,unc,mMP,i2n,mLP)  # generate B__struct (Benders master structure)
        S = gen_S2_type(B.data,mSP) # generate S__struct (Benders subproblem structure)
    end
    if a == 3 # (if Zaker_Bend)
        B = gen_B3_type(c,a,q,j,e,emb,stab,γs,ms,mp,ps,pp,unc,mMP,i2n)  # generate B__struct (Benders master structure)
        S = gen_S3_type(B.data,mSP) # generate S__struct (Benders subproblem structure)
    end

    return B,S
end

# */ --------------------------------------------------------------------------------------------- /* #
# */ --- lower level functions : subproblem ------------------------------------------------------ /* #
# */ --------------------------------------------------------------------------------------------- /* #

function set_Iex!(b::B2_type)::B2_type

    # function to update B2_type.temp.Ie (index of w subproblems solved exactly @ iter k)
    # for given B2_type (Benders master problem) -> algorithm 2 (Adapt_Bend)

    if b.data.ad==1 b.data.w=1 end



    for k in 1:b.data.w # for each k = 1,..,w (for each subset ℰₖ)
        #***** add switch

        iv,ie = findmax([b.data.mp.kappa[i]*b.data.mp.prob[i] for i in ms.I][b.data.kidx[k]].*(b.temp.θu[b.data.kidx[k]].-b.temp.θl[b.data.kidx[k]]))
        relGaps=([b.data.mp.kappa[i]*b.data.mp.prob[i] for i in ms.I][b.data.kidx[k]].*(b.temp.θu[b.data.kidx[k]].-b.temp.θl[b.data.kidx[k]]))./iv

        if b.data.ad==1
            #ie=rand((findall(relGaps.>=0.99999))); #RANDOMLY SELECT SUBPROBLEM IF THE GAPS ARE NUMERICALLY SIMILAR
            #ie=findall(relGaps.>=0.99999)[1]; #ALWAYS SELECT THE FIRST OF THE LIST WHEN TEH GAPS ARE TEH SAME
            ie=findall(relGaps.>=0.99999)[end]; #ALWAYS SELECT THE LAST OF THE LIST WHEN TEH GAPS ARE TEH SAME
        end
        b.temp.Ie[k] = b.data.iidx[k][ie] # store ie in B2_type.temp.Ie


    end

    return b
end

function solv_exact!(s::S1_type)::S1_type

    # function to update and solve S1_type.ex.m (exact subproblem) and store exact solution in S1_type.temp
    # for given S1_type (Benders subproblem) -> algorithm 1 (Stand_Bend)

    s.ex.objf =  s.ex.objfor + s.ex.vars'*s.temp.c[:,s.temp.i] # compute new objective of subproblem i (c₀ + cᵢᵀϕ)
    @objective(s.ex.m, Min, s.ex.objf) # update objective in JuMP model (c₀ + cᵢᵀϕ)
    fix.(s.ex.m[:x],s.temp.x[:,s.temp.i];force=true) # fix values of variables xᵢ in the model to decisions xᵢ set by the master problem
    
    optimize!(s.ex.m) # solve the JuMP model to optimality
    s.temp.θ[s.temp.i] = objective_value(s.ex.m) # store optimal objective θᵢ
    s.temp.λ[:,s.temp.i] .= dual.(FixRef.(s.ex.m[:x])) # store subgradient λᵢ wrt xᵢ

    return s
end


function solv_exact!(s::S2_type)::S2_type

    # function to update and solve S2_type.ex.m (exact subproblem) and store exact solution in S2_type.temp
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.ex.objf =  s.ex.objfor + s.ex.vars'*s.temp.c[:,s.temp.i] # compute new objective of subproblem i (c₀ + cᵢᵀϕ)
    #s.ex.objf =  s.ex.vars'*vcat([1.],s.temp.c[:,s.temp.i]) # compute new objective of subproblem i (c₀ + cᵢᵀϕ)
    @objective(s.ex.m, Min, s.ex.objf) # update objective in JuMP model (c₀ + cᵢᵀϕ)

    fix.(s.ex.m[:x],s.temp.x[:,s.temp.i];force=true) # fix values of variables xᵢ in the model to decisions xᵢ set by the master problem



    optimize!(s.ex.m) # solve the JuMP model to optimality
    s.temp.θl[s.temp.i] = objective_value(s.ex.m) # store optimal objective θᵢ
    s.temp.λl[:,s.temp.i] .= dual.(FixRef.(s.ex.m[:x])) # store subgradient λᵢ wrt xᵢ
    s.temp.ϕu[:,s.temp.i] .= value.(s.ex.m[:phi]) # store subgradient ϕᵢ wrt cᵢ
    ####
    # s.temp.θu[s.temp.i] = objective_value(s.ex.m)

    return s
end

function extended_sol!(s::S2_type)::S2_type

    # function to store new extended exact solution to S2_type.q (extended exact solutions stored)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    ϕ  = value.(s.ex.m[:phi] ) # extract the operational cost dependent on uncertain cost cᵢ
    c0 = value.(s.ex.m[:c0]) # extract the operational cost independent of uncertain costs c
    yG = value.(s.ex.m[:yG]) # extract the generation level of conventional generation
    yI = value.(s.ex.m[:yI]) # extract the charging power of storage generation
    yO = value.(s.ex.m[:yI]) # extract the discharging power of storage generation
    yL = value.(s.ex.m[:yL]) # extract the energy level of storage generation
    yS = value.(s.ex.m[:yS]) # extract the load shedding
    #x0 = value.(s.ex.m[:x0]) # extract the values of rhs parameters in the subproblem (eg, capacity)
    x0 = value.(s.ex.m[:x][1:s.data.unc.nx0]) # extract the values of rhs parameters in the subproblem (eg, capacity)
    #push!(s.q,Q__type(ϕ,c0,yG,yI,yO,yL,yS,x0)) # add new extended exact solution

    return s
end

function update_m_!(s::S2_type)::S2_type

    # function to add a new exact solution to S2_type.m (exact solutions stored)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.m.n += 1 # update number of solution stored (+1)
    s.m.θ[s.m.n] = s.temp.θl[s.temp.i] # store optimal solution θᵢ
    s.m.λ[s.m.n,:] .= s.temp.λl[:,s.temp.i] # store subgradient λᵢ wrt xᵢ
    s.m.ϕ[s.m.n,:] .= s.temp.ϕu[:,s.temp.i] # store subgradient ϕᵢ wrt cᵢ
    s.m.x[s.m.n,:] .= s.temp.x[:,s.temp.i] # store xᵢ
    s.m.c[s.m.n,:] .= s.temp.c[:,s.temp.i] # store cᵢ

    return s
end

function update_lb!(s::S2_type)::S2_type

    # function to add a new exact solution to S2_type.olb (oracle for lower bound)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.olb.cons[s.m.n] = @constraint(s.olb.m, s.olb.m[:phi] + s.olb.m[:gamma]'*s.m.c[s.m.n,:] >= 0. ) # add constraint (ϕ + γᵀcₙ >= 0) to the oracle
    s.olb.help[1:s.m.n] .= s.m.θ[1:s.m.n] .- sum!(ones(s.m.n),s.m.λ[1:s.m.n,:].*s.m.x[1:s.m.n,:]) # compute vector h used when solving the lb oracle (hₛ = θₛ - λₛᵀxₛ, ∀s=1,..,n)
    s.olb.n += 1 # update number of solution stored (+1)

    return s
end

function update_lbP!(s::S2_type)::S2_type


    set_normalized_coefficient.(s.olbP.m[:csr_lim],s.olbP.m[:μ][s.m.n],s.m.c[s.m.n,:])
    s.olbP.n += 1 # update number of solution stored (+1)

    unfix(s.olbP.m[:μ][s.olbP.n]);
    set_lower_bound(s.olbP.m[:μ][s.olbP.n],0.0);

    return s
end

function update_ub!(s::S2_type)::S2_type

    # function to add a new exact solution in S2_type.olb (oracle for upper bound)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.oub.cons[s.m.n] = @constraint(s.oub.m, s.oub.m[:phi] + s.oub.m[:gamma]'*s.m.x[s.m.n,:] <= 0. ) # add constraint (ϕ + γᵀxₙ <= 0) to the oracle
    s.oub.help[1:s.m.n] .= s.m.θ[1:s.m.n] .- sum!(ones(s.m.n),s.m.ϕ[1:s.m.n,:].*s.m.c[1:s.m.n,:]) # compute vector h used when solving the ub oracle (hₛ = θₛ - ϕₛᵀcₛ, ∀s)
    s.oub.n += 1 # update number of solution stored (+1)

    return s
end

function update_ubP!(s::S2_type)::S2_type

    #for x2 in 1:s.data.unc.nx
    #    set_normalized_coefficient(s.oubP.m[:csr_lim][x2],s.oubP.m[:ν][s.m.n],s.m.x[s.m.n,x2])
    #end

    set_normalized_coefficient.(s.oubP.m[:csr_lim],s.oubP.m[:ν][s.m.n],s.m.x[s.m.n,:])
    s.oubP.n += 1 # update number of solution stored (+1)

    unfix(s.oubP.m[:ν][s.oubP.n]);
    set_lower_bound(s.oubP.m[:ν][s.oubP.n],0.0);

    return s
end

function update_s!(s::S2_type)::S2_type

    # function to add a new exact solution to S2_type.m (exact solutions stored), S2_type.olb (oracle for lower bound), and S2_type.olb (oracle for upper bound)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    update_m_!(s) # store new exact solution in S2_type.m
    update_lb!(s) # store new exact solution in lower bound oracle
    update_ub!(s) # store new exact solution in upper bound oracle

    return s
end


function update_s!(b::B2_type,s::S2_type)::Tuple{B2_type,S2_type,Float64}
    # function to add a new exact solution to S2_type.m (exact solutions stored), S2_type.olb (oracle for lower bound), and S2_type.olb (oracle for upper bound)
    update_m_!(s) # store new exact solution in S2_type.m


    #***** add switch, for primal

    #t= @elapsed if b.data.emb==1 && s.m.n>1
    t= @elapsed if b.data.emb==1 && b.hist.k>=1
            for i in eachindex(b.data.ms.I)
                #@constraint(b.rmp.m, b.rmp.m[:phi][i] + b.rmp.m[:γ][1:b.data.unc.nc,i]'*s.m.c[s.m.n,1:b.data.unc.nc] >= s.m.θ[s.m.n]+s.m.λ[s.m.n,:]'*(b.rmp.m[:x][:,i]-s.m.x[s.m.n,:]))
                @constraint(b.rmp.m, b.rmp.m[:phi][i] + sum(b.rmp.m[:gamma][c,i]'*s.m.c[s.m.n,c] for c in 1:b.data.unc.nc) >= s.m.θ[s.m.n]+s.m.λ[s.m.n,:]'*(b.rmp.m[:x][:,i]-s.m.x[s.m.n,:]))
            #set_name(all_constraints(b.rmp.m,AffExpr,MOI.GreaterThan{Float64})[end],"csrCut_$(s.m.n)_$(i)")
            end
        end


    t+= @elapsed if b.data.emb==1 && b.hist.k>=1 && b.data.stab==1
        for i in eachindex(b.data.ms.I)
            #@constraint(b.rmp.m, b.rmp.m[:phi][i] + b.rmp.m[:γ][1:b.data.unc.nc,i]'*s.m.c[s.m.n,1:b.data.unc.nc] >= s.m.θ[s.m.n]+s.m.λ[s.m.n,:]'*(b.rmp.m[:x][:,i]-s.m.x[s.m.n,:]))
            @constraint(b.lmp.m, b.lmp.m[:phi][i] + sum(b.lmp.m[:gamma][c,i]'*s.m.c[s.m.n,c] for c in 1:b.data.unc.nc) >= s.m.θ[s.m.n]+s.m.λ[s.m.n,:]'*(b.lmp.m[:x][:,i]-s.m.x[s.m.n,:]))
        #set_name(all_constraints(b.rmp.m,AffExpr,MOI.GreaterThan{Float64})[end],"csrCut_$(s.m.n)_$(i)")
        end
    end


    update_lb!(s) # store new exact solution in lower bound oracle
    update_ub!(s) # store new exact solution in upper bound oracle

    # solve primal
    #update_lbP!(s) # calculate constraint coefficients
    #update_ubP!(s) # calculate constraint coefficients

    return b,s,t
end


function run_lb!(s::S2_type)::S2_type

    # function to solve S2_type.olb (oracle for lower bound) and store the solutions in S2_type.temp (temporary data)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.olb.objf = s.olb.vars'*vcat([1.],s.temp.c[:,s.temp.i]) # compute objective function for subproblem i (ϕ + γᵀcᵢ)
    @objective(s.olb.m, Min, s.olb.objf) # set objective function
    set_normalized_rhs.(s.olb.cons[1:s.olb.n], s.olb.help[1:s.olb.n].+ s.m.λ[1:s.olb.n,:]*s.temp.x[:,s.temp.i]) # update the rhs of the constraints (ϕ + γᵀcₛ >= hₛ + λₛᵀxᵢ, ∀s)
    optimize!(s.olb.m) # solve the JuMP model to optimality
    s.temp.θl[s.temp.i] = objective_value(s.olb.m) # store the valid lower bound θᵢ
    s.temp.λl[:,s.temp.i] .= s.m.λ[1:s.m.n,:]'*dual.(s.olb.cons[1:s.olb.n]) # store the valid subgradient λᵢ wrt xᵢ

    #println("Dual LBoracle: ",objective_value(s.olb.m))


    return s
end

function run_ub!(s::S2_type)::S2_type

    # function to solve S2_type.oub (oracle for upper bound) and store the solutions in S2_type.temp (temporary data)
    # for given S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.oub.objf = s.oub.vars'*vcat([1.],s.temp.x[:,s.temp.i]) # compute objective function for subproblem i (ϕ + γᵀxᵢ)
    #println("s.oub.objf: ",s.oub.objf)
    @objective(s.oub.m, Max, s.oub.objf) # set objective function
    # println(s.m.n)
    set_normalized_rhs.(s.oub.cons[1:s.oub.n], s.oub.help[1:s.oub.n].+ s.m.ϕ[1:s.oub.n,:]*s.temp.c[:,s.temp.i])  # update the rhs of the constraints (ϕ + γᵀxₛ <= hₛ + ϕₛᵀcᵢ, ∀s)
    #push!(N.UBOm,DataFrame(vcat(hcat([1.],s.temp.x[:,s.temp.i]',[1.]), hcat(Matrix(JuMP._standard_form_matrix(s.oub.m)[:A])[:,1:B.data.unc.nx+1],s.oub.help[1:s.oub.n].+ s.m.ϕ[1:s.oub.n,:]*s.temp.c[:,s.temp.i])),:auto))

    optimize!(s.oub.m) # solve the JuMP model to optimality
    #push!(N.ubo_sol, vcat(value(s.oub.m[:phi]),value.(s.oub.m[:gamma])))
    #push!(N.ubo_dual, dual.(all_constraints(s.oub.m,AffExpr,MOI.LessThan{Float64})))

    # println("Contribution of Special point ot UBO: $(dual(all_constraints(s.oub.m,AffExpr,MOI.LessThan{Float64})[1]))");

    if termination_status(s.oub.m)!=MOI.OPTIMAL
        println("Adjusting the algorithm for UBOracle")
        JuMP.set_optimizer_attribute(m,"method",2)
        JuMP.set_optimizer_attribute(m,"crossover",0)
        JuMP.set_optimizer_attribute(m,"Presolve",0)  #
        optimize!(s.oub.m)
        if termination_status(s.oub.m)!=MOI.OPTIMAL
            println("Readjusting the algorithm for UBOracle")
            JuMP.set_optimizer_attribute(m,"BarHomogeneous",1)
            JuMP.set_optimizer_attribute(m,"NumericFocus",3)
            JuMP.set_optimizer_attribute(m,"Presolve",0)  #
            optimize!(s.oub.m)
            if termination_status(s.oub.m)!=MOI.OPTIMAL
                throw("UB oracle has failed")
            end
        end
        JuMP.set_optimizer_attribute(m,"BarHomogeneous",0)
        JuMP.set_optimizer_attribute(m,"method",-1)
    end


    s.temp.θu[s.temp.i] = objective_value(s.oub.m) # store the valid upper bound θᵢ
    s.temp.ϕu[:,s.temp.i] .= s.m.ϕ[1:s.m.n,:]'*dual.(s.oub.cons[1:s.oub.n]) # store the valid subgradient ϕᵢ wrt cᵢ

    #println("Dual UBoracle: ",s.temp.θu[s.temp.i])

    return s
end


function run_lbP!(s::S2_type)::S2_type


    s.olbP.objf= sum(s.olbP.m[:μ][s2]*(s.m.θ[s2]+s.m.λ[s2,:]'*(s.temp.x[:,s.temp.i]-s.m.x[s2,:])) for s2 in 1:s.olbP.n);
    @objective(s.olbP.m,Max,s.olbP.objf)
    set_normalized_rhs.(s.olbP.m[:csr_lim],s.temp.c[:,s.temp.i])
    optimize!(s.olbP.m)
    s.temp.θl[s.temp.i] = objective_value(s.olbP.m) # store the valid upper bound θᵢ
    s.temp.λl[:,s.temp.i] .= s.m.λ[1:s.m.n,:]'*value.(s.olbP.m[:μ][1:s.olbP.n]);

    #println("Primal LBoracle: ",objective_value(s.olbP.m))

    return s
end


function run_ubP!(s::S2_type)::S2_type


    s.oubP.objf= sum(s.oubP.m[:ν][s2]*(s.m.θ[s2]+s.m.ϕ[s2,:]'*(s.temp.c[:,s.temp.i]-s.m.c[s2,:])) for s2 in 1:s.oubP.n);
    @objective(s.oubP.m,Min,s.oubP.objf)
    set_normalized_rhs.(s.oubP.m[:csr_lim],s.temp.x[:,s.temp.i])
    optimize!(s.oubP.m)
     s.temp.θu[s.temp.i] = objective_value(s.oubP.m) # store the valid upper bound θᵢ
     s.temp.ϕu[:,s.temp.i] .= s.m.ϕ[1:s.m.n,:]'*value.(s.oubP.m[:ν][1:s.oubP.n]);

    #println("Primal UBoracle: ",objective_value(s.oubP.m))

    return s
end




# */ --------------------------------------------------------------------------------------------- /* #


function step_a!(b::Union{B1_type,B2_type},s::Union{S1_type,S2_type},n::Union{N1_type,N2_type})

    # function to perform 'step_a' of the Std. and Adaptive Benders algorithm (solve the relaxed master problem and store decisions {xᵢ, ∀i})
    # for given B1_type, B2_type (Benders master problem) and S1_type, S2_type (Benders subproblem) -> algorithms 1,2 (Stand_Bend, Adaptive_Bend)
    optimize!(b.rmp.m) # solve the relaxed master problem to optimality
    b.temp.β .= value.(b.rmp.m[:beta]).data # store optimal values of {βᵢ, ∀i}
    b.hist.L[b.hist.k] = objective_value(b.rmp.m) # compute and store the lower bound Lₖ on the optimal objective @ iter k
    b.temp.xr = value.(b.rmp.m[:x]);

    return b,s,n
end

function step_b!(b::Union{B1_type,B2_type},n::Union{N1_type,N2_type})

    # function to perform 'step_b' of the Benders algorithm (compute lower bound on optimal objective)
    # for given B1_type,B2_type (Benders master problem) -> algorithm 1,2 (Std_Bend,Adapt_Bend)

    ##### LMP stabilisation
    if b.data.stab==1
        if b.hist.k==1
            b.temp.fr = value(b.rmp.m[:f]);
            b.temp.x = value.(b.rmp.m[:x]);
            b.temp.aux=val2dict(b.rmp.m[:pN]);

            #push!(n.LMP_target,b.hist.L[b.hist.k])
        else
            b.temp.lf = (b.hist.L[b.hist.k] + b.data.γs*(b.hist.U[b.hist.k-1] - b.hist.L[b.hist.k]))
            #push!(n.LMP_target,b.temp.lf)
            fix(b.lmp.m[:lf],b.temp.lf;force=true)

            @objective(b.lmp.m, Min, sum(b.data.mp.prob[i]*(b.lmp.m[:pN][g,i] - b.temp.aux[g][i])^2 for g in b.data.ms.M for i in b.data.ms.I0)) # update LMP objective (distance wrt decisions at previous iteration)

            optimize!(b.lmp.m);

            b.temp.x = value.(b.lmp.m[:x]);
            b.temp.fr = value(b.lmp.m[:f]);
            b.temp.aux=val2dict(b.lmp.m[:pN]);

        end

    else
        b.temp.x = value.(b.rmp.m[:x]);
        b.temp.fr = value(b.rmp.m[:f]);
    end

    #push!(n.xlmp,b.temp.x)

    return b,n
end



function step_c!(b::B1_type,s::S1_type,n::N1_type)::Tuple{B1_type,S1_type,N1_type}

    # function to perform 'step_c' of the Benders algorithm (solve each subproblem exactly and store optimal solutions)
    # for given B1_type (Benders master problem) and S1_type (Benders subproblem) -> algorithm 1 (Stand_Bend)
    s.temp.x .= b.temp.x

    adj=0.0
    for s.temp.i in 1:b.data.unc.ni
        solv_exact!(s) # solve subproblem i exactly and store solutions
    end
    b.temp.θ .= s.temp.θ # send optimal objectives {θᵢ, ∀i} to the Benders master problem
    b.temp.λ .= s.temp.λ # send subgradients {λᵢ, ∀i} to the Benders master problem

    #b.data.adj=adj

    return b,s,n
end


function step_c!(b::B2_type,s::S2_type,n::N2_type)::Tuple{B2_type,S2_type,N2_type}

    #***** add switch for solve all SP exactly, also a switch in the report function, maybe best to keep one reprot function but with switches instead of having multiple report functions.
    ### This needs to be updated to properly do more than 1 subproblem when not adaptive.

    for i in 1:b.data.unc.nx0 for j in 1:b.data.unc.ni
        if b.temp.x[i,j]<=exp10(-6)
            b.temp.x[i,j]=0
        end
    end end
    s.temp.x .= b.temp.x
    b.temp.ne =0; #Set exploration counter

    for s.temp.i in 1:b.data.unc.ni
        run_lb!(s)
        run_ub!(s)
        # run_lbP!(s) # call the oracle for lower bound at (xᵢ,cᵢ) (Primal version)
        # run_ubP!(s) # call the oracle for upper bound at (xᵢ,cᵢ) (Primal version)
    end

    b.temp.θl .= s.temp.θl #Pass the lower bound oracle for the s.temp.x point to the MP
    b.temp.θu .= s.temp.θu #Pass the upper bound oracle for the s.temp.x point to the MP
    b.temp.ne =0; #Set exploration counter

    #push!(n.LBOα,(b.temp.fr+sum(b.data.mp.df[i]*b.data.mp.kappa[i]*b.data.mp.prob[i]*b.temp.θl[b.data.unc.n2i[i]] for i in eachindex(b.data.unc.n2i))));
    #push!(n.UBOα,(b.temp.fr+sum(b.data.mp.df[i]*b.data.mp.kappa[i]*b.data.mp.prob[i]*b.temp.θu[b.data.unc.n2i[i]] for i in eachindex(b.data.unc.n2i))));
    #push!(n.LBOiα,[b.data.mp.df[i]*b.data.mp.prob[i] for i in ms.I][b.data.kidx[1]].*b.temp.θl[b.data.kidx[1]]);
    #push!(n.UBOiα,[b.data.mp.df[i]*b.data.mp.prob[i] for i in ms.I][b.data.kidx[1]].*b.temp.θu[b.data.kidx[1]]);
    #push!(n.gapi,[b.data.mp.df[i]*b.data.mp.prob[i] for i in ms.I][b.data.kidx[1]].*(b.temp.θu[b.data.kidx[1]].-b.temp.θl[b.data.kidx[1]]))     #store the gaps of all SP

    if b.hist.k==1
        b.temp.bub=b.hist.U[1]
        b.temp.Δ=b.temp.Δ
    else
        b.temp.bub=b.hist.U[b.hist.k-1]
        b.temp.Δ=(b.temp.bub-b.hist.L[b.hist.k])/b.temp.bub*100
    end


    for i_ne in 1:b.data.unc.ni
        set_Iex!(b)
        for s.temp.i in b.temp.Ie
            solv_exact!(s)
            push!(n.sSP,b.temp.Ie[1])
            b.temp.ne +=1
            if b.temp.ne==1
                push!(n.vSP,[b.temp.Ie[1]])
            else
                push!(n.vSP[end],b.temp.Ie[1])
            end
            update_s!(b,s)


        end

        a=time();
        for s.temp.i in 1:b.data.unc.ni
            run_lb!(s)
            run_ub!(s)
        end


        b.temp.θl .= s.temp.θl #Pass the lower bound oracle for the s.temp.x point to the MP
        b.temp.θu .= s.temp.θu #Pass the upper bound oracle for the s.temp.x point to the MP


        b.temp.lb=b.temp.fr+sum(b.data.mp.kappa[i]*b.data.mp.prob[i]*b.temp.θl[b.data.unc.n2i[i]] for i in eachindex(b.data.unc.n2i))
        b.temp.ub=b.temp.fr+sum(b.data.mp.kappa[i]*b.data.mp.prob[i]*b.temp.θu[b.data.unc.n2i[i]] for i in eachindex(b.data.unc.n2i))
        b.temp.δl=(b.temp.ub-b.temp.lb)/(b.temp.ub)*100.0;

        #push!(n.LBO,b.temp.lb);
        #push!(n.UBO,b.temp.ub);
        #push!(n.LBOi,[b.data.mp.df[i]*b.data.mp.prob[i] for i in ms.I][b.data.kidx[1]].*b.temp.θl[b.data.kidx[1]]);
        #push!(n.UBOi,[b.data.mp.df[i]*b.data.mp.prob[i] for i in ms.I][b.data.kidx[1]].*b.temp.θu[b.data.kidx[1]]);
        #push!(n.gapis,[b.data.mp.df[i]*b.data.mp.prob[i] for i in ms.I][b.data.kidx[1]].*(b.temp.θu[b.data.kidx[1]].-b.temp.θl[b.data.kidx[1]]))     #store the gaps of all SP


        if b.temp.lb>b.temp.bub
            break:nothing
        elseif b.temp.δl <= b.temp.Δ
            break:nothing
        elseif b.data.ad==0
            break:nothing
        end
    end

    #push!(n.nSP,b.temp.ne);


    b.temp.λl .= s.temp.λl # send valid subgradients {λᵢ, ∀i} to the Benders master problem
    b.temp.ϕu .= s.temp.ϕu # send valid subgradients {ϕᵢ, ∀i} to the Benders master problem

    return b,s,n
end


function step_e!(b::Union{B1_type,B2_type},n::Union{N1_type,N2_type})


    # function to perform 'step_e' of the Benders algorithm (compute upper bound on optimal objective)
    # for given B1_type, B2_type (Benders master problem) -> algorithms 1,2 (Stand_Bend,Adapt_Bend)
    θ::Array{Float64,1}=Array{Float64,1}(); U::Float64=Inf;

    b.data.algm == 1 ? θ=b.temp.θ : θ=b.temp.θu
    U=(b.temp.fr+sum(b.data.mp.kappa[i]*b.data.mp.prob[i]*θ[b.data.unc.n2i[i]] for i in eachindex(b.data.unc.n2i))) # compute and store the upper bound Uₖ on the optimal objective @ iter k = 1
    b.hist.k == 1 ? b.hist.U[b.hist.k]=U : b.hist.U[b.hist.k]=min(b.hist.U[b.hist.k-1],U)

    return b,n
end

function step_f!(b::Union{B1_type,B2_type},n::Union{N1_type,N2_type})

  # function to perform 'step_f' of the Benders algorithm (compute optimality gap Δ and update the relaxed master problem) and add Cuts
  # for given B1_type,B2_type (Benders master problem) -> algorithms 1,2 (Std_Bend,Adapt_Bend)

  b.temp.Δ = (b.hist.U[b.hist.k]-b.hist.L[b.hist.k])/b.hist.U[b.hist.k]*100 # compute optimality gap Δ = (Uₖ-Lₖ)/Uₖ @ iter k

  θ::Array{Float64,1}=Array{Float64,1}();
  λ::Array{Float64,2}=Array{Float64}(undef, 0, 0);

  b.data.algm == 1 ? θ=b.temp.θ : θ=b.temp.θl;
  b.data.algm == 1 ? λ=b.temp.λ : λ=b.temp.λl;

  avoid::Bool=false
  newC::Int64=0;

  #push!(n.δ,b.temp.Δ)
  #push!(n.θc,deepcopy(θ))

  if b.data.ac==0 # if avoid adding cuts
      for i in 1:b.data.unc.ni
          @constraint(b.rmp.m,b.rmp.m[:beta][b.data.unc.i2n[i]] >= θ[i]+λ[:,i]'*(b.rmp.m[:x][:,i].-b.temp.x[:,i])) # add constraints βᵢ >= θᵢ + λᵢᵀ(𝑥ᵢ-xᵢ) to the relaxed master problem
          if b.data.stab==1
              @constraint(b.lmp.m,b.lmp.m[:beta][b.data.unc.i2n[i]] >= θ[i]+λ[:,i]'*(b.lmp.m[:x][:,i].-b.temp.x[:,i])) # add constraints βᵢ >= θᵢ + λᵢᵀ(𝑥ᵢ-xᵢ) to the relaxed master problem
          end
      end
  else
      push!(n.λc,deepcopy(λ))
      push!(n.hc,deepcopy(θ.-[λ[:,i]'*b.temp.x[:,i] for i in 1:b.data.unc.ni]));
      push!(n.sc,[vcat(θ[i]-λ[:,i]'*b.temp.x[:,i],λ[:,i])'*n.rv for i in 1:b.data.unc.ni])

      newC=0;
      for i in 1:b.data.unc.ni
          avoid=false
          for k in 1:b.hist.k-1
              if abs(n.sc[b.hist.k][i] - n.sc[k][i])<=1e-8
                 if norm(vcat(θ[i]-λ[:,i]'*b.temp.x[:,i],λ[:,i])-vcat(n.hc[k][i],n.λc[k][:,i]))<=1e-8
                     avoid=true
                     break
                 end
             end
          end
          if avoid==false
              @constraint(b.rmp.m,b.rmp.m[:beta][b.data.unc.i2n[i]] >= θ[i]+λ[:,i]'*(b.rmp.m[:x][:,i].-b.temp.x[:,i])) # add constraints βᵢ >= θᵢ + λᵢᵀ(𝑥ᵢ-xᵢ) to the relaxed master problem
              #set_name(all_constraints(b.rmp.m,AffExpr,MOI.GreaterThan{Float64})[end],"Cut_$(b.hist.k)_$i");
              newC+=1;
              if b.data.stab==1
                  @constraint(b.lmp.m,b.lmp.m[:beta][b.data.unc.i2n[i]] >= θ[i]+λ[:,i]'*(b.lmp.m[:x][:,i].-b.temp.x[:,i])) # add constraints βᵢ >= θᵢ + λᵢᵀ(𝑥ᵢ-xᵢ) to the relaxed master problem
              end
          end
      end
  end
end




# */ --------------------------------------------------------------------------------------------- /* #

function step_0!(b::B2_type,s::S2_type,n::N2_type)::Tuple{B2_type,S2_type,N2_type}

    # function to perform 'step_0' of the Benders algorithm (solve subproblem exactly at special point (x_min,c_min) to make the adaptive oracles always feasible)
    # for given B2_type (Benders master problem) and S2_type (Benders subproblem) -> algorithm 2 (Adapt_Bend)

    s.temp.i = 1 # set i = 1
    # s.temp.x.=deepcopy(round.(b.temp.x;digits=4))
    s.temp.x.=deepcopy(b.temp.x)
    s.temp.c.=deepcopy(b.temp.c)
    # for s.temp.i in 1:b.data.unc.ni
    #     run_oracles!(s)
    # end
    #println("First time it solves the special point")
    solv_exact!(s) # solve subproblem 1 exactly (at point (x_min,c_min)) and store solutions
    # println(objective_value(s.ex.m))
    #println("SPECIAL POINT FINISHED")

    update_s!(b,s) # send the new exact solution to the adaptive oracles


    #### Special functions for mp.xh
    #mpxh=zeros(length(b.data.ms.M),length(b.data.ms.I))
    #mpxh2=zeros(length(b.data.ms.M),length(b.data.ms.I))
    #arr=[[b.data.mp.xh[i][n] for i in b.data.ms.M] for n in b.data.ms.I];
    #for j in eachindex(ms.I) mpxh[:,j]=arr[j,:][1] end
    #for j in eachindex(ms.M)
    #    mpxh2[j,:]=mpxh[j,:]
       # mpxh2[j,:]=(1+0.5rand())*mpxh[j,:]
    #end
    #mpxh=mpxh2;
    #s.temp.x .= vcat(mpxh,b.data.unc.h') # set decisions {xᵢ, ∀i} on historial data
    #s.temp.c .= b.data.unc.c' # set {cᵢ, ∀i}

    #***** add switch for evaluating more special points.


    #for s.temp.i in 1:b.data.unc.ni
       # solv_exact!(s)
       #run_lb!(s)
       #run_ub!(s) # solve the adaptive oracles for subproblem i and store the solutions
       #
       # println("dual co2 budget:,", dual(S.ex.m[:c21]))
    #end
    #b.temp.θl .= s.temp.θl # send valid lower bound {θᵢ, ∀i} to the Benders master problem
    #b.temp.θu .= s.temp.θu # send valid upper bound {θᵢ, ∀i} to the Benders master problem
    #b.temp.λl .= s.temp.λl # send valid subgradients {λᵢ, ∀i} to the Benders master problem
    #b.temp.ϕu .= s.temp.ϕu # send valid subgradients {ϕᵢ, ∀i} to the Benders master problem

    return b,s,n
end

function cutsignature(b::Union{B1_type,B2_type},n::Union{N1_type,N2_type})
    if b.data.ac==1
        for i in 1:b.data.unc.nx+1
            push!(n.rv, 0.5.+0.5*rand())
        end
    end
end


function step!(b::Union{B1_type,B2_type},s::Union{S1_type,S2_type},n::Union{N1_type,N2_type})

    # function to perform a Benders 'step'
    # for given B1_type or B2_type (Benders master problem) and S1_type or S2_type (Benders subproblem) -> algorithms 1,2 (Std,Adapt_Bend)

    b.hist.k += 1                                # update iter k number
    b.hist.T[b.hist.k,1] = @elapsed step_a!(b,s,n) # solve the relaxed master problem and store decisions {xᵢ, ∀i}
    b.hist.T[b.hist.k,2] = @elapsed step_b!(b,n)   # compute lower bound on optimal objective
    b.hist.T[b.hist.k,3] = @elapsed step_c!(b,s,n) # solve w subproblem exactly and store optimal solutions
    b.hist.T[b.hist.k,4] = 0.0;
    b.hist.T[b.hist.k,5] = @elapsed step_e!(b,n)   # compute upper bound on optimal objective
    b.hist.T[b.hist.k,6] = @elapsed step_f!(b,n)   # compute optimality gap Δ and update the relaxed master problem

    #push!(n.tRMP, b.hist.T[b.hist.k,1])
    #push!(n.tLMP, b.hist.T[b.hist.k,2])
    #push!(n.tSP, b.hist.T[b.hist.k,3])
    #push!(n.ttotal, sum(b.hist.T[b.hist.k,:]))

    return b,s,n
end

# */ --------------------------------------------------------------------------------------------- /* #


function solve_Benders!(b::Union{B1_type,B2_type},s::Union{S1_type,S2_type},n::Union{N1_type,N2_type};print_output::Int64=1)

    print_output == 1 ? print0(b) : nothing
    b.data.algm ==2 ? step_0!(b,s,n) : nothing # perform step0 to solve the special point (x_min,c_min)
    b.data.emb == 1 ? b=embedLBOracle!(b,s) : nothing # embed LBO
    b.data.ac == 1 ? cutsignature(b,n) : nothing #Create signature for cut duplication avoidance

    while (b.hist.k < b.data.J) # check the iteration number k is lower than the limit J
        step!(b,s,n) # perform a Benders step
        print_output == 1 ? try print_info(b,n) catch; end : nothing
        # println("---------------------------------------------------------------------------------------------------------------------------------------------------")
        (b.temp.Δ <= b.data.δ) ? break : nothing # check if the optimality gap Δₖ is lower than the targer tolerance δ
    end
    print_output == 1 ? print_summary(b,n) : nothing
    #print_output == 1 ? report(b,n) : nothing # Hongyu needs to check if this is consistent

    return b,s,n
end

# */ --------------------------------------------------------------------------------------------- /* #



function embedLBOracle!(b::Union{B1_type,B2_type},s::Union{S1_type,S2_type})

    #### RMP
    #@variable(b.rmp.m, phi[b.data.ms.I], container=Array)
    @variable(b.rmp.m, gamma[1:b.data.unc.nc,i in eachindex(b.data.ms.I)] >= 0.0)


    for i in eachindex(b.data.ms.I), nn in [s.m.n]
        #@constraint(b.rmp.m, b.rmp.m[:phi][i] + b.rmp.m[:γ][1:b.data.unc.nc,i]'*s.temp.c[1:b.data.unc.nc,s.m.n] >= s.m.θ[nn]+s.m.λ[nn,:]'*(b.rmp.m[:x][:,i]-s.m.x[nn,:]))
        @constraint(b.rmp.m, b.rmp.m[:phi][i] + sum(b.rmp.m[:gamma][c,i]'*s.m.c[s.m.n,c] for c in 1:b.data.unc.nc) >= s.m.θ[s.m.n]+s.m.λ[s.m.n,:]'*(b.rmp.m[:x][:,i]-s.m.x[s.m.n,:]))
        #set_name(all_constraints(b.rmp.m,AffExpr,MOI.GreaterThan{Float64})[end],"csrCut_$(s.m.n)_$(i)")
    end
    @constraint(b.rmp.m,cOR[i in eachindex(b.data.ms.I)],b.rmp.m[:beta][b.data.unc.i2n[i]]==b.rmp.m[:phi][i]+b.rmp.m[:gamma][:,i]'*b.data.unc.c[i,:])

    if b.data.stab==1
        @variable(b.lmp.m, gamma[1:b.data.unc.nc,i in eachindex(b.data.ms.I)] >= 0.0)
        for i in eachindex(b.data.ms.I), nn in [s.m.n]
            @constraint(b.lmp.m, b.lmp.m[:phi][i] + sum(b.lmp.m[:gamma][c,i]'*s.m.c[s.m.n,c] for c in 1:b.data.unc.nc) >= s.m.θ[s.m.n]+s.m.λ[s.m.n,:]'*(b.lmp.m[:x][:,i]-s.m.x[s.m.n,:]))
        end
        @constraint(b.lmp.m,cOR[i in eachindex(b.data.ms.I)],b.lmp.m[:beta][b.data.unc.i2n[i]]==b.lmp.m[:phi][i]+b.lmp.m[:gamma][:,i]'*b.data.unc.c[i,:])
    end


    return b

end

println("done in $(round(time()-a,digits=1))s ")
