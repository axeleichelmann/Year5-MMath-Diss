################# MODEL LINKING #########################
############################################################


if alg.al!=0
    print("Loading SP model...          "); a=time();
    m = Model(HiGHS.Optimizer)
    MOI.set(m, MOI.Silent(), true)
    #  m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env),"OutputFlag"=>0,"Method"=>2,"Crossover"=>0,"Presolve"=>1))
    SP!(m,ps,pp);                   # subproblem
    println("done in $(round(time()-a,digits=1))s ")

    print("Loading RMP model...         "); a=time();
    m2 = Model(HiGHS.Optimizer)
    MOI.set(m2, MOI.Silent(), true)
    #m2 = Model(optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env),"OutputFlag"=>0,"Method"=>1,"Crossover"=>0,"Presolve"=>1))
    # m2 = Model(optimizer_with_attributes(()->CPLEX.Optimizer(),"CPX_PARAM_QPMETHOD"=>4,"CPXPARAM_SolutionType"=>2))
    RMP!(m2,ms,mp);                 # master problem
    println("done in $(round(time()-a,digits=1))s ")

    if alg.stab==1
        print("Loading LMP model...         "); a=time();
        #m3 = Model(optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env),"OutputFlag"=>0,"Method"=>1,"Crossover"=>0,"Presolve"=>1))
        m3 = Model(HiGHS.Optimizer)
        LMP!(m3,ms,mp);                 # same model as the LMP
        println("done in $(round(time()-a,digits=1))s ")
        else
        m3= Model();
    end
end



############# LINKING MODELS (Models need to be already in memory with the names of the definitions file)
if !fileExists || updateData
    print("Creating linking table...    "); a=time();
    lT,i2n=linkingTable(data["Linkage"]); ### The full group of tables (The information has to be already loaded into memory)
    jldopen("$path/$jldFile", "a+") do f f["lT"]=lT; f["i2n"]=i2n; end
    println("done in $(round(time()-a,digits=1))s ");
else
    print("Restoring linking table...   "); a=time();
    lT,i2n=load.("$path/$jldFile",["lT","i2n"]);
    println("done in $(round(time()-a,digits=1))s ");
end


n2i=Dict(i2n[i]=>i for i in keys(i2n));

if alg.al!=0
    print("Linking models...            "); a=time();
    addLinking!(lT);
    if alg.stab==1 addLinkingLMP!(lT,"m3") end
    println("done in $(round(time()-a,digits=1))s ")

    #### UNC Type
    unc=writeUNC(lT,"u",n2i,i2n);
end


if alg.al==0
    print("Loading and optimising undecomposed model...         "); a=time();
    #mU = Model(optimizer_with_attributes(()->Gurobi.Optimizer(gurobi_env),"OutputFlag"=>1,"Method"=>2,"Crossover"=>1,"Presolve"=>1))
    mU = Model(HiGHS.Optimizer)
    Undecomposed!(mU,ms,mp,lT);                 # master problem
    optimize!(mU);
    println("done in $(round(time()-a,digits=1))s ")
end
