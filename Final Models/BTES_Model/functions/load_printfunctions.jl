# */ --------------------------------------------------------------------------------------------- /* #

function print0(b::JuMP.Model,c::Int64)

    println(" ")
    println(" */"*"-"^ 32*"/*")
    println(" case               : $(c)")
    nv = length(all_variables(b))
    nc = 0
    nc += num_constraints(b,GenericAffExpr{Float64,VariableRef},MOI.EqualTo{Float64})
    nc += num_constraints(b,GenericAffExpr{Float64,VariableRef},MOI.GreaterThan{Float64})
    nc += num_constraints(b,GenericAffExpr{Float64,VariableRef},MOI.LessThan{Float64})
    nc += num_constraints(b,VariableRef,MOI.EqualTo{Float64})
    nc += num_constraints(b,VariableRef,MOI.GreaterThan{Float64})
    nc += num_constraints(b,VariableRef,MOI.LessThan{Float64})
    nvl = length("$nv")
    ncl = length("$nc")
    println(" variables          : $(round(nv*exp10(1-nvl);digits=1)) x 10^$(nvl-1)")
    println(" constraints        : $(round(nc*exp10(1-ncl);digits=1)) x 10^$(ncl-1)")
    println(" */"*"-"^ 32*"/*")
    println(" ")

end

function print0(b::Union{B1_type,B2_type})

    alg::String="";
    b.data.algm == 1 ? alg="std_bend" : alg="adapt_bend"

    println(" ")
    println(" */"*"-"^ 32*"/*")
    println(" algorithm          : " * alg)
    println(" case               : $(b.data.case)")
    if b.data.algm==2
        print(" subs per iter      : ")
        b.data.ad == 1 ? println("adaptive") : println(b.data.w)
    end
    println(" investment  nodes  : $(length(b.data.ms.I0))")
    println(" operational nodes  : $(length(b.data.ms.I))")
    b.data.stab==1 ? println(" stab ON, Factor    : $(b.data.γs)") : nothing
    println(" */"*"-"^ 32*"/*")
end


# */ --------------------------------------------------------------------------------------------- /* #

function print_info(b::B1_type,n::N1_type)

    a1 = "$(round(b.temp.Δ;digits=3))" * "0" ^ (5 - length("$(round(b.temp.Δ%1;digits=3))"))
    a2 = "$(round(sum(b.hist.T[b.hist.k,:]);digits=2))"
    str = " k =" * (" " ^ (4-length("$(b.hist.k)"))) * "$(b.hist.k), "
    str *= "δ =" * (" " ^ (7-length(a1))) * "$(a1) %, "
    str *= "t =" * (" " ^ (8-length(a2))) * "$(a2) s"
    println(str)

end

function print_info(b::B2_type,n::N2_type)

    a1 = "$(round(b.temp.Δ;digits=3))" * "0" ^ (6 - length("$(round(b.temp.Δ%1;digits=3))"))
    a2 = "$(round(sum(b.hist.T[b.hist.k,:]);digits=3))"
    t1 = "$(round(sum(b.hist.T[b.hist.k,1]);digits=3))"
    t2 = "$(round(sum(b.hist.T[b.hist.k,2]);digits=3))"
    t3 = "$(round(sum(b.hist.T[b.hist.k,3]);digits=3))"
    t4 = "$(round(sum(b.hist.T[b.hist.k,4]);digits=3))"
    t5 = "$(round(sum(b.hist.T[b.hist.k,5]);digits=3))"
    t6 = "$(round(sum(b.hist.T[b.hist.k,6]);digits=3))"
    str = " k =" * (" " ^ (4-length("$(b.hist.k)"))) * "$(b.hist.k), "
    str *= "n =" * (" " ^ (4-length("$(b.temp.ne)"))) * "$(b.temp.ne), "
    str *= "δ =" * (" " ^ (8-length(a1))) * "$(a1) %, "
    str *= "t =" * (" " ^ (8-length(a2))) * "$(a2) s, "
    str *= "t1 =" * (" " ^ (8-length(t1))) * "$(t1) s, "
    str *= "t2 =" * (" " ^ (10-length(t2))) * "$(t2) s, "
    str *= "t3 =" * (" " ^ (11-length(t3))) * "$(t3) s, "
    str *= "t4 =" * (" " ^ (6-length(t4))) * "$(t4) s, "
    str *= "t5 =" * (" " ^ (8-length(t5))) * "$(t5) s, "
    str *= "t6 =" * (" " ^ (8-length(t6))) * "$(t6) s "
    #str *= "SP =" * (" ") * "$(n.vSP[b.hist.k])"
    # str *= "$(b.temp.ne) solved exactly"
    println(str)

end

# */ --------------------------------------------------------------------------------------------- /* #

function str_fcn1(x::Array{Float64,1},s::String)::String

    str = " "*s*" "^(9-length(s))
    for i in 1:length(x)
        str *= (" "^(7-length(string(x[i])))*string(x[i]))
    end

    return str
end

function print_summary(b,n)

    #print_summary_a(b)
    print_summary_b(b,n);

end

function print_summary_a(b)

    println(" */"*"-"^ 32*"/*")
    println(" ")
    println(" */"*"-"^ 68*"/*")

    #if b.data.algm==2 && b.data.stab==1
        println("As the problem has been modified, to retrieve similar results the RMP will be solved again")
        optimize!(b.rmp.m)
    #end

    x0 = value.(b.rmp.m[:pC]).data


    x1 = round.(x0[:,1];digits=1)
    x2 = round.(x0[:,2:end];digits=1)
    str_tech = b.data.ms.G
    int1 = " tech."*" "^9*"i1"
    println(" ")
    #println(" co2 emission limit : $( b.data.case <= 0 ? "known" : "uncertain" )")
    #println(" co2 emission cost  : $( b.data.case <= 1 ? "known" : "uncertain" )")
    #println(" uranium cost       : $( b.data.case <= 2 ? "known" : "uncertain" )")
    println(" investment nodes   : $( length(b.data.ms.I0) )")
    println(" operational nodes  : $( length(b.data.ms.I) )")
    println(" optimal objective  : $( round(b.hist.U[b.hist.k]*exp10(-11);digits=3) ) x 10^11 £")
    println(" ")
    println(" */"*"-"^ 68*"/*")
    println(" ")
    println(" optimal investments (GW) @ 0 years")
    println(" " * "-" ^ (length(int1)-1))
    println(int1)
    println(" " * "-" ^ (length(int1)-1))
    for p in b.data.ms.G
        println(str_fcn1(x1[p,:],str_tech[p]))
    end
    println(" "*"-"^(length(int1)-1))
    println(" ")
    println(" optimal investments (GW) @ 5 years")
    #if (b.data.case<=2)
        int2 = " tech.    "
        for n in 1:size(x2)[2]
            int2 *= " "^(7-length("i$n"))*"i$n"
        end
        println(" " * "-" ^ (length(int2)-1))
        println(int2)
        println(" " * "-" ^ (length(int2)-1))
        for p in b.data.ms.G
            println(str_fcn1(x2[p,:],str_tech[p]))
        end
        println(" "*"-"^(length(int2)-1))
    #=
    else
        int2 = [" tech.    "," tech.    "," tech.    "]
        for j in 1:3
            for i in 1:Int64(size(x2)[2]/3)
                int2[j] *= " "^(7-length("i$(9*(j-1)+i)"))*"i$(9*(j-1)+i)"
            end
        end
        for j in 1:3
            println(" " * "-" ^ (length(int2[j])-1))
            println(int2[j])
            println(" " * "-" ^ (length(int2[j])-1))
            for p in b.data.ms.P
                println(str_fcn1(x2[p,9*(j-1)+1:9*j],str_tech[p]))
            end
            println(" " * "-" ^ (length(int2[j])-1))
        end
    end
    =#
    println(" ")
    println(" */"*"-"^ 68*"/*")
end

function print_summary_a(b::JuMP.Model,c::Int64)

    println(" */"*"-"^ 32*"/*")
    println(" ")
    println(" */"*"-"^ 68*"/*")
    x0 = value.(b[:x0])
    x1 = round.(x0[:,1];digits=1)
    x2 = round.(x0[:,2:end];digits=1)
    str_tech = ["coal","coalccs","OCGT","CCGT","diesel","nuclear","pumpL","pumpH","lithium","onwind","offwind","solar"]
    int1 = " tech."*" "^9*"i1"
    println(" ")
    println(" co2 emission limit : $( c <= 0 ? "known" : "uncertain" )")
    println(" co2 emission cost  : $( c <= 1 ? "known" : "uncertain" )")
    println(" uranium cost       : $( c <= 2 ? "known" : "uncertain" )")
    println(" investment nodes   : $( size(b[:x0])[2] )")
    println(" operational nodes  : $( size(b[:x] )[2] )")
    println(" optimal objective  : $( round(objective_value(b)*exp10(-5);digits=3) ) x 10^11 £")
    println(" ")
    println(" */"*"-"^ 68*"/*")
    println(" ")
    println(" optimal investments (GW) @ 0 years")
    println(" " * "-" ^ (length(int1)-1))
    println(int1)
    println(" " * "-" ^ (length(int1)-1))
    for p in 1:size(b[:x0])[1]
        println(str_fcn1(x1[p,:],str_tech[p]))
    end
    println(" "*"-"^(length(int1)-1))
    println(" ")
    println(" optimal investments (GW) @ 5 years")
    if (c<=2)
        int2 = " tech.    "
        for n in 1:size(x2)[2]
            int2 *= " "^(7-length("i$n"))*"i$n"
        end
        println(" " * "-" ^ (length(int2)-1))
        println(int2)
        println(" " * "-" ^ (length(int2)-1))
        for p in 1:size(b[:x0])[1]
            println(str_fcn1(x2[p,:],str_tech[p]))
        end
        println(" "*"-"^(length(int2)-1))
    else
        int2 = [" tech.    "," tech.    "," tech.    "]
        for j in 1:3
            for i in 1:Int64(size(x2)[2]/3)
                int2[j] *= " "^(7-length("i$(9*(j-1)+i)"))*"i$(9*(j-1)+i)"
            end
        end
        for j in 1:3
            println(" " * "-" ^ (length(int2[j])-1))
            println(int2[j])
            println(" " * "-" ^ (length(int2[j])-1))
            for p in 1:size(b[:x0])[1]
                println(str_fcn1(x2[p,9*(j-1)+1:9*j],str_tech[p]))
            end
            println(" " * "-" ^ (length(int2[j])-1))
        end
    end
    println(" ")
    println(" */"*"-"^ 68*"/*")
end

function print_summary_b(b,n)

    achievedGap=(b.hist.U[b.hist.k]-b.hist.L[b.hist.k])/b.hist.U[b.hist.k]*100;

    println(" ")
    println(" computational results:")
    println(" " * "-" ^ 78)
    println(" ϵ-target (%)" * " | " * "iters     time (s)" * "  " * "step_a!(%) step_b!(%) step_c!(%) steps_e,f!(%)")
    println(" " * "-" ^ 78)
    δ = exp10(2)*(b.hist.U[1:b.hist.k].-b.hist.L[1:b.hist.k])./b.hist.U[1:b.hist.k]
    for e in [1.00, 0.10, 0.01, 0.001]
        n =findmin(max.(δ.-e,0))[2]
        t  = sum(b.hist.T[1:n,:])
        tm = sum(b.hist.T[1:n,1])
        tl = sum(b.hist.T[1:n,2])
        ts = sum(b.hist.T[1:n,3])
        to = sum(b.hist.T[1:n,5])+sum(b.hist.T[1:n,6])
        rm = round(exp10(2)*tm/t;digits=2)
        rl = round(exp10(2)*tl/t;digits=2)
        rs = round(exp10(2)*ts/t;digits=2)
        ro = round(exp10(2)*to/t;digits=2)
        str_e = " $(round(e;digits=3))"
        str_n = "$(n)"
        str_t = "$(Int64(round(t;digits=0)))"
        str_e *= "0" ^ (6-length(str_e))
        str_n *= " " ^ (10-length(str_n))
        str_t *= " " ^ (9-length(str_t))
        str_m = "$(rm)"
        str_l = "$(rl)"
        str_s = "$(rs)"
        str_o = "$(ro)"
        str_m = " " ^ (3-findfirst(isequal('.'),str_m)) * str_m
        str_l = " " ^ (3-findfirst(isequal('.'),str_l)) * str_l
        str_s = " " ^ (4-findfirst(isequal('.'),str_s)) * str_s
        str_o = " " ^ (3-findfirst(isequal('.'),str_o)) * str_o
        str_m *= "0" ^ (5-length(str_m))
        str_l *= "0" ^ (5-length(str_l))
        str_s *= "0" ^ (6-length(str_s))
        str_o *= "0" ^ (5-length(str_o))
        str_m *= " " ^ (10-length(str_m))
        str_l *= " " ^ (10-length(str_l))
        str_s *= " " ^ (9-length(str_s))
        if e>achievedGap
            println(str_e * "         | " * str_n * str_t * "  " * str_m * str_l * str_s * str_o)
        end
    end
    println(" " * "-" ^ 78)
    if typeof(b)==B2_type
        println("LB: ", b.hist.L[b.hist.k], ", ", "UB: ", b.hist.U[b.hist.k] )#, ", ", "total SP evaluations: " , length(n.LBO))
    else
        println("LB: ", b.hist.L[b.hist.k], ", ", "UB: ", b.hist.U[b.hist.k], ", ", "total SP evaluations: ", b.hist.k*length(b.data.ms.I))
    end
    println(" ");
    println(" */"*"-"^ 78*"/*");
end

function emptyTrace() GenericTrace("",Dict{Symbol,Any}()) end

function report(b,n)

    Count=sum(N.nSP[:])
    LB=[];LMP_target=[];UB=[];iteration=[];LBOα=[];UBOα=[];
    LB_marker=Vector{Union{Float64,Missing}}(missing, Count);UB_marker=Vector{Union{Float64,Missing}}(missing, Count);LMP_target_marker=Vector{Union{Float64,Missing}}(missing, Count);LBOα_marker=Vector{Union{Float64,Missing}}(missing, Count);UBOα_marker=Vector{Union{Float64,Missing}}(missing, Count);
    LB_marker[1]=N.LB[1];UB_marker[1]=N.UB[1];LMP_target_marker[1]=N.LMP_target[1];LBOα_marker[1]=N.LBOα[1];UBOα_marker[1]=N.UBOα[1]
    marker_index=Vector{Union{Int,Missing}}(missing, Count)
    marker_index[1]=1
    for i in 2:B.hist.k
        marker_index[i]=sum(N.nSP[i] for i in 1:(i-1))+1
            LB_marker[sum(N.nSP[i] for i in 1:(i-1))+1]=N.LB[i]
            UB_marker[sum(N.nSP[i] for i in 1:(i-1))+1]=N.UB[i]
            UBOα_marker[sum(N.nSP[i] for i in 1:(i-1))+1]=N.UBOα[i]
            LBOα_marker[sum(N.nSP[i] for i in 1:(i-1))+1]=N.LBOα[i]
            LMP_target_marker[sum(N.nSP[i] for i in 1:(i-1))+1]=N.LMP_target[i]
    end

    for i in 1:B.hist.k
        for j in 1:N.nSP[i]
            push!(iteration,i)
            push!(LB,N.LB[i])
            push!(LMP_target,N.LMP_target[i])
            push!(LBOα,N.LBOα[i])
            push!(UBOα,N.UBOα[i])
            push!(UB,N.UB[i])
        end
    end

    xrmp=DataFrame()
    for i in 1:B.hist.k
            append!(xrmp, DataFrame(hcat(repeat(i:i,length(B.data.ms.G)),B.data.ms.G,N.xrmp[i][1:length(B.data.ms.G),:]), :auto))
    end
    rename!(xrmp,vcat("iteration","technology",B.data.ms.I))

    xlmp=DataFrame()
    for i in 1:B.hist.k
            append!(xlmp, DataFrame(hcat(repeat(i:i,length(B.data.ms.G)),B.data.ms.G,N.xlmp[i][1:length(B.data.ms.G),:]), :auto))
    end
    rename!(xlmp,vcat("iteration","technology",B.data.ms.I))

    xrmp_agg=DataFrame()
    split_G=[split(i,"_") for i in B.data.ms.G]
    for i in 1:B.hist.k
            for j in String.(unique!([split.(B.data.ms.G,"_")[i][1] for i in 1:length(B.data.ms.G)]))
            append!(xrmp_agg, DataFrame(hcat(i,j,sum(N.xrmp[i][k,:] for k in 1:length(B.data.ms.G) if split_G[k][1]==j)'),:auto))
    end
    end
    rename!(xrmp_agg,vcat("iteration","technology",B.data.ms.I))

    xlmp_agg=DataFrame()
    for i in 1:B.hist.k
            for j in String.(unique!([split.(B.data.ms.G,"_")[i][1] for i in 1:length(B.data.ms.G)]))
            append!(xlmp_agg, DataFrame(hcat(i,j,sum(N.xlmp[i][k,:] for k in 1:length(B.data.ms.G) if split_G[k][1]==j)'),:auto))
    end
    end
    rename!(xlmp_agg,vcat("iteration","technology",B.data.ms.I))


    DB_marker=DataFrame(LB_marker=LB_marker,LMP_target_marker=LMP_target_marker,LBOα_marker=LBOα_marker,UBOα_marker=UBOα_marker,UB_marker=UB_marker)
    DF_marker=DB_marker[2:end,:];
    namesDF=["LB(marker)","LMP_target(marker)","LBO_alpha(marker)","UBO_alpha(marker)","UB(marker)"]
    indicesDF=[1,2,3,4,5];
    coloursDF=["blue","purple","gold","green","red"];
    markersDF=["square","diamond",5,5,"square"];
    Traces_marker=[emptyTrace() for i in eachindex(indicesDF)];

    for i in eachindex(indicesDF)
        Traces_marker[i].fields[:type]="scatter";
        Traces_marker[i].fields[:x]=collect(2:Count).+0.1*i;
        Traces_marker[i].fields[:y]=DF_marker[:,indicesDF[i]];
        Traces_marker[i].fields[:name]=namesDF[i];
        Traces_marker[i].fields[:hoverinfo]="x+y";
        Traces_marker[i].fields[:textposition]="top center";
        Traces_marker[i].fields[:textfont]=Dict(:family =>  "Raleway, sans-serif");
        Traces_marker[i].fields[:visible]=true;
        Traces_marker[i].fields[:mode]="markers"
        Traces_marker[i].fields[:marker]=Dict(:size=>7,:color=>coloursDF[i],:symbol=>markersDF[i]);
        Traces_marker[i].fields[:showlegend]=true;
        Traces_marker[i].fields[:yaxis]="y";
    end
    Traces_marker[end].fields[:x]=collect(2:Count).+0.7;

    # idx=collect(1:Count)
    # plot(idx[2:end], LB[2:end],ylim=(0, 1.8*exp10(6)),label="LB: lower bound for the whole problem",markershape = :utriangle, markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0,color="SkyBlue",left_margin = 5Plots.mm, right_margin = 15Plots.mm)
    # plot!(idx[2:end].+0.14, LMP_target[2:Count],xlabel="subproblem sampling", ylabel="value",markershape =:hexagon,markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0,label="LMP target", color="DarkViolet")
    # # plot!(twinx(),idx[2:end], yticks=2:10:idx[end], ylabel="iteration", [SP.iteration[2:Count]],seriestype = :scatter,markersize=2,legend=false,left_margin = 5Plots.mm, right_margin = 15Plots.mm)
    # # plot!(idx[2:end].+0.28, SP.real_value[2:Count], markershape =:circle, label="Real function value", markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0,color="Gold")
    # plot!(idx[2:end].+0.7, N.UBO[2:Count], markershape =:circle, label="UBO: oralce upper bound",color="Crimson",markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0)
    # plot!(idx[2:end].+0.84, N.LBO[2:Count], markershape =:circle, label="LBO: oracle lower bound",color="SkyBlue",markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0)
    # p1=plot!(idx[2:end].+0.9, UB[2:Count], markershape =:utriangle, label="UB: upper bound for the whole problem",color="Crimson",markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0)
    # savefig(p1, pwd()*"/Hongyu/results/log.pdf")

    DB=DataFrame(iteration=iteration,SP_sampling=collect(1:Count),LB=LB,LMP_target=LMP_target,LBOα=LBOα,UBOα=UBOα,LBO=N.LBO,UBO=N.UBO,UB=UB)
    DF=DB[2:end,:];
    namesDF=["iteration","LB","LMP_target","LBO_alpha","UBO_alpha","LBO","UBO","UB"]
    indicesDF=[1,3,4,5,6,7,8,9];
    coloursDF=["gray","blue","purple","gold","green","blue","red","red"];
    markersDF=["circle","square","diamond","circle","circle","circle","circle","square"];
    Traces=[emptyTrace() for i in eachindex(indicesDF)];

    for i in eachindex(indicesDF)
        Traces[i].fields[:type]="scatter";
        Traces[i].fields[:x]=collect(DF[:,:SP_sampling]).+0.1*(i-1);
        Traces[i].fields[:y]=DF[:,indicesDF[i]];
        Traces[i].fields[:name]=namesDF[i];
        Traces[i].fields[:hoverinfo]="x+y";
        Traces[i].fields[:textposition]="top center";
        Traces[i].fields[:textfont]=Dict(:family =>  "Raleway, sans-serif");
        Traces[i].fields[:visible]=true;
        Traces[i].fields[:mode]="lines"
        Traces[i].fields[:marker]=Dict(:size=>7,:color=>coloursDF[i],:symbol=>markersDF[i]);
        Traces[i].fields[:showlegend]=true;
        Traces[i].fields[:yaxis]="y";
    end
    Traces[6].fields[:mode]="markers+lines";
    Traces[7].fields[:mode]="markers+lines";

    #Traces[end-1].fields[:line]=Dict(:dash=>"dash")
    Traces[1].fields[:mode]="markers";
    Traces[1].fields[:marker][:size]=7;
    Traces[1].fields[:yaxis]="y2";
    Traces[1].fields[:showlegend]=false;

    LayoutGraph=  Layout(Dict{Symbol,Any}(
                    :plot_bgcolor=>"white",
                    #:title => "Objective",
                    :xaxis => Dict(
                            # :domain => [0.35,0.5],
                            :showgrid => true,
                            # :gridcolor => "black",
                            :linecolor=>"black",
                            # :nticks   => 40,
                            # :range    => [0.130,0.170],
                            # :range    => [0,40],
                            :title    => Dict(:text => "Evaluation")
                            ),
                    :yaxis => Dict(
                            #:domain => [0.35,0.95],
                            :showgrid => false,
                            :linecolor=>"black",
                            # :gridcolor => "black",
                            # :range    => [0,4e6],
                            :title    => "Cost",
                            ),
                    :yaxis2 => Dict(
                            # :domain => [0.1,0.3],
                            :showgrid => true,
                            :linecolor=>"black",
                            # :gridcolor => "black",
                            # :range    => [0,100],
                            :title    => "Iteration",
                            :overlaying => "y",
                            :side => "right"
                            ),
                    :legend => Dict(
                            :x=>0.75,
                            :y=>0.95)
                    ) );

    Graph1=PlotlyJS.Plot(vcat(Traces,Traces_marker), LayoutGraph)
    open(pwd()*"/$modelpath/results/0102log.html","w") do f
        PlotlyJS.PlotlyBase.to_html(f, Graph1; include_plotlyjs="cdn", full_html=true)
    end
    savefig(Graph1,pwd()*"/$modelpath/results/0102log.pdf")

    ##time
    DFt=DataFrame(iteration=collect(1:b.hist.k),total=n.ttotal,tRMP=n.tRMP,tLMP=n.tLMP,tSP=n.tSP)

    XLSX.writetable(pwd()*"/$modelpath/results/0305log.xlsx",overwrite=true,log=(collect(DataFrames.eachcol(DB)), DataFrames.names(DB)),
    time=(collect(DataFrames.eachcol(DFt)), DataFrames.names(DFt)),xrmp=(collect(DataFrames.eachcol(xrmp)), DataFrames.names(xrmp)),
    xlmp=(collect(DataFrames.eachcol(xlmp)), DataFrames.names(xlmp)),xrmp_agg=(collect(DataFrames.eachcol(xrmp_agg)), DataFrames.names(xrmp_agg)),
    xlmp_agg=(collect(DataFrames.eachcol(xlmp_agg)), DataFrames.names(xlmp_agg)))

end

function report_realval(b,n)

    Count=sum(N.nSP[:])
    LB=[];LMP_target=[];UB=[];iteration=[];LBOα=[];UBOα=[];Rval=[];
    LB_marker=Vector{Union{Float64,Missing}}(missing, Count);UB_marker=Vector{Union{Float64,Missing}}(missing, Count);LMP_target_marker=Vector{Union{Float64,Missing}}(missing, Count);LBOα_marker=Vector{Union{Float64,Missing}}(missing, Count);UBOα_marker=Vector{Union{Float64,Missing}}(missing, Count);Rval_marker=Vector{Union{Float64,Missing}}(missing, Count);
    LB_marker[1]=N.LB[1];UB_marker[1]=N.UB[1];LMP_target_marker[1]=N.LMP_target[1];LBOα_marker[1]=N.LBOα[1];UBOα_marker[1]=N.UBOα[1];Rval_marker[1]=N.Rval[1]
    marker_index=Vector{Union{Int,Missing}}(missing, Count)
    marker_index[1]=1
    for i in 2:B.hist.k
        marker_index[i]=sum(N.nSP[i] for i in 1:(i-1))+1
            LB_marker[sum(N.nSP[i] for i in 1:(i-1))+1]=N.LB[i]
            UB_marker[sum(N.nSP[i] for i in 1:(i-1))+1]=N.UB[i]
            UBOα_marker[sum(N.nSP[i] for i in 1:(i-1))+1]=N.UBOα[i]
            LBOα_marker[sum(N.nSP[i] for i in 1:(i-1))+1]=N.LBOα[i]
            Rval_marker[sum(N.nSP[i] for i in 1:(i-1))+1]=N.Rval[i]
            LMP_target_marker[sum(N.nSP[i] for i in 1:(i-1))+1]=N.LMP_target[i]
    end

    for i in 1:B.hist.k
        for j in 1:N.nSP[i]
            push!(iteration,i)
            push!(LB,N.LB[i])
            push!(LMP_target,N.LMP_target[i])
            push!(LBOα,N.LBOα[i])
            push!(UBOα,N.UBOα[i])
            push!(UB,N.UB[i])
            push!(Rval,N.Rval[i])
        end
    end


    DB_marker=DataFrame(LB_marker=LB_marker,LMP_target_marker=LMP_target_marker,LBOα_marker=LBOα_marker,UBOα_marker=UBOα_marker,UB_marker=UB_marker,Rval=Rval)
    DF_marker=DB_marker[2:end,:];
    namesDF=["LB(marker)","LMP_target(marker)","LBO_alpha(marker)","UBO_alpha(marker)","UB(marker)","Rval(marker)"]
    indicesDF=[1,2,3,4,5,6];
    coloursDF=["blue","purple","gold","green","red","yellow"];
    markersDF=["square","diamond",5,5,"square","diamond"];
    Traces_marker=[emptyTrace() for i in eachindex(indicesDF)];

    for i in eachindex(indicesDF)
        Traces_marker[i].fields[:type]="scatter";
        Traces_marker[i].fields[:x]=collect(2:Count).+0.1*i;
        Traces_marker[i].fields[:y]=DF_marker[:,indicesDF[i]];
        Traces_marker[i].fields[:name]=namesDF[i];
        Traces_marker[i].fields[:hoverinfo]="x+y";
        Traces_marker[i].fields[:textposition]="top center";
        Traces_marker[i].fields[:textfont]=Dict(:family =>  "Raleway, sans-serif");
        Traces_marker[i].fields[:visible]=true;
        Traces_marker[i].fields[:mode]="markers"
        Traces_marker[i].fields[:marker]=Dict(:size=>7,:color=>coloursDF[i],:symbol=>markersDF[i]);
        Traces_marker[i].fields[:showlegend]=true;
        Traces_marker[i].fields[:yaxis]="y";
    end
    Traces_marker[end].fields[:x]=collect(2:Count).+0.7;

    # idx=collect(1:Count)
    # plot(idx[2:end], LB[2:end],ylim=(0, 1.8*exp10(6)),label="LB: lower bound for the whole problem",markershape = :utriangle, markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0,color="SkyBlue",left_margin = 5Plots.mm, right_margin = 15Plots.mm)
    # plot!(idx[2:end].+0.14, LMP_target[2:Count],xlabel="subproblem sampling", ylabel="value",markershape =:hexagon,markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0,label="LMP target", color="DarkViolet")
    # # plot!(twinx(),idx[2:end], yticks=2:10:idx[end], ylabel="iteration", [SP.iteration[2:Count]],seriestype = :scatter,markersize=2,legend=false,left_margin = 5Plots.mm, right_margin = 15Plots.mm)
    # # plot!(idx[2:end].+0.28, SP.real_value[2:Count], markershape =:circle, label="Real function value", markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0,color="Gold")
    # plot!(idx[2:end].+0.7, N.UBO[2:Count], markershape =:circle, label="UBO: oralce upper bound",color="Crimson",markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0)
    # plot!(idx[2:end].+0.84, N.LBO[2:Count], markershape =:circle, label="LBO: oracle lower bound",color="SkyBlue",markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0)
    # p1=plot!(idx[2:end].+0.9, UB[2:Count], markershape =:utriangle, label="UB: upper bound for the whole problem",color="Crimson",markersize=2,dpi=400, grid=false, linewidth=0.5,framestyle=:box,foreground_color_legend = nothing,background_color_legend=nothing,markerstrokewidth = 0)
    # savefig(p1, pwd()*"/Hongyu/results/log.pdf")

    DB=DataFrame(iteration=iteration,SP_sampling=collect(1:Count),LB=LB,LMP_target=LMP_target,LBOα=LBOα,UBOα=UBOα,LBO=N.LBO,UBO=N.UBO,UB=UB,Rval=Rval)
    DF=DB[2:end,:];
    namesDF=["iteration","LB","LMP_target","LBO_alpha","UBO_alpha","LBO","UBO","UB","Rval"]
    indicesDF=[1,3,4,5,6,7,8,9,10];
    coloursDF=["gray","blue","purple","gold","green","blue","red","red","yellow"];
    markersDF=["circle","square","diamond","circle","circle","circle","circle","square","diamond"];
    Traces=[emptyTrace() for i in eachindex(indicesDF)];

    for i in eachindex(indicesDF)
        Traces[i].fields[:type]="scatter";
        Traces[i].fields[:x]=collect(DF[:,:SP_sampling]).+0.1*(i-1);
        Traces[i].fields[:y]=DF[:,indicesDF[i]];
        Traces[i].fields[:name]=namesDF[i];
        Traces[i].fields[:hoverinfo]="x+y";
        Traces[i].fields[:textposition]="top center";
        Traces[i].fields[:textfont]=Dict(:family =>  "Raleway, sans-serif");
        Traces[i].fields[:visible]=true;
        Traces[i].fields[:mode]="lines"
        Traces[i].fields[:marker]=Dict(:size=>7,:color=>coloursDF[i],:symbol=>markersDF[i]);
        Traces[i].fields[:showlegend]=true;
        Traces[i].fields[:yaxis]="y";
    end
    Traces[6].fields[:mode]="markers+lines";
    Traces[7].fields[:mode]="markers+lines";

    #Traces[end-1].fields[:line]=Dict(:dash=>"dash")
    Traces[1].fields[:mode]="markers";
    Traces[1].fields[:marker][:size]=7;
    Traces[1].fields[:yaxis]="y2";
    Traces[1].fields[:showlegend]=false;

    LayoutGraph=  Layout(Dict{Symbol,Any}(
                    :plot_bgcolor=>"white",
                    #:title => "Objective",
                    :xaxis => Dict(
                            # :domain => [0.35,0.5],
                            :showgrid => true,
                            # :gridcolor => "black",
                            :linecolor=>"black",
                            # :nticks   => 40,
                            # :range    => [0.130,0.170],
                            # :range    => [0,40],
                            :title    => Dict(:text => "Evaluation")
                            ),
                    :yaxis => Dict(
                            #:domain => [0.35,0.95],
                            :showgrid => false,
                            :linecolor=>"black",
                            # :gridcolor => "black",
                            # :range    => [0,4e6],
                            :title    => "Cost",
                            ),
                    :yaxis2 => Dict(
                            # :domain => [0.1,0.3],
                            :showgrid => true,
                            :linecolor=>"black",
                            # :gridcolor => "black",
                            # :range    => [0,100],
                            :title    => "Iteration",
                            :overlaying => "y",
                            :side => "right"
                            ),
                    :legend => Dict(
                            :x=>0.75,
                            :y=>0.95)
                    ) );

    Graph1=PlotlyJS.Plot(vcat(Traces,Traces_marker), LayoutGraph)
        open(pwd()*"/$modelpath/results/0105log.html","w") do f
            PlotlyJS.PlotlyBase.to_html(f, Graph1; include_plotlyjs="cdn", full_html=true)
        end
        savefig(Graph1,pwd()*"/$modelpath/results/0105log.pdf")

    ##time
    DFt=DataFrame(iteration=collect(1:b.hist.k),total=n.ttotal,tRMP=n.tRMP,tLMP=n.tLMP,tSP=n.tSP)

    XLSX.writetable(pwd()*"/$modelpath/results/0102log.xlsx",overwrite=true,log=(collect(DataFrames.eachcol(DB)), DataFrames.names(DB)),
    time=(collect(DataFrames.eachcol(DFt)), DataFrames.names(DFt)))

end

# report(B,N)
