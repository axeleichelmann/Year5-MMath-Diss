cd(dirname(@__FILE__)); 
flush(stdout);
println("*/"*repeat("-",40)*"/*");

modelpath="OREIA3_2017";                         # Path to the optimisation models
path="$modelpath/data";                             # Specifies the relative path where data is held

include("functions/load_packages.jl");              # Loads the required packages: CSV, DataFrames, Gurobi, JuMP, etc.
include("functions/load_support.jl");               # Loads the data support functions
include("$modelpath/load_models.jl");               # Contains the models (SP and RMP)


#### Data Loading & Restoring
mainFile="main.xlsx";                            # Name of the index file
jldFile="data.jld2";                                # Name of jld2 file to be saved or read
updateData=true;                                   # Flag to ask it to update the data, and generate a new jld2 file
#####

fileExists=jldFile in cd(readdir,path);
if !fileExists || updateData
    flush(stdout); print("Generating data tables...    "); a=time();
    data=generateTablesfromXLSXNew2(path,mainFile);   # Data loaded from XLSX files into dictionaries (new version, experimental)
    println("done in $(round(time()-a,digits=1))s ")
else
    flush(stdout); print("Restoring data tables...     "); a=time();
    data=load("$path/$jldFile","data");
    println("done in $(round(time()-a,digits=1))s ")
end

############# ALGORITHM GLOBAL PARAMETERS DEFINITION
algDict=loadParameters(data["Algorithm"],"parameters"); # Loads the tables for the general setup of the solution algorithm
alg=transformStructure(algDict,"alg");                  # Creates a structure for the algorithm parameters

############### DATA LOADING

include("$path/data_preparation.jl");             # Custom combination of functions, according to requirements of model (IMPORTANT)
include("functions/link_models.jl");              # Create optimisation models and link them for decomposition

############### PLANNER PART
if alg.al!=0
    include("functions/load_structures.jl");          # Add stochastic planner's structures
    include("functions/load_functions.jl");           # Add stochastic planner's functions
    include("functions/generate_structures.jl");      # Generates objects B and S, neccesary for the planner to work [with data loaded here] and set the special point
#######

############### OUTPUT PRINTING FUNCTIONS
    include("functions/load_printfunctions.jl");      # Has the functions to print out the resultss

    B, S, N = solve_Benders!(B,S,N);                            # Starts solving the problem

end;



#####  -------  CREATE EXCEL FILE CONTAINING RESULTS OF OPERATIONAL MODEL & PRODUCE GRAPHS FOR RESULTS  --------  #####

if alg.al == 0   ##### FOR UNDECOMPOSED MODEL (FOR 1 INVESTMENT NODE)

    ###### CHECK THAT THE SHEDS ARE ALL 0 AT INVESTMENT NODE 1
    sum(sum(value.(mU[:qHShedP]["N1","$(hs)",:]).data) for hs in union(ps.HnS,ps.HOUSES))  # Total House & Heat Store Positive Shed at Investment Node 1
    sum(sum(value.(mU[:qHShedN]["N1","$(hs)",:]).data) for hs in union(ps.HnS,ps.HOUSES))  # Total House & Heat Store Negative Shed at Investment Node 1
    sum(sum(value.(mU[:gShed]["N1","$r",:]).data) for r in ps.R)  # Total Renewables Shed at Investment Node 1
    sum(sum(value.(mU[:dShed]["N1","$d",:]).data) for d in ps.D)  # Total Demand Shed at Investment Node 1

    ##### CREATE EXCEL FILE CONTAINING RESULTS FOR UNDECOMPOSED MODEL
    include("OREIA3_2017/rp_invest_undecomp.jl")

    ######## CREATE PLOTS FOR UNDECOMPOSED MODEL
    include("OREIA3_2017/result_plots_undecomp.jl")

else ##### ---- FOR DECOMPOSED MODEL ---- ####

    ###### CHECK THAT THE SHEDS ARE ALL 0 
    sum(sum(value.(S.ex.m[:qHShedP]["$(hs)",:]).data) for hs in union(ps.HnS,ps.HOUSES))  # Total House & Heat Store Positive Shed
    sum(sum(value.(S.ex.m[:qHShedN]["$(hs)",:]).data) for hs in union(ps.HnS,ps.HOUSES))  # Total House & Heat Store Negative Shed
    sum(sum(value.(S.ex.m[:gShed]["$r",:]).data) for r in ps.R)  # Total Renewables Shed
    sum(sum(value.(S.ex.m[:dShed]["$d",:]).data) for d in ps.D)  # Total Demand Shed

    ##### CREATE EXCEL FILE CONTAINING RESULTS FOR DECOMPOSED MODEL
    include("OREIA3_2017/rp_invest.jl")

    ######## CREATE PLOTS FOR DECOMPOSED MODEL
    include("OREIA3_2017/result_plots.jl")

end

