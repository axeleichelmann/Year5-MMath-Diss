cd(dirname(@__FILE__)); flush(stdout); println("*/"*repeat("-",40)*"/*");

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
include("OREIA3_2017/rp_invest.jl")

sum(value.(S.ex.m[:pHShedP]["House1",:]).data)

B.temp.x # get final investment values --- can also use S.temp.x