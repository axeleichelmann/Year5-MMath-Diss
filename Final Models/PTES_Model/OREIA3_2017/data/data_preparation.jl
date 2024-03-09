################# DATA PREPARATION #########################
############################################################

if !fileExists || updateData

    #### 1. OPERATIONAL PART ##################################
    inheritData!(data["Operation"], "fuels" => "generators", "Name" => "Fuel")
    augmentSlices!(data["Operation"], "slices") # Add information about length and offset for each slice

    ### Time series composer
    genTimeSeries!(data["Operation"], "slices", "time_series_info", "time_series")

    ###### Add the timeseries to the demand, generation and areas
    inheritData!(data["Operation"], "time_series" => "generators", "Name" => "time_series", ["Series"])
    inheritData!(data["Operation"], "time_series" => "demands", "Name" => "time_series", ["Series"])
    inheritData!(data["Operation"], "time_series" => "buses", "Name" => "time_series", ["Series"])

    ##### Scale the "Series" by the "Scaling" factor
    scaleData!(data["Operation"], "demands", "Series", "Scaling")
    scaleData!(data["Operation"], "demands", "Series", "Scaling")

    #####Define Sets / Paramters from tables - Specific functions

    psDict = loadSets(data["Operation"], "sets")
    ppDict = loadParameters(data["Operation"], "parameters")

    #### 2. INVESTMENT PART ##################################
    augmentStructure(data["Investment"], "structure") # Creates stochastic path based on tree with probabilities
    inheritData!(data["Investment"], "stages" => "structure", "Stage" => "Level") # Inherits data about node length and type of node (operation or investment)

    addMapTable(data["Investment"], "structure", "Node") #(Optional) Once the information about stages has been inherited then sets of parents (for investment and operation) are calculated.
    addYrs(data["Investment"], "structure")

    msDict = loadSets(data["Investment"], "sets")
    mpDict = loadParameters(data["Investment"], "parameters")

    #### Store data in JLD2 file
    save("$path/$jldFile", "data", data)
    jldopen("$path/$jldFile", "a+") do f
        f["ppDict"] = ppDict
        f["psDict"] = psDict
        f["mpDict"] = mpDict
        f["msDict"] = msDict
    end

    println("done in $(round(time()-a,digits=1))s ")

else
    print("Restoring data preparation...")
    a = time()
    ppDict, psDict, mpDict, msDict = load.("$path/$jldFile", ["ppDict", "psDict", "mpDict", "msDict"])
    println("done in $(round(time()-a,digits=1))s ")

end


##### Convert dictonaries into tight structures

ps = transformStructure(psDict, "ps");
pp = transformStructure(ppDict, "pp");

ms = transformStructure(msDict, "ms");
mp = transformStructure(mpDict, "mp");


#### Other operations
##Quickly populate a set of nodes according to some table
for s in ms.S, n in ms.IL[s]
    for g in ms.M
        mp.imin[g][n] = data["Investment"]["dataHist"]["Cap_S$(s)"][findfirst(data["Investment"]["dataHist"]["Name"] .== g)]
    end
end
