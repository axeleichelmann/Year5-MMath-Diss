flush(stdout); print("Loading support functions... "); a=time();

function next(h::Int64,HSl::Union{Array{Int64,1},UnitRange{Int64}};periods=1)
    out=0;
    if !(h in HSl)
        throw(ArgumentError("period $h is not in the set"))
    else
        for p in 1:periods
            if h in HSl[end]
                out=HSl[1];
            else
                out=h+1
            end
            h=out;
        end
    end
    return out
end



function prev(h::Int64,HSl::Union{Array{Int64,1},UnitRange{Int64}};periods=1)
    out=0;
    if !(h in HSl)
        throw(ArgumentError("period $h is not in the set"))
    else
        for p in 1:periods
            if h in HSl[1]
                out=HSl[end];
            else
                out=h-1
            end
            h=out
        end
    end
    return out
end


function next(h::Int64,slst,slln;periods=1)
    # slst corresponds to the start of the slice and slln to the length, these are standard properties added when the augmentSlices is used.

    @assert keys(slst)==keys(slln) "The indexing elements of the slice start positions and lengths are not the same, please index by the same property"

    st=collect(values(slst));
    en=collect(values(slst))+collect(values(slln)).-1
    out=0;

    if h>maximum(en) || h<minimum(st)
        throw(ArgumentError("period $h is not in the set"))
    else

        for p in 1:periods
            if h in en
                out=st[findfirst(h.==en)];
            else
                out=h+1
            end
            h=out;
        end

    end
    return out

end

function prev(h::Int64,slst,slln;periods=1)
    # slst corresponds to the start of the slice and slln to the length, these are standard properties added when the augmentSlices is used.

    @assert keys(slst)==keys(slln) "The indexing elements of the slice start positions and lengths are not the same, please index by the same property"

    st=collect(values(slst));
    en=collect(values(slst))+collect(values(slln)).-1
    out=0;

    if h>maximum(en) || h<minimum(st)
        throw(ArgumentError("period $h is not in the set"))
    else

        for p in 1:periods
            if h in st
                out=en[findfirst(h.==st)];
            else
                out=h-1
            end
            h=out
        end
    end
    return out

end


function H2S(setH,slStart,slLength)

    outputDict=Dict();

    ks=keys(slLength);

    @assert ks==keys(slStart) "The indexing elements of the slice start positions and lengths are not the same, please index by the same property"

    for k in ks
        for h in slStart[k]:slStart[k]+slLength[k]-1
            outputDict[h]=k
        end
    end

    return outputDict


end



function iS(A)
    if !ismissing(A)
        !(typeof(A)<:Array)
            if typeof(A) <: InlineString
                A=string(A)
            end

    end
    return A
end



function generateTables(path::String,filename::String)
    readCSV(file)=CSV.read(file,DataFrame;comment="#",silencewarnings=true, missingstring=["","-"],ignoreemptyrows=true); #Shorthand function to read CSV files

    basepath=pwd();


    cd(path*"/");
    folders=readCSV(filename);
    indices=[readCSV(folders.Folder[i]*"/index.csv") for i in eachindex(folders.Folder)]; # Read indices on folders
    tables=[[readCSV(folders.Folder[i]*"/"*indices[i].filename[j]) for j in eachindex(indices[i].filename)] for i in eachindex(folders.Folder)];
    columns=[[names(tables[i][j]) for j in eachindex(indices[i].filename)] for i in eachindex(folders.Folder)];
    cd("..");

    #outputTables=Dict{Any,Any}(folders.Data[i]=>Dict{Any,Any}(indices[i].table_name[j]=>Dict{Any,Any}(columns[i][j][k]=>tables[i][j][:,columns[i][j][k]] for k in eachindex(columns[i][j])) for j in eachindex(indices[i].table_name)) for i in eachindex(folders.Folder));

    outputTables=Dict{Any,Any}(folders.Data[i]=>Dict{Any,Any}(indices[i].table_name[j]=>Dict{Any,Any}(columns[i][j][k]=>iS.(tables[i][j][:,columns[i][j][k]]) for k in eachindex(columns[i][j])) for j in eachindex(indices[i].table_name)) for i in eachindex(folders.Folder));

    cd(basepath)

    return outputTables

end

function dropComments(df)

    listComments=occursin.("#",string.(df[!,names(df)[1]]));
    validRows=setdiff(1:nrow(df),listComments)

    df=df[validRows,:]

    #df .= ifelse.(df .== "", missing, df);
    #df .= ifelse.(df .== "-", missing, df);

    return df

end

function generateTablesfromXLSXNew(path::String,filename::String)
    basepath=pwd();
    cd(path*"/");

    folders=DataFrame(XLSX.readtable(filename,"main")...)
    folders=dropComments(folders)

    indices=[DataFrame(XLSX.readtable(folders.Filename[i],"index")...) for i in eachindex(folders.Filename)]
    indices=dropComments.(indices)

    #@btime tables=[[readXLSX(folders.Filename[i],iS(indices[i].sheet_name[j])) for j in eachindex(indices[i].sheet_name)] for i in eachindex(folders.Filename)]; ### This is the slow one
    tables=[[dropComments(DataFrame(XLSX.readtable(folders.Filename[i],iS(indices[i].sheet_name[j]))...)) for j in eachindex(indices[i].sheet_name)] for i in eachindex(folders.Filename)]; ### This is the slow one
    #tables=dropComments.(tables)

    columns=[[names(tables[i][j]) for j in eachindex(indices[i].sheet_name)] for i in eachindex(folders.Filename)];

    cd(basepath)

    outputTables=Dict{Any,Any}(folders.Data[i]=>Dict{Any,Any}(indices[i].table_name[j]=>Dict{Any,Any}(columns[i][j][k]=>iS.(tables[i][j][:,columns[i][j][k]]) for k in eachindex(columns[i][j])) for j in eachindex(indices[i].table_name)) for i in eachindex(folders.Filename));


    return outputTables

end

function generateTablesfromXLSX(path::String,filename::String)
    basepath=pwd();
    cd(path*"/");
    function readXLSX(file,sheet)
        CSV.write("temp.csv",DataFrame(XLSX.readtable(file,sheet)...));
        CSV.read("temp.csv",DataFrame;comment="#",silencewarnings=true, missingstring=["","-"],ignoreemptyrows=true); #Shorthand function to read CSV files
    end
    folders=readXLSX(filename,"main");
    indices=[readXLSX(folders.Filename[i],"index") for i in eachindex(folders.Filename)]
    tables=[[readXLSX(folders.Filename[i],iS(indices[i].sheet_name[j])) for j in eachindex(indices[i].sheet_name)] for i in eachindex(folders.Filename)]; ### This is the slow one

    columns=[[names(tables[i][j]) for j in eachindex(indices[i].sheet_name)] for i in eachindex(folders.Filename)];
    rm("temp.csv");
    cd("..");

    cd(basepath)

    outputTables=Dict{Any,Any}(folders.Data[i]=>Dict{Any,Any}(indices[i].table_name[j]=>Dict{Any,Any}(columns[i][j][k]=>iS.(tables[i][j][:,columns[i][j][k]]) for k in eachindex(columns[i][j])) for j in eachindex(indices[i].table_name)) for i in eachindex(folders.Filename));


    return outputTables
end


function indexElement(data::Dict,path::NTuple{3,String})
    return Dict{Any,Any}(data[path[1]][path[2]][path[3]][i]=>Dict{Any,Any}(j=>data[path[1]][path[2]][j][i] for j in eachindex(data[path[1]][path[2]])) for i in eachindex(data[path[1]][path[2]][path[3]]))
end

function indexElement(data::Dict,path::Tuple{String,String,Array{String,1}})
    return Dict{Any,Any}((data[path[1]][path[2]][path[3][1]][i],data[path[1]][path[2]][path[3][2]][i])=>Dict{Any,Any}(j=>data[path[1]][path[2]][j][i] for j in eachindex(data[path[1]][path[2]])) for i in eachindex(data[path[1]][path[2]][path[3][1]]))
end

function indexElement(datawpath::Dict,path1::String,path2::String)
    return Dict{Any,Any}(datawpath[path1][path2][i]=>Dict{Any,Any}(j=>datawpath[path1][j][i] for j in eachindex(datawpath[path1])) for i in eachindex(datawpath[path1][path2]))
end


function inheritData!(obj1::Dict,obj1_field::String,obj2::Dict)
    for i in eachindex(obj1)
        if obj1_field in keys(obj1[i])
            for j in eachindex(obj2[obj1[i][obj1_field]])
                if !(j in keys(obj1[i]))
                    if !ismissing(obj2[obj1[i][obj1_field]][j])
                        obj1[i][j]=obj2[obj1[i][obj1_field]][j]
                    end
                end
            end
        end
    end
    return obj1
end


function inheritData!(obj1::Dict,obj1_field::String,obj2::Dict,obj2_field::String)

    mapProp=Dict(obj2[j][obj2_field]=>j for j in keys(obj2))


    for i in eachindex(obj1)
        if obj1_field in keys(obj1[i])
            for j in eachindex(obj2[mapProp[obj1[i][obj1_field]]])
                if !(j in keys(obj1[i]))
                    if !ismissing(obj2[mapProp[obj1[i][obj1_field]]])
                        obj1[i][j]=obj2[mapProp[obj1[i][obj1_field]]][j]
                    end
                end
            end
        end
    end
    return obj1
end




function inheritData!(data::Dict,obj::Pair,matchn::Pair;overwrite=false)

    matchPos=Any[];

    for j in eachindex(data[obj.second][matchn.second])
        if ismissing(data[obj.second][matchn.second][j])
            push!(matchPos,missing)
        else
            push!(matchPos,findfirst(data[obj.second][matchn.second][j].==data[obj.first][matchn.first]))
        end
    end

    for i in keys(data[obj.first])
        if !overwrite
            if !(i in keys(data[obj.second]))
                if !(i in keys(data[obj.second])) data[obj.second][i]=Any[]; end
                for j in eachindex(matchPos)
                    if ismissing(matchPos[j])
                        push!(data[obj.second][i],missing);
                    else
                        push!(data[obj.second][i],data[obj.first][i][matchPos[j]])
                    end
                end
            end
        else

                data[obj.second][i]=Any[];
                for j in eachindex(matchPos)
                    if ismissing(matchPos[j])
                        push!(data[obj.second][i],missing);
                    else
                        push!(data[obj.second][i],data[obj.first][i][matchPos[j]])
                    end
                end



        end
    end
end

####This
function inheritData!(data::Dict,obj::Pair,matchn::Pair,fields::Array{String,1};outputChange=[""],overwrite=false)


    matchPos=Any[];
    if outputChange==[""] outputChange=fields end

    for j in eachindex(data[obj.second][matchn.second])
        if ismissing(data[obj.second][matchn.second][j])
            push!(matchPos,missing)
        else
            push!(matchPos,findfirst(data[obj.second][matchn.second][j].==data[obj.first][matchn.first]))
        end
    end

    for i in keys(data[obj.first])
        if !overwrite
            #if !(i in keys(data[obj.second])) && i in fields
            if (!(i in keys(data[obj.second])) && i in fields) || (i in fields && !(outputChange[findfirst(i.==fields)] in keys(data[obj.second])))
                 #data[obj.second][i]=Any[];
                 data[obj.second][outputChange[findfirst(i.==fields)]]=Any[];
                for j in eachindex(matchPos)
                    #println("i=$i,j=$j,mp=$(matchPos[j])")
                    if ismissing(matchPos[j])
                        push!(data[obj.second][outputChange[findfirst(i.==fields)]],missing);
                        #println("isMissing");
                    else
                        push!(data[obj.second][outputChange[findfirst(i.==fields)]],data[obj.first][i][matchPos[j]])
                        #println(data[obj.first][i][matchPos[j]])
                    end
                end
            end
        else
            if i in fields
                data[obj.second][outputChange[findfirst(i.==fields)]]=Any[];
                for j in eachindex(matchPos)
                    if ismissing(matchPos[j])
                        push!(data[obj.second][outputChange[findfirst(i.==fields)]],missing);
                    else
                        push!(data[obj.second][outputChange[findfirst(i.==fields)]],data[obj.first][i][matchPos[j]])
                    end
                end
            end
        end
    end
end





function inheritSelectData!(obj1::Dict,obj1_field::String,obj2::Dict,field::String)
    for i in eachindex(obj1)
        if obj1_field in keys(obj1[i])
            obj1[i][field]=obj2[obj1[i][obj1_field]][field]
        end
    end
    return obj1
end


function dropmiss!(obj)
    for i in eachindex(obj)
        for j in eachindex(obj[i])
            ismissing(obj[i][j]) ? pop!(obj[i],j) : nothing
        end
    end
end


function matchData!(object,path,format,data,format_start,format_end)
    for k in keys(object)
        startRow=findfirst((data[path][object[k]["Table"]][format[1]].==object[k][format_start[1]]).*(data[path][object[k]["Table"]][format[2]].==object[k][format_start[2]]))
        endRow=findfirst((data[path][object[k]["Table"]][format[1]].==object[k][format_end[1]]).*(data[path][object[k]["Table"]][format[2]].==object[k][format_end[2]]))
        object[k]["Series"]=data[path][object[k]["Table"]][object[k]["Column"]][startRow:endRow]
    end
end


function collapseDataBackup(object::Dict{Any,Any},indexer,sublevel)
    karray=collect(keys(object));
    uIndex=unique([karray[i][indexer] for i in eachindex(karray)]);
    sIndex=unique([karray[i][sublevel] for i in eachindex(karray)]);

    object_out=Dict{Any,Any}(i=>Dict() for i in uIndex);
    nms=["Circular","Weight","slStart","slEnd","Series","slices"]

    for i in uIndex
        slCirc=[object[(i,j)]["Circular"] for j in sIndex];
        slWeight=[object[(i,j)]["Weight"] for j in sIndex];
        slLength=[length(object[(i,j)]["Series"]) for j in sIndex];
        slSeries=Float64[];
        for j in sIndex
            append!(slSeries,object[(i,j)]["Series"].*object[(i,j)]["Scaling"])
        end

        slStart=vcat(1,1 .+cumsum(slLength[j] for j in 1:length(sIndex)-1))
        slEnd=cumsum(slLength[j] for j in eachindex(sIndex));

        slices=[sIndex[j] for j in eachindex(sIndex)];

        object_out[i]=Dict(nms[1]=>slCirc,nms[2]=>slWeight,nms[3]=>slStart,nms[4]=>slEnd,nms[5]=>slSeries,nms[6]=>slices);


        #vls=[slCirc,slWeight,slStart,slEnd,slSeries];
        #object_out[i]=Dict(nms[k]=>vls[k] for k in eachindex(nms));
    end


    return object_out

end


function collapseData(object::Dict{Any,Any},indexer,sublevel)
    karray=collect(keys(object));
    uIndex=sort(unique([karray[i][indexer] for i in eachindex(karray)]));
    sIndex=sort(unique([karray[i][sublevel] for i in eachindex(karray)]));

    object_out=Dict{Any,Any}(i=>Dict() for i in uIndex);
    nms=["Series","slices"]

    for i in uIndex
        slSeries=Float64[];
        slices=String[];
        for j in sIndex
            append!(slSeries,object[(i,j)]["Series"].*object[(i,j)]["Scaling"])

        end

        slices=[sIndex[j] for j in eachindex(sIndex)];
        object_out[i]=Dict(nms[1]=>slSeries,nms[2]=>slices);


        #vls=[slCirc,slWeight,slStart,slEnd,slSeries];
        #object_out[i]=Dict(nms[k]=>vls[k] for k in eachindex(nms));
    end


    return object_out

end


function generateTimeSeriesBackup(tableData,dataPath,slices,ts_info,format)
    ts=indexElement(tableData,(dataPath,ts_info,["Name","Slice"]));
    inheritData!(ts,"Slice",slices)
    matchData!(ts,dataPath,("Date","Time"),tableData,(format[1][1],format[1][2]),(format[2][1],format[2][2]));
    ts=collapseData(ts,1,2);

    return ts
    ##
    #time_series=indexElement(tablesD,("Operation","time_series_info",["Name","Slice"]));
    #inheritData!(time_series,"Slice",slices);
    #matchData!(time_series,"Operation",("Date","Time"),tablesD,("Start_date","Start_time"),("End_date","End_time"));
    #time_series=collapseData(time_series,1,2);
    ##
end

function subindexElement(data,filt::Pair{String,String})
    d=filter(g->g.second[filt.first]==filt.second,data)
    return d
end




function scaleSelectData!(object,property,scalingval);
    for k in keys(object)
        object[k][property].*=object[k][scalingval];
    end
    return object
end

function augmentSlicesOld!(slices);
    offset::Int64=0;
    element::Int64=0;
    orderedSlices=sort(collect(keys(slices)));
    for i in orderedSlices
    #for i in keys(slices)
        dstart=DateTime(Date(slices[i]["Start_date"],"d/m/y"),slices[i]["Start_time"]);
        dend=DateTime(Date(slices[i]["End_date"],"d/m/y"),slices[i]["End_time"]);
        stepMin=Minute(slices[i]["Time_step"]*60); ### time step in minutes

        ln=length(dstart:stepMin:dend);

        offset=ln*element+1;
        element+=1;

        slices[i]["Length"]=ln;
        slices[i]["Offset"]=offset;
        slices[i]["Order"]=element;
        slices[i]["Time_step"]=stepMin;

    end
end




function augmentSlices!(datawpath,table_sl);
    offset::Int64=0;
    element::Int64=0;
    ln::Int64=0;

    offsetV::Array{Int64,1}=Int64[];
    elementV::Array{Int64,1}=Int64[];
    lnV::Array{Int64,1}=Int64[];

    #orderedSlices=sort(collect(keys(table_sl[uniqueField])));


    if typeof(datawpath[table_sl]["Start_date"][1])==Date
        datawpath[table_sl]["Start_date"]=[Dates.format(datawpath[table_sl]["Start_date"][i], "dd/mm/yyyy") for i in eachindex(datawpath[table_sl]["Start_date"])]
        datawpath[table_sl]["End_date"]=[Dates.format(datawpath[table_sl]["End_date"][i], "dd/mm/yyyy") for i in eachindex(datawpath[table_sl]["End_date"])]
    end


    for i in eachindex(datawpath[table_sl]["Start_date"])
    #for i in orderedSlices
    #for i in keys(slices)

            dstart=DateTime(Date(datawpath[table_sl]["Start_date"][i],"d/m/y"),datawpath[table_sl]["Start_time"][i]);
            dend=DateTime(Date(datawpath[table_sl]["End_date"][i],"d/m/y"),datawpath[table_sl]["End_time"][i]);
            stepMin=Minute(datawpath[table_sl]["Time_step"][i]*60); ### time step in minutes
            ln=length(dstart:stepMin:dend);

            #offset=ln*element+1;
            if length(offsetV)==0
                offset=1
            else
                offset=sum(lnV)+1
            end
            #offset=ln*element+1;
            element+=1;
            push!(offsetV,offset);
            push!(elementV,element);
            push!(lnV,ln);

    end




    datawpath[table_sl]["Length"]=lnV;
    datawpath[table_sl]["Offset"]=offsetV;
    datawpath[table_sl]["Order"]=elementV;
    datawpath[table_sl]["Offset2"]=offsetV.+lnV.-1;
    #tables_sl["Time_step"]=stepMin;

    off,pos=findmax(datawpath[table_sl]["Offset"]);
    tsLength=off+datawpath[table_sl]["Length"][pos]-1;
    datawpath[table_sl]["Set"]=1:tsLength

    datawpath[table_sl]["Subset"]=[offsetV[i]:offsetV[i].+lnV[i].-1 for i in eachindex(offsetV)]


end


function generateTimeSeries(tableData,dataPath,slices,ts_info,format)
    ts=indexElement(tableData,(dataPath,ts_info,["Name","Slice"]));
    #### Find the start of each slice in the corresponding table
    rearrangeSet=sort(collect(keys(ts)))

    #for t in keys(ts)
    for t in rearrangeSet
        stpos=findfirst((tableData[dataPath][ts[t]["Table"]]["Date"].==slices[t[2]][format[1]]) .* (tableData[dataPath][ts[t]["Table"]]["Time"].==slices[t[2]][format[2]]))
        ts[t]["Series"]=tableData[dataPath][ts[t]["Table"]][ts[t]["Column"]][stpos:stpos+slices[t[2]]["Length"]-1]
    end
    ts=collapseData(ts,1,2);

    return ts
end



function generateTimeSeries(tableData,dataPath,slices,ts_info,format)


    ts=indexElement(tableData,(dataPath,ts_info,["Name","Slice"]));
    #### Find the start of each slice in the corresponding table
    rearrangeSet=sort(collect(keys(ts)))

    #for t in keys(ts)
    for t in rearrangeSet
        stpos=findfirst((tableData[dataPath][ts[t]["Table"]]["Date"].==slices[t[2]][format[1]]) .* (tableData[dataPath][ts[t]["Table"]]["Time"].==slices[t[2]][format[2]]))
        ts[t]["Series"]=tableData[dataPath][ts[t]["Table"]][ts[t]["Column"]][stpos:stpos+slices[t[2]]["Length"]-1]
    end
    ts=collapseData(ts,1,2);

    return ts
end




function generateTimeSeries(tableData,dataPath,slices,ts_info,format,indexingProp)

    mapProp=Dict(slices[j][indexingProp]=>j for j in keys(slices));

    ts=indexElement(tableData,(dataPath,ts_info,["Name","Slice"]));
    #### Find the start of each slice in the corresponding table
    rearrangeSet=sort(collect(keys(ts)))

    #for t in keys(ts)
    for t in rearrangeSet
        stpos=findfirst((tableData[dataPath][ts[t]["Table"]]["Date"].==slices[mapProp[t[2]]][format[1]]) .* (tableData[dataPath][ts[t]["Table"]]["Time"].==slices[mapProp[t[2]]][format[2]]))
        ts[t]["Series"]=tableData[dataPath][ts[t]["Table"]][ts[t]["Column"]][stpos:stpos+slices[mapProp[t[2]]]["Length"]-1]
    end
    ts=collapseData(ts,1,2);

    return ts
end


function genTimeSeries!(datawpath,tableSlices,tableTS,outputName;tsCol="Name",sliceCol="Slice",format=("Start_date","Start_time"))

    off,pos=findmax(datawpath[tableSlices]["Offset"]);
    tsLength=off+datawpath[tableSlices]["Length"][pos]-1;

    ts_names=unique(datawpath[tableTS][tsCol]);

    datawpath[outputName]=Dict("Name"=>[ts_names[i] for i in eachindex(ts_names)],"Series"=>Any[zeros(tsLength) for i in eachindex(ts_names)]);
    #datawpath[outputName]["Index"]=[i for i in eachindex(ts_names)]

    for i in eachindex(datawpath[tableTS][tsCol])
        cf=datawpath[tableTS]["Name"][i],datawpath[tableTS]["Table"][i],datawpath[tableTS]["Column"][i],datawpath[tableTS]["Scaling"][i],datawpath[tableTS]["Slice"][i] #TS info for each position
        tsdata=datawpath[cf[2]][cf[3]]; ### time series
        matchSlice=findfirst(datawpath[tableSlices][tsCol].==cf[5])  #matching slice col

        if typeof(datawpath[cf[2]]["Date"])==Array{Date,1} && typeof(datawpath[tableSlices]["Start_date"][matchSlice])==String
            stpos=findfirst((datawpath[cf[2]]["Date"].==Date(datawpath[tableSlices]["Start_date"][matchSlice],"d/m/y")) .* (datawpath[cf[2]]["Time"].==datawpath[tableSlices]["Start_time"][matchSlice]))
            matchTS=findfirst(datawpath[outputName][tsCol].==cf[1])
        else
            stpos=findfirst((datawpath[cf[2]]["Date"].==datawpath[tableSlices]["Start_date"][matchSlice]) .* (datawpath[cf[2]]["Time"].==datawpath[tableSlices]["Start_time"][matchSlice]))
            matchTS=findfirst(datawpath[outputName][tsCol].==cf[1])
        end

        datawpath[outputName]["Series"][matchTS][datawpath[tableSlices]["Offset"][matchSlice]:datawpath[tableSlices]["Offset"][matchSlice]+datawpath[tableSlices]["Length"][matchSlice]-1]=cf[4].*tsdata[stpos:stpos+datawpath[tableSlices]["Length"][matchSlice]-1]
        #datawpath[outputName]["Series"][matchTS][datawpath[tableSlices]["Offset"][matchSlice]:datawpath[tableSlices]["Offset"][matchSlice]+datawpath[tableSlices]["Length"][matchSlice]-1].=tsdata[1:100]

        #println("[$(datawpath[tableSlices]["Offset"][matchSlice]):$(datawpath[tableSlices]["Offset"][matchSlice]+datawpath[tableSlices]["Length"][matchSlice]-1)]")
        #println("$(length(tsdata[(stpos):(stpos+datawpath[tableSlices]["Length"][matchSlice]-1)]))")
    end
end


function formSets(tabdata,path)
    setList=findall(tabdata[path]["properties"]["set_name"].!=="")
    setNames=tabdata[path]["properties"]["set_name"][setList]
    setIndexing=tabdata[path]["properties"]["property"][setList]
    setGroup=tabdata[path]["properties"]["table_name"][setList]

    sets=Dict{Any,Any}(setNames[i]=>collect(keys(indexElement(tabdata,(path,setGroup[i],setIndexing[i])))) for i in eachindex(setList))
    return sets
end

function formSets(tabdata,path,subpath::String)


    setList=findall(tabdata[path][subpath]["set_name"].!=="")
    setNames=tabdata[path][subpath]["set_name"][setList]
    setIndexing=tabdata[path][subpath]["indexing"][setList]
    setGroup=tabdata[path][subpath]["table_name"][setList]

    sets=Dict{Any,Any}(setNames[i]=>sort(collect(keys(indexElement(tabdata,(path,setGroup[i],setIndexing[i]))))) for i in eachindex(setList))
    return sets
end



function addVariables!(m::String,variables,sets)

    for v in keys(variables)
        if "lower_bound" in keys(variables[v])
            eval(Meta.parse("@variable($m,$(v)[$(sets[variables[v]["index1"]]),$(sets[variables[v]["index2"]])],lower_bound=$(variables[v]["lower_bound"]))"))
        else
            eval(Meta.parse("@variable($m,$(v)[$(sets[variables[v]["index1"]]),$(sets[variables[v]["index2"]])])"))
        end

    end

end


function addIntIndex!(tabdata,path,object,full_set)

    tabdata[path][object]["Index"]=[i for i in eachindex(tabdata[path][object][full_set])]

end


function addIntIndex!(tabwpath,object,full_set)

    tabwpath[object]["Index"]=[i for i in eachindex(tabwpath[object][full_set])]

end



function formSet(object,pair::Pair)
    return sort(collect(keys(subindexElement(object,pair.first=>pair.second))))
end

function formSet(object)
    return sort(collect(keys(object)))
end


function scaleData!(tabwpath,object,toscale,scaling)

    for i in eachindex(tabwpath[object][toscale])
        if !ismissing(tabwpath[object][toscale][i])
            tabwpath[object][toscale][i]*=tabwpath[object][scaling][i]
        end
    end


end




function combineTables(datawpath,outputTable,tabs,indexing;deleteSource=false)

    Names=[unique(datawpath[tabs[i]][indexing]) for i in eachindex(tabs)]
    Keys=[collect(keys(datawpath[tabs[i]])) for i in eachindex(tabs)]

    unNames=[];
    unKeys=[];
    for t in eachindex(tabs)
        for n in eachindex(Names[t]) push!(unNames,Names[t][n]) end;
        for n in eachindex(Keys[t]) push!(unKeys,Keys[t][n]) end;
    end

    sort!(unique!(unNames)); sort!(unique!(unKeys));

    outputDict=Dict{Any,Any}(unKeys[j]=>Any[missing for i in eachindex(unNames)] for j in eachindex(unKeys))
    outputDict[indexing]=unNames;

    for t in eachindex(tabs)
        for i in Keys[t]
            for j in Names[t]
                idxOut=findfirst(j.==outputDict[indexing])
                idxSrc=findfirst(j.==Names[t])
                outputDict[i][idxOut]=datawpath[tabs[t]][i][idxSrc]
            end
        end
        if deleteSource
            pop!(datawpath,tabs[t])
        end
    end

    datawpath[outputTable]=Dict{Any,Any}()
    datawpath[outputTable]=outputDict;

end


#ps=loadSets(data["Operation"],"sets");

function loadSetsOld(datawpath,setTable);

    modelNameCol="model_name";
    tableNameCol="table_name";
    indexingCol="indexing";
    filterCon="filtergroup_";
    filterVal="filterval_";

    output=Dict{Any,Any}();

    headers=collect(keys(datawpath[setTable]));
    maxFilters=sum(occursin.(filterCon,headers));

    filters=headers[findall(occursin.(filterCon,headers))]
    filtervals=replace.(filters,filterCon=>filterVal)

    for i in eachindex(datawpath[setTable][modelNameCol])
        pr=datawpath[setTable][modelNameCol][i],datawpath[setTable][tableNameCol][i],datawpath[setTable][indexingCol][i],filters;

        filttoapply=pr[4][findall(.!(ismissing.([datawpath[setTable][pr[4][j]][i] for j in eachindex(pr[4])])))];

        dataOutput=deepcopy(datawpath[pr[2]][pr[3]]);

        for k in filttoapply
            filteredData=datawpath[setTable][k][i]
            k2=replace.(k,filterCon=>filterVal)
            filteredVal=datawpath[setTable][k2][i]

            missingV=findall(ismissing.(datawpath[pr[2]][filteredData]));
            datawpath[pr[2]][filteredData][missingV].=""

            filtOutput=datawpath[pr[2]][pr[3]][findall(filteredVal.==datawpath[pr[2]][filteredData])]

            datawpath[pr[2]][filteredData][findall(ismissing.(datawpath[pr[2]][filteredData]))].=missing
            intersect!(dataOutput,filtOutput)
            datawpath[pr[2]][filteredData][missingV].=missing
        end

        output[datawpath[setTable][modelNameCol][i]]=dataOutput


    end

    return output

end





function joinSets!(setdict,output,tojoin::Array;deleteSource=false)

    combinedElements=Any[];
    for i in eachindex(tojoin)
        for j in eachindex(setdict[tojoin[i]])
            push!(combinedElements,setdict[tojoin[i]][j])
        end
        if deleteSource
            pop!(setdict,tojoin[i])
        end

    end

    setdict[output]=combinedElements;

end


function diffSets!(setdict,output,todiff::Array;deleteSource=false)

    diffElements=deepcopy(setdict[todiff[1]]);

    for i in 2:length(todiff)
        setdiff!(diffElements,setdict[todiff[i]])


        if deleteSource
            pop!(setdict,todiff[i])
        end

    end

    if deleteSource pop!(setdict,todiff[1])  end

    setdict[output]=diffElements;

end



function loadParameters(datawpath,paramTable);

    modelNameCol="model_name";
    tableNameCol="table_name";
    indexingCol="indexing";
    colNameCol="column_name";

    #### Multi dimension array
    tableData="table_data";
    rowData="row_data";
    columnData="column_data";
    ####################

    output=Dict{Any,Any}();

    headers=collect(keys(datawpath[paramTable]));
    maxDepth=sum(occursin.(indexingCol,headers));

    indxIndex=headers[findall(occursin.(indexingCol,headers))]
    indxTables=replace.(indxIndex,indexingCol=>tableNameCol)
    indxColumns=replace.(indxIndex,indexingCol=>colNameCol)
    
    for i in eachindex(datawpath[paramTable][modelNameCol])

        ## (part 1) the original had indxIndex instead of IndxTables... indxtables always exist whereas indexing not always, therefore this change is needed for some cases.
        #indx=indxIndex[findall(.!(ismissing.([datawpath[paramTable][indxIndex[j]][i] for j in eachindex(indxIndex)])))]; #name of indexing columns to apply
        indx=indxTables[findall(.!(ismissing.([datawpath[paramTable][indxTables[j]][i] for j in eachindex(indxTables)])))]; #name of indexing columns to apply

        ## (part 2) the original had indxIndex instead of IndxTables... indxtables always exist whereas indexing not always, therefore this change is needed for some cases.
        #idx=replace.(indx,indexingCol=>"")
        idx=replace.(indx,tableNameCol=>"")
        dp=length(idx);

        if dp<=1
            if dp==0
                pr=[datawpath[paramTable][modelNameCol][i],datawpath[paramTable][tableNameCol][i],datawpath[paramTable][indexingCol][i],datawpath[paramTable][colNameCol][i]];
            else
                pr=[datawpath[paramTable][modelNameCol][i],datawpath[paramTable][tableNameCol*idx[1]][i],datawpath[paramTable][indexingCol*idx[1]][i],datawpath[paramTable][colNameCol*idx[1]][i]];
            end

            if !ismissing(pr[3])

                    dataOutput=deepcopy(datawpath[pr[2]][pr[4]]);
                    dataKeys=deepcopy(datawpath[pr[2]][pr[3]]);

                    output[datawpath[paramTable][modelNameCol][i]]=Dict{Any,Any}(dataKeys.=>dataOutput)

            else
                    dataOutput=deepcopy(datawpath[pr[2]][pr[4]][1]);
                    output[datawpath[paramTable][modelNameCol][i]]=dataOutput

            end

        elseif dp>=2

            rowField=findfirst("_"*replace(datawpath[paramTable][rowData][i][findfirst("\$", datawpath[paramTable][rowData][i])[1]:end],"\$"=>"").==idx)
            columnField=findfirst("_"*replace(datawpath[paramTable][columnData][i][findfirst("\$", datawpath[paramTable][columnData][i])[1]:end],"\$"=>"").==idx)
            columnName=replace(datawpath[paramTable][columnData][i][1:findfirst("\$", datawpath[paramTable][columnData][i])[1]],"\$"=>"")
            ### row field and column field positions with respect to idx (list of all instances with numbering)

            #### dataTable Information ######
            colTable=datawpath[paramTable][tableNameCol*idx[columnField]][i]
            colColumn=datawpath[paramTable][indexingCol*idx[columnField]][i]

            rowTable=datawpath[paramTable][tableNameCol*idx[rowField]][i]
            rowColumn=datawpath[paramTable][indexingCol*idx[rowField]][i]

            dataTable=datawpath[paramTable][tableData][i];
            colNameTable=columnName;

            headersTable=collect(keys(datawpath[dataTable]))
            rowsTable=datawpath[dataTable][rowColumn];

            ####################

            ### Sets for the columns from the parameters table in the order of idx
            tabStruct=[datawpath[paramTable][indxTables[j]][i] for j in eachindex(idx)] ## Order of how the dictionary will be created e.g. generators -> structure
            colStruct=[datawpath[paramTable][indxIndex[j]][i] for j in eachindex(idx)]
            indStruct=[datawpath[tabStruct[j]][colStruct[j]] for j in eachindex(idx)]

            outputData=Dict{Any,Any}(indStruct[1][l]=>Dict{Any,Any}(indStruct[2][n]=>missing for n in eachindex(indStruct[2])) for l in eachindex(indStruct[1]))

            #### Populate Data

                #Check how's the indexing of the outputData dictionary to that of the data in the dataTable
                if tabStruct[1]==rowTable #This means that the table and the outer index of the outputdata are the same
                    for ii in indStruct[1]
                        for jj in indStruct[2]
                            if ii in rowsTable && colNameTable*string(jj) in headersTable
                                corrRow=findfirst(ii.==rowsTable);
                                outputData[ii][jj]=datawpath[dataTable][colNameTable*string(jj)][corrRow]
                            end
                        end
                    end
                elseif tabStruct[2]==rowTable #This means that the tables are crossed between outputdata and data tables
                    for ii in indStruct[1]
                        for jj in indStruct[2]
                            if colNameTable*string(ii) in headersTable && jj in rowsTable
                                corrRow=findfirst(jj.==rowsTable)
                                outputData[ii][jj]=datawpath[dataTable][colNameTable*string(ii)][corrRow]
                            end
                        end
                    end
                end



            output[datawpath[paramTable][modelNameCol][i]]=outputData;
        end


    end

    return output

end


function loadSets(datawpath,setTable);

    modelNameCol="model_name";
    typeCol="element_table";
    indexingCol="indexing_table";
    indexingProp="indexing_property";

    matchingCol="matching_property";
    innerIndex="inner_indexing";

    filterCon="filtergroup_";
    filterVal="filterval_";


    output=Dict{Any,Any}();

    headers=collect(keys(datawpath[setTable]));


    maxFilters=sum(occursin.(filterCon,headers)); #Number of filters
    filterKeys=replace.(headers[findall(occursin.(filterCon,headers))],filterCon=>"");


    for i in eachindex(datawpath[setTable][modelNameCol])

        #pr=outputName,element,indexing_element,indexing_property,matching_prop,filterKeys
        pr=datawpath[setTable][modelNameCol][i],datawpath[setTable][typeCol][i],datawpath[setTable][indexingCol][i],datawpath[setTable][indexingProp][i],datawpath[setTable][matchingCol][i],datawpath[setTable][innerIndex][i],filterKeys

        filttoapply=pr[7][findall(.!(ismissing.([datawpath[setTable][filterCon.*pr[7][j]][i] for j in eachindex(pr[7])])))];

        #if pr[2]==pr[3]

        if ismissing(pr[5]) ###Needs to remove the need of columns
        
            dataOutput=deepcopy(datawpath[pr[3]][pr[4]]);

            for k in filterCon.*filttoapply
                filteredData=datawpath[setTable][k][i]
                k2=replace.(k,filterCon=>filterVal)
                filteredVal=datawpath[setTable][k2][i]

                missingV=findall(ismissing.(datawpath[pr[3]][filteredData]));
                datawpath[pr[3]][filteredData][missingV].=""

                filtOutput=datawpath[pr[3]][pr[4]][findall(string(filteredVal).==string.(datawpath[pr[3]][filteredData]))]

                if sum(ismissing.(datawpath[pr[3]][filteredData]))>=1
                    datawpath[pr[3]][filteredData][findall(ismissing.(datawpath[pr[3]][filteredData]))].=missing
                end
                intersect!(dataOutput,filtOutput)
                if length(missingV)>0
                    datawpath[pr[3]][filteredData][missingV].=missing
                end
            end

                output[datawpath[setTable][modelNameCol][i]]=dataOutput

        else


            dataOutput=Dict{Any,Any}(datawpath[pr[3]][pr[4]].=>"")

            for l in keys(dataOutput)
                #dataOutput[l]=datawpath[pr[2]][pr[4]][findall(datawpath[pr[2]][pr[5]].==l)]
                dataOutput[l]=datawpath[pr[2]][pr[6]][findall(datawpath[pr[2]][pr[5]].==l)] #Inner indexing column correction
                if !ismissing(dataOutput[l])
                    if length(dataOutput[l])>0
                        if typeof(dataOutput[l][1])<:Array || typeof(dataOutput[l][1])<:UnitRange
                            dataOutput[l]=dataOutput[l][1]
                        end
                    end
                end

            end

            for k in filterCon.*filttoapply
                filteredData=datawpath[setTable][k][i]
                k2=replace.(k,filterCon=>filterVal)
                filteredVal=datawpath[setTable][k2][i]

                 missingV=findall(ismissing.(datawpath[pr[2]][filteredData]));
                 datawpath[pr[2]][filteredData][missingV].=""

                #filtOutput=datawpath[pr[2]][pr[4]][findall(filteredVal.==datawpath[pr[2]][filteredData])]
                filtOutput=datawpath[pr[2]][pr[6]][findall(filteredVal.==datawpath[pr[2]][filteredData])]

                #Add the try bypass
                try
                    datawpath[pr[2]][filteredData][findall(ismissing.(datawpath[pr[2]][filteredData]))].=missing
                catch end

                for l in keys(dataOutput)
                    intersect!(dataOutput[l],filtOutput)
                end

                #Add the try bypass
                try
                    datawpath[pr[2]][filteredData][missingV].=missing
                catch end
            end

            output[datawpath[setTable][modelNameCol][i]]=dataOutput

        end


    end

    return output

end



function structureDefinitionBak(obj,nameV)

    keyList=sort(collect(keys(obj)));
    objTypes=DataType[];

    for k in keyList
        if typeof(obj[k])<:Dict
            kL=collect(keys(obj[k]))

            kType=eval(Symbol(replace(replace(string(typeof([kL[i] for i in eachindex(kL)])),"Vector{"=>""),"}"=>"")));
            emptyVals=[];
            for m in eachindex(kL)
                if !ismissing(obj[k][kL[m]])
                    for j in eachindex(obj[k][kL[m]])
                        push!(emptyVals,obj[k][kL[m]][j])
                    end
                else
                    push!(emptyVals,missing)
                end
            end
            vType=typeof([emptyVals[i] for i in eachindex(emptyVals)]);


            push!(objTypes,Dict{kType,vType})
        else

            #if isa(obj[k],Array) || isa(obj[k],UnitRange) || isa(obj[k],PooledArray)
            if isa(obj[k],Array) || isa(obj[k],UnitRange)
                objPush=typeof([obj[k][i] for i in eachindex(obj[k])]);
                if objPush==Array{Int64,1}
                    rangeComp=obj[k][1]:obj[k][end];
                    if collect(rangeComp)==[obj[k][i] for i in eachindex(obj[k])]
                        objPush=typeof(rangeComp)
                    end
                end
                push!(objTypes,objPush)
            else
                push!(objTypes,typeof(obj[k]))

            end

        end
    end


    io = IOBuffer();
    write(io,"struct $(nameV)_type\n")
    for k in eachindex(keyList)
        write(io,"  $(keyList[k])::$(objTypes[k])\n")
    end
    write(io,"end")
    #println(String(take!(io)))
    eval(Meta.parse(String(take!(io))))

end

macro var2string(arg)
   string(arg)
end

function transformStructureBak(obj,nameV);

    structureDefinition(obj,nameV); # Form the structure

    fn=fieldnames(eval(Symbol("$(nameV)_type")));
    ft=[fieldtype(eval(Symbol("$(nameV)_type")),fn[i]) for i in eachindex(fn)]

    stpar="$(nameV)_type("
    for n in eachindex(fn)
        println("ft $(ft[n]) , fn $(fn[n])")
        stpar*="$(convert(ft[n],obj[String(fn[n])]))"
        if n!=length(fn) stpar*="," else stpar*=")" end
    end

    return eval(Meta.parse(stpar))

end



function transformStructure(obj,nameV);

    structureDefinition(obj,nameV); # Form the structure

    fn=fieldnames(eval(Symbol("$(nameV)_type")));
    ft=[fieldtype(eval(Symbol("$(nameV)_type")),fn[i]) for i in eachindex(fn)]

    stpar="$(nameV)_type("
    for n in eachindex(fn)
        if isa(obj[String(fn[n])],Dict)
            #stpar*="$(convert(ft[n],obj[String(fn[n])]))"
            stpar*="$(convertDictionary(obj[String(fn[n])]))"
        else
            #println("Entering Here\ns")
            #println(ft[n],",",obj[String(fn[n])])
            if ft[n]==UnitRange{Int64}
                stpar*="$(obj[String(fn[n])][1]:obj[String(fn[n])][end])"
            else
                stpar*="$(convert(ft[n],obj[String(fn[n])]))"
            end
        end
        if n!=length(fn) stpar*="," else stpar*=")" end

    end

    #println(stpar)
    #return stpar
    return eval(Meta.parse(stpar))

end


function transformStructure!(obj,nameV);

    structureDefinition(obj,nameV); # Form the structure

    fn=fieldnames(eval(Symbol("$(nameV)_type")));
    ft=[fieldtype(eval(Symbol("$(nameV)_type")),fn[i]) for i in eachindex(fn)]

    stpar="$(nameV)=$(nameV)_type("
    for n in eachindex(fn)
        if isa(obj[String(fn[n])],Dict)
            #stpar*="$(convert(ft[n],obj[String(fn[n])]))"
            stpar*="$(convertDictionary(obj[String(fn[n])]))"
            #println("$(convertDictionary(obj[String(fn[n])]))")
        else
            #println("Entering Here\ns")
            stpar*="$(convert(ft[n],obj[String(fn[n])]))"
            #println("$(convert(ft[n],obj[String(fn[n])]))")
        end
        if n!=length(fn) stpar*="," else stpar*=")" end

    end

    #println(stpar)
    #return stpar
    return eval(Meta.parse(stpar))

end


function convertDictionary(dict)
    ks=collect(keys(dict));
    kst=typeof([ks[i] for i in eachindex(ks)]);
    ks=convert(kst,ks);

    kv=collect(values(dict));
    kvt=typeof([kv[i] for i in eachindex(kv)]);
    kv=convert(kvt,kv);

    return Dict(ks.=>kv)

end


function setDictType(dict)
    ks=collect(keys(dict));
    kst=typeof([ks[i] for i in eachindex(ks)]);
    #ks=convert(kst,ks);

    kv=collect(values(dict));
    kvt=typeof([kv[i] for i in eachindex(kv)]);
    #kv=convert(kvt,kv);

    kst=kst.parameters[1]
    kvt=kvt.parameters[1]

    return kst,kvt

end




function structureDefinition(obj,nameV)

    keyList=sort(collect(keys(obj)));
    objTypes=DataType[];

    for k in keyList
        if isa(obj[k],Dict)
            kType,vType=setDictType(obj[k])
            ##if kType<:InlineString kType=String end
            ##if vType<:InlineString vType=String end
            #println(kType,vType) ####
            push!(objTypes,Dict{kType,vType})
        else
            #if isa(obj[k],Array) || isa(obj[k],UnitRange) || isa(obj[k],PooledArray)
            if isa(obj[k],Array) || isa(obj[k],UnitRange)
                objPush=typeof([obj[k][i] for i in eachindex(obj[k])]);
                ##objPush=typeof([(typeof(obj[k][i])<:InlineString ? String(obj[k][i]) : obj[k][i]) for i in eachindex(obj[k])]);
                if objPush==Array{Int64,1}
                    rangeComp=obj[k][1]:obj[k][end];
                    if collect(rangeComp)==[obj[k][i] for i in eachindex(obj[k])]
                        objPush=typeof(rangeComp)
                    end
                end
                push!(objTypes,objPush)
            else
                ##if typeof(obj[k])<:InlineString
                ##    push!(objTypes,String)
                ##else
                    push!(objTypes,typeof(obj[k]))
                    #println(typeof(obj[k]))
                ##end
            end

        end
    end


    io = IOBuffer();
    write(io,"mutable struct $(nameV)_type\n")
    for k in eachindex(keyList)
        write(io,"  $(keyList[k])::$(objTypes[k])\n")
    end
    write(io,"end")
    #println(String(take!(io)))
    eval(Meta.parse(String(take!(io))))

end


function addOperationalLinkingOld!(m,linkingTable,parameterVariable,setVariable)

    linkingTable=data["Linkage"]["elements"];
    parameterVariable=pp;
    setVariable=ps;

    lVar="linking_var";
    #lVar="link";
    oParam="operation_parameter";
    oSet="operation_set";

    uniquelVars=unique(linkingTable[lVar]);
    lengthVars=zeros(Int64,length(uniquelVars));

    namesS=fieldnames(typeof(setVariable));

    for i in eachindex(uniquelVars)
        pos=findall(uniquelVars[i].==linkingTable[lVar]);
        for j in pos
            if ismissing(linkingTable[oSet][j])
                lengthVars[i]+=1;
            else
                npos=findfirst(Symbol(linkingTable[oSet][j]).==namesS);
                lengthVars[i]+=length(collect(getfield(setVariable,npos)));
            end
        end
    end

    for i in eachindex(uniquelVars)
        eval(Meta.parse("@variable(m,$(uniquelVars[i])[1:$(lengthVars[i])] == 1.0)"))
    end


    return m

end



function addOperationalLinking!(m,linkingTable,setVariable;lVar="linking_var",oVar="operation_variable",oSet="operation_set")

    #linkingTable=data["Linkage"]["elements"];
    #setVariable=ps;


    uniquelVars=unique(linkingTable[lVar]);
    lengthVars=zeros(Int64,length(uniquelVars));

    namesS=fieldnames(typeof(setVariable));

    for i in eachindex(uniquelVars)
        pos=findall(uniquelVars[i].==linkingTable[lVar]);
        for j in pos
            if ismissing(linkingTable[oSet][j])
                lengthVars[i]+=1;
            else
                npos=findfirst(Symbol(linkingTable[oSet][j]).==namesS);
                lengthVars[i]+=length(collect(getfield(setVariable,npos)));
            end
        end
    end


    varEq=[Any[] for i in eachindex(uniquelVars)];

    for i in eachindex(uniquelVars)
        ii=0;
        pos=findall(uniquelVars[i].==linkingTable[lVar]);
        for j in pos
            npos=findfirst(Symbol(linkingTable[oSet][j]).==namesS);
            nvar=linkingTable[oVar][j];
            for j in collect(getfield(setVariable,npos))
                ii+=1;
                push!(varEq[i],(ii,nvar,j))
            end
        end


        #eval(Meta.parse("@variable($m,$(uniquelVars[i])[1:$ii])"))
        println("@variable($m,$(uniquelVars[i])[1:$ii])")
        for jj in 1:ii
            #eval(Meta.parse("@constraint($m,$(uniquelVars[i])[$jj]==$m[:$(varEq[i][jj][2])][\"$(varEq[i][jj][3])\"])"))
            println("@constraint($m,$(uniquelVars[i])[$jj]==$m[:$(varEq[i][jj][2])][\"$(varEq[i][jj][3])\"])")
        end

    end


    return m

end





function augmentStructure(datawpath,tableName)
    nodeCol="Node";
    parentCol="Parent_Node";
    probCol="Probability";

    investCol="Investment" # Name of investment column (if applicable)
    operCol="Operation" # Name of investment column (if applicable)

    ##### Name for new calculated columns;
    cpCol="CondProb"; #Calculated "conditional" probability
    levCol="Level";
    pathCol="Path";
    paternCol="Parents";

    pathInv="ParInv"; #Parent investment nodes
    pathOp="ParOp"; # #Parent operation nodes
    ###########



    nodeList=collect(datawpath[tableName][nodeCol]);
    rootNode=nodeList[findfirst(ismissing.(collect(datawpath[tableName][parentCol])))];
    parentList=collect(datawpath[tableName][parentCol]);
    probList=collect(datawpath[tableName][probCol]);

    dictParent=Dict(nodeList[i]=>parentList[i] for i in eachindex(nodeList))
    dictProb=Dict(nodeList[i]=>probList[i] for i in eachindex(nodeList))

    family=Dict{Any,Any}([i=>Any[] for i in nodeList]);
    famProb=Dict{Any,Any}([i=>Float64[] for i in nodeList]);
    nodeLev=Dict{Any,Int64}([i=>0 for i in nodeList]);
    famCProb=Dict{Any,Any}([i=>0.0 for i in nodeList]);

    orderR=Dict{Int64,Any}([i=>missing for i in eachindex(nodeList)]);

    for i in nodeList
        push!(family[i],i) #Adding its name to the family
        push!(famProb[i],dictProb[i]) #Adding its name to the family
        orderR[findfirst(i.==nodeList)]=i;

        if rootNode!=i
            push!(family[i],dictParent[i])             #Adding its own name to the family
            push!(famProb[i],dictProb[family[i][end]]) #Adding its own probability to the probability family
            while family[i][end]!=rootNode
                push!(family[i],dictParent[family[i][end]])
                push!(famProb[i],dictProb[family[i][end]])
            end
        end
        nodeLev[i]=length(famProb[i])
        famCProb[i]=prod(famProb[i])
    end

    cpVec=Any[famCProb[orderR[i]] for i in eachindex(nodeList)]
    famVec=Any[family[orderR[i]] for i in eachindex(nodeList)]
    levVec=Any[nodeLev[orderR[i]] for i in eachindex(nodeList)]

    parVec=Any[family[orderR[i]][2:end] for i in eachindex(nodeList)]

    datawpath[tableName][cpCol]=cpVec
    datawpath[tableName][levCol]=levVec
    datawpath[tableName][pathCol]=famVec
    datawpath[tableName][paternCol]=parVec

end



function addMapTable(datawpath,structureTable,NameIndex)

    colInvestment="Investment"
    colOperation="Operation"
    colParents="Parents"
    #####
    colParentInv="ParentsInv"
    colParentOper="ParentsOper"



    Nrofcols=eachindex(datawpath[structureTable][colParents]);

    aux=datawpath[structureTable][NameIndex][findall([datawpath[structureTable][colInvestment][i] for i in Nrofcols].==1)];




    InvL=[intersect(aux,[datawpath[structureTable][colParents][i] for i in eachindex(datawpath[structureTable][colParents])][j]) for j in Nrofcols]

    aux=datawpath[structureTable][NameIndex][findall([datawpath[structureTable][colOperation][i] for i in Nrofcols].==1)];
    OpeL=[intersect(aux,[datawpath[structureTable][colParents][i] for i in eachindex(datawpath[structureTable][colParents])][j]) for j in Nrofcols]

    datawpath[structureTable][colParentInv]=InvL;
    datawpath[structureTable][colParentOper]=OpeL;


end




function addInvestmentLinking!(m,linkingTable,setVariable)

    #linkingTable=data["Linkage"]["elements"];
    #setVariable=ps;

    lVar="linking_var";
    oVar="investment_variable";
    oSet="passed_set";
    iSet="investment_set";


    uniquelVars=unique(linkingTable[lVar]);
    lengthVars=zeros(Int64,length(uniquelVars));

    namesS=fieldnames(typeof(setVariable));

    for i in eachindex(uniquelVars)
        pos=findall(uniquelVars[i].==linkingTable[lVar]);
        for j in pos
            if ismissing(linkingTable[oSet][j])
                lengthVars[i]+=1;
            else
                npos=findfirst(Symbol(linkingTable[oSet][j]).==namesS);
                lengthVars[i]+=length(collect(getfield(setVariable,npos)));
            end
        end
    end


    varEq=[Any[] for i in eachindex(uniquelVars)];

    for i in eachindex(uniquelVars)
        ii=0;
        pos=findall(uniquelVars[i].==linkingTable[lVar]);
        for j in pos
            npos=findfirst(Symbol(linkingTable[oSet][j]).==namesS);
            npos2=findfirst(Symbol(linkingTable[iSet][j]).==namesS);
            nvar=linkingTable[oVar][j];
            for j in collect(getfield(setVariable,npos))
                for k in collect(getfield(setVariable,npos2))
                    ii+=1;
                    push!(varEq[i],(ii,nvar,j,k))
                end
            end
        end


        eval(Meta.parse("@variable($m,$(uniquelVars[i])[1:$ii])"))
        for jj in 1:ii
            eval(Meta.parse("@constraint($m,$(uniquelVars[i])[$jj]==$m[:$(varEq[i][jj][2])][\"$(varEq[i][jj][3])\",\"$(varEq[i][jj][4])\"])"))
            #println("@constraint($m,$(uniquelVars[i])[$jj]==$m[:$(varEq[i][jj][2])][\"$(varEq[i][jj][3])\",\"$(varEq[i][jj][4])\"])")
        end

    end


    return m

end




function linkingTable(linkFolder)

    ### Required fields
    namesOfTables=["definitions","variables","indices"]; #Information about models, name of variables to be matched, name of elements
    definitionFields=["model_name","set_container","model_type","parameter_container"];
    variableFields=["linking_variable","operation_variable","operation_set","investment_variable","investment_set","investment_nodeset","type","sp_minimum","sp_scale"];
    indicesFields=["operation_set","operation_index","investment_set","investment_index"];
    vectTables=[definitionFields,variableFields,indicesFields];
    dictTables=Dict(namesOfTables[i]=>vectTables[i] for i in eachindex(namesOfTables));

    #### Optional Field
    scalingN="scaling"


    ### Check for required fields
    tablesAvailable=collect(keys(linkFolder));
    innerTables=Dict{Any,Any}();

    for i in namesOfTables
        i in tablesAvailable ? innerTables[i]=dictTables[i] : throw(ArgumentError("Required table \"$i\" missing"))
    end


    innerAvailable=Dict(i=>collect(keys(linkFolder[i])) for i in tablesAvailable);

    for i in namesOfTables for j in innerTables[i]
            j in innerAvailable[i] ? nothing : throw(ArgumentError("Required column \"$j\" from table \"$i\" missing"))
    end end

    #### Table definitions
    modelNames=linkFolder[namesOfTables[1]][definitionFields[1]];
    modelSets=linkFolder[namesOfTables[1]][definitionFields[2]];
    modelTypes=linkFolder[namesOfTables[1]][definitionFields[3]];
    modelParameters=linkFolder[namesOfTables[1]][definitionFields[4]];

    #### Table variables
    lVars=linkFolder[namesOfTables[2]][variableFields[1]]
    operVars=linkFolder[namesOfTables[2]][variableFields[2]]
    operSet=linkFolder[namesOfTables[2]][variableFields[3]]
    invVars=linkFolder[namesOfTables[2]][variableFields[4]]
    invSet=linkFolder[namesOfTables[2]][variableFields[5]]
    invNodeSet=linkFolder[namesOfTables[2]][variableFields[6]]
    invType=linkFolder[namesOfTables[2]][variableFields[7]]
    spMin=linkFolder[namesOfTables[2]][variableFields[8]]
    spScal=linkFolder[namesOfTables[2]][variableFields[9]]

    scalingVar=[1.0 for i in eachindex(lVars)];
    if scalingN in innerAvailable[namesOfTables[2]]
        scalingVar=linkFolder[namesOfTables[2]][scalingN]
        for j in findall(ismissing.(linkFolder[namesOfTables[2]][scalingN])) scalingVar[j]=1.0 end
    end

    #### Table Indices
    idxOperSet=linkFolder[namesOfTables[3]][indicesFields[1]];
    idxOperIdx=linkFolder[namesOfTables[3]][indicesFields[2]];
    idxInvSet=linkFolder[namesOfTables[3]][indicesFields[3]];
    idxInvIdx=linkFolder[namesOfTables[3]][indicesFields[4]];

    ######## OPERATIONAL LINKING
    opModelRows=findfirst(modelTypes.=="operation"); ## Only for 1 sp model done right now
    invModelRows=findfirst(modelTypes.=="investment"); ## Only for 1 sp model done right now

    #addOperationalLinking!(eval(Meta.parse(modelNames[k])),linkFolder[namesOfTables[2]],lVars;lVar=variableFields[1],oVar=variableFields[2],oSet=variableFields[3])

    uniquelVars=unique(lVars);
    lengthVars=zeros(Int64,length(uniquelVars),2); # 2 because an operational length and an investment length


    dictMatch=Array{Any,1}();


#######
    for i in eachindex(uniquelVars)
        pos=findall(uniquelVars[i].==lVars);
        indexGeneral=0;
        for j in pos
            setVariableOp=eval(Meta.parse((modelSets[opModelRows])));
            setVariableInv=eval(Meta.parse((modelSets[invModelRows])));
            nposInvNode=findfirst(Symbol(invNodeSet[j]).==collect(fieldnames(typeof(setVariableInv))))
            invNodeList=getfield(setVariableInv,nposInvNode);
            #invNodeDict=Dict{Any,Int64}(invNodeList[i]=>i for i in eachindex(invNodeList));
            invNodeDict=Dict{Int64,Any}(i=>invNodeList[i] for i in eachindex(invNodeList));

            #spVariableVal=
            if ismissing(operSet[j])
                #lengthVars[i][1]+=1;
                #lengthVars[i][2]+=1;



                indexGeneral+=1
                v=eval(Meta.parse("minimum(["*modelParameters[invModelRows]*"."*spMin[j]*"[nn] for nn in "*modelSets[invModelRows]*"."*invNodeSet[j]*"])*"*string(spScal[j])));
                # if spScal[j]<=1
                #         v=floor(v,digits=4);
                #     else
                #         v=ceil(v,digits=4);
                #     end
                push!(dictMatch,[uniquelVars[i],indexGeneral,modelNames[opModelRows],operVars[j],missing,modelNames[invModelRows],invVars[j],missing,invNodeDict,scalingVar[j],invType[j],v])






            else


                nposOp=findfirst(Symbol(operSet[j]).==collect(fieldnames(typeof(setVariableOp)))); ##Operational set
                nposInv=findfirst(Symbol(invSet[j]).==collect(fieldnames(typeof(setVariableInv)))); ##Operational set
                #lengthVars[i][1]+=length(collect(getfield(setVariableOp,nposOp)))
                #lengthVars[i][2]+=length(collect(getfield(setVariableInv,nposInv)))

                for k in findall(operSet[j].==idxOperSet)
                    indexGeneral+=1
                    v=eval(Meta.parse("minimum(["*modelParameters[invModelRows]*"."*spMin[j]*"[\""*idxInvIdx[k]*"\"][nn] for nn in "*modelSets[invModelRows]*"."*invNodeSet[j]*"])*"*string(spScal[j])));
                    # if spScal[j]<=1
                    #     v=floor(v,digits=4);
                    # else
                    #     v=ceil(v,digits=4);
                    # end
                    push!(dictMatch,[uniquelVars[i],indexGeneral,modelNames[opModelRows],operVars[j],idxOperIdx[k],modelNames[invModelRows],invVars[j],idxInvIdx[k],invNodeDict,scalingVar[j],invType[j],v])





                end

            end
        end
    end

    rRows=findall([dictMatch[i][11] for i in 1:length(dictMatch)].=="r");
    hRows=findall([dictMatch[i][11] for i in 1:length(dictMatch)].=="h");
    cRows=findall([dictMatch[i][11] for i in 1:length(dictMatch)].=="c");

    dictMatch=vcat(dictMatch[rRows],dictMatch[hRows],dictMatch[cRows]);
    for i in eachindex(dictMatch) dictMatch[i][2]=i end

    ind2nodeDict=Dict{Any,Any}();
    for i in eachindex(dictMatch)
        if !ismissing(dictMatch[i][9])
            ind2nodeDict=dictMatch[i][9];#### This has to change
        end
    end

    return dictMatch,ind2nodeDict

end


function addLinking!(LinkingTable)

    tVars=[LinkingTable[i][1] for i in eachindex(LinkingTable)]
    uVars=unique(tVars);
    lnVars=zeros(Int64,length(uVars))

    for i in eachindex(uVars)
        pVars=findall(uVars[i].==tVars);
        lnVars[i]=length(pVars);
        eval(Meta.parse("@variable($(LinkingTable[pVars[1]][3]),$(uVars[i])[1:$(lnVars[i])])")) ### Add variables to operational Model
        eval(Meta.parse("@variable($(LinkingTable[pVars[1]][6]),$(uVars[i])[1:$(lnVars[i]),1:$(length(LinkingTable[1][9]))])")) ### Add variables to operational Model

        for jj in 1:lnVars[i]
            s="";
            if ismissing(LinkingTable[pVars[jj]][5])
                eval(Meta.parse("@constraint($(LinkingTable[pVars[1]][3]),$(LinkingTable[pVars[jj]][10])*$(LinkingTable[pVars[1]][3])[:$(uVars[i])][$jj]==$(LinkingTable[pVars[1]][3])[:$(LinkingTable[pVars[jj]][4])])"))
            else
                typeof(LinkingTable[pVars[jj]][5])==String ? s="\"" : s=""
                eval(Meta.parse("@constraint($(LinkingTable[pVars[1]][3]),$(LinkingTable[pVars[jj]][10])*$(LinkingTable[pVars[1]][3])[:$(uVars[i])][$jj]==$(LinkingTable[pVars[1]][3])[:$(LinkingTable[pVars[jj]][4])][$s$(LinkingTable[pVars[jj]][5])$s])"))
            end
            #eval(Meta.parse
            for kk in 1:length(LinkingTable[1][9])
                s2="";
                typeof(LinkingTable[1][9][kk])==String ? s2="\"" : s2=""
                if ismissing(LinkingTable[pVars[jj]][8])
                    eval(Meta.parse("@constraint($(LinkingTable[pVars[1]][6]),$(LinkingTable[pVars[1]][6])[:$(uVars[i])][$jj,$kk]==$(LinkingTable[pVars[1]][6])[:$(LinkingTable[pVars[jj]][7])][$s2$((LinkingTable[pVars[1]][9])[kk])$s2])"))
                else
                    s="";
                    typeof(LinkingTable[pVars[jj]][8])==String ? s="\"" : s=""
                    eval(Meta.parse("@constraint($(LinkingTable[pVars[1]][6]),$(LinkingTable[pVars[1]][6])[:$(uVars[i])][$jj,$kk]==$(LinkingTable[pVars[1]][6])[:$(LinkingTable[pVars[jj]][7])][$s$(LinkingTable[pVars[jj]][8])$s,$s2$((LinkingTable[pVars[1]][9])[kk])$s2])"))
                end
            end
        end
    end

end



function addLinkingLMP!(LinkingTable,nameLMP)

    tVars=[LinkingTable[i][1] for i in eachindex(LinkingTable)]
    uVars=unique(tVars);
    lnVars=zeros(Int64,length(uVars))

    lToriginal=deepcopy(LinkingTable);

    for i in eachindex(LinkingTable)
        LinkingTable[i][6]=nameLMP
    end

    for i in eachindex(uVars)
        pVars=findall(uVars[i].==tVars);
        lnVars[i]=length(pVars);

        eval(Meta.parse("@variable($(LinkingTable[pVars[1]][6]),$(uVars[i])[1:$(lnVars[i]),1:$(length(LinkingTable[1][9]))])")) ### Add variables to operational Model

        for jj in 1:lnVars[i]

            for kk in 1:length(LinkingTable[1][9])
                s2="";
                typeof(LinkingTable[1][9][kk])==String ? s2="\"" : s2=""
                if ismissing(LinkingTable[pVars[jj]][8])
                    eval(Meta.parse("@constraint($(LinkingTable[pVars[1]][6]),$(LinkingTable[pVars[1]][6])[:$(uVars[i])][$jj,$kk]==$(LinkingTable[pVars[1]][6])[:$(LinkingTable[pVars[jj]][7])][$s2$((LinkingTable[pVars[1]][9])[kk])$s2])"))
                else
                    s="";
                    typeof(LinkingTable[pVars[jj]][8])==String ? s="\"" : s=""
                    eval(Meta.parse("@constraint($(LinkingTable[pVars[1]][6]),$(LinkingTable[pVars[1]][6])[:$(uVars[i])][$jj,$kk]==$(LinkingTable[pVars[1]][6])[:$(LinkingTable[pVars[jj]][7])][$s$(LinkingTable[pVars[jj]][8])$s,$s2$((LinkingTable[pVars[1]][9])[kk])$s2])"))
                end
            end
        end
    end

    LinkingTable=lToriginal;

end

function writeUNC(linkingTable,nameV,n2i,i2n)

    io = IOBuffer();
    keyList=["ni","nx","nx0","nh","nc","h","c","n2i","i2n"];
    objTypes=["Int64","Int64","Int64","Int64","Int64","Array{Float64,2}","Array{Float64,2}","Dict{Any,Int64}","Dict{Int64,Any}"];

    write(io,"struct $(nameV)_type\n")
        for k in eachindex(keyList)
            write(io,"  $(keyList[k])::$(objTypes[k])\n")
        end
    write(io,"end")
    eval(Meta.parse(String(take!(io))))

    ########
    spNumber=length(linkingTable[1][9]); ### Number of subp to solve
    
    elM=unique([linkingTable[i][11] for i in 1:length(linkingTable)]);
    "r" in elM ? rRowLength=findlast([linkingTable[i][11] for i in 1:length(linkingTable)].=="r") : rRowLength=0; 
    "h" in elM ? hRowLength=findlast([linkingTable[i][11] for i in 1:length(linkingTable)].=="h")-rRowLength : hRowLength=0;
    "c" in elM ? cRowLength=findlast([linkingTable[i][11] for i in 1:length(linkingTable)].=="c")-(rRowLength+hRowLength) : cRowLength=0;
    ###

    keyList=["ni","nx","nx0","nh","nc","h","c"];

    hvec=zeros(spNumber,hRowLength);
    cvec=zeros(spNumber,cRowLength);

    for j in rRowLength+1:rRowLength+hRowLength
        nn=1
        for n in keys(linkingTable[j][9])
            s2="";
            typeof(linkingTable[j][9][n])==String ? s2="\"" : s2=""
            if ismissing(linkingTable[j][5])
                hvec[n,j-rRowLength]=fix_value(eval(Meta.parse("$(linkingTable[j][6])[:$(linkingTable[j][7])][$s2$(linkingTable[j][9][n])$s2]")))
                #println("hvec[$(n),$(j-rRowLength)]=fix_value($(linkingTable[j][6])[:$(linkingTable[j][7])][$s2$(linkingTable[j][9][n])$s2])")
            else
                s="";
                typeof(linkingTable[j][8])==String ? s="\"" : s=""
                hvec[n,j-rRowLength]=fix_value(eval(Meta.parse("$(linkingTable[j][6])[:$(linkingTable[j][7])][$s$(linkingTable[j][8])$s,$s2$(linkingTable[j][9][n])$s2]")));
                #println("hvec[$n,$(j-rRowLength)]=fix_value($(linkingTable[j][6])[:$(linkingTable[j][7])][$s$(linkingTable[j][8])$s,$s2$(linkingTable[j][9][n])$s2])")
            end
            #hvec[nn,n2i[n]]=fix_value(eval(Meta.parse("$(linkingTable[j][6])[:$(linkingTable[j][7])][$s2$(linkingTable[j][9][n])$s2]")))

            nn+=1
        end
    end

    for j in rRowLength+hRowLength+1:rRowLength+hRowLength+cRowLength
        for n in keys(linkingTable[j][9])
            s2="";
            typeof(linkingTable[j][9][n])==String ? s2="\"" : s2=""
            if ismissing(linkingTable[j][5])
                cvec[n,j-rRowLength-hRowLength]=fix_value(eval(Meta.parse("$(linkingTable[j][6])[:$(linkingTable[j][7])][$s2$(linkingTable[j][9][n])$s2]")))
                #println("cvec[$(n),$(j-rRowLength-hRowLength)]=fix_value($(linkingTable[j][6])[:$(linkingTable[j][7])][$s2$(linkingTable[j][9][n])$s2])")
            else
                s="";
                typeof(linkingTable[j][8])==String ? s="\"" : s=""
                cvec[n,j-rRowLength-hRowLength]=fix_value(eval(Meta.parse("$(linkingTable[j][6])[:$(linkingTable[j][7])][$s$(linkingTable[j][8])$s,$s2$(linkingTable[j][9][n])$s2]")))
                #println("cvec[$(n),$(j-rRowLength-hRowLength)]=fix_value($(linkingTable[j][6])[:$(linkingTable[j][7])][$s$(linkingTable[j][8])$s,$s2$(linkingTable[j][9][n])$s2])")
            end


        end
    end




    write(io,"$(nameV)_type($spNumber,$(rRowLength+hRowLength),$rRowLength,$(hRowLength),$(cRowLength),$(hvec),$(cvec),$(n2i),$(i2n))\n")
    eval(Meta.parse(String(take!(io))))


end

function addYrs(datawpath,tableName)

    outputcol="YearOff";

    nodecol="Node";
    lengthcol="Length";
    parentscol="Parents";

    nodeLen=Dict(datawpath[tableName][nodecol][i]=>datawpath[tableName][lengthcol][i] for i in eachindex(datawpath[tableName][nodecol]))
    lengthV=[0 for i in eachindex(datawpath[tableName][nodecol])]

    for i in eachindex(datawpath[tableName][nodecol])
        for j in datawpath[tableName][parentscol][i]
            lengthV[i]+=nodeLen[j]
        end
    end

    yrDict=Dict(datawpath[tableName][nodecol][i]=>lengthV[i] for i in eachindex(datawpath[tableName][nodecol]))

    datawpath[tableName][outputcol]=yrDict
end


function setSP!(b,lt)


    rh=vcat(findall([lT[i][11]=="r" for i in eachindex(lT)]),findall([lT[i][11]=="h" for i in eachindex(lT)]));
    c=findfirst([lT[i][11]=="c" for i in eachindex(lT)]):findlast([lT[i][11]=="c" for i in eachindex(lT)]);
    coff=findfirst([lT[i][11]=="c" for i in eachindex(lT)])-1

    for i in rh
        b.temp.x[lT[i][2],1] = lT[i][12];
    end

    for i in c
        b.temp.c[i-coff,1] = lT[i][12];
    end

    #for i in 1:size(b.data.unc.c)[1]
    #    b.temp.c[:,i]=b.temp.c[:,1]
    #end

end

function val2dict(obj)

    @assert typeof(obj)<:JuMP.Containers.DenseAxisArray

    axL=obj.axes;
    axLV=length(axL);

    dictRet=Dict{Any,Any}();

    if axLV==1
        dictRet=Dict(i => value(obj[i]) for i in  axL[1])
    elseif axLV==2
        dictRet=Dict(i => Dict(j => value(obj[i,j]) for j in axL[2]) for i in  axL[1])
    elseif axLV==3
        dictRet=Dict(i => Dict(j => Dict(k => value(obj[i,j,k]) for k in axL[3]) for j in axL[2]) for i in  axL[1])
    else
        throw("The function has not been defined for this type of object")
    end

    return dictRet

end


"""
    compactPeriods!(SetData)

Given a dictionary (SetData), the function compacts a list of periods if the datatype is arrays of array, into an array, without touching the indexing keys of the dictionary.
Using the append function; then if the resulting set of integers is contiguous, it will return a dictionary with UnitRanges, if not then it will be of Array{Int64,1}.
Bespoke functions for periods

# Example

Given a dictionary where the values are arrays of arrays (in this example an array of UnitRange)

```julia-repl
julia> psDict["H_Sc"]
Dict{Any, Any} with 3 entries:
  2 => UnitRange{Int64}[5:5, 6:6, 7:7, 8:8]
  3 => UnitRange{Int64}[9:9, 10:10, 11:11, 12:13]
  1 => UnitRange{Int64}[1:1, 2:2, 3:3, 4:4]

```

Compacts it, dropping the second level

```julia-repl
julia> compactPeriods!(psDict["H_Sc"])
Dict{Any, Any} with 3 entries:
  2 => 5:8
  3 => 9:13
  1 => 1:4
```

"""
function compactPeriods!(SetData)
#This function helps to append arrays of unit ranges into a single array or unitrange. (Applicable to periods)
  for k in keys(SetData)
      output=Int64[];
      for i in eachindex(SetData[k])
        append!(output,collect(SetData[k][i]))
      end
      sort!(output)
        if collect(output[1]:output[end]) == output
          output=output[1]:output[end]
        end
      SetData[k]=output
  end
  return SetData
end


"""
    createPeriodSubset!(outputDict,outGroup,outerIndex,innerIndex)

Creates a subset of periods, where there is substitution between elements and adds the dictionary to the `outputDict` in the key `outGroup`
See the examples for better understanding. Bespoke functions for periods

# Examples

Given a dictionary "Sce" (in this example a list of slices for by scenario)
```julia-repl
julia> psDict["Sce"]
Dict{Any, Any} with 3 entries:
  2 => ["season1_2", "season2_2", "season3_2", "season4_2"]
  3 => ["season1_3", "season2_3", "season3_3", "season4_3"]
  1 => ["season1_1", "season2_1", "season3_1", "season4_1"]
```

And a list of periods per slice:

```julia-repl
julia> psDict["H_Sl"]
Dict{Any, Any} with 12 entries:
    "season2_3" => 10:10
    "season3_3" => 11:11
    "season4_3" => 12:13
    "season3_2" => 7:7
    "season4_1" => 4:4
    "season2_2" => 6:6
    "season1_3" => 9:9
    "season4_2" => 8:8
    "season3_1" => 3:3
    "season2_1" => 2:2
    "season1_2" => 5:5
    "season1_1" => 1:1
```

Substitute the periods instead of the slices in the scenario list:

```julia-repl
julia> createPeriodSubset!(psDict,"H_Sc","Sce","H_Sl")
Dict{Any, Any} with 3 entries:
  2 => UnitRange{Int64}[5:5, 6:6, 7:7, 8:8]
  3 => UnitRange{Int64}[9:9, 10:10, 11:11, 12:13]
  1 => UnitRange{Int64}[1:1, 2:2, 3:3, 4:4]
```

And store it in the desired key:

```julia-repl
julia> psDict["H_Sc"]
Dict{Any, Any} with 3 entries:
  2 => UnitRange{Int64}[5:5, 6:6, 7:7, 8:8]
  3 => UnitRange{Int64}[9:9, 10:10, 11:11, 12:13]
  1 => UnitRange{Int64}[1:1, 2:2, 3:3, 4:4]

```
"""
function createPeriodSubset!(outputDict,outGroup,outerIndex,middleIndex)
    outputDict[outGroup]=Dict{Any,Any}(i=>[outputDict[middleIndex][outputDict[outerIndex][i][j]] for j in eachindex(outputDict[outerIndex][i])] for i in keys(outputDict[outerIndex]));
end





function fillEmptyStructure(nameV);

    fn=fieldnames(eval(Symbol("$(nameV)_type")));
    ft=[fieldtype(eval(Symbol("$(nameV)_type")),fn[i]) for i in eachindex(fn)]

    stpar="$(nameV)_type("
    for n in eachindex(fn)
        stpar*="$(ft[n])()"
        if n!=length(fn) stpar*="," else stpar*=")" end
    end
    return eval(Meta.parse(stpar))
end



function generateTablesfromXLSXNew2(path::String,filename::String)
    basepath=pwd();
    cd(path*"/");

    folders=DataFrame(XLSX.readtable(filename,"main"))
    folders=dropComments(folders)

    indices=[DataFrame(XLSX.readtable(folders.Filename[i],"index")) for i in eachindex(folders.Filename)]
    indices=dropComments.(indices)

    tables=[];
    for i in eachindex(folders.Filename)
        subtables=Array{Union{DataFrame,Missing}}(missing,length(indices[i].sheet_name));
        XLSX.openxlsx(folders.Filename[i]) do xf
            for j in eachindex(indices[i].sheet_name)
                subtables[j]=dropComments(DataFrame(XLSX.gettable(xf[iS(indices[i].sheet_name[j])])))
            end
        end
        push!(tables,subtables)
    end    
    
    columns=[[names(tables[i][j]) for j in eachindex(indices[i].sheet_name)] for i in eachindex(folders.Filename)];

    cd(basepath)

    outputTables=Dict{Any,Any}(folders.Data[i]=>Dict{Any,Any}(indices[i].table_name[j]=>Dict{Any,Any}(columns[i][j][k]=>iS.(tables[i][j][:,columns[i][j][k]]) for k in eachindex(columns[i][j])) for j in eachindex(indices[i].table_name)) for i in eachindex(folders.Filename));

    return outputTables

end


println("done in $(round(time()-a,digits=1))s ");

