## Display Results of Investment Model
header=reshape(vcat(["Node"],mU[:pN].axes[2]),1,:);
cols=hcat(mU[:pN].axes[1],value.(mU[:pN]).data);
Full=vcat(header,cols)
CSV.write("$path/Optimal_investments.csv",Tables.table(Full),header=false)
Full

#=

## Display Results of Investment Model for Node 1
header=reshape(vcat(["Investment Index","Investment Name"],mU[:pN].axes[2]),1,:);
cols=hcat(data["Linkage"]["indices"]["investment_index"],data["Linkage"]["indices"]["operation_index"],value.(mU[:pN]).data)
Full=vcat(header,cols)
Full

## Display results of Operational Model for Thermal Technology
# display results for "Diesel" at Node 1 for each time step, Repeat for "OCGT", "CCGT", "Coal"
time_steps_diesel=vcat(["Time_step"],mU[:pG].axes[3])
power_output_diesel=vcat(["Power Output (MWh)"],value.(mU[:pG]["N1","Diesel",:]).data)  ## IS IT (MWh)???
diesel_full = hcat(time_steps_diesel,power_output_diesel)
diesel_full

## Display Results of Operational Model for Renewable technology
# display results for "Wind" at Node 1 for each time step, Repeat for "Solar"
time_steps_wind=vcat(["Time_step"],mU[:pR].axes[3])
power_output_wind=vcat(["Power Output (MWh)"],value.(mU[:pR]["N1","Wind",:]).data)
wind_full = hcat(time_steps_wind,power_output_wind)
wind_full

=#

# A=LowerBoundRef.(mU[:pC]);

# open("test.csv", "w") do file
#     for i in ms.I 
#         i==ms.I[1] ? write(file,"Node/Tech,") : nothing
#         write(file,i)
#         i==ms.I[end] ? write(file,"\n") : write(file,",")
#     end
#     for m in ms.M
#         for n in ms.I
#             n==ms.I[1] ? write(file,m,",") : nothing
#             write(file, string(dual(A[m,n])))
#             n == ms.I[end] ? write(file,"\n") : write(file,",") 
#         end
#     end
# end




# for i in ms.I
#     open("renewables$i.csv", "w") do file

#     for r in ps.R 
#         r==ps.R[1] ? write(file,"Period/Tech,") : nothing
#         write(file,r)
#         r==ps.R[end] ? write(file,"\n") : write(file,",")
#     end
#     for t in ps.H
#         write(file, "$t,")
#         for r in ps.R
#             write(file, string(round(value(mU[:pR][i,r,t]),digits=5)))
#             r == ps.R[end] ? write(file,"\n") : write(file,",") 
#         end
#     end
# end
# end

# DF=CSV.read("renewablesN1.csv",DataFrame)

# #### OPERATIONAL DICTIONARY
# d2e=Dict("RenProd"=>("pR","value"),
#         "ThermalProd"=>("pG","value"),
#         "HeatStoCharging"=>("pH_in","value"),
#         "HeatStoDischarging"=>("pH_out","value"),
#         "ExampleDual"=>("c01","dual"),
#         )


# XLSX.openxlsx("OutputsOper.xlsx",mode="w") do xf
#     sheet=xf[1]
#     XLSX.rename!(sheet, "blank")
#     for k in collect(keys(d2e))
#         for (ii,i) in enumerate(mU[Symbol(d2e[k][1])].axes[1])
#             XLSX.sheetnames(xf)[1]=="blank" ? XLSX.rename!(sheet,"$(k)_$i") : XLSX.addsheet!(xf,"$(k)_$i")
#             sheet=xf[XLSX.sheetcount(xf)]
#             if d2e[k][2]=="value"
#                 mtx=DataFrame(value.(mU[Symbol(d2e[k][1])]).data[ii,:,:]'[:,:],:auto)
#             elseif d2e[k][2]=="dual"
#                 mtx=DataFrame(dual.(mU[Symbol(d2e[k][1])]).data[ii,:,:]'[:,:],:auto)
#             end
#             tdf=DataFrame("Period"=>mU[Symbol(d2e[k][1])].axes[3])
#             rename!(mtx,mU[Symbol(d2e[k][1])].axes[2])
#             XLSX.writetable!(sheet,hcat(tdf,mtx),anchor_cell=XLSX.CellRef("A1"))
#         end
#     end
# end


# #### INVESTMENT DICTIONARY
# d2e=Dict("NewCap"=>("pN","value"),
#         "TotCap"=>("pC","value"),
#         "CapDual"=>("c02i","dual"),
#         )

# XLSX.openxlsx("OutputsInv.xlsx",mode="w") do xf
#             sheet=xf[1]
#             XLSX.rename!(sheet, "blank")
#             for k in collect(keys(d2e))
#                 XLSX.sheetnames(xf)[1]=="blank" ? XLSX.rename!(sheet,k) : XLSX.addsheet!(xf,k)
#                 sheet=xf[XLSX.sheetcount(xf)]
#                 if d2e[k][2]=="value"
#                     mtx=DataFrame(value.(mU[Symbol(d2e[k][1])]).data,:auto)
#                 elseif d2e[k][2]=="dual"
#                     mtx=DataFrame(dual.(mU[Symbol(d2e[k][1])]).data,:auto)
#                 end
#                 tdf=DataFrame("Technology"=>mU[Symbol(d2e[k][1])].axes[1])
#                 rename!(mtx,mU[Symbol(d2e[k][1])].axes[2])
#                 XLSX.writetable!(sheet,hcat(tdf,mtx),anchor_cell=XLSX.CellRef("A1"))
#                 #XLSX.writetable!(sheet,mtx,anchor_cell=XLSX.CellRef("B1"))
#             end
#     end