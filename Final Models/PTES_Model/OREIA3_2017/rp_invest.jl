# ## Display Results of Investment Model
# header=reshape(vcat(["Node"],mU[:pN].axes[2]),1,:);
# cols=hcat(mU[:pN].axes[1],value.(mU[:pN]).data);
# Full=vcat(header,cols)
# CSV.write("$path/Optimal_investments.csv",Tables.table(Full),header=false)
# Full


######## ----------  CODE TO DISPLAY RESULTS FOR ALG == 1  ----------- #########

### ---- Display Results of Operational Model for 1 investment node ---- ###
# N1 Newly invested capacity
avg_power_outputs = DataFrame(
                        hcat(["Diesel","CCGT","OCGT","Coal","Wind","Onshore_Wind","Solar","Nuclear","BioEnergy", "PumpHydro_Charging","PumpHydro_Discharging", "PumpHydro_Storage","Heat_Pump","HL1","HL2","HL3","L1"],
                             [sum(value.(S.ex.m[:pG]["Diesel",:]).data)/length(S.ex.m[:pG]["Diesel",:]), sum(value.(S.ex.m[:pG]["CCGT",:]).data)/length(S.ex.m[:pG]["CCGT",:]), 
                              sum(value.(S.ex.m[:pG]["OCGT",:]).data)/length(S.ex.m[:pG]["OCGT",:]), sum(value.(S.ex.m[:pG]["Coal",:]).data)/length(S.ex.m[:pG]["Coal",:]), 
                              sum(value.(S.ex.m[:pR]["Wind",:]).data)/length(S.ex.m[:pR]["Wind",:]), sum(value.(S.ex.m[:pR]["Onshore_Wind",:]).data)/length(S.ex.m[:pR]["Onshore_Wind",:]), sum(value.(S.ex.m[:pR]["Solar",:]).data)/length(S.ex.m[:pR]["Solar",:]), 
                              sum(value.(S.ex.m[:pNuc]["Nuclear",:]).data)/length(S.ex.m[:pNuc]["Nuclear",:]), sum(value.(S.ex.m[:pG]["BioEnergy",:]).data)/length(S.ex.m[:pG]["BioEnergy",:]), 
                              sum(value.(S.ex.m[:pS_in]["PumpedHydro",:]).data)/length(S.ex.m[:pS_in]["PumpedHydro",:]), sum(value.(S.ex.m[:pS_out]["PumpedHydro",:]).data)/length(S.ex.m[:pS_out]["PumpedHydro",:]), sum(value.(S.ex.m[:qB]["PumpedHydro_Store",:]).data)/length(S.ex.m[:qB]["PumpedHydro_Store",:]),
                              sum(value.(S.ex.m[:pH_heat]["HPump",:]).data)/length(S.ex.m[:pH_heat]["HPump",:]),
                              sum(value.(S.ex.m[:pHL]["HL1",:]).data)/length(S.ex.m[:pHL]["HL1",:]), sum(value.(S.ex.m[:pHL]["HL2",:]).data)/length(S.ex.m[:pHL]["HL2",:]), sum(value.(S.ex.m[:pHL]["HL3",:]).data)/length(S.ex.m[:pHL]["HL3",:]),
                              sum(value.(S.ex.m[:pL]["L1",:]).data)/length(S.ex.m[:pL]["L1",:])],
                             [findmax(value.(S.ex.m[:pG]["Diesel",:]).data)[1], findmax(value.(S.ex.m[:pG]["CCGT",:]).data)[1], findmax(value.(S.ex.m[:pG]["OCGT",:]).data)[1], findmax(value.(S.ex.m[:pG]["Coal",:]).data)[1], 
                              findmax(value.(S.ex.m[:pR]["Wind",:]).data)[1], findmax(value.(S.ex.m[:pR]["Onshore_Wind",:]).data)[1], findmax(value.(S.ex.m[:pR]["Solar",:]).data)[1], 
                              findmax(value.(S.ex.m[:pNuc]["Nuclear",:]).data)[1], findmax(value.(S.ex.m[:pG]["BioEnergy",:]).data)[1], 
                              findmax(value.(S.ex.m[:pS_in]["PumpedHydro",:]).data)[1], findmax(value.(S.ex.m[:pS_out]["PumpedHydro",:]).data)[1], findmax(value.(S.ex.m[:qB]["PumpedHydro_Store",:]).data)[1],
                              findmax(value.(S.ex.m[:pH_heat]["HPump",:]).data)[1],
                              findmax(value.(S.ex.m[:pHL]["HL1",:]).data)[1], findmax(value.(S.ex.m[:pHL]["HL2",:]).data)[1], findmax(value.(S.ex.m[:pHL]["HL3",:]).data)[1],
                              findmax(value.(S.ex.m[:pL]["L1",:]).data)[1]],
                             [findmax(value.(S.ex.m[:pG]["Diesel",:]).data)[2], findmax(value.(S.ex.m[:pG]["CCGT",:]).data)[2], findmax(value.(S.ex.m[:pG]["OCGT",:]).data)[2], findmax(value.(S.ex.m[:pG]["Coal",:]).data)[2], 
                             findmax(value.(S.ex.m[:pR]["Wind",:]).data)[2], findmax(value.(S.ex.m[:pR]["Onshore_Wind",:]).data)[2], findmax(value.(S.ex.m[:pR]["Solar",:]).data)[2], 
                             findmax(value.(S.ex.m[:pNuc]["Nuclear",:]).data)[2], findmax(value.(S.ex.m[:pG]["BioEnergy",:]).data)[2], 
                             findmax(value.(S.ex.m[:pS_in]["PumpedHydro",:]).data)[2], findmax(value.(S.ex.m[:pS_out]["PumpedHydro",:]).data)[2], findmax(value.(S.ex.m[:qB]["PumpedHydro_Store",:]).data)[2],
                             findmax(value.(S.ex.m[:pH_heat]["HPump",:]).data)[2],
                             findmax(value.(S.ex.m[:pHL]["HL1",:]).data)[2], findmax(value.(S.ex.m[:pHL]["HL2",:]).data)[2], findmax(value.(S.ex.m[:pHL]["HL3",:]).data)[2],
                             findmax(value.(S.ex.m[:pL]["L1",:]).data)[2]]
                            ),
                            ["Generator", "Avg. Power Output (GW)", "Max. Power Output (GW)","Index of Max Power Output value"])



#### ----- Calculate dual value of pN lower_bound constraint ----- ####
# wind_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["Wind","N1"]))
# onshore_wind_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["Onshore_Wind","N1"]))
# solar_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["Solar","N1"]))
# pumped_hydro_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["PumpedHydro","N1"]))
# pumped_hydro_store_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["PumpedHydro_Store","N1"]))

# PHES_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["PHES","N1"]))
# PHES_E_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["PHES_E","N1"]))

# #hpump_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["HPump","N1"]))
# TES_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["TTES","N1"]))

# diesel_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["Diesel","N1"]))
# OCGT_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["OCGT","N1"]))
# CCGT_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["CCGT","N1"]))
# coal_dual = dual.(LowerBoundRef(B.rmp.m[:pN]["Coal","N1"]))

# # Calculate Upper Bound Dual for non-investable technologies
# nuclear_dual = dual.(B.rmp.m[:c03i])
# biomass_dual = dual.(B.rmp.m[:c04i])





pN_vals = DataFrame(hcat(
                    axes(B.rmp.m[:pN])[1], 
                    value.(B.rmp.m[:pN]).data, 
                    value.(B.rmp.m[:pC]).data,
                    [mp.capex[i] for i in axes(B.rmp.m[:pN])[1]], 
                    [mp.cf[i] for i in axes(B.rmp.m[:pN])[1]],
                    [1/mp.yr*((mp.capex[i]/mp.lf[i])+(mp.cf[i])) for i in axes(B.rmp.m[:pN])[1]],
                    [dual.(LowerBoundRef(B.rmp.m[:pN][i,"N1"])) for i in axes(B.rmp.m[:pN])[1]]
                    ),
                    ["Investment","pN (GW(h)(/*C))","pC (GW(h)(/*C))","CAPEX (GBP/MW(h)(/*C))", "FixOM","Price/hour","Lower Bound Dual"])

## OPERATIONAL DICTIONARY
d2e_op = Dict("Heating Network" => ("pH_heat","value", "pHL","value","qH","value","qHS","value","tInt","value","qHShedP","value","qHShedN","value"),
              "Electrical Network" => ("pG","value","pR","value","pNuc","value","pS_in","value","pS_out","value","qB","value","pL","value","c19","dual"))

mU = S.ex.m

# Create Spreadsheet
XLSX.openxlsx("Operational Results Plots/ResultsAlg1.xlsx",mode="w") do xf
    sheet=xf[1]
    XLSX.rename!(sheet, "Investments")
    sheet["A1"] = "Investment Model Final Objective Function Value:"
    sheet["B1"] = "LB: "
    sheet["C1"] = B.hist.L[B.hist.k]
    sheet["D1"] = "UB: "
    sheet["E1"] = B.hist.U[B.hist.k]
    sheet["A2"] = "Total SP Evaluations"
    sheet["B2"] = B.hist.k*length(B.data.ms.I)
    sheet["C2"] = "Total Time Taken"
    sheet["D2"] = sum(B.hist.T[:,:])
	sheet["H1"] = "Lambda Value"
	sheet["I1"] = pp.lambda["PTES"]
	sheet["H2"] = "CO2 Tax"
	sheet["I2"] = mp.cco2["N1"]
    XLSX.writetable!(sheet, pN_vals, anchor_cell=XLSX.CellRef("A4"))
    XLSX.writetable!(sheet, avg_power_outputs, anchor_cell=XLSX.CellRef("J4"))
    for k in collect(keys(d2e_op))
        XLSX.sheetnames(xf)[1]=="blank" ? XLSX.rename!(sheet,"$k") : XLSX.addsheet!(xf,"$k")
        sheet=xf[XLSX.sheetcount(xf)]
        mtx_all = DataFrame()
        for i in 2:2:length(d2e_op[k])
            if d2e_op[k][i] == "value"
                mtx_new = DataFrame(value.(mU[Symbol(d2e_op[k][i-1])]).data[:,:]',:auto)
                names = ["$(d)_$(d2e_op[k][i-1])_value" for d in mU[Symbol(d2e_op[k][i-1])].axes[1]]
                rename!(mtx_new,names)
            elseif d2e_op[k][i] == "dual"
                mtx_new = DataFrame(dual.(mU[Symbol(d2e_op[k][i-1])]).data[:,:]',:auto)
                names = ["$(d)_$(d2e_op[k][i-1])_dual" for d in mU[Symbol(d2e_op[k][i-1])].axes[1]]
                rename!(mtx_new,names)
            end
            mtx_all = hcat(mtx_new,mtx_all)
        end
        tdf=DataFrame("Period"=>1:length(ps.H))
        XLSX.writetable!(sheet,hcat(tdf,mtx_all),anchor_cell=XLSX.CellRef("A1"))
    end
end




#### ----- CALCULATE PRICES/HOUR/MW OF EACH TECHNOLOGY ----- ####
wind_pph = 1/mp.yr*((mp.capex["Wind"]/mp.lf["Wind"])+(mp.cf["Wind"]))
onshore_wind_pph = 1/mp.yr*((mp.capex["Onshore_Wind"]/mp.lf["Onshore_Wind"])+(mp.cf["Onshore_Wind"]))
solar_pph = 1/mp.yr*((mp.capex["Solar"]/mp.lf["Solar"])+(mp.cf["Solar"]))

diesel_pph = 1/mp.yr*((mp.capex["Diesel"]/mp.lf["Diesel"])+(mp.cf["Diesel"]))
OCGT_pph = 1/mp.yr*((mp.capex["OCGT"]/mp.lf["OCGT"])+(mp.cf["OCGT"]))
CCGT_pph =  1/mp.yr*((mp.capex["CCGT"]/mp.lf["CCGT"])+(mp.cf["CCGT"]))
coal_pph = 1/mp.yr*((mp.capex["Coal"]/mp.lf["Coal"])+(mp.cf["Coal"]))
biomass_pph = 1/mp.yr*((mp.capex["BioEnergy"]/mp.lf["BioEnergy"])+(mp.cf["BioEnergy"]))

nuclear_pph = 1/mp.yr*((mp.capex["Nuclear"]/mp.lf["Nuclear"])+(mp.cf["Nuclear"]))









