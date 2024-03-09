####### ----------- CODE TO DISPLAY RESULTS IN EXCEL SHEET FOR UNDECOMPOSED MODEL  -------------- #######


###  ------   CODE TO PRODUCE OUTPUT EXCEL FILE WHEN THERE IS ONLY 1 INVESTMENT NODE  ------- ###

pN_vals = DataFrame(hcat(axes(mU[:pN])[1], 
                        value.(mU[:pN]).data, 
                        [mp.capex["M$i"] for i in 1:length(collect(keys(mp.capex)))], 
                        [mp.cf["M$i"] for i in 1:length(collect(keys(mp.cf)))]),
                        ["Investment","pN","CAPEX (GBP/MW(h)(/*C))", "FixOM"])

## OPERATIONAL DICTIONARY
d2e_op = Dict("House & Heat Store Shed" => ("qHShedP","value","qHShedN","value"),
              "Heat Store Energy" => ("qH","value"),
              "House Temp" => ("tInt","value"),
              "Heat line Transfer" => ("pHL", "value"),
              "Bus LMPs" => ("c18","dual"))

#### --------- CREATE SPREADSHEET --------- ####
XLSX.openxlsx("ResultsAlg0.xlsx",mode="w") do xf
    sheet=xf[1]
    XLSX.rename!(sheet, "Investments")
    sheet["A1"] = "Final Investment Model Objective Value:"
    sheet["B1"] = value.(mU[:f]) + sum(mp.kappa[i]*mp.prob[i]*value.(mU[:beta])[i] for i in ms.I)
    XLSX.writetable!(sheet, pN_vals, anchor_cell=XLSX.CellRef("A4"))
    for k in collect(keys(d2e_op))
        XLSX.sheetnames(xf)[1]=="blank" ? XLSX.rename!(sheet,"$k") : XLSX.addsheet!(xf,"$k")
        sheet=xf[XLSX.sheetcount(xf)]
        mtx_all = DataFrame()
        for i in 2:2:length(d2e_op[k])
            if d2e_op[k][i] == "value"
                mtx_new = DataFrame(value.(mU[Symbol(d2e_op[k][i-1])]).data[1,:,:]',:auto)
                names = ["$(d)_$(d2e_op[k][i-1])_value" for d in mU[Symbol(d2e_op[k][i-1])].axes[2]]
                rename!(mtx_new,names)
            elseif d2e_op[k][i] == "dual"
                mtx_new = DataFrame(dual.(mU[Symbol(d2e_op[k][i-1])]).data[1,:,:]',:auto)
                names = ["$(d)_$(d2e_op[k][i-1])_dual" for d in mU[Symbol(d2e_op[k][i-1])].axes[2]]
                rename!(mtx_new,names)
            end
            mtx_all = hcat(mtx_new,mtx_all)
        end
        tdf=DataFrame("Period"=>1:length(ps.H))
        XLSX.writetable!(sheet,hcat(tdf,mtx_all),anchor_cell=XLSX.CellRef("A1"))
    end
end

