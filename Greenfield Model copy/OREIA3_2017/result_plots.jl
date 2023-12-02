### DEFINE FUNCTION TO EXRACT THE RESULTS THAT WE WANT TO PLOT FROM S.ex.m
function heat_line(hl::String)
    return value.(S.ex.m[:pHL][hl, :]).data
end

function house_temp(house::String)
    return value.(S.ex.m[:tInt][house, :]).data
end

function heat_store_temp(hns::String)
    return value.(S.ex.m[:qH][hns, : ]).data/value.(S.ex.m[:QMass][hns])
end

function heat_pump_lmp(hp::String)
    return dual.(S.ex.m[:c19][pp.hSupp(hp), :]).data
end

time_periods = ps.H

heatline_plots = Plots.plot(time_periods, heat_line("HL1"), label="HL1", xlabel="Time Period", ylabel="Power in Heat Line (MW)");
for hl in ps.HL[Not("HL1")]
    heatline_plots.plot!(time_periods, heat_line(hl), label=hl, xlabel="Time Period", ylabel="Power in Heat Line (MW)");
end
Plots.savefig("HL_plots.pdf")

house_temp_plots = Plots.plot(time_periods, [house_temp(house) for house in ps.HOUSES], xlabel="Time Period", ylabel="Internal House Temperature (C)");


heat_store_temp_plots = Plots.plot(time_periods, [heat_store_temp(hns) for hns in ps.HnS], xlabel="Time Period", ylabel="Internal Heat Store Temperature (C)");


Plots.savefig("operational_results.pdf")