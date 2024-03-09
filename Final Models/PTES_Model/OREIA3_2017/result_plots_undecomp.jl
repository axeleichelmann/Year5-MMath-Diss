### DEFINE FUNCTION TO EXRACT THE RESULTS THAT WE WANT TO PLOT FROM S.ex.m
function heat_line(hl::String)
    return value.(mU[:pHL]["N1", hl, :]).data
end

function house_temp(house::String)
    return value.(mU[:tInt]["N1", house, :]).data
end

function heat_store_temp(hns::String)
    return value.(mU[:qH]["N1", hns, : ]).data/value.(mU[:QMass]["N1", hns])
end

function heat_pump_lmp(hp::String)
    return dual.(mU[:c18]["N1", pp.hSupp[hp], :]).data
end

function thermal_generator(t::String)
    return value.(mU[:pG]["N1",t,:]).data
end

function renewable_generator(r::String)
    return value.(mU[:pR]["N1",r,:]).data
end

function electric_pumps(es::String)
    return value.(mU[:pS_in]["N1", es,:]).data - value.(mU[:pS_out]["N1",es,:]).data
end

function electric_storage(ens::String)
    return value.(mU[:qB]["N1", ens,:]).data
end

time_periods = ps.H

heatlines = Plots.plot(time_periods, heat_line(ps.HL[1]), label=ps.HL[1], xlabel="Time Period", ylabel="Power in Heat Line (MW)");
if length(ps.HL)>1
    for hl in ps.HL[2:end]
        heatlines = Plots.plot!(time_periods, heat_line(hl), label=hl, xlabel="Time Period", ylabel="Power in Heat Line (MW)");
    end
end
Plots.savefig("Operational Results Plots/HL_plots.pdf")

## Create graphs with house/heat store temperature against time period
temp_plots = Plots.plot(time_periods, house_temp(ps.HOUSES[1]), label=ps.HOUSES[1], xlabel="Time Period", ylabel="Internal Temp. (C)");
if length(ps.HOUSES)>1
    for house in ps.HOUSES[2:end]
        temp_plots = Plots.plot!(time_periods, house_temp(house), label=house, xlabel="Time Period", ylabel="Internal Temp. (C)");
    end
end
for hns in ps.HnS
    temp_plots = Plots.plot!(time_periods, heat_store_temp(hns), label=hns, xlabel="Time Period", ylabel="Internal Temp. (C)");
end
Plots.savefig("Operational Results Plots/temp_plots.pdf")

heat_pump_lmp_plots = Plots.plot(time_periods, heat_pump_lmp(ps.HP[1]), label=ps.HP[1], xlabel="Time Period", ylabel="LMP");
if length(ps.HP)>1
    for hp in ps.HP[2:end]
        heat_pump_lmp_plots = Plots.plot!(time_periods, heat_pump_lmp(hp), label=hp, xlabel="Time Period", ylabel="LMP");
    end
end
Plots.savefig("Operational Results Plots/HP_LMPs_plot.pdf")

# generators_plots = Plots.plot(time_periods, )


### CODE TO GET GENERATOR POWER, HOUSE & HEAT STORE TEMP, BUS 2 LMPS PLOTS ALL IN ONE GRAPH

## Add Heat Line data to the everything plot
# everything_plot = Plots.plot(time_periods, heat_line(ps.HL[1]), 
#                              label=ps.HL[1], xlabel="Time Period");
# if length(ps.HL)>1
#     for hl in ps.HL[2:end]
#         heatlines = Plots.plot!(time_periods, heat_line(hl), label=hl, xlabel="Time Period");
#     end
# end

# for house in ps.HOUSES    ## Add House Temperatures to the everything plot
#     everything_plot = Plots.plot!(time_periods, house_temp(house), label=house, xlabel="Time Period");
# end
 
# for hns in ps.HnS     ## Add Heat Store temperatures to the everything plot
#     everything_plot = Plots.plot!(time_periods, heat_store_temp(hns), label=hns, xlabel = "Time Period");
# end

# everything_plot = Plots.plot(time_periods, heat_pump_lmp(ps.HP[1]), 
#                              label="$(ps.HP[1]) Bus LMP", xlabel="Time Period");
# if length(ps.HP)>1
#     for hp in ps.HP[2:end]   ## Add Heat Pump LMP to the everything_plot
#         everything_plot = Plots.plot!(time_periods, heat_pump_lmp(hp), label="$(hp) Bus LMP", xlabel="Time Period");
#     end
# end


generator_plot = Plots.plot(time_periods, renewable_generator(ps.R[1]), label = ps.R[1])
if length(ps.R)>1
    for r in ps.R    ## Add renewable generators to the everything_plot
        generator_plot = Plots.plot!(time_periods, renewable_generator(r), label=r);
    end
end

for t in ps.T    ## Add thermal generators to the everything plot
    generator_plot = Plots.plot!(time_periods, thermal_generator(t), label=t);
end

for es in ps.ES   ## Add electric pump charge/discharge to everything plot
    generator_plot = Plots.plot!(time_periods, electric_pumps(es), label=es, ylabel="Power Generated/Discharged (MW)");
end

# for ens in ps.EnS   ## Add electric storages to everything plot
#     everything_plot = Plots.plot!(time_periods, electric_storage(ens), label=ens);
# end

Plots.savefig("Operational Results Plots/generator_plot.pdf")

