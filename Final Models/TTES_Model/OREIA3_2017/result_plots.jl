### DEFINE FUNCTION TO EXRACT THE RESULTS THAT WE WANT TO PLOT FROM S.ex.m
function heat_line(hl::String)
    return value.(S.ex.m[:pHL][hl, :]).data
end

function house_temp(house::String)
    return value.(S.ex.m[:tInt][house, :]).data
end

function heat_store_temp(hns::String)
    return value.(S.ex.m[:qH][hns,:]).data/value.(S.ex.m[:QMass][hns])
end

function house_energy(house::String)
    return value.(S.ex.m[:qHS][house,:]).data .- minimum(value.(S.ex.m[:qHS][house,:]).data)
end

function heat_store_energy(hns::String)
    return value.(S.ex.m[:qH][hns,:]).data .- minimum(value.(S.ex.m[:qH][hns,:]).data)
end

function heat_pump_lmp(hp::String)
    return dual.(S.ex.m[:c19][pp.hSupp[hp], :]).data
end

function thermal_generator(t::String)
    return value.(S.ex.m[:pG][t,:]).data
end

function renewable_generator(r::String)
    return value.(S.ex.m[:pR][r,:]).data
end

function nuclear_generator(n::String)
    return value.(S.ex.m[:pNuc][n,:]).data
end


function electric_pumps(es::String)
    return value.(S.ex.m[:pS_in][es,:]).data - value.(S.ex.m[:pS_out][es,:]).data
end

function electric_storage(ens::String)
    return value.(S.ex.m[:qB][ens,:]).data
end

time_periods = ps.H




####### ---------- CREATE PLOT SHOWING THE ENERGY STORED IN THE ELECTRIC STORES AND THE POWER CHARGED AND DISCHARGED AT EACH TIME PERIOD ------- ########
electric_stores = Plots.plot(time_periods, electric_storage(ps.EnS[1]), label=ps.EnS[1], xlabel="Time Period", ylabel="Energy in Storage (GWh)");
if length(ps.EnS)>1
    for ens in ps.EnS[2:end]
        electric_stores = Plots.plot!(time_periods, electric_storage(ens), label=ens, xlabel="Time Period", ylabel="Energy in Storage (GWh)");
    end
end
Plots.savefig("Operational Results Plots/EnS_plots.pdf")

electric_stores_charging = Plots.plot(time_periods, electric_pumps(ps.ES[1]), label=ps.ES[1], xlabel="Time Period", ylabel="Power Charged into Storage (GW)");
if length(ps.ES)>1
    for es in ps.ES[2:end]
        electric_stores_charging = Plots.plot!(time_periods, electric_pumps(es), label=es, xlabel="Time Period", ylabel="Power Charged into Storage (GW)");
    end
end
Plots.savefig("Operational Results Plots/ES_plots.pdf")


####### ---------- CREATE PLOT SHOWING POWER FLOW THROUGH EACH HEAT LINE AT EACH TIME PERIOD ---------- #######

heatlines = Plots.plot(time_periods, heat_line(ps.HL[1]), label=ps.HL[1], xlabel="Time Period", ylabel="Power in Heat Line (GW)");
if length(ps.HL)>1
    for hl in ps.HL[2:end]
        heatlines = Plots.plot!(time_periods, heat_line(hl), label=hl, xlabel="Time Period", ylabel="Power in Heat Line (GW)");
    end
end
Plots.savefig("Operational Results Plots/HL_plots.pdf")




####### ----------- CREATE PLOT SHOWING INTERNAL TEMPERATURE OF HOUSE/HEAT STORE AT EACH TIME PERIOD ---------- ########

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


######   ------------   CREATE PLOT WITH TOTAL ENERGY STORED IN HOUSE/HEAT STORE AT EACH TIME PERIOD  ----------   ########

stored_energy_plots = Plots.plot(time_periods, house_energy(ps.HOUSES[1]), label = ps.HOUSES[1], xlabel="Time Period", ylabel="Energy Stored (GWh)");
if length(ps.HOUSES)>1
    for house in ps.HOUSES[2:end]
        stored_energy_plots = Plots.plot!(time_periods, house_energy(house), label=house, xlabel="Time Period", ylabel="Energy Stored (GWh)");
    end
end
for hns in ps.HnS
    stored_energy_plots = Plots.plot!(time_periods, heat_store_energy(hns), label=hns, xlabel="Time Period", ylabel="Energy Stored (GWh)");
end
Plots.savefig("Operational Results Plots/stored_energy_plots.pdf")



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




#####    ---------   CREATE PLOT WITH POWER OUTPUTS OF EACH GENERATOR FOR ALL TIME PERIODS    ------------    #######

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
 

for n in ps.NUC     ## Add Nuclear generators to everything plot
    generator_plot = Plots.plot!(time_periods, nuclear_generator(n), label=n, ylabel="Power Generated/Discharged (MW)")    
end


Plots.savefig("Operational Results Plots/generator_plot.pdf")

