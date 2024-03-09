function SP!(m::JuMP.Model,ps,pp)::JuMP.Model

    # */ -- Auxiliary info ----------------------------------------- /* #
    HS = H2S(ps.H,pp.st,pp.ln); # map of period (h) to slice
    sc = Dict(h => pp.ts[HS[h]]*pp.wg[HS[h]]/(pp.ln[HS[h]]*pp.ts[HS[h]]) for h in ps.H); # Auxiliary dictionary (optional) that contains the scaling factors given a period h

    # */ -- VARIABLES ---------------------------------------------- /* #
    # - Thermal generation variables
    @variable(m, pG[ps.T,ps.H] >= 0);      #power output of thermal technology for each period

    # - Renewable generation variables
    @variable(m, pR[ps.R,ps.H] >= 0);      #power output of renewable technology for each period

    # - Nuclear generation variables
    @variable(m, pNuc[ps.NUC,ps.H] >= 0);

    # - Electric store variables
    @variable(m, pS_in[ps.ES,ps.H]>= 0);      #power charging of storage
    @variable(m, pS_out[ps.ES,ps.H]>= 0);      #power discharging of storage
    @variable(m, qB[ps.EnS,ps.H]>=0);          #energy in store

    # - Heat store variables
    @variable(m,  pH_heat[ps.HP,ps.H]>=0);         #POWER INTO OF HEAT STORE (Increasing) [GW]
    @variable(m,  qH[ps.HnS,ps.H]>= 0);        #ENERGY IN THE HEAT STORE [GWh]
    @variable(m,  qHS[ps.HOUSES,ps.H]>= 0);      #ENERGY IN THE HOUSES
    @variable(m,  tInt[ps.HOUSES,ps.H]);        #TEMP IN THE HOUSES [GWh]
    @variable(m, qHShedP[union(ps.HOUSES,ps.HnS),ps.H]>=0);      #Power flexibility for heatstores positive (to keep it always feasible)
    @variable(m, qHShedN[union(ps.HOUSES,ps.HnS),ps.H]>=0);       #Power flexibility for heatstores neg (to keep it always feasible)

    # - Line / Bus variables
    @variable(m, pL[ps.L,ps.H]);           #power flow into line l [GW]
    @variable(m, pHL[ps.HL,ps.H]>=0);         #power flow into heat line hl
    @variable(m, delta[ps.B,ps.H]);        #voltage angle [rad]


    # - Shed/Curtail variables
    @variable(m, gShed[ps.R,ps.H]>= 0);    #renewable shed at period t
    @variable(m, dShed[ps.D,ps.H]>= 0);    #demand shed at period t

    
    #### Investment set variables
    @variable(m,pub[ps.T]);             # upper bounds for thermal
    @variable(m,rub[ps.R]);             # renewable scaling
    @variable(m,sub[ps.ES]);            # upper bounds for stores (electric) - power
    @variable(m,eub[ps.EnS]);           # upper bounds for stores (electric) - energy
    @variable(m,hub[ps.HP]);            # upper bounds for heat pumps - power
    @variable(m,QMass[ps.HnS]);        # QMass of TES
    @variable(m,lub[ps.L]);             # thermal limits on lines
    @variable(m,hlub[ps.HL]>= 0);           # thermal upper bound on Heat Lines
    @variable(m,nub[ps.NUC]);           # upper bounds for nuclear power output

    @variable(m,co2l);                  #CO2 level (passed from MP)
    @variable(m,sdm);                  #Demand scaling (passed from MP)

    @variable(m,co2gen);                  #CO2 cost (passed from MP)
    
    # */ -- EXPRESSIONS ---------------------------------------------- /* #
    # - Cost gathering expression -

    ##      VARIABLE COSTS (MOSTLY OPERATION) - PER PERIOD = VARIABLE O&M + Fuel Costs without CO2 Tax [T£/period]
    @expression(m,Var_Cost[h in ps.H], sum(pG[g,h]*(pp.cvg[g]+pp.fc[g]/pp.eta[g]) for g in ps.T) + 
                                       sum(pR[r,h]*pp.cvg[r] for r in ps.R) +
                                       sum(pS_in[s,h]*pp.cvg[s] for s in ps.ES) +
                                       sum(pH_heat[s,h]*pp.cvg[s] for s in ps.HP) +
                                       sum(pNuc[n,h]*(pp.cvg[n]+pp.fc[n]/pp.eta[n]) for n in ps.NUC)
                                       );


    ##      SHED / CURTAIL COSTS  - PER PERIOD = Load shed + Generation curtailment [T£/period]
    @expression(m,Sh_Cost[h in ps.H], sum(gShed[r,h]*(pp.ccurt[r]) for r in ps.R) + 
                                       sum(dShed[d,h]*pp.cs[d] for d in ps.D) +
                                       sum((qHShedP[hs,h]+qHShedN[hs,h])*pp.ccurt[hs] for hs in union(ps.HOUSES,ps.HnS))
                                       );

    @expression(m, c0, sum(sc[h]*(Var_Cost[h]+Sh_Cost[h]) for h in ps.H));


    # - Proper system constraints -
    # */ ------------ conventional generator limits + ramping limits ------- /* #
    @constraint(m, c01[g=ps.T, h=ps.H],  pG[g,h] <= pub[g]);    # maximum capacity limitation on conventional generation

    # */ --------------- electric and thermal storage technology limits ------------ /* #
    @constraint(m, c02[b=ps.EnS, h=ps.H], qB[b,h]<=eub[b]);   # energy upper limit in electric stores
    @constraint(m, c03[b=ps.ES, h=ps.H], pS_in[b,h]<=sub[b]);   # charging power limits
    @constraint(m, c04[b=ps.ES, h=ps.H], pS_out[b,h]<=sub[b]);   # discharging power limits for electric store

    @constraint(m, c08[b=ps.EnS, h=ps.H], qB[b,next(h,pp.st,pp.ln)] - qB[b,h] == pp.ts[HS[h]]*(pp.eta[pp.lp[b]]*pS_in[pp.lp[b],h] - pS_out[pp.lp[b],h])); # storage energy inventory for electric

    @constraint(m,c09[r=ps.R, h=ps.H], pR[r,h]==rub[r]/pp.hc[r]*pp.pr[r][h]);  # generation limit on renewables

    # */ ---------- Nuclear Limits ----------- /* #
    @constraint(m, c26[n=ps.NUC, h=ps.H], pNuc[n,h] <= nub[n]);


    # */ ------------- heat store & house constraints ---------------------- /* #
    @constraint(m, c13[hs=ps.HOUSES, h=ps.H], tInt[hs,h] == qHS[hs,h]/pp.qmass[hs]);
    @constraint(m, c14[hs=ps.HOUSES, h=ps.H], tInt[hs,h] >= pp.tmin[hs]);
    @constraint(m, c15[hs=ps.HOUSES, h=ps.H], tInt[hs,h] <= pp.tmax[hs]);

    @constraint(m, c05[b=ps.HnS,h=ps.H], qH[b,h] <= QMass[b]*(pp.tmax[b]));   # Upper bound for energy in TES
    @constraint(m, c10[b=ps.HnS,h=ps.H], qH[b,h] >= QMass[b]*(pp.tmin[b]));   # Lower bound for energy in TES
    @constraint(m, c06[b=ps.HP, h=ps.H], pH_heat[b,h]<=hub[b]);   # charging power limits

    # # */ ------ Heat Node Balance Constraint ------ /* #
    @constraint(m, c16[hn=ps.HN, h=ps.H], 
        sum(pH_heat[hp,h]*pp.eta[hp] for hp in ps.HP if pp.hDeliv[hp]==hn) == sum(pHL[hl,h] for hl in ps.HL if pp.hfm[hl]==hn));

    # */ ------ Heat Store & House Energy balance Constraint ------ /* #
    # Heat balance for Heat Stores
    @expression(m, PLoss[hns in ps.HnS], pp.lambda[hns]*QMass[hns]);
    @constraint(m, c17[hns=ps.HnS,h=ps.H], qH[hns,next(h,pp.st,pp.ln)] == qH[hns,h] + qHShedP[hns,h] - qHShedN[hns,h] +  pp.ts[HS[h]]*(pp.text[hns][h]*PLoss[hns] - pp.lambda[hns]*qH[hns,h]
                                            + (sum(pHL[hl,h] for hl in ps.HL if pp.hto[hl]==hns)) - sum(pHL[hl,h] for hl in ps.HL if pp.hfm[hl]==hns)));
    
    # Heat balance for Houses
    @constraint(m, c18[hs=ps.HOUSES,h=ps.H], qHS[hs,next(h,pp.st,pp.ln)]== qHS[hs,h] + qHShedP[hs,h] - qHShedN[hs,h] + pp.ts[HS[h]]*(pp.ploss[hs]*(pp.text[hs][h]-tInt[hs,h])
                                            + sum(pHL[hl,h] for hl in ps.HL if pp.hto[hl]==hs) - sum(pHL[hl,h] for hl in ps.HL if pp.hfm[hl]==hs)));


    # */ -------------------- KCL ----------------------- /* #
    @constraint(m, c19[b=ps.B,h=ps.H], sum(pG[g,h] for g in ps.T_B[b]) + sum(pR[r,h] for r in ps.R_B[b])  + sum(pNuc[n,h] for n in ps.NUC_B[b]) +
                sum(pL[l,h] for l in ps.L_to[b]) - sum(pL[l,h] for l in ps.L_fm[b]) + sum(dShed[d,h] for d in ps.D_B[b]) + sum(pS_out[b,h] for b in ps.ES_B[b])
                == sdm*sum(pp.pd[d][h] for d in ps.D_B[b]) + sum(gShed[r,h] for r in ps.R_B[b]) + sum(pS_in[s,h] for s in ps.ES_B[b]) + sum(pH_heat[hp,h] for hp in ps.HP_B[b]));

    # */ -------------------- KVL ----------------------- /* #
    @constraint(m,c20[l=ps.L,h=ps.H],pL[l,h]==1/pp.imp[l]*(delta[pp.to[l],h]-delta[pp.fm[l],h]));

 
    # */ -------------- LINE LIMITS --------------------- /* #
    @constraint(m, c21[l=ps.L,h=ps.H], pL[l,h] <=  lub[l]);
    @constraint(m, c22[l=ps.L,h=ps.H], pL[l,h] >= -lub[l]);

    # */ -------------- HEAT LINE LIMITS ---------------- /* #
    @constraint(m, c23[hl=ps.HL,h=ps.H], pHL[hl,h] <=  hlub[hl]);   

#	@constraint(m, c48[h=ps.H[1417:15313]], pHL["HL3",h] == 0);
#	@constraint(m, c48[h=ps.H[1:1417]], pHL["HL3",h] == 0);    # force tank to hold heat during Jan & Feb
#	@constraint(m, c49[h=ps.H[15313:end]], pHL["HL3",h] == 0);   # force tank to hold heat during october-december

    # */ ----------------------- CO2 budget ------------------------ /* #
    @constraint(m, c30, co2gen==sum(sum(pp.eg[g]/pp.eta[g]*pG[g,h] for g in ps.T)*sc[h] for h in ps.H));  # compute cost dependent on CO₂ level
    @constraint(m, c25, co2gen<=co2l);


    # */ -- Objective function (will later be rewritten to include uncertainty in cost coeff, depending on node) ---------------------------------------------- /* #
    @objective(m,Min,c0);


    return m

end


function RMP!(m,ms,mp)

    # problem structure
    # min   f + ∑ᵢ πᵢ βᵢ)
    # s.t.  f = q₀ᵀx₀ + ∑ᵢ qᵢᵀxᵢ
    #      (x₀,x₁,..,xᵢ,..,xₙ) ∈ 𝕏
    #       βᵢ ∈ Θᵢ,  ∀i


    # */ -- variables ---------------------------------------------- /* #
    @variable(m, f) # investment-only cost
    @variable(m, pN[ms.M,ms.I0] >= 0.0001) # newly installed capacity in generators
    @variable(m, pC[ms.M,ms.I] >= 0) # cummulative capacity of generators

    @variable(m, beta[ms.I] >= 0) # operational cost of node i

    # */ -- auxiliary variables ------------------------------------- /* #
    @variable(m, dL[ms.I])  # Demand scaling
    @variable(m, dC[ms.I]) # CO2 level 
    @variable(m, cco[ms.I]) # CO2 cost

    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Min, f + sum(mp.kappa[i]*mp.prob[i]*beta[i] for i in ms.I) ) # investment plus operational cost (10⁶£)

    # */ -- constraints -------------------------------------------- /* #
    @constraint(m, c01i, f == 1/mp.yr*(sum(mp.prob[i0]*sum(mp.capex[g]*pN[g,i0]/mp.lf[g] for g in ms.M) for i0 in ms.I0) + sum(mp.kappa[i]*mp.prob[i]*sum(mp.cf[g]*pC[g,i] for g in ms.M) for i in ms.I))) # compute investment-only cost
    @constraint(m, c02i[g=ms.M,i=ms.I], pC[g,i] == mp.imin[g][i] + sum(pN[g,i])) # compute accumulated capacity at node i
    @constraint(m, c03i, pN["Nuclear","N1"] <= 0.0002)
    @constraint(m, c04i, pN["BioEnergy", "N1"] <= 0.0002)
    @constraint(m, c05i, pN["PumpedHydro_Store", "N1"] <= 0.002)


#	@constraint(m, c06i, pN["Onshore_Wind", "N1"] <= 30.00)
#	@constraint(m, c07i, pN["Wind", "N1"] <= 30.00)
#	@constraint(m, c07i, pN["Solar", "N1"] >= 50.00)

    for i in ms.I
        fix(cco[i],mp.cco2[i];force=true) # CO2 Cost (Exogenus)
        fix(dC[i],mp.co2lim[i];force=true)
        fix(dL[i],mp.dsc[i];force=true)
    end

    return m
end



function LMP!(m,ms,mp)

 
    # problem structure
    # min   f + ∑ᵢ πᵢ βᵢ)
    # s.t.  f = q₀ᵀx₀ + ∑ᵢ qᵢᵀxᵢ
    #      (x₀,x₁,..,xᵢ,..,xₙ) ∈ 𝕏
    #       βᵢ ∈ Θᵢ,  ∀i


    # */ -- variables ---------------------------------------------- /* #
    @variable(m, f) # investment-only cost
    @variable(m, pN[ms.M,ms.I0] >= 0) # newly installed capacity in generators
    @variable(m, pC[ms.M,ms.I] >= 0) # cummulative capacity of generators

    @variable(m, beta[ms.I] >= 0) # operational cost of node i

    # */ -- auxiliary variables ------------------------------------- /* #
    @variable(m, dL[ms.I])  # Demand scaling
    @variable(m, dC[ms.I]) # CO2 level 
    @variable(m, cco[ms.I]) # CO2 cost

    @variable(m, lf, container=Array)   # level method factor


    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Min, sum( (pN[g,i] - 0)^2  for g in ms.M for i in ms.I0)) # investment plus operational cost (10⁶£)

    # */ -- constraints -------------------------------------------- /* #
    @constraint(m, c00i, f + sum(mp.kappa[i]*mp.prob[i]*beta[i] for i in ms.I) <= lf )
    @constraint(m, c01i, f == 1/mp.yr*(sum(mp.prob[i0]*sum(mp.capex[g]*pN[g,i0]/mp.lf[g] for g in ms.M) for i0 in ms.I0) + sum(mp.kappa[i]*mp.prob[i]*sum(mp.cf[g]*pC[g,i] for g in ms.M) for i in ms.I))) # compute investment-only cost
    @constraint(m, c02i[g=ms.M,i=ms.I], pC[g,i] == mp.imin[g][i] + sum(pN[g,i])) # compute accumulated capacity at node i
   

    for i in ms.I
        fix(cco[i],mp.cco2[i];force=true) # CO2 Cost (Exogenus)
        fix(dC[i],mp.co2lim[i];force=true)
        fix(dL[i],mp.dsc[i];force=true)
    end

    return m
end


function Undecomposed!(m,ms,mp,lT)

    #INVESTMENT PART
    # */ -- variables ---------------------------------------------- /* #
    @variable(m, f) # investment-only cost
    @variable(m, pN[ms.M,ms.I0] >= 0) # newly installed capacity in generators
    @variable(m, pC[ms.M,ms.I] >= 0) # cummulative capacity of generators

    @variable(m, beta[ms.I] >= 0) # operational cost of node i

    # */ -- obj function ------------------------------------------- /* #
    @objective(m, Min, f + sum(mp.kappa[i]*mp.prob[i]*beta[i] for i in ms.I) ) # investment plus operational cost (10⁶£)

    # */ -- constraints -------------------------------------------- /* #
    @constraint(m, c01i, f >= 1/mp.yr*(sum(mp.prob[i0]*sum(mp.capex[g]*pN[g,i0]/mp.lf[g] for g in ms.M) for i0 in ms.I0) + sum(mp.kappa[i]*mp.prob[i]*sum(mp.cf[g]*pC[g,i] for g in ms.M) for i in ms.I))) # compute investment-only cost
    
    #Green Field capacity
    @constraint(m, c02i[g=ms.M,i=ms.I], pC[g,i] == mp.imin[g][i] + sum(pN[g,i])) # compute accumulated capacity at node i
    
    #Cumulative capacity on tree
    #@constraint(m, c02i[g=ms.G,i=ms.I], pC[g,i] == mp.xh[g][i] + sum(pN[g,i0] for i0 in ms.map[i])) # compute accumulated capacity at node i

    


    #OPERATIONAL PART

  
    # */ -- Auxiliary info ----------------------------------------- /* #
    HS = H2S(ps.H,pp.st,pp.ln); # map of period (h) to slice
    sc = Dict(h => pp.ts[HS[h]]*pp.wg[HS[h]]/(pp.ln[HS[h]]*pp.ts[HS[h]]) for h in ps.H); # Auxiliary dictionary (optional) that contains the scaling factors given a period h

    # */ -- VARIABLES ---------------------------------------------- /* #
    # - Thermal generation variables
    @variable(m, pG[ms.I,ps.T,ps.H] >= 0)      #power output of thermal technology for each period

    # - Renewable generation variables
    @variable(m, pR[ms.I,ps.R,ps.H] >= 0)      #power output of renewable technology for each period

    # - Electric store variables
    @variable(m, pS_in[ms.I,ps.ES,ps.H]>= 0);      #power charging of storage
    @variable(m, pS_out[ms.I,ps.ES,ps.H]>= 0);      #power discharging of storage
    @variable(m, qB[ms.I,ps.EnS,ps.H]>=0);          #energy in store

    # - Heat store variables
    @variable(m,  pH_heat[ms.I,ps.HP,ps.H]>=0);         #POWER INTO OF HEAT STORE (Increasing) [GW]
    @variable(m,  qH[ms.I,ps.HnS,ps.H]>= 0);        #ENERGY IN THE HEAT STORE [GWh]
    @variable(m,  qHS[ms.I,ps.HOUSES,ps.H]>= 0);      #ENERGY IN THE HOUSES
    @variable(m,  tInt[ms.I,ps.HOUSES,ps.H]);          #TEMP IN THE HEAT STORE [GWh]
    @variable(m, qHShedP[ms.I,union(ps.HOUSES,ps.HnS),ps.H]>=0);      #Power flexibility for heatstores positive (to keep it always feasible)
    @variable(m, qHShedN[ms.I,union(ps.HOUSES,ps.HnS),ps.H]>=0);       #Power flexibility for heatstores neg (to keep it always feasible)


    # - Line / Bus variables
    @variable(m, pL[ms.I,ps.L,ps.H]);           #power flow into line l [GW]
    @variable(m, delta[ms.I,ps.B,ps.H]);        #voltage angle [rad]
    @variable(m, pHL[ms.I,ps.HL,ps.H]>=0);         #power flow into heat line hl

    # - Shed/Curtail variables
    @variable(m, gShed[ms.I,ps.R,ps.H]>= 0);    #renewable shed at period t
    @variable(m, dShed[ms.I,ps.D,ps.H]>= 0);    #demand shed at period t


    #### Investment set variables
    @variable(m,pub[ms.I,ps.T]);             #upper bounds for thermal
    @variable(m,rub[ms.I,ps.R]);             #renewable scaling
    @variable(m,sub[ms.I,ps.ES]);            #upper bounds for stores (electric) - power
    @variable(m,eub[ms.I,ps.EnS]);            #upper bounds for stores (electric) - energy
    @variable(m,hub[ms.I,ps.HP]);            #upper bounds for heat pumps - power
    @variable(m, QMass[ms.I,ps.HnS]);           #QMass of TES
    @variable(m,lub[ms.I,ps.L]);             #thermal limits on lines
    @variable(m,hlub[ms.I,ps.HL]>= 0);           #thermal upper bound on Heat Lines


    @variable(m,co2l[ms.I]);                  #CO2 level (passed from MP)
    @variable(m,sdm[ms.I]);                  #Demand scaling (passed from MP)

    @variable(m,co2c[ms.I]);                  #CO2 cost (passed from MP)

    # */ -- EXPRESSIONS ---------------------------------------------- /* #
    # - Cost gathering expression -

    ##      VARIABLE COSTS (MOSTLY OPERATION) - PER PERIOD = VARIABLE O&M + Fuel Costs without CO2 Tax [T£/period]
    @expression(m,Var_Cost[i in ms.I, h in ps.H], sum(pG[i,g,h]*(pp.cvg[g]+pp.fc[g]/pp.eta[g]) for g in ps.T) + 
                                                sum(pR[i,r,h]*pp.cvg[r] for r in ps.R) +
                                                sum(pS_in[i,s,h]*pp.cvg[s] for s in ps.ES) +
                                                sum(pH_heat[i,s,h]*pp.cvg[s] for s in ps.HP)
                                                );


    ##      SHED / CURTAIL COSTS  - PER PERIOD = Load shed + Generation curtailment [T£/period]
    @expression(m,Sh_Cost[i in ms.I, h in ps.H], sum(gShed[i,r,h]*(pp.ccurt[r]) for r in ps.R) + 
                                                sum(dShed[i,d,h]*pp.cs[d] for d in ps.D) +
                                                sum((qHShedP[i,hs,h]+qHShedN[i,hs,h])*pp.ccurt[hs] for hs in union(ps.HOUSES,ps.HnS))
                                                );    

    @expression(m, c0[i in ms.I],sum(sc[h]*(Var_Cost[i,h]+Sh_Cost[i,h]) for h in ps.H))

    # */ -- Objective function (will later be rewritten) ---------------------------------------------- /* #

    # - Proper system constraints -
    # */ ------------ conventional generator limits+ramping limits ------- /* #
    @constraint(m, c01[i=ms.I,g=ps.T, h=ps.H],  pG[i,g,h] <= pub[i,g]);    # maximum capacity limitation on conventional generation

    # */ --------------- electric and thermal storage technology limits ------------ /* #
    @constraint(m, c02[i=ms.I, b=ps.EnS, h=ps.H], qB[i,b,h]<=eub[i,b]); # energy upper limit in electric stores
    @constraint(m, c03[i=ms.I, b=ps.ES, h=ps.H], pS_in[i,b,h]<=sub[i,b]);   # charging power limits
    @constraint(m, c04[i=ms.I, b=ps.ES, h=ps.H], pS_out[i,b,h]<=sub[i,b]);   # discharging power limits for electric store

    @constraint(m, c05[i=ms.I,b=ps.HP, h=ps.H], pH_heat[i,b,h]<=hub[i,b]);   # charging power limits

    @constraint(m, c06[i=ms.I,b=ps.EnS, h=ps.H], qB[i,b,next(h,pp.st,pp.ln)] - qB[i,b,h] == pp.ts[HS[h]]*(pp.eta[pp.lp[b]]*pS_in[i,pp.lp[b],h] - pS_out[i,pp.lp[b],h])) # storage energy inventory for electric

    @constraint(m,c07[i=ms.I,r=ps.R, h=ps.H],pR[i,r,h]==rub[i,r]/pp.hc[r]*pp.pr[r][h]);


    # */ ------------- heat store & house constraints ---------------------- /* #
    @constraint(m, c08[i=ms.I, hs=ps.HOUSES, h=ps.H], tInt[i,hs,h] == qHS[i,hs,h]/pp.qmass[hs]);
    @constraint(m, c09[i=ms.I, hs=ps.HOUSES, h=ps.H], tInt[i,hs,h] >= pp.tmin[hs]);
    @constraint(m, c10[i=ms.I, hs=ps.HOUSES, h=ps.H], tInt[i,hs,h] <= pp.tmax[hs]);

    @constraint(m, c11[i=ms.I, b=ps.HnS,h=ps.H], qH[i,b,h] <= QMass[i,b]*(pp.tmax[b]));   # Upper bound for energy in TES
    @constraint(m, c12[i=ms.I, b=ps.HnS,h=ps.H], qH[i,b,h] >= QMass[i,b]*(pp.tmin[b]));   # Lower bound for energy in TES
    @constraint(m, c13[i=ms.I, b=ps.HP, h=ps.H], pH_heat[i,b,h]<=hub[i,b]);   # charging power limits


    # # */ ------ Heat Node Balance Constraint ------ /* #
    @constraint(m, c14[i=ms.I, hn=ps.HN, h=ps.H], 
    sum(pH_heat[i,hp,h]*pp.eta[hp] for hp in ps.HP if pp.hDeliv[hp]==hn) == sum(pHL[i,hl,h] for hl in ps.HL if pp.hfm[hl]==hn));

    # */ ------ Heat Store & House Energy balance Constraint ------ /* #
    # Heat balance for Heat Stores
    @expression(m, PLoss[i in ms.I, hns in ps.HnS], pp.lambda[hns]*QMass[i,hns]);
    @constraint(m, c15[i=ms.I,hns=ps.HnS], qH[i,hns,1] == QMass[i,hns]*pp.tmin[hns]);  # Set initial temperature of heat store to min temp
    @constraint(m, c16[i=ms.I,hns=ps.HnS,h=ps.H], qH[i,hns,next(h,pp.st,pp.ln)] == qH[i,hns,h] + qHShedP[i,hns,h] - qHShedN[i,hns,h] +  pp.ts[HS[h]]*(pp.text[hns][h]*PLoss[i,hns] - pp.lambda[hns]*qH[i,hns,h]
                                        + sum(pHL[i,hl,h] for hl in ps.HL if pp.hto[hl]==hns) - sum(pHL[i,hl,h] for hl in ps.HL if pp.hfm[hl]==hns)));

    # Heat balance for Houses
    @constraint(m, c17[i=ms.I,hs=ps.HOUSES,h=ps.H], qHS[i,hs,next(h,pp.st,pp.ln)]== qHS[i,hs,h] + qHShedP[i,hs,h] - qHShedN[i,hs,h] + pp.ts[HS[h]]*(pp.ploss[hs]*(pp.text[hs][h]-tInt[i,hs,h])
                                        + sum(pHL[i,hl,h] for hl in ps.HL if pp.hto[hl]==hs) - sum(pHL[i,hl,h] for hl in ps.HL if pp.hfm[hl]==hs)));



    # */ -------------------- KCL ----------------------- /* #
    @constraint(m, c18[i=ms.I,b=ps.B,h=ps.H], sum(pG[i,g,h] for g in ps.T_B[b]) + sum(pR[i,r,h] for r in ps.R_B[b])  +
    sum(pL[i,l,h] for l in ps.L_to[b]) - sum(pL[i,l,h] for l in ps.L_fm[b]) + sum(dShed[i,d,h] for d in ps.D_B[b]) + sum(pS_out[i,b,h] for b in ps.ES_B[b])
     == sdm[i]*sum(pp.pd[d][h] for d in ps.D_B[b]) + sum(gShed[i,r,h] for r in ps.R_B[b]) + sum(pS_in[i,s,h] for s in ps.ES_B[b]) + sum(pH_heat[i,hp,h] for hp in ps.HP_B[b]))

    # */ -------------------- KVL ----------------------- /* #
    @constraint(m,c19[i=ms.I,l=ps.L,h=ps.H],pL[i,l,h]==1/pp.imp[l]*(delta[i,pp.to[l],h]-delta[i,pp.fm[l],h]));
    
    # */ -------------- LINE LIMITS --------------------- /* #
    @constraint(m, c20[i=ms.I,l=ps.L,h=ps.H], pL[i,l,h] <=  lub[i,l]);
    @constraint(m, c21[i=ms.I,l=ps.L,h=ps.H], pL[i,l,h] >= -lub[i,l]);

    # */ -------------- HEAT LINE LIMITS ---------------- /* #
    @constraint(m, c22[i=ms.I,hl=ps.HL,h=ps.H], pHL[i,hl,h] <=  hlub[i,hl]);

    # */ ----------------------- CO2 budget ------------------------ /* #
    @expression(m, cco[i=ms.I], sum(sum(pp.eg[g]/pp.eta[g]*pG[i,g,h] for g in ps.T)*sc[h] for h in ps.H))      # compute cost dependent on CO₂ level
    @constraint(m, c23[i=ms.I], cco[i]<=co2l[i]);


    #LINKING PART
    #Objective linking
    @constraint(m,csr_ObjLink[i in ms.I],beta[i]==(c0[i]+mp.cco2[i]*cco[i]));

    #### Investment set variables
    op2inv=Dict(lT[i][5]=>lT[i][8] for i in eachindex(lT))
    @constraint(m,csr_pub[i in ms.I,g in ps.T], pC[op2inv[g],i]==pub[i,g]);
    @constraint(m,csr_sub[i in ms.I,s in ps.ES], pC[op2inv[s],i]==sub[i,s]);
    @constraint(m,csr_eub[i in ms.I,s in ps.EnS], pC[op2inv[s],i]==eub[i,s]);
    @constraint(m,csr_hub[i in ms.I,s in ps.HP], pC[op2inv[s],i]==hub[i,s]);
    @constraint(m,csr_lub[i in ms.I,l in ps.L], pC[op2inv[l],i]==lub[i,l]);
    @constraint(m,csr_rub[i in ms.I,r in ps.R], pC[op2inv[r],i]==rub[i,r]);
    @constraint(m,csr_QMass[i in ms.I,hns in ps.HnS], pC[op2inv[hns],i]==QMass[i,hns]);
    @constraint(m,csr_hlub[i in ms.I,hl in ps.HL], pC[op2inv[hl],i]==hlub[i,hl]);
    @constraint(m,csr_co2l[i in ms.I], co2l[i]==mp.co2lim[i]);
    @constraint(m,csr_dsc[i in ms.I], sdm[i]==-1.0*mp.dsc[i]);

    return m
end


