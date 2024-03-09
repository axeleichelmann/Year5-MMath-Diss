#Generate structures

print("Generating B,S structures... "); a=time();

if alg.al==2
    N=fillEmptyStructure("N2");
elseif alg.al==1
    N=fillEmptyStructure("N1");
end

B,S= generate_structures(path,alg.al,0.,alg.w,alg.J,alg.gap,alg.emb,alg.stab,alg.gammas,alg.ad,alg.ac,alg.pg,ms,mp,ps,pp,unc,m2,m,i2n;mLP=m3);

println("done in $(round(time()-a,digits=1))s ");

#Set special point
setSP!(B,lT);  ## Needed the B structure and linking table lT
