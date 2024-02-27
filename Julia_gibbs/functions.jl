using Pkg;

using Statistics, Random, NamedArrays, DataStructures, CSV, DataFrames, GLM, StatsBase, Distributions;

isControl(sampleID::String) = (sampleID in ctrlID)
isControl(sampleID::String3) = (sampleID in ctrlID)
isNotControl(sampleID::String) = !(sampleID in ctrlID)
isNotControl(sampleID::String3) = !(sampleID in ctrlID);

cat!(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)

function mySaver(X::Dict; fileRoot="") 
    for key in keys(X)
        df = DataFrame(X[key], names(X[key], 2))
        CSV.write(fileRoot*string(key)*".csv", df)
    end
end

function getData_mats(data; mitochan::String, chan::String,
    ctrlID::Vector{String} = String[], pts::Vector{String} = String[], ctrlOnly::Bool = false, getIndex::Bool = true) 
    if size(ctrlID,1)==0 & ("sbjType" in names(data))
        ctrlID = unique( data[data.sbjType.=="control", :sampleID] )
    elseif size(ctrlID,1)==0 & !("sbjType" in names(data))
        error("Must pass `ctrlID` or have a column in the data called `sbjType` which indicates which sample IDs are `control` or `patient`.")
    end

    sbj = sort(unique(data[:,:sampleID]) )
    if size(pts,1)==0 
        pts = sbj[ [!(x in ctrlID) for x in sbj] ]
    end

    isControl(sampleID) = (sampleID in ctrlID)
    isNotControl(sampleID) = !(sampleID in ctrlID)
    ctrlData = filter(:sampleID => isControl, data)
    
    xCtrl = [] 
    yCtrl = []
    indexCtrl = []
    ind = 1
    for ctrl in ctrlID
        crlData = filter(:sampleID => id->id==ctrl, ctrlData)
        n = sum( crlData[:,channel].==mitochan ) 
        xCtrl = [xCtrl; crlData[crlData[:,channel].==mitochan, "Value"]]
        yCtrl = [yCtrl; crlData[crlData[:,channel].==chan, "Value"]]
        indexCtrl = [indexCtrl; repeat([ind], n)]
        ind += 1
    end
    ctrlMat = hcat(xCtrl, yCtrl)
    if ctrlOnly 
        if getIndex
            return Dict([ ("ctrl",ctrlMat), ("index",indexCtrl)])
        else 
            return ctrlMat
        end
    else 
        ptsData = filter(:sampleID=>isNotControl, data)
        xPat = []
        yPat = []
        indexPat = []
        for pat in pts
            patData = filter(:sampleID => id->id==pat, ptsData)
            n = sum( patData[:,channel].==mitochan ) 
            xPat = [xPat; patData[patData[:,channel].==mitochan, "Value"]]
            yPat = [yPat; patData[patData[:,channel].==chan, "Value"]]
            indexPat = [indexPat; repeat([ind], n)]
            ind = ind + 1
        end
    end 
    ptsMat = hcat(xPat, yPat);
    
    if getIndex  
        return Dict([ ("ctrl",ctrlMat), ("ctrlIndex",indexCtrl), ("pts",ptsMat), ("patIndex",indexPat) ])
    else
        return Dict( [ ("ctrl",ctrlMat), ("pts",ptsMat) ] )
    end
end



function gibbs_sampler(data, mitochan, chan; warmup=20000, iter=1000, thin=1)

    iterTotal = (warmup + iter) * thin
    warmIter = warmup*thin

    allSbjs = unique( data[:,"sampleID"] )
    nSbjs = length( allSbjs )

    ctrlData = filter(:sbjType => type-> type=="control", data)
    ctrlIDs = unique( ctrlData[:,:sampleID] )
    nCtrl = length(ctrlIDs)
    
    ptsData = filer(:sbjType => type-> type=="pateint", data)
    ptsIDs = unique( ptsData[:,:sampleID] )
    nPts = length(ptsIDs)
    
    mlab = [string("m[", sampleID, "]") for sampleID in allSbjs]
    clab = [string("c[", sampleID, "]") for sampleID in allSbjs]
    ctrlMat = ctrlData[:,[mitochan, chan]]
    patMat = dataMats[:, [mitochan, chan]]
    indices = [unique(ctrlIndex); unique(patIndex) ]

    
    nFib = [ sum(data[:,"sampleID"] .== sampleID) for sampleID in allSbjs ]
    fibreCount = NamedArray(nFib, allSbjs)

    ######
    # Calculate hyper-parameter values
    ######
    
    grad = Vector{Float64}(undef, nCtrl)
    inter = Vector{Float64}(undef, nCtrl)
    prec = Vector{Float64}(undef, nCtrl)
    for ctrl_id in ctrlIDs
        crl_data = filter(:sampleID => id->id==ctrl_id, ctrlData)
        
        xCtrl = ctrlMat[ctrlIndex.==i, 1]
        yCtrl = ctrlMat[ctrlIndex.==i, 2]
        mod = lm(@formula(chan ~ mitochan), crl_data)
        grad[i] = coef(mod)[2]
        inter[i] = coef(mod)[1]
        prec[i] = 1 / dispersion(mod.model)^2
    end
    slope_mean = mean(grad)
    inter_mean = mean(inter)
    tau_mean = mean(prec)
    
    #######
    # Set hyper-parameter values
    ######
    nSyn = 1000
    synMin = minimum([ctrlMat[:,1]; xPat] )-2.0
    synMax = maximum([ctrlMat[:,1]; xPat] )+2.0
    xSyn = collect(range(synMin, synMax, length=nSyn))

    hyperParams = ["mean_mu_m", "prec_mu_m", "mean_mu_c", "prec_mu_c", 
     "shape_tau_m", "rate_tau_m", "shape_tau_c", "rate_tau_c",
     "shape_tau", "rate_tau", "alpha_pi", "beta_pi", "tau_def"]
    
    hyperTheta = NamedArray(zeros(13), hyperParams)
    
    hyperTheta["mean_mu_m"] = slope_mean
    hyperTheta["prec_mu_m"] = 1 / 0.25^2
    hyperTheta["mean_mu_c"] = inter_mean
    hyperTheta["prec_mu_c"] = 1 / 0.25^2
    
    tau_m_var = 50
    tau_m_mean = 50
    hyperTheta["shape_tau_m"] = tau_m_mean^2 / tau_m_var
    hyperTheta["rate_tau_m"] = tau_m_mean / tau_m_var 
    
    tau_c_var = 50
    tau_c_mean = 50
    hyperTheta["shape_tau_c"] = tau_c_mean^2 / tau_c_var
    hyperTheta["rate_tau_c"] = tau_c_mean / tau_c_var

    tau_var = 10
    hyperTheta["shape_tau"] = tau_mean^2 / tau_var
    hyperTheta["rate_tau"] = tau_mean / tau_var

    hyperTheta["alpha_pi"] = 1
    hyperTheta["beta_pi"] = 1
    hyperTheta["tau_def"] = 0.00001

    ######
    # data summaries - hard coded for three controls
    ######
    data_summs = Dict{String, Vector{Float64}}()
    data_summs["sum_x"] = Vector{Float64}(undef, nCtrl)
    data_summs["sum_y"] = Vector{Float64}(undef, nCtrl)
    data_summs["sum_xsq"] = Vector{Float64}(undef, nCtrl)
    data_summs["sum_ysq"] = Vector{Float64}(undef, nCtrl)
    data_summs["sum_xy"] = Vector{Float64}(undef, nCtrl)
    for i in 1:nCtrl
        data_summs["sum_x"][i] = sum(ctrlMat[ctrlIndex.==i,1])
        data_summs["sum_y"][i] = sum(ctrlMat[ctrlIndex.==i,2])
        data_summs["sum_xsq"][i] = sum(ctrlMat[ctrlIndex.==i,1].^2)
        data_summs["sum_ysq"][i] = sum(ctrlMat[ctrlIndex.==i,2].^2)
        data_summs["sum_xy"][i] = sum(ctrlMat[ctrlIndex.==i,1].*ctrlMat[ctrlIndex.==i,2])
    end
    
    # data summaries for patient - sums must be calculated during inference
    # based on current classification
    data_summs["xsq_pat"] = xPat.*xPat
    data_summs["ysq_pat"] = yPat.*yPat
    data_summs["xy_pat"] = xPat.*yPat

    ######
    # PRIOR DRAWS
    ######
    
    # draw from the prior distributions
    paramNames = reduce(vcat, (["tau_c", "tau_m", "mu_c", "mu_m"], clab, mlab, ["tau_norm", "probdiff", "m_pred", "c_pred"]))
    nParam = length(paramNames)
    prior = NamedArray(zeros(iter, nParam))
    
    setnames!(prior, paramNames, 2)
    
    prior[:,"tau_c"] = [ rand( Gamma( hyperTheta["shape_tau_c"], 1/hyperTheta["rate_tau_c"] ) )  for _ in 1:iter]
    prior[:,"tau_m"] = [ rand(Gamma(hyperTheta["shape_tau_m"], 1/hyperTheta["rate_tau_m"])) for _ in 1:iter]
    prior[:,"mu_c"] = [ rand(Normal(hyperTheta["mean_mu_c"], 1/sqrt(hyperTheta["prec_mu_c"]))) for _ in 1:iter]
    prior[:,"mu_m"] = [ rand(Normal(hyperTheta["mean_mu_m"], 1/sqrt(hyperTheta["prec_mu_m"]))) for _ in 1:iter]
    for i in 1:nSbj 
      prior[:, clab[i]] = rand.(Normal.(prior[:,"mu_c"], 1 ./sqrt.(prior[:,"tau_c"]))) 
      prior[:, mlab[i]] = rand.(Normal.(prior[:,"mu_m"], 1 ./sqrt.(prior[:,"tau_m"])))
    end

    prior[:,"tau_norm"] = [rand(Gamma(hyperTheta["shape_tau"], 1/hyperTheta["rate_tau"])) for _ in 1:iter]
    prior[:,"probdiff"] = [rand(Beta(hyperTheta["alpha_pi"], hyperTheta["beta_pi"])) for _ in 1:iter]

    prior[:, "c_pred"] = rand.(Normal.(prior[:,"mu_c"], 1 ./sqrt.(prior[:,"tau_c"]))) 
    prior[:, "m_pred"] = rand.(Normal.(prior[:,"mu_m"], 1 ./sqrt.(prior[:,"tau_m"]))) 

    
    ############
    ######
    # POSTERIOR DRAWS
    ######
    ############
    
    ######
    ### set-up
    ######
    # matrices to store output
    post = NamedArray(zeros(iter, nParam))
    setnames!(post, paramNames, 2)
    
    fibreLab = ["classif["*string(i)*"]" for i in 1:nPat]
    classifs_mat = NamedArray(zeros(iter, nPat))
    setnames!(classifs_mat, fibreLab, 2)
    
    # create vector of parameters and sample from prior to initialise chain
    theta = NamedArray(zeros(nParam), paramNames)
    ### Initial theta 
    # draw from (hyper-)priors for the first element in our Markov chain
    theta["tau_c"] = rand(Gamma(hyperTheta["shape_tau_c"], 1/hyperTheta["rate_tau_c"]) )
    theta["tau_m"] = rand(Gamma(hyperTheta["shape_tau_m"], 1/hyperTheta["rate_tau_m"]) )
    theta["mu_c"] = rand(Normal(hyperTheta["mean_mu_c"], 1/sqrt(hyperTheta["prec_mu_c"])))
    theta["mu_m"] = rand(Normal(hyperTheta["mean_mu_m"], 1/sqrt(hyperTheta["prec_mu_m"])))
    for i in 1:nSbj 
      theta[clab[i]] = rand(Normal(theta["mu_c"], 1/sqrt(theta["tau_c"]))) 
      theta[mlab[i]] = rand(Normal(theta["mu_m"], 1/sqrt(theta["tau_m"])))
    end
    theta["tau_norm"] = rand(Gamma(hyperTheta["shape_tau"], 1/hyperTheta["rate_tau"]))
    theta["probdiff"] = rand(Beta(hyperTheta["alpha_pi"], hyperTheta["beta_pi"]))
    
    ### Initial classification
    denDef = theta["probdiff"] .* pdf.(Normal.(theta[mlab[nSbj]].*xPat .+ theta[clab[nSbj]], 1/sqrt(hyperTheta["tau_def"]) ), yPat )
    denNorm = (1 - theta["probdiff"]) .* pdf.(Normal.(theta[mlab[nSbj]] .*xPat .+ theta[clab[nSbj]], 1/sqrt(theta["tau_norm"]) ), yPat )
    u = rand.(Uniform.(0.0, denDef .+denNorm))
    classifs = u .< denDef
    
    likeCtrl = .!classifs
    nlikeCtrl = classifs

    # population densities
    theta["m_pred"] = rand(Normal(theta["mu_m"], 1/sqrt(theta["tau_m"])))
    theta["c_pred"] = rand(Normal(theta["mu_c"], 1/sqrt(theta["tau_c"])))
    
    ######
    # MCMC loop
    ######
    outCount = 1
    for tt in 1:iterTotal 
        sum_likeCtrl = sum(likeCtrl)
        sum_nlikeCtrl = nPat - sum_likeCtrl

        # update tau_c - intercept precision
        theta["tau_c"] = rand(Gamma(hyperTheta["shape_tau_c"]+nSbj/2, 1/(hyperTheta["rate_tau_c"] + 0.5*sum( (theta["mu_c"] .- theta[clab]).^2)) ) )
        # update tau_m - slope precision
        theta["tau_m"] = rand(Gamma(hyperTheta["shape_tau_m"]+nSbj/2, 1/(hyperTheta["rate_tau_m"] + 0.5*sum((theta["mu_m"] .- theta[mlab]).^2)) ))

        # update mu_c - expected intercept
        muc_prec = hyperTheta["prec_mu_c"] + theta["tau_c"]*nSbj
        theta["mu_c"] = rand(Normal( (hyperTheta["mean_mu_c"]*hyperTheta["prec_mu_c"] + sum(theta[clab])*theta["tau_c"] )/muc_prec, 1/sqrt(muc_prec)))

        # update mu_m - expected slope
        mum_prec = hyperTheta["prec_mu_m"] + theta["tau_m"]*nSbj
        theta["mu_m"] = rand(Normal( (hyperTheta["mean_mu_m"]*hyperTheta["prec_mu_m"] + sum(theta[mlab])*theta["tau_m"] )/mum_prec, 1/sqrt(mum_prec)) )
        
        # update c_i - the intercept for the control subjects
        for i in 1:nCtrl 
            tau_cCrl = theta["tau_norm"]*nFib_ctrl[i]
            mu_cCrl = (theta["tau_norm"]*(data_summs["sum_y"][i] - theta[mlab[i]]*data_summs["sum_x"][i]) ) / tau_cCrl
            theta[clab[i]] = rand(Normal( (theta["mu_c"]*theta["tau_c"] + mu_cCrl*tau_cCrl) / (theta["tau_c"] + tau_cCrl), 1/sqrt(theta["tau_c"] + tau_cCrl) ))
        end

        # update c_nSbj - the intercept for the patient subject
        tau_cPat = theta["tau_norm"]*sum_likeCtrl + hyperTheta["tau_def"]*sum_nlikeCtrl
        mu_cPat = ( theta["tau_norm"]*sum(yPat[likeCtrl] .- theta[mlab[nSbj]] .*xPat[likeCtrl])
                  + hyperTheta["tau_def"]*sum(yPat[nlikeCtrl] .- theta[mlab[nSbj]] .*xPat[nlikeCtrl]) ) / tau_cPat
        
        theta[clab[nSbj]] = rand(Normal( (theta["mu_c"]*theta["tau_c"] + mu_cPat*tau_cPat)/(theta["tau_c"]+tau_cPat), 1/sqrt(theta["tau_c"] + tau_cPat)))
        
        # update m_i - the slopes for the control subjects
        for i in 1:nCtrl
            tau_mCrl = theta["tau_norm"] *data_summs["sum_xsq"][i]
            mu_mCrl = (theta["tau_norm"] *(data_summs["sum_xy"][i] - theta[clab[i]]*data_summs["sum_x"][i] ) ) / tau_mCrl
            theta[mlab[i]] = rand(Normal( (theta["mu_m"]*theta["tau_m"] + mu_mCrl*tau_mCrl)/(theta["tau_m"] + tau_mCrl), 1/sqrt(theta["tau_m"]+tau_mCrl) ))
        end
        # update mu_nSbj - the slope of the patient subject
        tau_mPat = theta["tau_norm"]*sum(data_summs["xsq_pat"][likeCtrl]) + hyperTheta["tau_def"]*sum(data_summs["xsq_pat"][nlikeCtrl])
        mu_mPat = ( theta["tau_norm"]*sum(data_summs["xy_pat"][likeCtrl] .- theta[clab[nSbj]] .*xPat[likeCtrl])
                  + hyperTheta["tau_def"]*sum(data_summs["xy_pat"][nlikeCtrl] .- theta[clab[nSbj]] .*xPat[nlikeCtrl]) ) / tau_mPat
        theta[mlab[nSbj]] = rand(Normal( (theta["mu_m"]*theta["tau_m"] + mu_mPat*tau_mPat)/(theta["tau_m"] + tau_mPat), 1/sqrt(theta["tau_m"] + tau_mPat) ) )

        # update tau - the error in the control and like-control patient fibres
        sq_diff = 0.0 # calculate the squared difference between the expected expression level and observed expression level
        for i in 1:nCtrl
            sq_diff += sum( (theta[mlab[i]] .*ctrlMat[ctrlIndex.==i,1] .+theta[clab[i]] .- ctrlMat[ctrlIndex.==i,2]) .^2 )
        end
        sq_diff += sum( (theta[mlab[nSbj]] .*xPat[likeCtrl] .+ theta[clab[nSbj]] .- yPat[likeCtrl]) .^2)
        # update tau - the model error for like-control patients
        theta["tau_norm"] = rand(Gamma( hyperTheta["shape_tau"] + 0.5*(sum(nFib_ctrl) + sum_likeCtrl), 1 /(hyperTheta["rate_tau"] + 0.5*sq_diff) ) )

        # update proportion of deficiency
        theta["probdiff"] = rand(Beta( hyperTheta["alpha_pi"] + sum_nlikeCtrl, hyperTheta["beta_pi"] + sum_likeCtrl ) )
        
        # densities for the classification
        denDef = theta["probdiff"] .* pdf.(Normal.(theta[mlab[nSbj]] .*xPat .+ theta[clab[nSbj]], 1/sqrt(hyperTheta["tau_def"]) ), yPat)
        denNorm = (1 - theta["probdiff"]) .* pdf.(Normal.(theta[mlab[nSbj]] .*xPat .+ theta[clab[nSbj]], 1/sqrt(theta["tau_norm"]) ), yPat)
        # draw nPat uniform random draws between 0 and dens_total
        u = rand.(Uniform.(0.0, (denDef .+ denNorm)))
        # if the random draw for fibre i is greater than the density of the like control model
        # fibre i is classified as deficient
        classifs = u .<denDef
        
        likeCtrl = .!classifs
        nlikeCtrl = classifs

        # population densities
        theta["m_pred"] = rand(Normal(theta["mu_m"], 1/sqrt(theta["tau_m"])))
        theta["c_pred"] = rand(Normal(theta["mu_c"], 1/sqrt(theta["tau_c"])))
    
        if (tt>warmIter) && mod(tt-warmIter,thin)==0 
            post[outCount, paramNames] = theta[paramNames]
            classifs_mat[outCount,:] = nlikeCtrl
            outCount += 1
        end
    end
    
    ######
    ### PREDICTIVE INTERVALS
    ######
    
    # only "care" about predictive interval of patient 
    p = (0.025, 0.5, 0.975)
    
    postpred = Matrix{Float64}(undef, nSyn, 4)
    m_post = post[:,mlab[nSbj]]
    c_post = post[:,clab[nSbj]]
    tau_post = post[:,"tau_norm"]
    
    priorpred = Matrix{Float64}(undef, nSyn, 4)
    m_prior = prior[:,mlab[nSbj]]
    c_prior = prior[:,clab[nSbj]]
    tau_prior = prior[:,"tau_norm"]

    
    for j in 1:nSyn
        postpred_x = rand.(Normal.(m_post .*xSyn[j] .+ c_post, 1 ./sqrt.(tau_post) ))
        priorpred_x = rand.(Normal.(m_prior .*xSyn[j] .+ c_prior, 1 ./sqrt.(tau_prior) ))
        postpred[j,:] = [xSyn[j]; collect(quantile(postpred_x, p))]
        priorpred[j,:] = [xSyn[j];  collect(quantile(priorpred_x, p))]
    end

    postpred_named = NamedArray( postpred )
    priorpred_named = NamedArray( priorpred )
    
    setnames!(postpred_named, ["mitochan", "lwr", "med", "upr"], 2)
    setnames!(priorpred_named, ["mitochan", "lwr", "med", "upr"], 2)
    
    output = Dict()
    output["POST"] = post
    output["PRIOR"] = prior
    output["POSTPRED"] = postpred_named
    output["PRIORPRED"] = priorpred_named
    output["CLASSIF"] = classifs_mat
    
    return output
end

function gibbs_sampler_bigHierarchy(data, mitochan, chan; warmup=20000, iter=1000, thin=1)

    iterTotal = (warmup + iter) * thin
    warmIter = warmup*thin
    
    ## get some data summaries ready to go
    data[:,"mitochan_sq"] = data[:,mitochan] .* data[:,mitochan]
    data[:,"chan_sq"] = data[:,chan] .* data[:,chan]
    data[:,"mitochan_chan"] = data[:,mitochan] .* data[:,chan] 

    
    allSbjs = unique( data[:,"sampleID"] )
    nSbjs = length( allSbjs )
    nSbjs_half = nSbjs/2
    
    ctrlData = filter(:sbjType => tt-> tt=="control", data)
    ctrlIDs = unique( ctrlData[:,:sampleID] )
    nCtrl = length(ctrlIDs)
    
    ptsData = filter(:sbjType => tt-> tt=="patient", data)
    ptsIDs = unique( ptsData[:,:sampleID] )
    nPts = length(ptsIDs)
    uID = ptsData[:,:uniqueID]
    sID = ptsData[:,:sampleID]
    
    mlab = [string("m[", sampleID, "]") for sampleID in allSbjs]
    clab = [string("c[", sampleID, "]") for sampleID in allSbjs]
    plab = [string("pi[", sampleID, "]") for sampleID in ptsIDs]
    
    nFib = [ sum(data[:,"sampleID"] .== sampleID) for sampleID in allSbjs ]
    fibreCount = NamedArray(nFib, allSbjs)
    nFib_ctrl = sum( fibreCount[ctrlIDs] )
    
    ######
    # Calculate hyper-parameter values
    ######

    
    grad = NamedArray(zeros(nCtrl), ctrlIDs)
    inter = NamedArray(zeros(nCtrl), ctrlIDs)
    prec = NamedArray(zeros(nCtrl), ctrlIDs)
    for ctrl_id in ctrlIDs
        crl_data = filter(:sampleID => id->id==ctrl_id, ctrlData)
        crl_data[:,:x] = crl_data[:, mitochan]
        crl_data[:,:y] = crl_data[:, chan]
        mod = lm(@formula(y ~ x), crl_data)
        grad[ctrl_id] = coef(mod)[2]
        inter[ctrl_id] = coef(mod)[1]
        prec[ctrl_id] = 1 / dispersion(mod.model)^2
    end
    slope_mean = mean(grad)
    inter_mean = mean(inter)
    tau_mean = mean(prec)
    
    #######
    # Set hyper-parameter values
    ######
    
    nSyn = 1000
    synMin = minimum( data[:, mitochan]) - 2.0
    synMax = maximum( data[:, chan] ) + 2.0
    xSyn = collect(range(synMin, synMax, length=nSyn))
    
    hyperParams = ["mean_mu_m", "prec_mu_m", "mean_mu_c", "prec_mu_c", 
     "shape_tau_m", "rate_tau_m", "shape_tau_c", "rate_tau_c",
     "shape_tau", "rate_tau", "alpha_pi", "beta_pi", "tau_def"]
    
    hyperTheta = NamedArray(zeros(13), hyperParams)
    
    hyperTheta["mean_mu_m"] = slope_mean
    hyperTheta["prec_mu_m"] = 1 / 0.25^2
    hyperTheta["mean_mu_c"] = inter_mean
    hyperTheta["prec_mu_c"] = 1 / 0.25^2
    
    tau_m_var = 50
    tau_m_mean = 50
    hyperTheta["shape_tau_m"] = tau_m_mean^2 / tau_m_var
    hyperTheta["rate_tau_m"] = tau_m_mean / tau_m_var 
    
    tau_c_var = 50
    tau_c_mean = 50
    hyperTheta["shape_tau_c"] = tau_c_mean^2 / tau_c_var
    hyperTheta["rate_tau_c"] = tau_c_mean / tau_c_var
    
    tau_var = 10
    hyperTheta["shape_tau"] = tau_mean^2 / tau_var
    hyperTheta["rate_tau"] = tau_mean / tau_var
    
    hyperTheta["alpha_pi"] = 1
    hyperTheta["beta_pi"] = 1
    hyperTheta["tau_def"] = 0.0001
    
    ######
    # data summaries - hard coded for three controls
    ######
    
    ctrlData_summs = Dict{String, NamedArray}()
    ctrlData_summs["sum_mitochan"] = NamedArray(zeros(nCtrl), ctrlIDs)
    ctrlData_summs["sum_chan"] = NamedArray(zeros(nCtrl), ctrlIDs)
    ctrlData_summs["sum_mitochan_sq"] = NamedArray(zeros(nCtrl), ctrlIDs)
    ctrlData_summs["sum_chan_sq"] = NamedArray(zeros(nCtrl), ctrlIDs)
    ctrlData_summs["sum_mitochan_chan"] = NamedArray(zeros(nCtrl), ctrlIDs)
    for ctrl_id in ctrlIDs
        crl_data = filter(:sampleID => id->id==ctrl_id, ctrlData)
        ctrlData_summs["sum_mitochan"][ctrl_id] = sum( crl_data[:,mitochan] )
        ctrlData_summs["sum_chan"][ctrl_id] = sum( crl_data[:,chan] )
        ctrlData_summs["sum_mitochan_sq"][ctrl_id] = sum( crl_data[:,:mitochan_sq] )
        ctrlData_summs["sum_chan_sq"][ctrl_id] = sum( crl_data[:, :chan_sq] )
        ctrlData_summs["sum_mitochan_chan"][ctrl_id] = sum( crl_data[:, :mitochan_chan] )
    end
        
    ######
    # PRIOR DRAWS
    ######
    
    # draw from the prior distributions
    paramNames = reduce(vcat, (["tau_c", "tau_m", "mu_c", "mu_m"], clab, mlab, plab, ["tau_norm", "m_pred", "c_pred"]))
    nParam = length(paramNames)
    prior = NamedArray(zeros(iter, nParam))
    
    setnames!(prior, paramNames, 2)
    
    prior[:,"tau_c"] = [ rand( Gamma( hyperTheta["shape_tau_c"], 1/hyperTheta["rate_tau_c"] ) )  for _ in 1:iter]
    prior[:,"tau_m"] = [ rand(Gamma( hyperTheta["shape_tau_m"], 1/hyperTheta["rate_tau_m"])) for _ in 1:iter]
    prior[:,"mu_c"] = [ rand(Normal(hyperTheta["mean_mu_c"], 1/sqrt(hyperTheta["prec_mu_c"]))) for _ in 1:iter]
    prior[:,"mu_m"] = [ rand(Normal(hyperTheta["mean_mu_m"], 1/sqrt(hyperTheta["prec_mu_m"]))) for _ in 1:iter]
    for i in 1:nSbjs
      prior[:, clab[i]] = rand.(Normal.(prior[:,"mu_c"], 1 ./sqrt.(prior[:,"tau_c"]))) 
      prior[:, mlab[i]] = rand.(Normal.(prior[:,"mu_m"], 1 ./sqrt.(prior[:,"tau_m"])))
    end
    for pil in plab
        prior[:,pil] = [rand(Beta(hyperTheta["alpha_pi"], hyperTheta["beta_pi"])) for _ in 1:iter]
    end
    
    prior[:,"tau_norm"] = [rand(Gamma(hyperTheta["shape_tau"], 1/hyperTheta["rate_tau"])) for _ in 1:iter]
    
    prior[:, "c_pred"] = rand.(Normal.(prior[:,"mu_c"], 1 ./sqrt.(prior[:,"tau_c"]))) 
    prior[:, "m_pred"] = rand.(Normal.(prior[:,"mu_m"], 1 ./sqrt.(prior[:,"tau_m"]))) ;

    
    ############
    ######
    # POSTERIOR DRAWS
    ######
    ############
    
    ######
    ### set-up
    ######
    
    # matrices to store output
    post = NamedArray(zeros(nParam, iter))
    setnames!(post, paramNames, 1)
    
    classifs_mat = NamedArray( zeros(nrow(ptsData), iter) )
    setnames!(classifs_mat, uID, 1)
    
    ptsData[:,:denDef] .= 0.0
    ptsData[:,:denNorm] .= 0.0
    ptsData[:,:classif] .= true
    
    # create vector of parameters and sample from prior to initialise chain
    theta = NamedArray(zeros(nParam), paramNames)
    ### Initial theta 
    # draw from (hyper-)priors for the first element in our Markov chain
    theta["tau_c"] = rand(Gamma(hyperTheta["shape_tau_c"], 1/hyperTheta["rate_tau_c"]) )
    theta["tau_m"] = rand(Gamma(hyperTheta["shape_tau_m"], 1/hyperTheta["rate_tau_m"]) )
    theta["mu_c"] = rand(Normal(hyperTheta["mean_mu_c"], 1/sqrt(hyperTheta["prec_mu_c"])))
    theta["mu_m"] = rand(Normal(hyperTheta["mean_mu_m"], 1/sqrt(hyperTheta["prec_mu_m"])))
    for i in 1:nSbjs 
      theta[clab[i]] = rand(Normal(theta["mu_c"], 1/sqrt(theta["tau_c"]))) 
      theta[mlab[i]] = rand(Normal(theta["mu_m"], 1/sqrt(theta["tau_m"])))
    end
    theta["tau_norm"] = rand(Gamma(hyperTheta["shape_tau"], 1/hyperTheta["rate_tau"]))
    for pil in plab
        theta[pil] = rand(Beta(hyperTheta["alpha_pi"], hyperTheta["beta_pi"]))
    end
    
    ### Initial classification
    for pat in ptsIDs
        pil = "pi["*pat*"]"
        ml = "m["*pat*"]"
        cl = "c["*pat*"]"
        pat_sID = sID.==pat
        ptsData[pat_sID,:denDef] .= theta[pil] .* pdf.( Normal.(theta[ml].*ptsData[pat_sID, mitochan] .+ theta[cl], 1/sqrt(hyperTheta["tau_def"])), ptsData[pat_sID, chan] )
        ptsData[pat_sID,:denNorm] .= (1 - theta[pil]) .* pdf.( Normal.(theta[ml].*ptsData[pat_sID, mitochan] .+ theta[cl], 1/sqrt(theta["tau_norm"]) ), ptsData[pat_sID, chan] )
    end
    
    u = rand.(Uniform.(0.0, ptsData[:,:denDef] .+ ptsData[:,:denNorm]))
    ptsData[:,:classif] .= (u .< ptsData[:,:denDef])
    
    likeCtrl = .!ptsData[:,:classif]
    nlikeCtrl =  ptsData[:,:classif]
    
    # population densities
    theta["m_pred"] = rand(Normal(theta["mu_m"], 1/sqrt(theta["tau_m"])))
    theta["c_pred"] = rand(Normal(theta["mu_c"], 1/sqrt(theta["tau_c"])))
    
    sum_likeCtrl = NamedArray( zeros(nPts), ptsIDs )
    sum_nlikeCtrl = NamedArray( zeros(nPts), ptsIDs )
    for pat in ptsIDs
        sum_nlikeCtrl[pat] = sum( ptsData[sID.==pat, :classif] )
        sum_likeCtrl[pat] = fibreCount[pat] - sum_nlikeCtrl[pat]
    end
    
    ######
    # MCMC loop
    ######
    
    for tt in 1:iterTotal 
        for pat in ptsIDs
            sum_nlikeCtrl[pat] = sum( ptsData[sID.==pat, "classif"] )
            sum_likeCtrl[pat] = fibreCount[pat] - sum_nlikeCtrl[pat]
        end
        
        # println(theta)
        # println("tau_m, tau_c")
        # update tau_c - intercept precision
        theta["tau_c"] = rand(Gamma(hyperTheta["shape_tau_c"]+nSbjs_half, 1/(hyperTheta["rate_tau_c"] + 0.5*sum( (theta["mu_c"] .- theta[clab]).^2)) ) )
        # update tau_m - slope precision
        theta["tau_m"] = rand(Gamma(hyperTheta["shape_tau_m"]+nSbjs_half, 1/(hyperTheta["rate_tau_m"] + 0.5*sum( (theta["mu_m"] .- theta[mlab]).^2)) ) )

        # println("mu_m, mu_c")
        # update mu_c - expected intercept
        muc_prec = hyperTheta["prec_mu_c"] + theta["tau_c"]*nSbjs
        theta["mu_c"] = rand(Normal( (hyperTheta["mean_mu_c"]*hyperTheta["prec_mu_c"] + sum(theta[clab])*theta["tau_c"] )/muc_prec, 1/sqrt(muc_prec)))
        # update mu_m - expected slope
        mum_prec = hyperTheta["prec_mu_m"] + theta["tau_m"]*nSbjs
        theta["mu_m"] = rand(Normal( (hyperTheta["mean_mu_m"]*hyperTheta["prec_mu_m"] + sum(theta[mlab])*theta["tau_m"] )/mum_prec, 1/sqrt(mum_prec)) )

        # println("c_ctrl")
        # update c_i - the intercept for the control subjects
        for id in ctrlIDs 
            cl = "c["*id*"]"
            ml = "m["*id*"]"
            tau_cCrl = theta["tau_norm"]*fibreCount[id]
            mu_cCrl = (theta["tau_norm"]*(ctrlData_summs["sum_chan"][id] - theta[ml]*ctrlData_summs["sum_mitochan"][id]) ) / tau_cCrl
            theta[cl] = rand(Normal( (theta["mu_c"]*theta["tau_c"] + mu_cCrl*tau_cCrl) / (theta["tau_c"] + tau_cCrl), 1/sqrt(theta["tau_c"] + tau_cCrl) ))
        end
        
        # println("c_pts")
        # update c_nSbj - the intercept for the patient subject
        for pat in ptsIDs
            cl = "c["*pat*"]"
            ml = "m["*pat*"]"
            pat_likeCtrl = (sID.==pat) .& likeCtrl
            pat_nlikeCtrl = (sID.==pat) .& nlikeCtrl
        
            tau_cPat = theta["tau_norm"]*sum_likeCtrl[pat] + hyperTheta["tau_def"]*sum_nlikeCtrl[pat]
            mu_cPat = ( theta["tau_norm"]*sum(ptsData[pat_likeCtrl, chan] .- theta[ml] .*ptsData[pat_likeCtrl, mitochan] )
                      + hyperTheta["tau_def"]*sum( ptsData[pat_nlikeCtrl, chan] .- theta[ml] .*ptsData[pat_nlikeCtrl, mitochan]) ) / tau_cPat
            
            theta[cl] = rand(Normal( (theta["mu_c"]*theta["tau_c"] + mu_cPat*tau_cPat)/(theta["tau_c"]+tau_cPat), 1/sqrt(theta["tau_c"] + tau_cPat)))
        end

        # println("m_ctrl")
         # update m_i - the slopes for the control subjects
        for id in ctrlIDs
            cl = "c["*id*"]"
            ml = "m["*id*"]"
            tau_mCrl = theta["tau_norm"] *ctrlData_summs["sum_mitochan_sq"][id]
            mu_mCrl = ( theta["tau_norm"] *(ctrlData_summs["sum_mitochan_chan"][id] - theta[cl]*ctrlData_summs["sum_mitochan"][id] ) ) / tau_mCrl
            theta[ml] = rand(Normal( (theta["mu_m"]*theta["tau_m"] + mu_mCrl*tau_mCrl)/(theta["tau_m"] + tau_mCrl), 1/sqrt(theta["tau_m"]+tau_mCrl) ))
        end

        # println("m_pts")
        # update mu_i - the slope of the patient subjects
        for pat in ptsIDs
            cl = "c["*pat*"]"
            ml = "m["*pat*"]"
            pat_likeCtrl = (sID.==pat) .& likeCtrl
            pat_nlikeCtrl = (sID.==pat) .& nlikeCtrl
            
            tau_mPat = theta["tau_norm"]*sum( ptsData[pat_likeCtrl, "mitochan_sq"] ) + hyperTheta["tau_def"]*sum( ptsData[pat_nlikeCtrl, "mitochan_sq"] )
            
            mu_mPat = ( theta["tau_norm"]*sum( ptsData[pat_likeCtrl, "mitochan_chan"] .- theta[cl] .*ptsData[pat_likeCtrl, mitochan])
                      + hyperTheta["tau_def"]*sum( ptsData[pat_nlikeCtrl, "mitochan_chan"] .- theta[cl] .*ptsData[pat_nlikeCtrl, mitochan]) ) / tau_mPat
            theta[ml] = rand(Normal( (theta["mu_m"]*theta["tau_m"] + mu_mPat*tau_mPat)/(theta["tau_m"] + tau_mPat), 1/sqrt(theta["tau_m"] + tau_mPat) ) )
        end

        # println("squared error: ctrl")
        # update tau - the error in the control and like-control patient fibres
        sq_diff = 0.0 # calculate the squared difference between the expected expression level and observed expression level
        for id in ctrlIDs
            ml = "m["*id*"]"
            cl = "c["*id*"]"
            ctrl_ind = ctrlData[:,"sampleID"] .== id
            sq_diff += sum( (theta[ml] .*ctrlData[ctrl_ind, mitochan] .+theta[cl] .- ctrlData[ctrl_ind, chan] ) .^2 )
        end
        # println("squared error: pat")
        for pat in ptsIDs
            ml = "m["*pat*"]"
            cl = "c["*pat*"]"
            pat_sID_likeCtrl = (sID .== pat) .& likeCtrl
            sq_diff += sum( (theta[ml] .*ptsData[pat_sID_likeCtrl, mitochan] .+ theta[cl] .- ptsData[pat_sID_likeCtrl, chan]) .^2 )
        end
        # update tau - the model error for like-control patients
        theta["tau_norm"] = rand(Gamma( hyperTheta["shape_tau"] + 0.5*(nFib_ctrl + sum(sum_likeCtrl) ), 1 /(hyperTheta["rate_tau"] + 0.5*sq_diff) ) )

        # println("prop def")
        # update proportion of deficiency and densities
        for pat in ptsIDs
            pil = "pi["*pat*"]"
            ml = "m["*pat*"]"
            cl = "c["*pat*"]"
            pat_sID = sID .== pat
            theta[pil] = rand(Beta( hyperTheta["alpha_pi"] + sum_nlikeCtrl[pat], hyperTheta["beta_pi"] + sum_likeCtrl[pat] ) )
            
            ptsData[pat_sID,"denDef"] .= theta[pil] .* pdf.(Normal.(theta[ml].*ptsData[pat_sID, mitochan] .+ theta[cl], 1/sqrt(hyperTheta["tau_def"])), ptsData[pat_sID, chan] )
            ptsData[pat_sID,"denNorm"] .= (1 - theta[pil]) .* pdf.(Normal.(theta[ml].*ptsData[pat_sID, mitochan] .+ theta[cl], 1/sqrt(theta["tau_norm"]) ), ptsData[pat_sID, chan] )
        end
        
        u = rand.(Uniform.(0.0, ptsData[:,"denDef"] .+ptsData[:,"denNorm"]))
        ptsData[:,"classif"] .= (u .< ptsData[:,"denDef"])
        
        likeCtrl = .!ptsData[:,"classif"]
        nlikeCtrl =  ptsData[:,"classif"]
    
        # population densities
        theta["m_pred"] = rand(Normal(theta["mu_m"], 1/sqrt(theta["tau_m"])))
        theta["c_pred"] = rand(Normal(theta["mu_c"], 1/sqrt(theta["tau_c"])))

        outIter = tt - warmIter
        if outIter > 0 && mod(outIter,thin)==0 
            post[paramNames, outIter] = theta[paramNames]
            classifs_mat[:, outIter] = ptsData[:, "classif"]
        end
    end
    
    """
    ######
    ### PREDICTIVE INTERVALS
    ######
    
    # only "care" about predictive interval of patient 
    p = (0.025, 0.5, 0.975)
    
    postpred = Matrix{Float64}(undef, nSyn, 4)
    m_post = post[:,mlab[nSbj]]
    c_post = post[:,clab[nSbj]]
    tau_post = post[:,"tau_norm"]
    
    priorpred = Matrix{Float64}(undef, nSyn, 4)
    m_prior = prior[:,mlab[nSbj]]
    c_prior = prior[:,clab[nSbj]]
    tau_prior = prior[:,"tau_norm"]

    
    for j in 1:nSyn
        postpred_x = rand.(Normal.(m_post .*xSyn[j] .+ c_post, 1 ./sqrt.(tau_post) ))
        priorpred_x = rand.(Normal.(m_prior .*xSyn[j] .+ c_prior, 1 ./sqrt.(tau_prior) ))
        postpred[j,:] = [xSyn[j]; collect(quantile(postpred_x, p))]
        priorpred[j,:] = [xSyn[j];  collect(quantile(priorpred_x, p))]
    end

    postpred_named = NamedArray( postpred )
    priorpred_named = NamedArray( priorpred )
    
    setnames!(postpred_named, ["mitochan", "lwr", "med", "upr"], 2)
    setnames!(priorpred_named, ["mitochan", "lwr", "med", "upr"], 2)
    """
    output = Dict()
    output["POST"] = post
    output["PRIOR"] = prior
    output["CLASSIF"] = classifs_mat
    
    return output
end