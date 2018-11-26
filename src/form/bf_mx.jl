export
SDPUBFPowerModel, SDPUBFForm,
LPUBFForm, LPfullUBFPowerModel, LPfullUBFForm, LPdiagUBFPowerModel, LPdiagUBFForm,
SOCUBFForm, SOCConicUBFPowerModel, SOCNLPUBFPowerModel, SOCConicUBFForm, SOCNLPUBFForm

""
abstract type AbstractNLPUBFForm <: PMs.AbstractBFQPForm end

""
abstract type AbstractConicUBFForm <: PMs.AbstractBFConicForm end

AbstractUBFForm = Union{AbstractNLPUBFForm, AbstractConicUBFForm}


"SDP BFM per Gan and Low 2014, PSCC"
abstract type SDPUBFForm <: AbstractConicUBFForm end


"SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as an QCP"
abstract type SOCNLPUBFForm <: AbstractNLPUBFForm end

"SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as a SOC"
abstract type SOCConicUBFForm <: AbstractConicUBFForm end

SOCUBFForm = Union{SOCNLPUBFForm, SOCConicUBFForm}


"Abstract form for linear unbalanced power flow models"
abstract type AbstractLPUBFForm <: AbstractNLPUBFForm end

"Simplified BFM per Gan and Low 2014, PSCC, using matrix variables for power, voltage and current"
abstract type LPfullUBFForm <: AbstractLPUBFForm end

"LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current"
abstract type LPdiagUBFForm <: AbstractLPUBFForm end


""
const SDPUBFPowerModel = GenericPowerModel{SDPUBFForm}

"default SDP unbalanced DistFlow constructor"
SDPUBFPowerModel(data::Dict{String,Any}; kwargs...) =
GenericPowerModel(data, SDPUBFForm; kwargs...)

""
const SOCNLPUBFPowerModel = GenericPowerModel{SOCNLPUBFForm}

"default SOC unbalanced DistFlow constructor"
SOCNLPUBFPowerModel(data::Dict{String,Any}; kwargs...) =
GenericPowerModel(data, SOCNLPUBFForm; kwargs...)

""
const SOCConicUBFPowerModel = GenericPowerModel{SOCConicUBFForm}

"default SOC unbalanced DistFlow constructor"
SOCConicUBFPowerModel(data::Dict{String,Any}; kwargs...) =
GenericPowerModel(data, SOCConicUBFForm; kwargs...)

""
const LPfullUBFPowerModel = GenericPowerModel{LPfullUBFForm}

"default LP unbalanced DistFlow constructor"
LPfullUBFPowerModel(data::Dict{String,Any}; kwargs...) =
GenericPowerModel(data, LPfullUBFForm; kwargs...)

""
const LPdiagUBFPowerModel = GenericPowerModel{LPdiagUBFForm}

"default LP unbalanced DistFlow constructor"
LPdiagUBFPowerModel(data::Dict{String,Any}; kwargs...) =
GenericPowerModel(data, LPdiagUBFForm; kwargs...)

function variable_tp_branch_current(pm::GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_branch_series_current_prod_hermitian(pm; kwargs...)
end

function variable_tp_voltage(pm::GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_voltage_prod_hermitian(pm; kwargs...)
end

""
function variable_tp_voltage_prod_hermitian(pm::GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    buses = ref(pm, nw, :bus)

    wmaxdict = Dict{Int64, Any}()
    for i in ids(pm, nw, :bus)
        wmaxdict[i] = ref(pm, nw, :bus, i, "vmax").values*ref(pm, nw, :bus, i, "vmax").values'
    end
    wi_re = Matrix(n_cond,n_cond)
    wi_im = Matrix(n_cond,n_cond)

    for c1 in 1:n_cond
        for c2 in 1:n_cond
            if bounded
                if c1<=c2
                    wi_re[c1,c2] = @variable(pm.model,
                    [i in ids(pm, nw, :bus)], basename="$(nw)_$(c1)_$(c2)_wr",
                    lowerbound = -wmaxdict[i][c1,c2],
                    upperbound =  wmaxdict[i][c1,c2],
                    # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                    )
                end
                if c1<c2
                    wi_im[c1,c2] = @variable(pm.model,
                    [i in ids(pm, nw, :bus)], basename="$(nw)_$(c1)_$(c2)_wi",
                    lowerbound = -wmaxdict[i][c1,c2],
                    upperbound =  wmaxdict[i][c1,c2],
                    # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                    )
                end
            else
                if c1<=c2
                    wi_re[c1,c2] = @variable(pm.model,
                    [i in ids(pm, nw, :bus)], basename="$(nw)_$(c1)_$(c2)_wr",
                    # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                    )
                end
                if c1<c2
                    wi_im[c1,c2] = @variable(pm.model,
                    [i in ids(pm, nw, :bus)], basename="$(nw)_$(c1)_$(c2)_wi",
                    # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                    )
                end
            end
        end
    end

    #Store dictionary with matrix variables by bus
    W_re_dict = Dict{Int64, Any}()
    W_im_dict = Dict{Int64, Any}()

    # w is going to be the vector of diagonal elements of W
    for c in conductor_ids(pm)
        var(pm, nw, c)[:w] = Dict()
    end

    #reshape vectors to matrices
    for i in ids(pm, nw, :bus)
        if n_cond  == 3
            W_re_dict[i] =
            [wi_re[1,1][i] wi_re[1,2][i] wi_re[1,3][i];
            wi_re[1,2][i] wi_re[2,2][i] wi_re[2,3][i];
            wi_re[1,3][i] wi_re[2,3][i] wi_re[3,3][i]]

            W_im_dict[i] =
            [0              wi_im[1,2][i]   wi_im[1,3][i];
            -wi_im[1,2][i]   0               wi_im[2,3][i];
            -wi_im[1,3][i]   -wi_im[2,3][i]   0]
        elseif n_cond == 4
            W_re_dict[i] =
            [wi_re[1,1][i] wi_re[1,2][i] wi_re[1,3][i] wi_re[1,4][i];
            wi_re[1,2][i] wi_re[2,2][i] wi_re[2,3][i] wi_re[2,4][i];
            wi_re[1,3][i] wi_re[2,3][i] wi_re[3,3][i] wi_re[3,4][i];
            wi_re[1,4][i] wi_re[2,4][i] wi_re[3,4][i] wi_re[4,4][i]]

            W_im_dict[i] =
            [0              wi_im[1,2][i]   wi_im[1,3][i]   wi_im[1,4][i];
            -wi_im[1,2][i]   0               wi_im[2,3][i]   wi_im[2,4][i];
            -wi_im[1,3][i]   -wi_im[2,3][i]   0               wi_im[3,4][i];
            -wi_im[1,4][i]   -wi_im[2,4][i]   -wi_im[3,4][i]   0]
        else
            error("this number of conductors is not supported")
        end

        for c in conductor_ids(pm)
            var(pm, nw, c)[:w][i] = wi_re[c,c][i]
            setlowerbound(wi_re[c,c][i], (buses[i]["vmin"][c])^2)
            setupperbound(wi_re[c,c][i], (buses[i]["vmax"][c])^2)
        end
    end
    var(pm, nw)[:W_re] = W_re_dict
    var(pm, nw)[:W_im] = W_im_dict
end

function variable_tp_branch_series_current_prod_hermitian(pm::GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    branches = ref(pm, nw, :branch)
    buses = ref(pm, nw, :bus)

    cmax = Dict([(key, zeros(n_cond)) for key in keys(branches)])

    for (key, branch) in branches
        bus_fr = buses[branch["f_bus"]]
        bus_to = buses[branch["t_bus"]]

        vmin_fr = bus_fr["vmin"].values
        vmin_to = bus_to["vmin"].values

        vmax_fr = bus_fr["vmax"].values
        vmax_to = bus_to["vmax"].values

        y_fr_mag = abs.(branch["g_fr"].values + im* branch["b_fr"].values)
        y_to_mag = abs.(branch["g_to"].values + im* branch["b_to"].values)

        smax = branch["rate_a"].values
        cmaxfr = smax./vmin_fr + vmax_fr.*y_fr_mag
        cmaxto = smax./vmin_to + vmax_to.*y_to_mag

        cmax[key] = max.(cmaxfr, cmaxto)
    end

    ccm_re = Matrix(n_cond,n_cond)
    ccm_im = Matrix(n_cond,n_cond)

    for c1 in 1:n_cond
        for c2 in 1:n_cond
            if bounded
                if c1<=c2
                    ccm_re[c1,c2] = @variable(pm.model,
                    [l in ids(pm, nw, :branch)], basename="$(nw)_$(c1)_$(c2)_ccmr",
                    lowerbound = -(cmax[l]*cmax[l]')[c1,c2],
                    upperbound =  (cmax[l]*cmax[l]')[c1,c2],
                    # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                    )
                end
                if c1<c2
                    ccm_im[c1,c2] = @variable(pm.model,
                    [l in ids(pm, nw, :branch)], basename="$(nw)_$(c1)_$(c2)_ccmi",
                    lowerbound = -(cmax[l]*cmax[l]')[c1,c2],
                    upperbound =  (cmax[l]*cmax[l]')[c1,c2],
                    # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                    )
                end
            else
                if c1<=c2
                    ccm_re[c1,c2] = @variable(pm.model,
                    [l in ids(pm, nw, :branch)], basename="$(nw)_$(c1)_$(c2)_ccmr",
                    # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                    )
                end
                if c1<c2
                    ccm_im[c1,c2] = @variable(pm.model,
                    [l in ids(pm, nw, :branch)], basename="$(nw)_$(c1)_$(c2)_ccmi",
                    # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                    )
                end
            end
        end
    end

    #Store dictionary with matrix variables by branch
    ccm_re_dict = Dict{Int64, Any}()
    ccm_im_dict = Dict{Int64, Any}()

    for c in conductor_ids(pm)
        var(pm, nw, c)[:cm] = Dict()
    end

    #reshape vectors to matrices
    for i in ids(pm, nw, :branch)
        if n_cond  == 3
            ccm_re_dict[i] =
            [ccm_re[1,1][i] ccm_re[1,2][i] ccm_re[1,3][i];
            ccm_re[1,2][i] ccm_re[2,2][i] ccm_re[2,3][i];
            ccm_re[1,3][i] ccm_re[2,3][i] ccm_re[3,3][i]]

            ccm_im_dict[i] =
            [0              ccm_im[1,2][i]   ccm_im[1,3][i];
            -ccm_im[1,2][i]   0               ccm_im[2,3][i];
            -ccm_im[1,3][i]   -ccm_im[2,3][i]   0]
        elseif n_cond == 4
            ccm_re_dict[i] =
            [ccm_re[1,1][i] ccm_re[1,2][i] ccm_re[1,3][i] ccm_re[1,4][i];
            ccm_re[1,2][i] ccm_re[2,2][i] ccm_re[2,3][i] ccm_re[2,4][i];
            ccm_re[1,3][i] ccm_re[2,3][i] ccm_re[3,3][i] ccm_re[3,4][i];
            ccm_re[1,4][i] ccm_re[2,4][i] ccm_re[3,4][i] ccm_re[4,4][i]]

            ccm_im_dict[i] =
            [0              ccm_im[1,2][i]   ccm_im[1,3][i]   ccm_im[1,4][i];
            -ccm_im[1,2][i]   0               ccm_im[2,3][i]   ccm_im[2,4][i];
            -ccm_im[1,3][i]   -ccm_im[2,3][i]   0               ccm_im[3,4][i];
            -ccm_im[1,4][i]   -ccm_im[2,4][i]   -ccm_im[3,4][i]   0]
        else
            error("this number of conductors is not supported")
        end

        for c in conductor_ids(pm)
            var(pm, nw, c)[:cm][i] = ccm_re[c,c][i]
        end
    end

    var(pm, nw)[:CC_re] = ccm_re_dict
    var(pm, nw)[:CC_im] = ccm_im_dict
end


""
function variable_tp_branch_flow(pm::GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    smaxdict = Dict{Tuple{Int64, Int64, Int64}, Any}()

    for (l,i,j) in ref(pm, nw, :arcs)
        cmax = ref(pm, nw, :branch, l, "rate_a").values./ref(pm, nw, :bus, i, "vmin").values
        smax = ref(pm, nw, :bus, i, "vmax").values.*cmax'
        smaxdict[(l,i,j)] = smax
    end

    P_mx = Matrix(n_cond,n_cond)
    Q_mx = Matrix(n_cond,n_cond)

    for c1 in 1:n_cond
        for c2 in 1:n_cond
            if bounded
                P_mx[c1,c2] = @variable(pm.model,
                [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(c1)_$(c2)_p",
                lowerbound = -smaxdict[(l,i,j)][c1,c2],
                upperbound =  smaxdict[(l,i,j)][c1,c2],
                # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                )

                Q_mx[c1,c2] = @variable(pm.model,
                [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(c1)_$(c2)_q",
                lowerbound = -smaxdict[(l,i,j)][c1,c2],
                upperbound =  smaxdict[(l,i,j)][c1,c2],
                # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                )
            else
                P_mx[c1,c2] = @variable(pm.model,
                [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(c1)_$(c2)_p",
                # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                )

                Q_mx[c1,c2] = @variable(pm.model,
                [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(c1)_$(c2)_q",
                # start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
                )
            end
        end
    end

    p_mx_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()
    q_mx_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()

    for c in conductor_ids(pm)
        var(pm, nw, c)[:p] = Dict()
        var(pm, nw, c)[:q] = Dict()
    end

    for i in ref(pm, nw, :arcs)
        if n_cond  == 3
            P =
            [P_mx[1,1][i] P_mx[1,2][i] P_mx[1,3][i];
            P_mx[2,1][i] P_mx[2,2][i] P_mx[2,3][i];
            P_mx[3,1][i] P_mx[3,2][i] P_mx[3,3][i]]

            Q =
            [Q_mx[1,1][i]   Q_mx[1,2][i]   Q_mx[1,3][i];
            Q_mx[2,1][i]   Q_mx[2,2][i]   Q_mx[2,3][i];
            Q_mx[3,1][i]   Q_mx[3,2][i]   Q_mx[3,3][i]]
        elseif n_cond == 4
            P =
            [P_mx[1,1][i] P_mx[1,2][i] P_mx[1,3][i] P_mx[1,4][i];
            P_mx[2,1][i] P_mx[2,2][i] P_mx[2,3][i] P_mx[2,4][i];
            P_mx[3,1][i] P_mx[3,2][i] P_mx[3,3][i] P_mx[3,4][i];
            P_mx[4,1][i] P_mx[4,2][i] P_mx[4,3][i] P_mx[4,4][i]]

            Q =
            [Q_mx[1,1][i] Q_mx[1,2][i] Q_mx[1,3][i] Q_mx[1,4][i];
            Q_mx[2,1][i] Q_mx[2,2][i] Q_mx[2,3][i] Q_mx[2,4][i];
            Q_mx[3,1][i] Q_mx[3,2][i] Q_mx[3,3][i] Q_mx[3,4][i];
            Q_mx[4,1][i] Q_mx[4,2][i] Q_mx[4,3][i] Q_mx[4,4][i]]
        else
            error("this number of conductors is not supported")
        end

        p_mx_dict[i] = P
        q_mx_dict[i] = Q

        for c in conductor_ids(pm)
            var(pm, nw, c)[:p][i] = P[c,c]
            var(pm, nw, c)[:q][i] = Q[c,c]
        end
    end
    var(pm, nw)[:P_mx] = p_mx_dict
    var(pm, nw)[:Q_mx] = q_mx_dict
end


""
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, i::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    constraint_tp_theta_ref(pm, nw, i)
end


"""
Defines branch flow model power flow equations
"""
function constraint_tp_flow_losses(pm::GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to) where T <: AbstractUBFForm
    p_to = var(pm, n, :P_mx)[t_idx]
    q_to = var(pm, n, :Q_mx)[t_idx]

    p_fr = var(pm, n, :P_mx)[f_idx]
    q_fr = var(pm, n, :Q_mx)[f_idx]

    w_to_re = var(pm, n, :W_re)[t_bus]
    w_fr_re = var(pm, n, :W_re)[f_bus]

    w_to_im = var(pm, n, :W_im)[t_bus]
    w_fr_im = var(pm, n, :W_im)[f_bus]

    ccm_re =  var(pm, n, :CC_re)[i]
    ccm_im =  var(pm, n, :CC_im)[i]

    @constraint(pm.model, p_fr + p_to .==  w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)' + r*ccm_re - x*ccm_im +  w_to_re*(g_sh_to)'  + w_to_im*(b_sh_to)')
    @constraint(pm.model, q_fr + q_to .==  w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)' + x*ccm_re + r*ccm_im +  w_to_im*(g_sh_to)'  - w_to_re*(b_sh_to)')
end

""
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, i) where T <: AbstractUBFForm
    nconductors = length(PMs.conductor_ids(pm))

    w_re = var(pm, n, :W_re)[i]
    w_im = var(pm, n, :W_im)[i]

    alpha = exp(-im*ThreePhasePowerModels.wraptopi(2 * pi / nconductors ))
    beta = (alpha*ones(nconductors)).^(0:nconductors-1)
    gamma = beta*beta'

    w_re_ref = real(gamma).*w_re[1,1]
    w_im_ref = imag(gamma).*w_re[1,1]
    @constraint(pm.model, diag(w_re)[2:nconductors]        .== diag(w_re_ref)[2:nconductors]) # first equality is implied
    @constraint(pm.model, mat2utrivec(w_re) .== mat2utrivec(w_re_ref))
    @constraint(pm.model, mat2utrivec(w_im) .== mat2utrivec(w_im_ref))
end


"""
Defines voltage drop over a branch, linking from and to side voltage
"""
function constraint_tp_voltage_magnitude_difference(pm::GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm) where T <: AbstractUBFForm
    w_fr_re = var(pm, n, :W_re)[f_bus]
    w_fr_im = var(pm, n, :W_im)[f_bus]

    w_to_re = var(pm, n, :W_re)[t_bus]
    w_to_im = var(pm, n, :W_im)[t_bus]

    p_fr = var(pm, n, :P_mx)[f_idx]
    q_fr = var(pm, n, :Q_mx)[f_idx]

    p_s_fr = p_fr - (w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)')
    q_s_fr = q_fr - (w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)')

    ccm_re =  var(pm, n, :CC_re)[i]
    ccm_im =  var(pm, n, :CC_im)[i]

    #KVL over the line:
    @constraint(pm.model, diag(w_to_re) .== diag(
    w_fr_re   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'
    + r*ccm_re*r' - x     *ccm_im*r' + x*ccm_re *x' + r*ccm_im *x'))
    @constraint(pm.model, mat2utrivec(w_to_re) .== mat2utrivec(
    w_fr_re   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'
    + r*ccm_re*r' - x     *ccm_im*r' + x*ccm_re *x' + r*ccm_im *x'))
    @constraint(pm.model, mat2utrivec(w_to_im) .== mat2utrivec(
    w_fr_im   - q_s_fr  *r' + p_s_fr*x'        - x*p_s_fr'    + r*q_s_fr'
    + x*ccm_re*r' + r     *ccm_im*r' - r*ccm_re *x' + x*ccm_im *x'))
end

""
function get_solution_tp(pm::GenericPowerModel, sol::Dict{String,Any})
    add_bus_voltage_setpoint(sol, pm)
    PMs.add_generator_power_setpoint(sol, pm)
    add_branch_flow_setpoint(sol, pm)
    add_branch_current_setpoint(sol, pm)
    PMs.add_dcline_flow_setpoint(sol, pm)

    PMs.add_kcl_duals(sol, pm)
    PMs.add_sm_duals(sol, pm) # Adds the duals of the transmission lines' thermal limits.

    if haskey(pm.setting, "output") && haskey(pm.setting["output"], "original_variables") && pm.setting["output"]["original_variables"] == true
        add_rank(sol, pm)
        add_is_ac_feasible(sol, pm)
        add_original_variables(sol, pm)
    end
end


""
function add_bus_voltage_setpoint(sol, pm::GenericPowerModel)
    PMs.add_setpoint(sol, pm, "bus", "vm", :w; scale = (x,item) -> sqrt(x))
    PMs.add_setpoint(sol, pm, "bus", "w",  :w)

    add_setpoint_mx(sol, pm, "bus", "W_re",  :W_re)
    add_setpoint_mx(sol, pm, "bus", "W_im",  :W_im)

end

""
function add_branch_flow_setpoint(sol, pm::GenericPowerModel)
    if haskey(pm.setting, "output") && haskey(pm.setting["output"], "branch_flows") && pm.setting["output"]["branch_flows"] == true
        PMs.add_setpoint(sol, pm, "branch", "pf", :p; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "qf", :q; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "pt", :p; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "qt", :q; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])])
        add_setpoint_mx(sol, pm, "branch", "Pf", :P_mx; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])])
        add_setpoint_mx(sol, pm, "branch", "Qf", :Q_mx; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])])
        add_setpoint_mx(sol, pm, "branch", "Pt", :P_mx; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])])
        add_setpoint_mx(sol, pm, "branch", "Qt", :Q_mx; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])])
    end
end


""
function add_branch_current_setpoint(sol, pm::GenericPowerModel)
    if haskey(pm.setting, "output") && haskey(pm.setting["output"], "branch_flows") && pm.setting["output"]["branch_flows"] == true
        PMs.add_setpoint(sol, pm, "branch", "cm", :cm; scale = (x,item) -> sqrt(x))
        PMs.add_setpoint(sol, pm, "branch", "cc", :cm)
        add_setpoint_mx(sol, pm, "branch", "CC_re", :CC_re)
        add_setpoint_mx(sol, pm, "branch", "CC_im", :CC_im)
    end
end


function add_rank(sol, pm::GenericPowerModel; tol = 1e-6)
    add_voltage_variable_rank(sol, pm; tol=tol)
    add_current_variable_rank(sol, pm; tol=tol)
    add_branch_flow_rank(sol, pm; tol=tol)
end

function add_voltage_variable_rank(sol, pm::GenericPowerModel; tol = 1e-6)
end

function add_current_variable_rank(sol, pm::GenericPowerModel; tol = 1e-6)
end

function add_branch_flow_rank(sol, pm::GenericPowerModel; tol = 1e-6)
end

function add_is_ac_feasible(sol, pm::GenericPowerModel)
end

function add_is_ac_feasible(sol, pm::GenericPowerModel{T}) where T <: AbstractUBFForm
    for (nw, network) in pm.ref[:nw]
        branch_feasibilities = [branch["rank"] <= 1 for (b, branch) in sol["branch"]]
        branch_feasbility = all(branch_feasibilities)
        bus_feasibilities = [bus["rank"] <= 1 for (b, bus) in sol["bus"]]
        bus_feasibility = all(bus_feasibilities)

        ac_feasibility = bus_feasibility && branch_feasbility
        sol["is_ac_feasible_if_radial"] = ac_feasibility
    end
end

function add_voltage_variable_rank(sol, pm::GenericPowerModel{T}; tol = 1e-6) where T <: AbstractUBFForm
    for (nw, network) in pm.ref[:nw]
        buses       = ref(pm, nw, :bus)
        for (b, bus) in buses
            # Wre, Wim = make_hermitian_matrix_variable(sol["bus"]["$b"]["w"].values, sol["bus"]["$b"]["wr"].values, sol["bus"]["$b"]["wi"].values)
            Wre = sol["bus"]["$b"]["W_re"]
            Wim = sol["bus"]["$b"]["W_im"]
            W = Wre + im*Wim
            sol["bus"]["$b"]["rank"] = rank(W.values, tol)
        end
    end
end

function add_current_variable_rank(sol, pm::GenericPowerModel{T}; tol = 1e-6) where T <: AbstractUBFForm
    for (nw, network) in pm.ref[:nw]
        branches       = ref(pm, nw, :branch)
        for (b, branch) in branches
            # CCre, CCim = make_hermitian_matrix_variable(sol["branch"]["$b"]["cc"].values, sol["branch"]["$b"]["ccr"].values, sol["branch"]["$b"]["cci"].values)
            CCre = sol["branch"]["$b"]["CC_re"]
            CCim = sol["branch"]["$b"]["CC_im"]
            CC = CCre + im*CCim
            sol["branch"]["$b"]["CC_rank"] = rank(CC.values, tol)
        end
    end
end

function add_branch_flow_rank(sol, pm::GenericPowerModel{T}; tol = 1e-6) where T <: AbstractUBFForm
    for (nw, network) in pm.ref[:nw]
        buses       = ref(pm, nw, :bus)
        for (b, branch) in ref(pm, nw, :branch)
            g_fr = diagm(branch["g_fr"].values)
            b_fr = diagm(branch["b_fr"].values)
            y_fr = g_fr + im* b_fr

            fbus = branch["f_bus"]
            Wre = sol["bus"]["$b"]["W_re"].values
            Wim = sol["bus"]["$b"]["W_im"].values
            W = Wre + im*Wim

            CCre = sol["branch"]["$b"]["CC_re"].values
            CCim = sol["branch"]["$b"]["CC_im"].values
            CC = CCre + im*CCim

            Pij = sol["branch"]["$b"]["Pf"].values
            Qij = sol["branch"]["$b"]["Qf"].values

            Sij = Pij + im*Qij
            Ssij = Sij - W'*y_fr'

            f = [W Ssij; Ssij' CC]
            sol["branch"]["$b"]["rank"] = rank(f, tol)
        end
    end
end



""
function add_original_variables(sol, pm::GenericPowerModel)
    if !haskey(pm.setting, "output") || !haskey(pm.setting["output"], "branch_flows") || pm.setting["output"]["branch_flows"] == false
        error(LOGGER, "deriving the original variables requires setting: branch_flows => true")
    end

    for (nw, network) in pm.ref[:nw]
        #find rank-1 starting points
        ref_buses   = find_ref_buses(pm, nw)
        #TODO develop code to start with any rank-1 W variable
        buses       = ref(pm, nw, :bus)
        arcs        = ref(pm, nw, :arcs)
        branches    = ref(pm, nw, :branch)
        #define sets to explore
        all_bus_ids             = Set([b for (b, bus)      in ref(pm, nw, :bus)])
        all_arc_from_ids        = Set([(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_from)])
        all_arc_to_ids          = Set([(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_to)])
        all_branch_ids          = Set([l for (l,i,j)       in ref(pm, nw, :arcs_from)])
        visited_arc_from_ids    = Set()
        visited_arc_to_ids      = Set()
        visited_bus_ids         = Set()
        visited_branch_ids      = Set()

        for b in ref_buses
            sol["bus"]["$b"]["va"] = [0, -2*pi/3, 2*pi/3] #TODO support arbitrary angles at the reference bus
            sol["bus"]["$b"]["vm"] = ref(pm, nw, :bus, b)["vm"].values
            push!(visited_bus_ids, b)
        end

        tt = 0
        while visited_branch_ids != all_branch_ids && visited_bus_ids != all_bus_ids
            tt = tt+1
            if tt >10000
                break
            end

            remaining_arc_from_ids = setdiff(all_arc_from_ids, visited_arc_from_ids)
            remaining_arc_to_ids = setdiff(all_arc_to_ids, visited_arc_to_ids)

            candidate_arcs_from = [(l,i,j) for (l,i,j) in remaining_arc_from_ids if i in visited_bus_ids && !(j in visited_bus_ids)]
            candidate_arcs_to   = [(l,i,j) for (l,i,j) in remaining_arc_to_ids   if i in visited_bus_ids && !(j in visited_bus_ids)]

            if !isempty(candidate_arcs_from)
                (l,i,j) = arc = candidate_arcs_from[1]
                g_fr = diagm(branches[l]["g_fr"].values)
                b_fr = diagm(branches[l]["b_fr"].values)
                y_fr = g_fr + im* b_fr
                g_to = diagm(branches[l]["g_to"].values)
                b_to = diagm(branches[l]["b_to"].values)
                y_to = g_to + im* b_to
                r = branches[l]["br_r"].values
                x = branches[l]["br_x"].values
                z = (r + im*x)
                Ui = sol["bus"]["$i"]["vm"].*exp.(im*sol["bus"]["$i"]["va"])

                Pij = sol["branch"]["$l"]["Pf"].values
                Qij = sol["branch"]["$l"]["Qf"].values
                Sij = Pij + im*Qij

                Ssij = Sij - Ui*Ui'*y_fr'
                Isij = (1/trace(Ui*Ui'))*(Ssij')*Ui
                Uj = Ui - z*Isij
                Iij = Isij + y_fr*Ui

                Isji = -Isij
                Iji = Isji + y_to*Uj

                sol["bus"]["$j"]["vm"] = abs.(Uj)
                sol["bus"]["$j"]["va"] = wraptopi(angle.(Uj))

                sol["branch"]["$l"]["cfm"] = abs.(Iij)
                sol["branch"]["$l"]["cfa"] = wraptopi(angle.(Iij))
                sol["branch"]["$l"]["ctm"] = abs.(Iji)
                sol["branch"]["$l"]["cta"] = wraptopi(angle.(Iji))
                #
                push!(visited_arc_from_ids, arc)
                push!(visited_branch_ids, l)
                push!(visited_bus_ids, j)

            elseif !isempty(candidate_arcs_to)
                (l,i,j) = arc = candidate_arcs_to[1]
                g_fr = diagm(branches[l]["g_to"].values)
                b_fr = diagm(branches[l]["b_to"].values)
                y_fr = g_fr + im* b_fr
                g_to = diagm(branches[l]["g_fr"].values)
                b_to = diagm(branches[l]["b_fr"].values)
                y_to = g_to + im* b_to
                r = branches[l]["br_r"].values
                x = branches[l]["br_x"].values
                z = (r + im*x)
                Ui = sol["bus"]["$i"]["vm"].*exp.(im*sol["bus"]["$i"]["va"])

                Pij = sol["branch"]["$l"]["Pf"].values
                Qij = sol["branch"]["$l"]["Qf"].values
                Sij = Pij + im*Qij

                Ssij = Sij - Ui*Ui'*y_fr'
                Isij = (1/trace(Ui*Ui'))*(Ssij')*Ui
                Uj = Ui - z*Isij
                Iij = Isij + y_fr*Ui

                Isji = -Isij
                Iji = Isji + y_to*Uj


                sol["bus"]["$j"]["vm"] = abs.(Uj)
                sol["bus"]["$j"]["va"] = wraptopi(angle.(Uj))

                sol["branch"]["$l"]["ctm"] = abs.(Iij)
                sol["branch"]["$l"]["cta"] = wraptopi(angle.(Iij))
                sol["branch"]["$l"]["cfm"] = abs.(Iji)
                sol["branch"]["$l"]["cfa"] = wraptopi(angle.(Iji))
                #
                push!(visited_arc_to_ids, arc)
                push!(visited_branch_ids, l)
                push!(visited_bus_ids, j)

            else #in case you have loops or multiple reference buses
                candidate_arcs = [(l,i,j) for (l,i,j) in remaining_arc_from_ids if i in visited_bus_ids && j in visited_bus_ids]
                (l,i,j) = arc = candidate_arcs[1]
                Sij = sol["branch"]["$l"]["pf"] + im* sol["branch"]["$l"]["qf"]
                Sji = sol["branch"]["$l"]["pt"] + im* sol["branch"]["$l"]["qt"]
                Ui = sol["bus"]["$i"]["vm"].*exp.(im*sol["bus"]["$i"]["va"])
                Uj = sol["bus"]["$j"]["vm"].*exp.(im*sol["bus"]["$j"]["va"])

                Iij = Sij./Ui
                Iji = Sji./Uj
                sol["branch"]["$l"]["cfm"] = abs.(Iij)
                sol["branch"]["$l"]["cfa"] = wraptopi(angle.(Iij))
                sol["branch"]["$l"]["ctm"] = abs.(Iji)
                sol["branch"]["$l"]["cta"] = wraptopi(angle.(Iji))
                push!(visited_arc_from_ids, arc)
                push!(visited_arc_to_ids, (l,j,i))
                push!(visited_branch_ids, l)
            end
        end
    end
end



"adds matrix values based on JuMP variables"
function add_setpoint_mx(
    sol,
    pm::GenericPowerModel,
    dict_name,
    param_name,
    variable_symbol;
    index_name = "index",
    default_value = (item) -> NaN,
    scale = (x,item,cnd) -> x,
    extract_var = (var,idx,item) -> var[idx],
    sol_dict = get(sol, dict_name, Dict{String,Any}()),
    conductorless = false
)

    if InfrastructureModels.ismultinetwork(pm.data)
        data_dict = pm.data["nw"]["$(pm.cnw)"][dict_name]
    else
        data_dict = pm.data[dict_name]
    end

    if length(data_dict) > 0
        sol[dict_name] = sol_dict
    end
    for (i,item) in data_dict
        idx = Int(item[index_name])
        sol_item = sol_dict[i] = get(sol_dict, i, Dict{String,Any}())
        if conductorless

        else
            num_conductors = length(conductor_ids(pm))
            sol_item[param_name] = MultiConductorMatrix{Real}([default_value(item) for i in 1:num_conductors, j in 1:num_conductors])
            try
                variable = extract_var(var(pm, pm.cnw, variable_symbol), idx, item)

                #sol_item[param_name] = scale(getvalue(variable), item) #TODO: what is the purpose of this line?
                sol_item[param_name] = MultiConductorMatrix{Real}(getvalue(variable))
            catch
            end
        end
    end
end
