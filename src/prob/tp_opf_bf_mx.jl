export run_tp_opf_bf_mx

""
function run_tp_opf_bf_mx(data::Dict{String,Any}, model_constructor, solver; kwargs...)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bf_mx; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function run_tp_opf_bf_mx(file::String, model_constructor, solver; kwargs...)
    data = ThreePhasePowerModels.parse_file(file)
    return PMs.run_generic_model(data, model_constructor, solver, post_tp_opf_bf_mx; solution_builder=get_solution_tp, multiconductor=true, kwargs...)
end


""
function post_tp_opf_bf_mx(pm::GenericPowerModel)
    # Variables
    variable_tp_voltage(pm)
    variable_tp_branch_current(pm)
    variable_tp_branch_flow(pm)
    variable_tp_generation(pm)
    variable_tp_load(pm)


    for c in PMs.conductor_ids(pm)
        # PMs.variable_generation(pm, cnd=c)
        PMs.variable_dcline_flow(pm, cnd=c)
    end

    # Constraints
    for i in ids(pm, :ref_buses)
        constraint_tp_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_tp_kcl_shunt_mx(pm, i)
    end


    for i in ids(pm, :gen)
        constraint_tp_kcl_gen(pm, i)
    end

    for i in ids(pm, :load)
        constraint_tp_kcl_load(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_tp_flow_losses(pm, i)

        constraint_tp_voltage_magnitude_difference(pm, i)

        constraint_tp_branch_current(pm, i)


        for c in PMs.conductor_ids(pm)
            PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

            PMs.constraint_thermal_limit_from(pm, i, cnd=c)
            PMs.constraint_thermal_limit_to(pm, i, cnd=c)
        end
    end

    for i in ids(pm, :dcline), c in PMs.conductor_ids(pm)
        PMs.constraint_dcline(pm, i, cnd=c)
    end

    PMs.objective_min_fuel_cost(pm)
end
