# Problem Specifications

## Optimal Power Flow (OPF) with On-Load Tap Changers (OLTC)

This problem is identical to `mc_opf`, except that all transformers are now modelled as on-load tap changers (OLTCs). Each phase has an individual tap ratio, which can be either variable or fixed, as specified in the data model.

### Objective

```julia
objective_min_fuel_cost(pm)
```

### Variables

```julia
variable_mc_voltage(pm)
variable_mc_branch_flow(pm)

for c in PMs.conductor_ids(pm)
    PMs.variable_generation(pm, cnd=c)
    PMs.variable_dcline_flow(pm, cnd=c)
end
variable_mc_transformer_flow(pm)
variable_mc_oltc_tap(pm)
```

### Constraints

```julia
constraint_mc_model_voltage(pm)

for i in PMs.ids(pm, :ref_buses)
    constraint_mc_theta_ref(pm, i)
end

for i in PMs.ids(pm, :bus), c in PMs.conductor_ids(pm)
    constraint_mc_power_balance(pm, i, cnd=c)
end

for i in PMs.ids(pm, :branch)
    constraint_mc_ohms_yt_from(pm, i)
    constraint_mc_ohms_yt_to(pm, i)

    for c in PMs.conductor_ids(pm)
        PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

        PMs.constraint_thermal_limit_from(pm, i, cnd=c)
        PMs.constraint_thermal_limit_to(pm, i, cnd=c)
    end
end

for i in PMs.ids(pm, :dcline), c in PMs.conductor_ids(pm)
    PMs.constraint_dcline(pm, i, cnd=c)
end

for i in PMs.ids(pm, :transformer)
    constraint_mc_oltc(pm, i)
end
```

## Optimal Power Flow (OPF) with Load Models (LM)

Unlike `mc_opf`, which models all loads as constant power loads, this problem specification additionally supports loads proportional to the voltage magnitude (a.k.a. constant current) and the square of the voltage magnitude (a.k.a. constant impedance). Each load now has associated active and reactive power variables. In `mc_opf`, loads are directly added as parameters in KCL.

### Objective

```julia
objective_min_fuel_cost(pm)
```

### Variables

```julia
variable_mc_voltage(pm)
variable_mc_branch_flow(pm)

for c in PMs.conductor_ids(pm)
    PMs.variable_generation(pm, cnd=c)
    PMs.variable_dcline_flow(pm, cnd=c)
end
variable_mc_transformer_flow(pm)
variable_mc_oltc_tap(pm)
```

### Constraints

```julia
constraint_mc_model_voltage(pm)

for i in PMs.ids(pm, :ref_buses)
    constraint_mc_theta_ref(pm, i)
end

for i in PMs.ids(pm, :bus)
    constraint_mc_power_balance_load(pm, i)
end

for id in PMs.ids(pm, :load)
    model = PMs.ref(pm, pm.cnw, :load, id, "model")
    if model=="constant_power"
        constraint_mc_load_power_setpoint(pm, id)
    elseif model=="proportional_vm"
        constraint_mc_load_power_prop_vm(pm, id)
    elseif model=="proportional_vmsqr"
        constraint_mc_load_power_prop_vmsqr(pm, id)
    else
        Memento.@error(LOGGER, "Unknown model $model for load $id.")
    end
end

for i in PMs.ids(pm, :branch)
    constraint_mc_ohms_yt_from(pm, i)
    constraint_mc_ohms_yt_to(pm, i)

    for c in PMs.conductor_ids(pm)
        PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

        PMs.constraint_thermal_limit_from(pm, i, cnd=c)
        PMs.constraint_thermal_limit_to(pm, i, cnd=c)
    end
end

for i in PMs.ids(pm, :dcline), c in PMs.conductor_ids(pm)
    PMs.constraint_dcline(pm, i, cnd=c)
end

for i in PMs.ids(pm, :transformer)
    constraint_mc_transformer(pm, i)
end
```

## Power Flow (PF) with Load Models (LM)

Unlike `mc_pf`, which models all loads as constant power loads, this problem specification additionally supports loads proportional to the voltage magnitude (a.k.a. constant current) and the square of the voltage magnitude (a.k.a. constant impedance). Each load now has associated active and reactive power variables. In `mc_pf`, loads are directly added as parameters in KCL.

### Variables

```julia
variable_mc_voltage(pm, bounded=false)
variable_mc_branch_flow(pm, bounded=false)
variable_mc_transformer_flow(pm, bounded=false)
variable_mc_load(pm)

for c in PMs.conductor_ids(pm)
    PMs.variable_generation(pm, bounded=false, cnd=c)
    PMs.variable_dcline_flow(pm, bounded=false, cnd=c)
end
```

### Constraints

```julia
constraint_mc_model_voltage(pm, bounded=false)

for (i,bus) in PMs.ref(pm, :ref_buses)
    constraint_mc_theta_ref(pm, i)

    for c in PMs.conductor_ids(pm)
        @assert bus["bus_type"] == 3
        PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c)
    end
end

for (i,bus) in PMs.ref(pm, :bus)
    constraint_mc_power_balance_load(pm, i)

    for c in PM.conductor_ids(pm)
        # PV Bus Constraints
        if length(PMs.ref(pm, :bus_gens, i)) > 0 && !(i in PMs.ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c)
            for j in PMs.ref(pm, :bus_gens, i)
                PMs.constraint_active_gen_setpoint(pm, j, cnd=c)
            end
        end
    end
end

for id in PMs.ids(pm, :load)
    model = PMs.ref(pm, pm.cnw, :load, id, "model")
    if model=="constant_power"
        constraint_mc_load_power_setpoint(pm, id)
    elseif model=="proportional_vm"
        constraint_mc_load_power_prop_vm(pm, id)
    elseif model=="proportional_vmsqr"
        constraint_mc_load_power_prop_vmsqr(pm, id)
    else
        Memento.@error(LOGGER, "Unknown model $model for load $id.")
    end
end

for i in PMs.ids(pm, :branch)
    constraint_mc_ohms_yt_from(pm, i)
    constraint_mc_ohms_yt_to(pm, i)
end

for (i,dcline) in PMs.ref(pm, :dcline), c in PMs.conductor_ids(pm)
    PMs.constraint_active_dcline_setpoint(pm, i, cnd=c)

    f_bus = PMs.ref(pm, :bus)[dcline["f_bus"]]
    if f_bus["bus_type"] == 1
        PMs.constraint_voltage_magnitude_setpoint(pm, f_bus["index"], cnd=c)
    end

    t_bus = PMs.ref(pm, :bus)[dcline["t_bus"]]
    if t_bus["bus_type"] == 1
        PMs.constraint_voltage_magnitude_setpoint(pm, t_bus["index"], cnd=c)
    end
end

for i in PMs.ids(pm, :transformer)
    constraint_mc_transformer(pm, i)
end
```

## Minimal Load Delta (MLD) Problem Specification

Load shed (continuous) problem. See "Relaxations of AC Maximal Load Delivery for Severe Contingency Analysis" by C. Coffrin _et al._ (DOI: [10.1109/TPWRS.2018.2876507](https://ieeexplore.ieee.org/document/8494809)) for single-phase case.

### Variables

```math
\begin{align}
\mbox{variables: } & \nonumber \\
& z^v_i \in \{0,1\}\ \ \forall i \in N \mbox{ - bus voltage on/off variable} \\
& z^g_i \in \{0,1\}\ \ \forall i \in G \mbox{ - generator on/off variable} \\
& z^{b}_i \in \{0,1\}\ \ \forall i \in B\mbox{ - storage on/off variable} \\
& z^d_i \in (0,1)\ \ \forall i \in L \mbox{ - continuous load shedding variable} \\
& z^s_i \in (0,1)\ \ \forall i \in H \mbox{ - continuous shunt shedding variable}
\end{align}
```

### Objective

```math
\begin{align}
\mbox{minimize: }\left (
\sum_{\substack{i\in N,c\in C}}{10 \left (1-z^v_i \right )} + \sum_{\substack{i\in L,c\in C}}{10 \omega_{i,c}\left |\Re{\left (S^d_i\right )}\right |\left ( 1-z^d_i \right ) } + \sum_{\substack{i\in H,c\in C}}{\left | \Re{\left (S^s_i \right )}\right | \left (1-z^s_i \right ) } + \sum_{\substack{i\in G,c\in C}}{\Delta^g_i } + \sum_{\substack{i\in B,c\in C}}{\Delta^b_i} \right )
\end{align}
```
where

```math
\begin{align}
\Delta^g_i &>= \left [\Re{\left (S^g_{i}(0) \right )} - \Re{\left (S^g_i \right )} \right ] \\
\Delta^g_i &>= -\left [\Re{\left (S^g_{i}(0) \right )} - \Re{\left (S^g_i \right )} \right ] \\
\Delta^b_i &>= \left [\Re{\left (S^b_{i}(0) \right )} - \Re{\left (S^b_i \right )} \right ] \\
\Delta^b_i &>= -\left [\Re{\left (S^b_{i}(0) \right )} - \Re{\left (S^b_i \right )} \right ]
\end{align}
```

### Constraints

```math
\begin{align}
\mbox{subject to: } & \nonumber \\
& z^v_i v^l_{i,c} \leq \left | V_{i,c} \right | \leq z_i^v v^u_{i,c}\ \ \forall i \in N,\forall c \in C \\
& z^g_i S^{gl}_{i,c} \leq S^g_{i,c} \leq z^g_i S^{gu}_{i,c}\ \ \forall i \in G,\forall c \in C \\
& \sum_{\substack{k\in G_i,c\in C}} S^g_{k,c} - \sum_{\substack{k\in L_i,c\in C}} z^d_k S^d_{k,c}- \sum_{\substack{k\in H_i,c\in C}} z^s_k Y^s_{k,c}\left | V_{i,c} \right |^2 \nonumber \\
& = \sum_{\substack{(i,j)\in E_i\cup E_i^R,c\in C}} S_{ij,c}\ \forall i \in N
\end{align}
```
