# Grab AeroSandbox
import aerosandbox as asb
import aerosandbox.library.aerodynamics as aero
import aerosandbox.library.atmosphere as atmo
from aerosandbox.library import mass_structural as lib_mass_struct
from aerosandbox.library import propulsion_electric as lib_prop_elec
from aerosandbox.library import propulsion_propeller as lib_prop_prop
from aerosandbox.library import power_human as lib_power_human
from aerosandbox.library import winds as lib_winds
from aerosandbox.library.airfoils import *
import copy
import matplotlib.pyplot as plt
import seaborn as sns
import json

sns.set(font_scale=1)

# region Setup
##### Initialize Optimization
opti = cas.Opti()
des_vars = {}

##### Caching
# Uncomment these lines to do ANALYSIS and OPTIMIZATION.
file_to_load_from = None
file_to_save_to = "des_vars.json"

# Uncomment these lines to do ANALYSIS on a FROZEN DESIGN.
# file_to_load_from = "des_vars.json"
# file_to_save_to = None

minimize = "-1*range/50"  # any "eval-able" expression

def des_var(  # design variable
        name,
        initial_guess,
        scale_factor,
        n_variables=1
):
    if file_to_load_from is None:
        var = scale_factor * opti.variable(n_variables)
        opti.set_initial(var, initial_guess)
        des_vars[name] = var
        return var
    else:
        var = scale_factor * opti.variable(n_variables)
        opti.set_initial(var, initial_guess)
        des_vars[name] = var
        with open(file_to_load_from, "r") as f:
            solved_des_vars = json.load(f)
        val = solved_des_vars[name]
        opti.set_initial(var, val)
        opti.subject_to(var == val)
        return var


def ops_var(  # operational variable
        initial_guess,
        scale_factor,
        n_variables=1
):
    var = scale_factor * opti.variable(n_variables)
    opti.set_initial(var, initial_guess)
    return var


##### Operating Parameters
structural_load_factor = 20  # over static with 10 deg deflection
mass_pilot = opti.parameter()
# opti.set_value(mass_pilot, 49.895) # 110 lbs, a light woman
opti.set_value(mass_pilot, 68.039) # 150 lbs, a light man
competition_altitude = 100
wind_speed_func = lambda alt: 0 # TODO placeholder
launch_altitude = opti.parameter()
opti.set_value(launch_altitude, 35 * 0.3048)
launch_airspeed = opti.parameter()
opti.set_value(launch_airspeed, 15 / 2.237)

##### Margins
structural_mass_margin_multiplier = opti.parameter()
opti.set_value(structural_mass_margin_multiplier, 1.25)

##### Simulation Parameters
n_timesteps = 100  # Only relevant if allow_trajectory_optimization is True.
# Quick convergence testing indicates you can get bad analyses below 150 or so...

##### Optimization bounds
min_speed = 1  # Specify a minimum speed - keeps the speed-gamma velocity parameterization from NaNing

# endregion

# region Trajectory Optimization Variables
##### Initialize trajectory optimization variables

x = ops_var(initial_guess=cas.linspace(0, 50, n_timesteps), scale_factor=5, n_variables=n_timesteps)

y = ops_var(initial_guess=cas.linspace(10, 0, n_timesteps), scale_factor=1, n_variables=n_timesteps)

airspeed = ops_var(initial_guess=4, scale_factor=4, n_variables=n_timesteps)
opti.subject_to([
    airspeed / min_speed > 0.1
])

flight_path_angle = ops_var(initial_guess=-12, scale_factor=1, n_variables=n_timesteps)
opti.subject_to([
    flight_path_angle / 90 < 1,
    flight_path_angle / 90 > -1,
])

alpha = ops_var(initial_guess=5, scale_factor=4, n_variables=n_timesteps)
opti.subject_to([
    alpha > -8,
    alpha < 12
])

thrust_force = ops_var(initial_guess=0, scale_factor=200, n_variables=n_timesteps)

net_accel_parallel = ops_var(initial_guess=0, scale_factor=1e-1, n_variables=n_timesteps)
net_accel_perpendicular = ops_var(initial_guess=0, scale_factor=1, n_variables=n_timesteps)

##### Set up time
time_nondim = cas.linspace(0, 1, n_timesteps)
endurance = ops_var(initial_guess=12, scale_factor=10)
opti.subject_to([
    endurance > 0
])
time = time_nondim * endurance
# endregion

# region Design Optimization Variables
##### Initialize design optimization variables (all units in base SI or derived units)

mass_total = des_var(name="mass_total", initial_guess=120, scale_factor=100)
opti.subject_to([
    mass_total > 1,
    mass_total <= 181.437, # from rules
])
max_mass_total = mass_total

### Initialize geometric variables
# wing
wing_span = des_var(name="wing_span", initial_guess=7.3152, scale_factor=5)
opti.subject_to([
    wing_span > 1,
    wing_span < 7.3152, # from rules
])

wing_root_chord = des_var(name="wing_root_chord", initial_guess=3, scale_factor=2)
opti.subject_to([
    wing_root_chord > 0.1,
    wing_root_chord < 6.096,
])

wing_taper_ratio = 0.5

# # hstab
# hstab_span = des_var(name="hstab_span", initial_guess=12, scale_factor=15)
# opti.subject_to(hstab_span > 0.1)
#
# hstab_chord = des_var(name="hstab_chord", initial_guess=3, scale_factor=2)
# opti.subject_to([hstab_chord > 0.1])
#
# hstab_twist_angle = ops_var(initial_guess=-7, scale_factor=2, n_variables=n_timesteps)
#
# # vstab
# vstab_span = des_var(name="vstab_span", initial_guess=7, scale_factor=8)
# opti.subject_to(vstab_span > 0.1)
#
# vstab_chord = des_var(name="vstab_chord", initial_guess=2.5, scale_factor=2)
# opti.subject_to([vstab_chord > 0.1])
#
# # fuselage
# boom_length = des_var(name="boom_length", initial_guess=23, scale_factor=2)  # TODO add scale factor
# opti.subject_to([
#     boom_length - vstab_chord - hstab_chord > wing_x_quarter_chord + wing_root_chord * 3 / 4
# ])

import dill as pickle

try:
    with open("wing_airfoil.cache", "rb") as f:
        wing_airfoil = pickle.load(f)
except:
    wing_airfoil = Airfoil(name="HALE_03", coordinates=r"C:\Projects\GitHub\Airfoils\HALE_03.dat")
    wing_airfoil.populate_sectional_functions_from_xfoil_fits(parallel=False)
    with open("wing_airfoil.cache", "wb+") as f:
        pickle.dump(wing_airfoil, f)
# wing_airfoil = e216

wing = asb.Wing(
    name="Main Wing",
    x_le=0,  # Coordinates of the wing's leading edge
    y_le=0,  # Coordinates of the wing's leading edge
    z_le=0,  # Coordinates of the wing's leading edge
    symmetric=True,
    xsecs=[  # The wing's cross ("X") sections
        asb.WingXSec(  # Root
            x_le=-wing_root_chord / 4,  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
            y_le=0,  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
            z_le=0,  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
            chord=wing_root_chord,
            twist=0,  # degrees
            airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
            control_surface_type='symmetric',
            # Flap # Control surfaces are applied between a given XSec and the next one.
            control_surface_deflection=0,  # degrees
        ),
        asb.WingXSec(  # Tip
            x_le=-wing_root_chord * wing_taper_ratio / 4,
            y_le=wing_span / 2,
            z_le=0,  # wing_span / 2 * cas.pi / 180 * 5,
            chord=wing_root_chord * wing_taper_ratio,
            twist=0,
            airfoil=wing_airfoil,
        ),
    ]
)

airplane = asb.Airplane(
    name="FlugTag1",
    x_ref=0,
    y_ref=0,
    z_ref=0,
    wings=[wing],
)

# endregion

# region Atmosphere
##### Atmosphere
P = atmo.get_pressure_at_altitude(y+competition_altitude)
rho = atmo.get_density_at_altitude(y+competition_altitude)
T = atmo.get_temperature_at_altitude(y+competition_altitude)
mu = atmo.get_viscosity_from_temperature(T)
a = atmo.get_speed_of_sound_from_temperature(T)
mach = airspeed / a
g = 9.81  # gravitational acceleration, m/s^2
q = 1 / 2 * rho * airspeed ** 2
# endregion

# region Aerodynamics
##### Aerodynamics

# Fuselage
CLA_fuse = 0
CDA_fuse = 0.25 # m^2, a guess for pilot drag in recumbent position

lift_fuse = CLA_fuse * q
drag_fuse = CDA_fuse * q

# wing
wing_Re = rho / mu * airspeed * wing.mean_geometric_chord()
wing_airfoil = wing.xsecs[0].airfoil  # type: asb.Airfoil
wing_Cl_inc = wing_airfoil.CL_function(alpha + wing.mean_twist_angle(), wing_Re, 0,
                                       0)  # Incompressible 2D lift coefficient
wing_CL = wing_Cl_inc * aero.CL_over_Cl(wing.aspect_ratio(), mach=mach,
                                        sweep=wing.mean_sweep_angle())  # Compressible 3D lift coefficient # TODO add in ground effect
lift_wing = wing_CL * q * wing.area()

wing_Cd_profile = wing_airfoil.CDp_function(alpha + wing.mean_twist_angle(), wing_Re, mach, 0)
drag_wing_profile = wing_Cd_profile * q * wing.area()

wing_oswalds_efficiency = 0.95  # TODO make this a function of taper ratio
drag_wing_induced = lift_wing ** 2 / (q * np.pi * wing.span() ** 2 * wing_oswalds_efficiency)

drag_wing = drag_wing_profile + drag_wing_induced

wing_Cm_inc = wing_airfoil.Cm_function(alpha + wing.mean_twist_angle(), wing_Re, 0,
                                       0)  # Incompressible 2D moment coefficient
wing_CM = wing_Cm_inc * aero.CL_over_Cl(wing.aspect_ratio(), mach=mach,
                                        sweep=wing.mean_sweep_angle())  # Compressible 3D moment coefficient
moment_wing = wing_CM * q * wing.area() * wing.mean_geometric_chord()

# Force totals
lift_force = lift_wing + lift_fuse
drag_force = drag_wing + drag_fuse

# endregion


# # region Stability
# ### Estimate aerodynamic center
# x_ac = (
#                wing.approximate_center_of_pressure()[0] * wing.area() +
#                hstab.approximate_center_of_pressure()[0] * hstab.area() * n_booms
#        ) / (
#                wing.area() + hstab.area() * n_booms
#        )
# static_margin_fraction = (x_ac - airplane.xyz_ref[0]) / wing.mean_geometric_chord()
# opti.subject_to([
#     static_margin_fraction == 0.2
# ])
#
# ### Trim
# net_pitching_moment = (
#         -wing.approximate_center_of_pressure()[0] * lift_wing + moment_wing
#         - hstab.approximate_center_of_pressure()[0] * lift_hstab * n_booms
# )
# opti.subject_to([
#     net_pitching_moment / 1e4 == 0  # Trim condition
# ])
#
# ### Size the tails off of tail volume coefficients
# Vh = hstab.area() * n_booms * (
#         hstab.approximate_center_of_pressure()[0] - wing.approximate_center_of_pressure()[0]
# ) / (wing.area() * wing.mean_geometric_chord())
# Vv = vstab.area() * n_booms * (
#         vstab.approximate_center_of_pressure()[0] - wing.approximate_center_of_pressure()[0]
# ) / (wing.area() * wing.span())
#
# hstab_effectiveness_factor = aero.CL_over_Cl(hstab.aspect_ratio()) / aero.CL_over_Cl(wing.aspect_ratio())
# vstab_effectiveness_factor = aero.CL_over_Cl(vstab.aspect_ratio()) / aero.CL_over_Cl(wing.aspect_ratio())
#
# opti.subject_to([
#     # Vh * hstab_effectiveness_factor > 0.3,
#     # Vh * hstab_effectiveness_factor < 0.6,
#     # Vh * hstab_effectiveness_factor == 0.45,
#     # Vv * vstab_effectiveness_factor > 0.02,
#     # Vv * vstab_effectiveness_factor < 0.05,
#     # Vv * vstab_effectiveness_factor == 0.035,
#     # Vh > 0.3,
#     # Vh < 0.6,
#     # Vh == 0.45,
#     Vv > 0.02,
#     # Vv < 0.05,
#     # Vv == 0.035,
# ])
#
# # endregion

# region Propulsion

### Propeller calculations
propeller_diameter = des_var(name="propeller_diameter", initial_guess=3, scale_factor=1)  # TODO scale factor
opti.subject_to([
    propeller_diameter / 0.1 > 1,
    propeller_diameter / 3.048 < 1
])

n_propellers = 1

# propeller_tip_mach = 0.36  # From Dongjoon, 4/30/20
# propeller_rads_per_sec = propeller_tip_mach * atmo.get_speed_of_sound_from_temperature(
#     atmo.get_temperature_at_altitude(20000)
# ) / (propeller_diameter / 2)
# propeller_rpm = propeller_rads_per_sec * 30 / cas.pi

area_propulsive = cas.pi / 4 * propeller_diameter ** 2 * n_propellers
propeller_coefficient_of_performance = 0.90  # a total WAG
powertrain_efficiency = 0.97 # A guess from bike data

power_out_propulsion_shaft = lib_prop_prop.propeller_shaft_power_from_thrust(
    thrust_force=thrust_force,
    area_propulsive=area_propulsive,
    airspeed=airspeed,
    rho=rho,
    propeller_coefficient_of_performance=propeller_coefficient_of_performance
)

power_out_propulsion = power_out_propulsion_shaft / powertrain_efficiency

power_out_max = des_var(name="power_out_max", initial_guess=4e2, scale_factor=2e2)
opti.subject_to([
    power_out_propulsion < power_out_max,
    power_out_max > 0
])


mass_propellers = n_propellers * lib_prop_prop.mass_hpa_propeller(
    diameter=propeller_diameter,
    max_power=power_out_max,
    include_variable_pitch_mechanism=False
)

mass_drivetrain = 3 # a total guess

# Total propulsion mass
mass_propulsion = mass_propellers + mass_drivetrain

### Power accounting
power_out = power_out_propulsion

power_in = lib_power_human.power_human(
        duration=endurance,
        dataset="Healthy Men"
        # dataset="First-Class Athletes"
        # dataset="World-Class Athletes"
    )
opti.subject_to([
    # thrust_force == 0 # No thrust
    thrust_force > 0,
    power_out < power_in
])

# endregion



### Structural mass

# Wing
n_ribs_wing = des_var(name="n_ribs_wing", initial_guess=15, scale_factor=10)
opti.subject_to([
    n_ribs_wing > 0,
])

mass_wing_primary = lib_mass_struct.mass_wing_spar(
    span=wing.span(),
    mass_supported=max_mass_total,
    # technically the spar doesn't really have to support its own weight (since it's roughly spanloaded), so this is conservative
    ultimate_load_factor=structural_load_factor,
) * 4
mass_wing_secondary = lib_mass_struct.mass_hpa_wing(
    span=wing.span(),
    chord=wing.mean_geometric_chord(),
    vehicle_mass=max_mass_total,
    n_ribs=n_ribs_wing,
    t_over_c=0.14,
    include_spar=False,
) * 4 # TODO review this number! Mark says 1.5! 4/30/2020
mass_wing = mass_wing_primary + mass_wing_secondary

# Fuselage

# The following taken from Daedalus:  # taken from Daedalus, http://journals.sfu.ca/ts/index.php/ts/article/viewFile/760/718
mass_fairings = 2.067
mass_landing_gear = 0.728

mass_fuse = mass_fairings + mass_landing_gear  # per fuselage

mass_structural = mass_wing + mass_fuse
mass_structural *= structural_mass_margin_multiplier

opti.subject_to([
    mass_total / 10 >= ( # TODO inequality only valid for fixed-mass aircraft!
            mass_pilot + mass_structural + mass_propulsion
    ) / 10
])

gravity_force = g * mass_total

# endregion

# region Dynamics

net_force_parallel_calc = (
        thrust_force * cosd(alpha) -
        drag_force -
        gravity_force * sind(flight_path_angle)
)
net_force_perpendicular_calc = (
        thrust_force * sind(alpha) +
        lift_force -
        gravity_force * cosd(flight_path_angle)
)

opti.subject_to([
    net_accel_parallel / 1e-2 == net_force_parallel_calc / mass_total / 1e-2,
    net_accel_perpendicular / 1e-1 == net_force_perpendicular_calc / mass_total / 1e-1,
])

speeddot = net_accel_parallel
gammadot = (net_accel_perpendicular / airspeed) * 180 / np.pi

trapz = lambda x: (x[1:] + x[:-1]) / 2

dt = cas.diff(time)
dx = cas.diff(x)
dy = cas.diff(y)
dspeed = cas.diff(airspeed)
dgamma = cas.diff(flight_path_angle)

xdot_trapz = trapz(airspeed * cosd(flight_path_angle))
ydot_trapz = trapz(airspeed * sind(flight_path_angle))
speeddot_trapz = trapz(speeddot)
gammadot_trapz = trapz(gammadot)

##### Winds

wind_speed = wind_speed_func(y)
wind_speed_midpoints = wind_speed_func(trapz(y))

# Total
opti.subject_to([
    dx / 1e4 == (xdot_trapz - wind_speed_midpoints) * dt / 1e4,
    dy / 1e2 == ydot_trapz * dt / 1e2,
    dspeed / 1e-1 == speeddot_trapz * dt / 1e-1,
    dgamma / 1e-2 == gammadot_trapz * dt / 1e-2,
])

# region Finalize Optimization Problem

##### Add initial state constraints
opti.subject_to([  # Push Launch
    x[0] == 0,
    y[0] == launch_altitude,
    airspeed[0] == launch_airspeed,
    flight_path_angle[0] == 0,
])

##### Add path constraints
opti.subject_to([
    y >= 0
])

##### Add final state constraints
opti.subject_to([  # Air Launch
    y[-1] == 0,
])


##### Useful metrics
wing_loading = 9.81 * max_mass_total / wing.area()
wing_loading_psf = wing_loading / 47.880258888889
empty_wing_loading = 9.81 * mass_structural / wing.area()
empty_wing_loading_psf = empty_wing_loading / 47.880258888889
propeller_efficiency = thrust_force * airspeed / power_out_propulsion_shaft
range = x[-1]

##### Add objective
objective = eval(minimize)

##### Add tippers
things_to_slightly_minimize = (
        wing_span / 80
        # + n_propellers / 1
        # + propeller_diameter / 2
        # + battery_capacity_watt_hours / 30000
        # + solar_area_fraction / 0.5
)

# Dewiggle
penalty = 0
penalty_denominator = n_timesteps
penalty += cas.sum1(cas.diff(thrust_force / 100) ** 2) / penalty_denominator
# penalty += cas.sum1(cas.diff(net_accel_parallel / 3) ** 2) / penalty_denominator
# penalty += cas.sum1(cas.diff(net_accel_perpendicular / 3) ** 2) / penalty_denominator
# penalty += cas.sum1(cas.diff(airspeed / 8) ** 2) / penalty_denominator
# penalty += cas.sum1(cas.diff(flight_path_angle / 10) ** 2) / penalty_denominator
penalty += cas.sum1(cas.diff(alpha / 5) ** 2) / penalty_denominator

opti.minimize(
    objective
    + penalty
    + 1e-6 * things_to_slightly_minimize
)
# endregion

# region Solve
p_opts = {}
s_opts = {}
s_opts["max_iter"] = 1e6  # If you need to interrupt, just use ctrl+c
s_opts["mu_strategy"] = "adaptive"
opti.solver('ipopt', p_opts, s_opts)

# endregion

if __name__ == "__main__":
    try:
        sol = opti.solve()
        # If successful, save a cache of the design variables
        solved_des_vars = {k: sol.value(des_vars[k]) for k in des_vars}
        if file_to_save_to is not None:
            with open(file_to_save_to, "w+") as f:
                json.dump(solved_des_vars, f)

    except:
        sol = opti.debug

    if np.abs(sol.value(penalty / objective)) > 0.01:
        print("\nWARNING: high penalty term! P/O = %.3f\n" % sol.value(penalty / objective))

    # # region Text Postprocessing & Utilities
    # ##### Text output
    o = lambda x: print(
        "%s: %f" % (x, sol.value(eval(x))))  # A function to Output a scalar variable. Input a variable name as a string
    outs = lambda xs: [o(x) for x in xs] and None  # input a list of variable names as strings
    print_title = lambda s: print("\n********** %s **********" % s.upper())

    print_title("Key Results")
    outs([
        "range",
        "endurance",
        "max_mass_total",
        "wing_span",
        "wing_root_chord"
    ])


    def qp(var_name):
        # QuickPlot a variable.
        fig = px.scatter(y=sol.value(eval(var_name)), title=var_name, labels={'y': var_name})
        fig.data[0].update(mode='markers+lines')
        fig.show()


    def qp2(x_name, y_name):
        # QuickPlot two variables.
        fig = px.scatter(
            x=sol.value(eval(x_name)),
            y=sol.value(eval(y_name)),
            title="%s vs. %s" % (x_name, y_name),
            labels={'x': x_name, 'y': y_name}
        )
        fig.data[0].update(mode='markers+lines')
        fig.show()


    def qp3(x_name, y_name, z_name):
        # QuickPlot two variables.
        fig = px.scatter_3d(
            x=sol.value(eval(x_name)),
            y=sol.value(eval(y_name)),
            z=sol.value(eval(z_name)),
            title="%s vs. %s" % (x_name, y_name),
            labels={'x': x_name, 'y': y_name},
            size_max=18
        )
        fig.data[0].update(mode='markers+lines')
        fig.show()


    s = lambda x: sol.value(x)

    draw = lambda: airplane.substitute_solution(sol).draw()


    # endregion

    # # Draw plots
    # def plot(x, y):
    #     # plt.plot(s(hour), s(y), ".-")
    #     plt.plot(s(x)[:dusk], s(y)[:dusk], '.-', color=(103 / 255, 155 / 255, 240 / 255), label="Day")
    #     plt.plot(s(x)[dawn:], s(y)[dawn:], '.-', color=(103 / 255, 155 / 255, 240 / 255))
    #     plt.plot(s(x)[dusk - 1:dawn + 1], s(y)[dusk - 1:dawn + 1], '.-', color=(7 / 255, 36 / 255, 84 / 255),
    #              label="Night")
    #     plt.legend()
    #
    #
    # fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8), dpi=200)
    # plot(hour, y / 1000)
    # ax.ticklabel_format(useOffset=False)
    # plt.xlabel("Hours after Solar Noon")
    # plt.ylabel("Altitude [km]")
    # plt.title("Altitude over a Day (Aug. 31)")
    # plt.tight_layout()
    # plt.savefig("outputs/altitude.png")
    # plt.show() if show_plots else plt.close(fig)
    #
    # fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8), dpi=200)
    # plot(hour, airspeed)
    # ax.ticklabel_format(useOffset=False)
    # plt.xlabel("Hours after Solar Noon")
    # plt.ylabel("True Airspeed [m/s]")
    # plt.title("True Airspeed over a Day (Aug. 31)")
    # plt.tight_layout()
    # plt.savefig("outputs/airspeed.png")
    # plt.show() if show_plots else plt.close(fig)
    #
    # fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8), dpi=200)
    # plot(hour, wing_CL)
    # ax.ticklabel_format(useOffset=False)
    # plt.xlabel("Hours after Solar Noon")
    # plt.ylabel("Lift Coefficient")
    # plt.title("Lift Coefficient over a Day (Aug. 31)")
    # plt.tight_layout()
    # plt.savefig("outputs/CL.png")
    # plt.close(fig)
    #
    # fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8), dpi=200)
    # plot(hour, net_power)
    # ax.ticklabel_format(useOffset=False)
    # plt.xlabel("Hours after Solar Noon")
    # plt.ylabel("Net Power [W] (positive is charging)")
    # plt.title("Net Power over a Day (Aug. 31)")
    # plt.tight_layout()
    # plt.savefig("outputs/net_power.png")
    # plt.show() if show_plots else plt.close(fig)
    #
    # fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8), dpi=200)
    # plot(hour, 100 * (battery_stored_energy_nondim + (1 - allowable_battery_depth_of_discharge)))
    # ax.ticklabel_format(useOffset=False)
    # plt.xlabel("Hours after Solar Noon")
    # plt.ylabel("State of Charge [%]")
    # plt.title("Battery Charge State over a Day")
    # plt.tight_layout()
    # plt.savefig("outputs/battery_charge.png")
    # plt.close(fig)
    #
    # fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8), dpi=200)
    # plot(hour, wing_Re)
    # ax.ticklabel_format(useOffset=False)
    # plt.xlabel("Hours after Solar Noon")
    # plt.ylabel("Wing Reynolds Number")
    # plt.title("Wing Reynolds Number over a Day (Aug. 31)")
    # plt.tight_layout()
    # plt.savefig("outputs/wing_Re.png")
    # plt.close(fig)
    #
    # fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8), dpi=200)
    # plot(x / 1000, y / 1000)
    # ax.ticklabel_format(useOffset=False)
    # plt.xlabel("Downrange Distance [km]")
    # plt.ylabel("Altitude [km]")
    # plt.title("Optimal Trajectory (Aug. 31)")
    # plt.tight_layout()
    # plt.savefig("outputs/trajectory.png")
    # plt.close(fig)
    #
    # fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8), dpi=200)
    # plot(hour, power_in)
    # ax.ticklabel_format(useOffset=False)
    # plt.xlabel("Hours after Solar Noon")
    # plt.ylabel("Power Generated [W]")
    # plt.title("Power Generated over a Day (Aug. 31)")
    # plt.tight_layout()
    # plt.savefig("outputs/power_in.png")
    # plt.close(fig)
    #
    # fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8), dpi=200)
    # plot(hour, power_out)
    # ax.ticklabel_format(useOffset=False)
    # plt.xlabel("Hours after Solar Noon")
    # plt.ylabel("Power Consumed [W]")
    # plt.title("Power Consumed over a Day (Aug. 31)")
    # plt.tight_layout()
    # plt.savefig("outputs/power_out.png")
    # plt.close(fig)
    #
    # # Draw mass breakdown
    # fig = plt.figure(figsize=(10, 8), dpi=200)
    # plt.suptitle("Mass Budget")
    #
    # ax_main = fig.add_axes([0.2, 0.3, 0.6, 0.6], aspect=1)
    # pie_labels = [
    #     "Payload",
    #     "Structural",
    #     "Propulsion",
    #     "Power Systems",
    #     "Avionics"
    # ]
    # pie_values = [
    #     s(mass_pilot),
    #     s(mass_structural),
    #     s(mass_propulsion),
    #     s(cas.mmax(mass_power_systems)),
    #     s(mass_avionics),
    # ]
    # colors = plt.cm.Set2(np.arange(5))
    # pie_format = lambda x: "%.1f kg\n(%.1f%%)" % (x * s(max_mass_total) / 100, x)
    # ax_main.pie(
    #     pie_values,
    #     labels=pie_labels,
    #     autopct=pie_format,
    #     pctdistance=0.7,
    #     colors=colors,
    #     startangle=120
    # )
    # plt.title("Overall Mass")
    #
    # ax_structural = fig.add_axes([0.05, 0.05, 0.3, 0.3], aspect=1)
    # pie_labels = [
    #     "Wing",
    #     "Stabilizers",
    #     "Fuses & Booms",
    #     "Margin"
    # ]
    # pie_values = [
    #     s(mass_wing),
    #     s(mass_hstab * n_booms + mass_vstab * n_booms),
    #     s(mass_fuse * n_booms),
    #     s(mass_structural - (mass_wing + n_booms * (mass_hstab + mass_vstab + mass_fuse))),
    # ]
    # colors = plt.cm.Set2(np.arange(5))
    # colors = np.clip(
    #     colors[1, :3] + np.expand_dims(
    #         np.linspace(-0.1, 0.2, len(pie_labels)),
    #         1),
    #     0, 1
    # )
    # pie_format = lambda x: "%.1f kg\n(%.1f%%)" % (x * s(mass_structural) / 100, x * s(mass_structural / max_mass_total))
    # ax_structural.pie(
    #     pie_values,
    #     labels=pie_labels,
    #     autopct=pie_format,
    #     pctdistance=0.7,
    #     colors=colors,
    #     startangle=60,
    # )
    # plt.title("Structural Mass*")
    #
    # ax_power_systems = fig.add_axes([0.65, 0.05, 0.3, 0.3], aspect=1)
    # pie_labels = [
    #     "Batt. Pack (Cells)",
    #     "Batt. Pack (Non-cell)",
    #     "Solar Cells",
    #     "Misc. & Wires"
    # ]
    # pie_values = [
    #     s(mass_battery_cells),
    #     s(mass_battery_pack - mass_battery_cells),
    #     s(mass_solar_cells),
    #     s(mass_power_systems - mass_battery_pack - mass_solar_cells),
    # ]
    # colors = plt.cm.Set2(np.arange(5))
    # colors = np.clip(
    #     colors[3, :3] + np.expand_dims(
    #         np.linspace(-0.1, 0.2, len(pie_labels)),
    #         1),
    #     0, 1
    # )[::-1]
    # pie_format = lambda x: "%.1f kg\n(%.1f%%)" % (
    #     x * s(mass_power_systems) / 100, x * s(mass_power_systems / max_mass_total))
    # ax_power_systems.pie(
    #     pie_values,
    #     labels=pie_labels,
    #     autopct=pie_format,
    #     pctdistance=0.7,
    #     colors=colors,
    #     startangle=15,
    # )
    # plt.title("Power Systems Mass*")
    #
    # plt.annotate(
    #     s="* percentages referenced to total aircraft mass",
    #     xy=(0.01, 0.01),
    #     # xytext=(0.03, 0.03),
    #     xycoords="figure fraction",
    #     # arrowprops={
    #     #     "color"     : "k",
    #     #     "width"     : 0.25,
    #     #     "headwidth" : 4,
    #     #     "headlength": 6,
    #     # }
    # )
    # plt.annotate(
    #     s="""
    #     Total mass: %.1f kg
    #     Wing span: %.2f m
    #     """ % (s(max_mass_total), s(wing.span())),
    #     xy=(0.03, 0.70),
    #     # xytext=(0.03, 0.03),
    #     xycoords="figure fraction",
    #     # arrowprops={
    #     #     "color"     : "k",
    #     #     "width"     : 0.25,
    #     #     "headwidth" : 4,
    #     #     "headlength": 6,
    #     # }
    # )
    #
    # plt.savefig("outputs/mass_pie_chart.png")
    # plt.show() if show_plots else plt.close(fig)
    #
    # # Write a mass budget
    # with open("outputs/mass_budget.csv", "w+") as f:
    #     from types import ModuleType
    #
    #     var_names = dir()
    #     f.write("Object or Collection of Objects, Mass [kg],\n")
    #     for var_name in var_names:
    #         if "mass" in var_name and not type(eval(var_name)) == ModuleType and not callable(eval(var_name)):
    #             f.write("%s, %f,\n" % (var_name, s(eval(var_name))))
    #
    # # Write a geometry spreadsheet
    # with open("outputs/geometry.csv", "w+") as f:
    #
    #     f.write("Design Variable, Value (all in base SI units or derived units thereof),\n")
    #     geometry_vars = [
    #         'wing.span()',
    #         'wing_root_chord',
    #         'wing_taper_ratio',
    #         'wing.xsecs[0].airfoil.name',
    #         '',
    #         'hstab.span()',
    #         'hstab.mean_geometric_chord()',
    #         'hstab.xsecs[0].airfoil.name',
    #         '',
    #         'vstab.span()',
    #         'vstab.mean_geometric_chord()',
    #         'vstab.xsecs[0].airfoil.name',
    #         '',
    #         'max_mass_total',
    #         '',
    #         'solar_area_fraction',
    #         'battery_capacity_watt_hours',
    #         'n_propellers',
    #         'propeller_diameter',
    #         '',
    #         'boom_length'
    #     ]
    #     for var_name in geometry_vars:
    #         if var_name == '':
    #             f.write(",,\n")
    #             continue
    #         try:
    #             value = s(eval(var_name))
    #         except:
    #             value = eval(var_name)
    #         try:
    #             f.write("%s, %f,\n" % (var_name, value))
    #         except:
    #             f.write("%s, %s,\n" % (var_name, value))
