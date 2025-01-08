
# Importing the header files
import pterasoftware as ps
import numpy as np
import matplotlib.pyplot as plt
from SelfFunctions import extract_forces, extract_second_cycle
import os 
import csv

# Function to write and save data onto the csv files
def Mk4SaveDataToCSV(Airfoil, FlappingPeriod, Va, AoA, mw_wingspan, 
                     mw_root_chord, taper_ratio, Lift,InducedDrag, file_name):
    
    # Ensure that Lift  and InducedDrag are numpy arrays
    Lift = np.asarray(Lift)
    InducedDrag = np.asarray(InducedDrag)

    # Total number of points for one flapping period
    num_points = Lift.size
    
    # Ensure that all arrays have the same length
    if not all(arr.size == num_points for arr in [InducedDrag]):
        raise ValueError(f"Lift, , and InducedDrag must all have the same length of {num_points} points.")
    
    # Create a normalized time array (range between 0 and 1)
    normalized_time = np.linspace(0, 1, num_points, endpoint=False)
    
    # Write the CSV file in the current directory
    with open(file_name, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        # # Write the header (column names)
        # writer.writerow(['Flapping Frequency', 'AirSpeed', 'Angle of Attack', 'Normalised Time', 'Lift', 'Induced Drag', 'Pitching Moment', 
        #                  'MW Root Airfoil', 'MW Root Chord', 'MW Wingspan', 'MW Tip Chord', 'Tail BPosition', 'Tail Type', 'Tail Root Chord', 
        #                  'Tail Tip Chord', 'Tail Wingspan', 'Tail Airfoil'])
        
        # Write the data for each time step
        for i in range(num_points):
            writer.writerow([ Airfoil, mw_wingspan, mw_root_chord, taper_ratio, FlappingPeriod, Va, AoA, normalized_time[i], Lift[i], InducedDrag[i]])
    print(f"Data has been successfully saved to {file_name}")



    
# Function to run the unsteady aerodynamic model with conceptual design models and (va, aoa, and flapping period) as the input parameters

def simulation(mw_airfoil = "naca8304", fp = 0.2, va = 5, aoa = 5, mw_wingspan = 1.4, mw_root_chord = 0.3, taper_ratio = 0.6):
    mw_tip_chord = taper_ratio * mw_root_chord
    example_airplane=ps.geometry.Airplane(
        name="naca8304",
        x_ref=0.11,
        y_ref=0.0,
        z_ref=0.0,
        s_ref=None,
        b_ref=None,
        c_ref=None,
        wings=[
            ps.geometry.Wing(
                name="Main Wing",
                x_le=0.0,
                y_le=0.0,
                z_le=0.0,
                symmetric=True,
                num_chordwise_panels=6,
                chordwise_spacing="cosine",
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        x_le=0.0,
                        y_le=0.0,
                        z_le=0.0,
                        twist=0.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.0,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=8,
                        spanwise_spacing="cosine",
                        chord=mw_root_chord,
                        airfoil=ps.geometry.Airfoil(
                            name= mw_airfoil,
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=400,
                        ),
                    ),
                    ps.geometry.WingCrossSection(
                        x_le=0.0,
                        y_le=mw_wingspan/2,
                        z_le=0.0,
                        chord=mw_tip_chord,
                        twist=0.0,
                        airfoil=ps.geometry.Airfoil(
                            name="naca8304",
                        ),
                    ),
                ],
            ),
            ps.geometry.Wing(
                name="V-Tail",
                x_le=0.45,
                y_le=0.0,
                z_le=0.0,
                num_chordwise_panels=6,
                chordwise_spacing="cosine",
                symmetric=True,
                wing_cross_sections=[
                    ps.geometry.WingCrossSection(
                        chord=0.2,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.1,
                        control_surface_deflection=0.0,
                        airfoil=ps.geometry.Airfoil(
                            name="naca0004",
                        ),
                        twist=0.0,
                    ),
                    ps.geometry.WingCrossSection(
                        x_le=0.19,
                        y_le=0.2,
                        z_le=0.003,
                        chord=0.01,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.1,
                        control_surface_deflection=0.0,
                        twist=0.0,
                        airfoil=ps.geometry.Airfoil(
                            name="naca0004",
                        ),
                    ),
                ],
            ),
        ],
    )
    main_wing_root_wing_cross_section_movement=ps.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[0],
    )
    main_wing_tip_wing_cross_section_movement=ps.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[1],
        sweeping_amplitude=30.0,
        #____________________________________________________________________
        sweeping_period=fp,
        #____________________________________________________________________
        sweeping_spacing="sine",
        pitching_amplitude=0.0,
        pitching_period=0.0,
        pitching_spacing="sine",
        heaving_amplitude=0.0,
        heaving_period=0.0,
        heaving_spacing="sine",
    )
    v_tail_root_wing_cross_section_movement=ps.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[1].wing_cross_sections[0],
    )
    v_tail_tip_wing_cross_section_movement=ps.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[1].wing_cross_sections[1],
    )
    main_wing_movement=ps.movement.WingMovement(
        base_wing=example_airplane.wings[0],
        wing_cross_sections_movements=[
            main_wing_root_wing_cross_section_movement,
            main_wing_tip_wing_cross_section_movement,
        ],
    )
    del main_wing_root_wing_cross_section_movement
    del main_wing_tip_wing_cross_section_movement
    v_tail_movement=ps.movement.WingMovement(
        base_wing=example_airplane.wings[1],
        wing_cross_sections_movements=[
            v_tail_root_wing_cross_section_movement,
            v_tail_tip_wing_cross_section_movement,
        ],
    )
    del v_tail_root_wing_cross_section_movement
    del v_tail_tip_wing_cross_section_movement
    airplane_movement=ps.movement.AirplaneMovement(
        base_airplane=example_airplane,
        wing_movements=[main_wing_movement,v_tail_movement],
    )
    del main_wing_movement
    del v_tail_movement
    example_operating_point=ps.operating_point.OperatingPoint(

        #____________________________________________________________________
        density=1.225,
        beta=0.0,
        velocity=va,
        alpha=aoa,
        nu=15.06e-6,
        external_thrust=0.0,
        #____________________________________________________________________

    )
    operating_point_movement=ps.movement.OperatingPointMovement(
        base_operating_point=example_operating_point,
        velocity_amplitude=0.0,
        velocity_period=0.0,
        velocity_spacing="sine",
    )
    movement=ps.movement.Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        num_steps=None,
        delta_time=None,
    )
    del airplane_movement
    del operating_point_movement
    example_problem=ps.problems.UnsteadyProblem(
        movement=movement,
    )
    example_solver=ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=example_problem,
    )
    del example_problem
    example_solver.run(
        logging_level="Warning",
        prescribed_wake=True,
    )

    lift, induced_drag, side_force, pitching_moment = extract_forces(example_solver)  

    return lift, induced_drag







