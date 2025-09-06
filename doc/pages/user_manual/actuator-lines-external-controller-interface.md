# Wind turbine controller interface for the Actuator Line (AL) model

The controller interface module for the actuator line model included in HORSES3D relies on the use of an external, precompiled controller dynamic library file (*.so) designed following the Bladed-style DISCON interface. When active, it shares the outputs of the Actuator Line (AL) module with the external controller and receives from it the demanded control variables, namely the collective blade pitch angle and generator torque. A simple 1DOF rotor model considering angular momentum conservation in the rotor is used to explicitly compute rotor dynamics. The following scheme depicts the interaction between the different HORSES3D modules and the controller interface,

<p align="center">
   <img src="../../images/controller_interface_scheme.svg" alt="Controller interface scheme">
</p>

For a given time step, the current state is evaluated and the AL module computes the source terms representing the rotating blades, as it would if the controller were not being used. These source terms are used to compute the flow field (via CFD solver) and some of the quantities required by the controller module (e.g. the total rotor aerodynamic torque). The exchange array `avrSWAP` is used to share information between the external controller and HORSES3D. Thus, it is accordingly updated before calling the DISCON routine defined in the external controller, and retrieved at the end of said call. External controller logics are not discussed in this documentation, since these are considered user-defined. However, two main outputs are expected by the controller interface to perform variable speed control; the demanded collective blade pitch angle $\theta$ and the demanded generator torque $T_g$. 

>[!NOTE] 
>The current implementation of the controller interface does not consider individual pitch control (IPC) nor yaw control.

These demanded values will be updated for their use in the next time step. On the other hand, the rotor dynamic model computes the rotor speed for the next time step. Hence, technically these 3 tasks (solving the flow field, calling the controller and computing the rotor speed) could be perfomed simultaneously once the AL source terms are computed, although in the implementation it is done sequentially.

>[!IMPORTANT] 
>The demanded values returned by the external controller are directly applied in the next time step, regardless of the time step size. This means that if any maximum/minimum pitch and generator torque rates are to be considered to saturate the demanded values, they must be applied within the external controller.

>[!WARNING] 
>In the current implementation of the controller interface, the time step between calls to the external controller is the same that the one used used to solve the flow field and defined in the driver `controlFile`. Moreover, this same time step is also used to solve the 1DOF rotor dynamic model. This value will generally be small enough that no control or rotor dynamic instability issues arise, however it must be considered when designing external controller logics.

# Usage
To activate the wind turbine controller module in your simulation, both input flags `use actuatorline` and `actuator use controller` must be defined and set to `.true.` in the simulation driver `controlFile`. The first option will require the definition of a file named `Act_ActuatorDef.dat` inside a directory named `ActuatorDef` located at the case root directory, containing the definition of all the AL model parameters for one or multiple wind turbines. Additionally, the second (controller) flag will require the definition of a series of controller-related inputs in this same file. These have to be appended as follows,

```
############################
# CONTROLLER PARAMETERS
############################
# Path to external DISCON controller *.so file(s) (absolute or relative to this file)
controller/DISCON_glin64.so
# Path to external DISCON controller input parameter file(s) (absolute or relative to this file)
controller/DISCON_DTU10MW.IN
# Moment of intertia of wind turbine(s) drivetrain (rotor + hub + shaft...) around rotation axis (kg m2)
2.0E8
# Gear ratio from LSS to HSS (>1)
50.0
# Generator efficiency
0.944
# Wind speed/direction probe position x,y,z
-15.0 0.0 119.0
# Controller startup time (s)
0.0
# Controller data file for restart (only read if "actuator controller restart" is set to ".true.") absolute or relative to this file
../RESULTS/D10_b_vr_ctrl_0000065279_Actuator_Line_turb_001.ctrl
```
>[!WARNING]
> No comment lines (#) in the snippet above must be removed nor added. The input `Act_ActuatorDef.dat` file parsing stage is line position-based and expects to find both the header and subsequent comment lines between the definition of each input parameter value(s). Thus, altering this structure will result in an error. This also implies that no blank lines can be used.

The controller parameters, in order, are
1. Path to the external controller(s) dynamic library file to use, relative to the `Act_ActuatorDef.dat` file or absolute.
2. Path to the external controller(s) dynamic library inputs file(s), relative to the `Act_ActuatorDef.dat` file or absolute. This file typically contains a set of control parameters specific to the used dynamic library. If your custom library does not require one, simply write a dummy value.
3. Moment of intertia of wind turbine(s) drivetrain (rotor + hub + shaft...) around rotation axis in the low speed side (LSS), expressed in kg·m2.
4. Gear ratio(s) of the drivetrain from low speed side (LSS, rotor) to high speed side (HSS, generator). Note that this parameter will necessarily have to be greater or equal (in the case of assuming direct drive) to 1.
5. Generator efficiency, between 0 and 1. If your external controller already includes this parameter in its inputs and takes it into account when computing the demanded generator torque, simply set this value to 1.
6. Wind speed probe position, written as $x$, $y$ and $z$ values separated by a single space. This probe position will be used to obtain the inflow wind velocity measurements (from CFD solution) at each time step that will be passed to the external controller as an input.
7. Time before controller startup from initial simulation time, in seconds. Before this time is reacher, neither the controller module be called nor the rotor 1DOF model will be active.
8. External controller data restart `.ctrl` file path, relative to the `Act_ActuatorDef.dat` file or absolute. This input will only be read if the `actuator controller restart` parameter is set to `.true.` in the main simulation `controlFile`.

>[!NOTE] 
>In the case of simulating more than one wind turbine, each input will be specified (in separate lines) as many times as required. The first one will correspond to turbine '001', the second one to turbine '002' and so on, in an analogous manner to how the AL parameters for different turbines are specified.

```
############################
# CONTROLLER PARAMETERS
############################
# Path to external DISCON controller *.so file(s) (absolute or relative to this file)
controller/DISCON_glin64.so
controller/DISCON_glin64.so
# Path to external DISCON controller input parameter file(s) (absolute or relative to this file)
controller/DISCON_DTU10MW.IN
controller/DISCON_DTU10MW.IN
# Moment of intertia of wind turbine(s) drivetrain (rotor + hub + shaft...) around rotation axis (kg m2)
2.0E8
2.0E8
# Gear ratio from LSS to HSS (>1)
50.0
50.0
# Generator efficiency
0.944
0.944
# Wind speed/direction probe position x,y,z
-15.0 0.0 119.0
600.0 0.0 119.0
# Controller startup time (s)
5.0
5.0
# Controller data file for restart (only read if "actuator controller restart" is set to ".true.") absolute or relative to this file
../RESULTS/D10_b_vr_ctrl_0000065279_Actuator_Line_turb_001.ctrl
../RESULTS/D10_b_vr_ctrl_0000065279_Actuator_Line_turb_002.ctrl
```
>[!NOTE] 
>The same precompiled controller library file (.so) and its external input file (.IN) can be used for as many wind turbines as necessary. The controller interface will dynamically load an independent copy of the controller in memory, avoiding segmentation issues.

# Custom controller status
The first channel in the `avrSWAP` array exchanged by the external controller DISCON routine and the main program is the simulation status indicator. The status considered in Bladed v4.3 and above are
<div style="margin-left: auto;
            margin-right: auto;
            width: 50%">

| Simulation status   | Definition                               |
| :-----------------: | :--------------------------------------- |
| 0                   | Initialisation: first call at time zero  |
| 1                   | Simulation discrete time step            |
| -1                  | Final call at the end of simulation      |

</div>

Additionally, HORSES3D's controller interface includes custom status values to address certain specific tasks like autosaving the external controller status or restarting it,

<div style="margin-left: auto;
            margin-right: auto;
            width: 50%">

| Custom simulation status  | Definition                                                        |
| :-----------------------: | :---------------------------------------------------------------- |
| 2                         | Restart: first call at restart of the simulation                  |
| 3                         | Write status, sent during autosave and/or final iteration saving  |

</div>

These states are set by HORSES3D controller interface module in the in each corresponding simulation stage and can be used by the external controller to define particular restarting and saving logics.

# Other custom channels
For the considered DISCON interface (Bladed v4.3 and above), channels 120 to 129 in the `avrSWAP` array used for data exchange between the main solver and the external controller are reserved for custom, user-defined variables. In HORSES3D's controller interface, the following custom channels are used:
* Channel 120: integer number (even though the array is strictly a float array) describing the wind turbine ID, starting at 1. This value is filled by HORSES3D and matches that of the wind turbine indexing within the AL module. It can be used in your external controller to properly save and store data concerning each wind turbine when simulating a farm.
* Channel 121: integer number (even though the array is strictly a float array) corresponding to the simulation interation. This value is filled/uptated by HORSES3D only when performing an autosave operation. It can be used in your external controller to save its data to disk (if necessary) at autosave time using the same identifier as that used in the CFD solution. For the end-of-simulation saving step, this value is set to `-1`.
* Channel 129: size of avrSWAP array.

# Resources

An extension to the open-source Delft Research Controller (DRC) prepared to have full compatibility with the external controller module for Actuator Lines in HORSES3D is available [here](https://github.com/EduardoJane/DRC_Fortran_H3D). This variable-speed wind turbine controller includes the ability to write its state to disk for later restart using the described [custom controller status](#custom-controller-status) and [custom channels](#other-custom-channels) at autosave and final simulation times as specified by the user in the main `controlFile`.