# Notes on development

## Goal

Develop a two-layer model with Skin, Core, Vein, and Artery temperatures.
Additionaly the parameters for this model will be located in a simple text file
format. The model will contain 5 segments:

1. Head
2. Upper torso
3. Lower torso/Abdomen
4. Arm
5. Leg

This model will be a somewhat simplified version of the [Salloum model](https://doi.org/10.1016/j.ijthermalsci.2006.06.017)

## Additional components

For the simulation runs I may want to have a simple input  file. 

### Simulation file

The simulation file should contain 3 parts 

1. Simulation meta parameters (e.g.) timestep
2. Initial temperatures for the segments
3. List of "run commands" that indicate enviroment temperatures and length of time

---

## Model

### Segment

Each segment is going to have 4 temperatures: a core, skin, vein, and artery temperature.
For each temperature, there will be a heat balance equation that is used to calculate
the temperature at each timestep 

#### Core

The core heat balanace equation should contain a few heat components:

1. Basal metabolic rate
2. Shivering rate
3. BAT*
4. Respiration*
5. Core-Skin Conduction
6. Artery heat exchange
7. Vein heat exchange
8. Core perfusion (i.e. capillary exchange)
9. Skin perfusion (i.e. capillary exchange to skin)

\*These components are going to be only for the upper torso, and the location of the
bat might move to another position in the body for the 

#### Skin

1. Basal metabolic rate
2. Core-Skin Conduction
3. Radiation, Evaporation, Convection
4. Skin perfusion

#### Artery

1. Artery core exchange
2. From mass flow entering

#### Vein

1. Vein core exchange
2. From mass flow entering
3. From perfusion

### Other considerations

I may change how the blood flow through the model a bit, but I think that this works pretty well.

## Simulation Structure

The simulation should have 3 main parts: 1) setup, 2) simulation loop, 3) output

In the setup section, the parameters and simulation files are parsed and the parameters are
used to build the simulation matrix/matrix. Additioaly during the setup step, the simulation
parameters should be read in. 

The simulation loop should have the following steps:

1. The thermoregulatory controls are calculated
2. The temperature / heatflow matrix is built, and multiplied
3. Bloodflow? (this could possibly come earlier)
4. The temperatures are updated

---

# Input formats

## Paramers file (*.params)

General idea for the initial parameter file

```
 # Comment Line  
 Head Parameters
 Upper Torso Parameters
 Lower Torso Parameters
 Arm Parameters
 Leg Parameters
 Whole body parameters?
``` 

## Simulation file (*.run)

General idea for the initial simulation file 

```
# Comment Line
timestep 

head_core_temp head_skin_temp head_vein_temp head_artery_temp
upper_core_temp upper_skin_temp upper_vein_temp upper_artery_temp
lower_core_temp lower_skin_temp lower_vein_temp lower_artery_temp
arm_core_temp arm_skin_temp arm_vein_temp arm_artery_temp
leg_core_temp leg_skin_temp leg_vein_temp leg_artery_temp

r1_steps temp_conv temp_rad humidity
r2_steps temp_conv temp_rad humidity
...
```

# Code/Data Structures

## Temperatures

The temperatures for the model will be stored as a vector, with the following structure:

`[{head}, {upper torso}, {lower torso}, {arm}, {leg}]`

with each of the `{segment}` components being formed with the flowing structure:

`..., core_temp, skin_temp, vein_temp, artery_temp, ...`

## Transport mechanisms

The main code for transport mechanisms is going to consist of a series of matrix
multiplications and additionas that will determine the heat flow into/out of each of the
temperature nodes. The matrixes should give the