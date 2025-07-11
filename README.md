READ LICENSE.

# RotorWhirl
Whirl prediction in rotors with an aerostatic thrust bearing at one end.

Theory developed and validated by Riju Chatterjee and Ashutosh Patel. The general idea behind (an earlier version of) the whirl prediction method is here: [Conference Paper](https://doi.org/10.1115/GT2022-82632), [Video Presentation](https://youtu.be/lfDOsH-XRDQ)

Some of these .m files do whirl prediction. Others fit and plot curves to model aerostatic bearing behavior.

## Example Usage for Homogenous Stepped Shaft

<img width="1237" height="983" alt="image" src="https://github.com/user-attachments/assets/02e034d8-4a6a-4e67-8ea9-27cf66c658d7" />

Define rotor, thrust bearing, material properties and simulation settings.
```
rotor_data = example_rotor()
bearing_data = example_bearing()
material_properties = default_material_properties()
simulation_settings = default_simulation_settings()
```
Make any changes required
e.g. change step sizes
```
rotor_data.shaft_step_values = [10 50 -50 -10]*1e-3
```
e.g. change number of processor cores
```
simulation_settings.num_processor_cores = 8
```
Run analysis
```
RunRotorWhirlSimulations(rotor_data,bearing_data,material_properties,simulation_settings)
```

<img width="1368" height="952" alt="image" src="https://github.com/user-attachments/assets/109f89d5-4a41-4172-941e-7008923fb4a5" />

<img width="1351" height="952" alt="image" src="https://github.com/user-attachments/assets/ae6020f6-1ff2-49e9-aaf5-425a5bf68064" />

<img width="875" height="656" alt="kbval1Frame1" src="https://github.com/user-attachments/assets/31cb5ae4-b00e-4aa4-ad42-80a29b30f52f" />

[^1]: https://ecommons.udayton.edu/graduate_theses/767 
