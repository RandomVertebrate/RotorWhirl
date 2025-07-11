READ LICENSE.

# RotorWhirl
Whirl prediction in rotors with an aerostatic thrust bearing at one end.

Theory developed and validated by Riju Chatterjee and Ashutosh Patel. The general idea behind (an earlier version of) the whirl prediction method is here: [Conference Paper](https://doi.org/10.1115/GT2022-82632), [Video Presentation](https://youtu.be/lfDOsH-XRDQ)

Some of these .m files do whirl prediction. Others fit and plot curves to model aerostatic bearing behavior.

## Examples
### Whirl prediction in homogenous partially hollow stepped shaft
<p align="center">
  <img src="figs/demo_rotor.svg" height="400" alt="Rotor geometry">
</p>

From the code directory, run:

```
rotor_data = example_rotor()
bearing_data = example_bearing()
material_properties = default_material_properties()
simulation_settings = default_simulation_settings()
RunRotorWhirlSimulations(rotor_data,bearing_data,material_properties,simulation_settings)
plot_crits
```
Results:
<p align="center">
  <img src="figs/demo_rotor_profile.svg" height="400" alt="Effective profile for bending stiffness">
</p><p align="center">
  <img src="figs/demoCampbellFrame1.png" width="300" alt="Campbell diagram 1"> <img src="figs/demoCampbellFrame6.png" width="300" alt="Campbell diagram 2">
</p><p align="center">
  <img src="figs/demo_modes1.svg" width="300" alt="Mode shapes 1"> <img src="figs/demo_modes2.svg" width="300" alt="Mode shapes 2"> <img src="figs/demo_modes3.svg" width="300" alt="Mode shapes 3"> <img src="figs/demo_modes4.svg" width="300" alt="Mode shapes 4">
</p><p align="center">
  <img src="figs/demo_rotor_critmap.svg" width="300" alt="Critical speed map">
</p>

[^1]: https://ecommons.udayton.edu/graduate_theses/767 
