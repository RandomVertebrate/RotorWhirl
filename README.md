Unauthorized use/duplication NOT ALLOWED.

# RotorWhirl
Whirl prediction in rotors with an aerostatic thrust bearing at one end.

Some of these .m files do whirl prediction. Others fit and plot curves to model aerostatic bearing behavior.

The general idea behind (an earlier version of) the whirl prediction method: https://doi.org/10.1115/GT2022-82632

## StiffShaftProfile
`StiffShaftProfile()` finds the effective profile as per the "45-degree rule"[^1] of a stepped shaft. For example:

![ArbitStiffShaft](https://github.com/RandomVertebrate/RotorWhirl/assets/54997017/ee3aff71-387d-41d7-bfe6-9c8e14cd018d)

[^1]: https://ecommons.udayton.edu/graduate_theses/767 
