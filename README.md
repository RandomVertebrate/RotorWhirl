Unauthorized use/duplication NOT ALLOWED.

# RotorWhirl
Whirl prediction in rotors with an aerostatic thrust bearing at one end.

Some of these .m files do whirl prediction. Others fit and plot curves to model aerostatic bearing behavior.

The general idea behind (an earlier version of) the whirl prediction method: [Conference Paper](https://doi.org/10.1115/GT2022-82632), [Video Presentation](https://youtu.be/lfDOsH-XRDQ)

## StiffShaftProfile
`StiffShaftProfile()` finds the effective profile as per the "45-degree rule"[^1] of a stepped shaft. For example:

![ArbitStiffShaft](https://github.com/RandomVertebrate/RotorWhirl/assets/54997017/76b25389-63f6-4f6c-8f54-819ed0a685fc)

[^1]: https://ecommons.udayton.edu/graduate_theses/767 
