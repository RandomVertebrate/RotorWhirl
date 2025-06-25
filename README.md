Unauthorized use/duplication of repository contents NOT PERMITTED.

# RotorWhirl
Whirl prediction in rotors with an aerostatic thrust bearing at one end.

Theory developed and validated by Riju Chatterjee and Ashutosh Patel. The general idea behind (an earlier version of) the whirl prediction method is here: [Conference Paper](https://doi.org/10.1115/GT2022-82632), [Video Presentation](https://youtu.be/lfDOsH-XRDQ)

Some of these .m files do whirl prediction. Others fit and plot curves to model aerostatic bearing behavior.

## Example
<p align = "center">
</br></br>
<img src = "https://github.com/user-attachments/assets/9d32956c-2113-4c61-a16a-eea292fbc7a7" width = "500">
</br>
<img src = "https://github.com/user-attachments/assets/985b157c-fa7f-4fab-b983-3c7411a2edfd" width = "500">
</br></br>
<img src = "https://github.com/user-attachments/assets/22eaf91c-0947-474a-bd1b-fdcdf25159e9" height = "250">
<img src = "https://github.com/user-attachments/assets/3f177b25-0eed-424e-8f15-c225bda81802" height = "250">
<img src = "https://github.com/user-attachments/assets/83b91bfd-c40b-4e75-8f0d-5c33010d21c7" height = "250">
</br></br>
</p>

## StiffShaftProfile
`StiffShaftProfile()` finds the effective profile as per the "45-degree rule"[^1] of a stepped shaft. For example:

![ArbitStiffShaft](https://github.com/RandomVertebrate/RotorWhirl/assets/54997017/76b25389-63f6-4f6c-8f54-819ed0a685fc)

[^1]: https://ecommons.udayton.edu/graduate_theses/767 
