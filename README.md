# Linear modulation and reception

A lab report with all answers, figures, motivations, and comments must be handed in. The lab report must be handed in by groups of 2 students. An otherwise excellent analysis and assessment would be of no practical use if we are not able to illustrate our findings. Matlab offers a multitude of tools to illustrate data. See some examples by typing

```
>> help graph2d
>> help graphics
```

The purpose of the given exercises is to create various signals, perform simple processing tasks and illustrate the results in appropriate ways. When handing in your report, remember to carefully label the figures. The following Matlab functions might be interesting for this purpose:

```
>> title
>> xlabel
>> ylabel
>> legend
```

## Pulse-amplitude modulation (PAM)
In this exercise, we simulate baseband, PAM-based, digital communications waveforms with Nyquist pulse-shaping. Furthermore, we investigate the effect of additive Gaussian noise at a receiver. The basic signal model is shown in Figure 1, and we study how to generate a symbol sequence ak, how to generate a modulation pulse shape p(t), how to perform the modulation with this pulse, and how to plot the eye-diagram in Matlab.

## Phase-shift keying (PSK)
In this second section, we consider two complex-valued symbol alphabets, phase shift keying (PSK) and quadrature amplitude modulation (QAM). We assume that the oversampling factor FsT = 1, so there is no sense to generate a pulse-shaping lter in this case (it would just act like a delay element).
