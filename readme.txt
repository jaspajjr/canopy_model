Repository for a model describing the growth and decline of a canopy using the function:
f(x) = (c / (1 + np.exp(-b1 * (x - m1)))) * (a + (c * np.exp((-1 * (np.exp(-b2 * (x - m2)))))))

The desired outputs if the project are:

Basics

- TT to anthesis
- TT Anthesis to senescence
- TT to senescence
- TT to GS31
- TT GS31 - anthesis
- TT Sowing to GS31
- PAR sowing - GS31
- PAR GS31 - anthesis
- PAR anthesis
- PAR anthesis - senescence
- PAR senescence

Absolute R:FR

- Max R:FR
- TT of max R:FR
- R:FR n days after sowing x 3
- R:FR n days before anthesis x 3
- R:FR n days after anthesis x 3

Rates of change

- Derivate at x
- Maximum rate of increase f'(x) = 0
- Maximum rate of decreas

Relative R:FR

- % change in R:FR n days after sowing x 3
- % change in R:FR n days before anthesis x 3
- % change in R:FR n days after anthesis x 3

Sum total of canopy duration / PAR interception
- integral between a & b
- time that canopy was above an objective threshold
- time that canopy was above a relative threshold

The model is now finished, integration above threshold, maximum rate of change were not deemed to be biologically relavent so were abandoned as objectives.
