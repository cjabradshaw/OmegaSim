# Simulation to compare weighted (ecological) impact magnitude (omega) of invasive species
<img align="right" src="www/InvaPact logo.jpg" alt="InvaPact logo" width="200" style="margin-top: 20px">

InvaPact Workshop<br>
Vers-Pont-du-Gard, Provence, France<br>
October 2023<br>

<a href="https://github.com/cjabradshaw">Corey Bradshaw</a><br>
Flinders University<br>

## Scripts
- <code>InvaPactSim.R</code>: script setting up simulation for 3 regions of set assessment characteristics, and the calculation of omega (weighted impact magnitude), mean omega per invasive species, and InvaPacts (affected species-weighted omega per invasive species)
- <code>InvaPact conf sensitivity.R</code>: script testing changes in the confidence of one region's assessments to determine affect on mean omega differences among regions 
- <code>InvaPact impact sensitivity.R</code>: script testing changes in the impact magnitude of one region's assessments to determine affect on mean omega differences among regions
- <code>InvaPact assessment sensitivity.R</code>: script testing changes in the maximum number of assessments per invasive species in one region to determine affect on mean omega differences among regions 

## Required R libraries
- <code>jtools</code>
- <code>ggplot2</code>
- <code>ggpubr</code>

<a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" width="170" style="margin-top: 20px"></a>
<a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="80" style="margin-top: 20px"></a>
