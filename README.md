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

<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University" width="140" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL" width="75" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.universite-paris-saclay.fr/"><img align="bottom-left" src="www/UPSlogo.png" alt="UPS" width="140" style="margin-top: 20px"></a>  &nbsp; &nbsp; <a href="https://www.cnrs.fr"><img align="bottom-left" src="www/CNRSlogo.png" alt="CNRS" width="60" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.axa.com/"><img align="bottom-left" src="www/AXAlogo.png" alt="AXA" width="50" style="margin-top: 20px"></a> &nbsp; &nbsp; &nbsp; <a href="https://anr.fr"><img align="bottom-left" src="www/ANRlogo.png" alt="ANR" width="150" style="margin-top: 20px"></a></p>
