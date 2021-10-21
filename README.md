[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# Infrastructure and long-duration Storage Integration using Green Hydrogen Tool (InSIGHT)
Simulating the application of hydrogen across multiple sectors, especially hard-to-decarbonize industries. Initially, a literature review of hydrogen's efficacy and historical and future uses was completed. Using a variety of sources, including Department of Energy simulation tools, data, and reports, a model for the generation of "green hydrogen" is being developed to determine when, how, and why hydrogen should (or should not) be used in certain applications, as seen from a economic and environmental standpoint.

# Solar Photovoltaic (PV) System Modeling
Modeling non-ideal photovoltaic (PV) cells in realistic environments. Graphics include IV and PV curves for sensitivity analyses. Additionally, the PV model is integrated with a battery storage system, Locational Marginal Pricing (LMP) data, etc, for optimizing an energy arbitrage scenario.

## Process
Using data from PV datasheets and weather databases (ex. DNI, DHI, GHI) we aim to simulate the power output of PV cells in differing conditions. Utilizing mathematical equations that describe the behavior of PV cells, including aspects such as internal resistance, we can create current-voltage (I-V) and power-voltage (P-V) curves and calculate the maximum power output.

## Code
Python scripts that contain classes for the PV model provide the modularity needed to switch between files and datasets easily. Jupyter Notebooks allow for quick and simple visualizations to be made from the output of the model, including I-V curve. Currently, we are using Matplotlib for plotting.

## Sources
The article authored by Marcelo Villalva, [Comprehensive Approach to Modeling and Simulation of Photovoltaic Arrays](10.1109/TPEL.2009.2013862), introduces how PV cells work and the fundamental equations and processes that desbrice ideal and non-ideal PV behavior.
Villalva also created an online [resource](https://www.dropbox.com/s/1nst7cbws6jsv5q/PV_MODEL_PACK_version_4_FEBRUARY_2015.zip?dl=0&file_subpath=%2FPV_MODEL_PACK_version_4_FEBRUARY_2015%2FModeling+Algorithms%2FPV_Modeling_Method_1_Iteractive_Villalva_Algorithm.m) as an example of calculating values such as optimal equivalent resistor values (for simulating a realistic PV device), the resulting output current for varying input voltage, and the diode ideality constant `a` (for accounting for lost current from electron reabsorbtion within the PV cell).
