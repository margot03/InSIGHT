# Real_modeling
Currently modeling photovoltaic (PV) cells in python. Eventually we will move on to model energy storage methods in combination with PV arrays to simulate a cost-minimizing consumer with a behind-the-meter PV array and battery. Our future goal is to model other forms of energy storage as well, including green hydrogen!

# Process
Using data from PV datasheets and weather databases (ex. DNI, DHI, GHI) we aim to simulate the power output of PV cells in differing conditions. Utilizing mathematical equations that describe the behavior of PV cells, including aspects such as internal resistance, we can create IV and PV curves and calculate the maximum power output.

# Code
Python scripts that contain classes for the PV model provide the modularity needed to switch between files and datasets easily. Jupyter Notebooks allow for quick and simple visualizations to be made from the output of the model, including IV curve. Currently, we are using Matplotlib for plotting.


# Sources
The Villalva paper introduces how PV cells work and the fundamental equations and processes that desbrice PV behavior (doi): 10.1109/TPEL.2009.2013862
The author of the Villalva paper also created this resource as an example of finding resistor values, current, and `a` iteratively: https://www.dropbox.com/s/1nst7cbws6jsv5q/PV_MODEL_PACK_version_4_FEBRUARY_2015.zip?dl=0&file_subpath=%2FPV_MODEL_PACK_version_4_FEBRUARY_2015%2FModeling+Algorithms%2FPV_Modeling_Method_1_Iteractive_Villalva_Algorithm.m
