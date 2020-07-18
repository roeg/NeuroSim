# NeuroSim

## Interface library providing support for simulation of synaptic input to biophysically and morphologically detailed single neuron models

This is a python library that enables simulation of processing of tens of thousands of synaptic inputs by detailed single neuron models. 
The backend for simulation of biophysical processes uses the NEURON software package.
An example model implementation used in my research can be found [here](https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=239145#tabs-2).

## Software requirements:
 
The model requires a Linux or Mac OS X installation.

Software packages that have to be installed:

Python 2.7
    
matplotlib
    
numpy
    
NEURON >= 7.2
    
Needs to be installed such that it can be imported as python module.
        For installation instructions see <https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3489>.
        Make sure to add the folder containing the executable files such as nrnivmodl to
        your PATH environment variable. For a typical installation, try:
        
`$ export PATH=[/path/to/neuron]/[processor architecture]/bin:$PATH`

sumatra and configparser (used for reading and writing parameter files)
        
Can be installed with pip:
        
`pip install sumatra`
            
`pip install configparser`
            
            
