# Chemotaxis

A simulation of bacterial/signal movement based on the idea of chemotaxis, that is, the movement of organisms in response to chemical stimulus.

This work was completed for my Master's thesis at the University of Utrecht, entitled 'A neuron based model for data storage and retrieval'. The model is used as inspiration while exploring mechanisms to control activity levels in the neuronal network.

## Usage

Clone submodules with
```
git submodule update --init --recursive
```

Compile on Windows with
```
nmake /f NMAKEFILE
```
(Linux makefile still under development).

System properties can be set in
```
config/config.json
```

Binary is located at
```
bin/chemotaxis.exe
```
