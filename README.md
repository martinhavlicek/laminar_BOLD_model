# laminar_BOLD_model
Matlab code of laminar BOLD signal model as described in Havlicek, M. &amp; Uludag, K. A dynamical model of the laminar BOLD response, , BioRxiv, doi: https://doi.org/10.1101/609099. Tested on Matlab R2015b.

Author: Martin Havlicek, Maastricht Unversity (m.havlicek@maastrichtuniversity.nl)

Main files:

LBR_model.m - main function containing laminar BOLD signal model (including two simple examples how to generate steady-state and dynamic              response)
LBR_parameters.m - supporting function for the LBR model, containg structure with all parameteres of the model. Output of this function is the first input to the LBR_model

Additional files:

neuronal_NVC_model.m - model to simulated neuronal and blood flow response, which can be used as a driving input for LBR_model in case of simulating dynamic responses. 
neuronal_NVC_parameters.m - supporting function for the neuronal_NVC_model, containg structure with all parameteres of the model. Output of this function is the first input to the LBR_model

Examples:


Updated: 10 August 2019 
