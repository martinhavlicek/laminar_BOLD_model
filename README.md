# laminar_BOLD_model
Matlab code of the laminar BOLD response (LBR) model as described in Havlicek, M. &amp; Uludag, K. (2019) A dynamical model of the laminar BOLD response, NeuroImage, doi: 10.1016/j.neuroimage.2019.116209. Tested with Matlab R2015b.

Author: Martin Havlicek, Maastricht Unversity (m.havlicek@maastrichtuniversity.nl)

### Main files:
* LBR_model.m - main function containing the laminar BOLD signal model (including two simple examples how to generate steady-state and dynamic              response).
* LBR_parameters.m - supporting function for the LBR model, containg structure with all parameteres of the model. Output of this function is the first input to the LBR_model.

### Additional files:
* neuronal_NVC_model.m - model to simulate neuronal and blood flow responses, which can be used as a driving input for LBR_model in case of simulating dynamic responses.
* neuronal_NVC_parameters.m - supporting function for the neuronal_NVC_model, containg structure with all parameteres of the model. Output of this function is the first input to the LBR_model.

### Examples:
#### *Steady-state:*
* Example_ss_variable_CBF_with_different_baselineCBV.m - showing an impact of changing increase of baseline CBV in the ascending vein on LBR profile (assuming variable laminar relative CBF profile as an input).
* Example_ss_linear_scaling.m - demonstrating linear scaling of BOLD response amplitude across depths with increasing relative CBF (or stimulus strenght), which is constant across depths.
* Example_ss_nonlinear_scaling.m - demonstrating non-linear scaling of BOLD response amplitude across depths with increasing relative CBF (or stimulus strenght), which is variable across depths.
#### *Dynamic response:*
* Example_dynamic_LBR.m - showing dynamic laminar BOLD response for short and long stimulus durations.

Updated: 28 November 2019 
