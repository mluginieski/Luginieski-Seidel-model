# Luginieski-Seidel-model

# About and Contact
The codes presented here are based on following studies:<br />
(1) [Effective thickness dependence on charge carrier transport in electrolyte-gated organic field-effect transistors](http://repositorio.utfpr.edu.br/jspui/handle/1/27497)<br />
(2) [General model for charge carriers transport in electrolyte-gated transistors](http://repositorio.utfpr.edu.br/jspui/handle/1/27497)<br />
<br />
The main code is the IdFLSm.py (Id current Fit of Luginieski-Seidel model) and is responsible for fitting the experimental data with Luginieski-Seidel model. It takes use of: IdcLSm.py (Id current of Luginieski-Seidel model), which is the code where the drain current is calculated; contants.py, which contains some important physical constants; unit_converter.py, which is responsible to convert input values into SI units.
CONTACT: mluginieski@ifsc.usp.br<br />
DATE: Nov 25, 2022 
<br /><br />
# Inputs:
IdFLSm needs the experimental data. By default it is considered a two collumn (Vd and Id) .csv file with ; as column separator. For a different data file configuration, code should be fixed. Additionaly, user shoud inform the guesses and bounds of the three fitted parameters: Dc, alpha and gamma.
<br /> IdcLSm needs an input file called as IdcLSm.input. This input file have the parameters for the drain current calculation. Every row contains: parameter name, parameter value and parameter unity, where every entry is separated by a simple tabular. Please see the IdcLSm.input file as an example.

# Outputs:
The output of all codes include: list of test-fitted parameters (optional), fitted parameters with standard errors, Root-Mean-Square Error (RMSE), R-squared (R²), graph with exp and fit data, DATA_CalculatedData.csv with input and calculated parameters and constants, CalculatedData.csv with the data calculated from IdcLSm.py (Vd, Id_sat, Id_lin, l, Vg, mu).

