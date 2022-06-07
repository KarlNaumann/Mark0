"""
Interface to the mark0_covid model that was originally developed
in C++ by Stanislao Gualdi. The corresponding C++ code /mark0_covid

This version of the model is from the  latest paper of Sharma et al., 
which includes optional shocks relating to the COVID-19 crisis.

"""
__author__ = "Karl Naumann-Woleske, Max Sina Knicker"
__credits__ = ["Max Sina Knicker", "Karl Naumann-Woleske"]
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = "Karl Naumann-Woleske"

import os
import yaml

import numpy as np
import pandas as pd

class Mark0_COVID(object):
    """ This class is the Python interface to the Mark0 COVID model
    dicussed in Sharma et al. (2021) "V –, U –, L – or W–shaped economic
    recovery after COVID-19: Insights from an Agent Based Model".

    The underlying methods are based on the authors' original code, written
    in C++ by Stanislao Gualdi. There have been minor modifications to this
    code, which can be found in the mark0_covid/ directory.

    Attributes
    ----------
    parameters : dict
        Dictionary of id:value for the model parameters
    hyper_parameters : dict
        Dictionary of id:value for the hyperparameters
    output_series : dict
        Dictionary of id:LaTeX name for each output to the Mark-0 model
    parameter_names : dict
        Dictionary of id:LaTeX symbol for each of the parameters in the model
    call_path : str 
        path to the mark0 executable file, e.g. ./models/mark0_covid/Debug/mark0_covid
    model_name : str, default: 'mark0_covid'
        Name of the model. This is used in the output filename generation
    base_call_func : str, default: './mark0_covid'
        commandline call to run the model 
    output : pd.DataFrame
        dataframe of the time-series of output once the model has been simulated

    Methods
    -------
    simulate(t_end: int, t_print: int, num_firms: int, seed: int, 
             seed_init: int, run_in_period: int, save: str,
             output_folder: str)
    """

    def __init__(self, parameters: dict = None, hyper_parameters: dict = None,
                 variables: list = None,
                 base_call: str = 'mark0_covid', call_path: str = ''):
        """ Initialisation of the Mark0_COVID class. Sets the parameters
        as well as the call function that executes the model.

        Parameters
        ----------
        parameters : dict, default: None
            Dictionary of the parameters used to call the Mark0-COVID
            model.
        hyper_parameters : dict, default: None
            Hyper-parameters for the model incl. runtime, number of firms, seed,
        variables: list[str], optional
            List of output series recorded from the model
        base_call : str, default: 'mark0_covid'
            The baseline function to call in the terminal to run the mark0 model
        call_path : str, default: ''
            Location of the of the mark0_covid.exe file. By default it assumes
            the executable is located in a folder 'mark0_covid/' in the same 
            directory as this python file.
        """

        # Mark-0 COVID Output series
        self.output_series = {
            't': r"Time",
            'u': r"Unemployment Rate (\%)",
            'bust': r"Bankruptcy Rate (\%)",
            'Pavg': r"Average Price $\bar{p}$",
            'Wavg': r"Average Wage $\bar{w}$",
            'S': r"Household Savings $S$",
            'Atot': r"Firm Balance $\sum_i\mathcal{E}_i$",
            'firm-savings': r"Firm deposits $\mathcal{E}^+$",
            'debt-tot': r"Firm liabilities $\mathcal{E}^-$",
            'inflation': r'Inflation rate $\pi$',
            'pi-avg': r"EMA Avg. Inflation $\tilde{\pi}$",
            'propensity': r"Consumption Rate $c$",
            'k': r"Debt-to-Equity Ratio",
            'Dtot': r"Total Demand $\sum_i D_i$",
            'rhom': r"Loan interest rate $\rho^l$",
            'rho': r"Central Bank Rate $\rho_0$",
            'rhop': r"Deposit interest rate $\rho^d$",
            'pi-used': r"Household inflation expectation $\hat{\pi}$",
            'tau_tar': r"Expect $\pi$ weight of CB target inflation",
            'tau_meas': r"Expect $\pi$ weight of EWMA inflation",
            'R': r"Hiring/firing rate",
        }

        # Mark-0 COVID Parameter Names in LaTeX format (for graphing)
        self.parameter_names = {
            'R': r'Hiring/firing rate, $R$',
            'alpha_gamma': r'Loan rate effect on $R$, $\alpha_\Gamma$',
            'gamma_0': r'Baseline gamma, $\Gamma_0$',
            'r': r'Adjustment ratio $\frac{\gamma_w}{\gamma_p}$',
            'gammap': r'Price-adjustment size, $\gamma_p$',
            'eta0m': r'Baseline firing propensity, $\eta_0$',
            'wage_factor': r'Wage adjustment to inflation',
            'rho_star': r'Baseline interest rate, $\rho^\star$',
            'phi_pi': r'CB reaction to inflation, $\phi_\pi$',
            'phi_e': r'CB reaction to employment $\phi_\varepsilon$',
            'pi_star': r'CB inflation target, $\pi^\star$',
            'e_star': r'CB unemployment target, $\hat{\varepsilon}^\star$',
            'theta': r'Default threshold, $\Theta$',
            'f': r'Bankruptcy effect on bank interest rates, $f$',
            'alpha_c': r'Real rate influence on consumption, $\alpha_c$',
            'c_0': r'Baseline propensity to consume $c_0$',
            'delta': r'Dividend share, $\delta$',
            'beta': r'Household intensity of choice, $\beta$',
            'tau_meas': r'Weight $\pi^{ema}$ on $\hat{\pi}$, $\tau^R$',
            'tau_tar': r'Weight $\pi^\star$ on $\hat{\pi}$, $\tau^T$',
            'phi': r'Revival frequency, $\phi$',
            'taupi': r'EWMA memory, $\omega$',
            'y0': r'Initial Production, $y_0$'
        }
        
        # If no variables given, use all available
        if variables == None:
            self.variables = list(self.output_series.keys())
            self.variables.remove('t')
        else:
            self.variables = variables

        self.phase_space = ()

        # Set the call-function passed to the command line
        if call_path == '':
            folder = os.path.dirname(os.path.abspath(__file__))
            self.call_path = folder + '/mark0_covid/' + base_call
        else:
            self.call_path = call_path

        self.model_name = 'mark0_covid'
        self.base_call_func = base_call

        # Check if we have all parameters and augment if necessary
        if parameters is None:
            parameters = self._read_parameters()
        else:
            for k, v in self._read_parameters().items():
                if k not in parameters:
                    parameters[k] = v

        # Check if we have all hyper parameters and augment if necessary
        if hyper_parameters is None:
            hyper_parameters = self._read_parameters(hyper=True)
        else:
            for k, v in self._read_parameters(hyper=True).items():
                if k not in hyper_parameters:
                    hyper_parameters[k] = v

        self.parameters = parameters
        self.hyper_parameters = hyper_parameters

    def simulate(self, t_end: int = None, t_print: int = None,
                 num_firms: int = None, seed: int = None, seed_init: int = None,
                 run_in_period: int = None, save: str = None,
                 output_folder: str = 'output') -> pd.DataFrame:
        """ Run the Mark-0 COVID Model. Optionally save its output as a txt file (tab
        separated), otherwise just return a DataFrame.

        Parameters
        ----------
        t_end : int, optional
            Total runtime of the model
        t_print : int, optional
            At which interval to save outputs
        num_firms : int, optional
            Number of firms in the simulation
        seed : int, optional
            Random seed of the simulation
        seed_init : int, optional
            Random seed for the initial values of the simulation
        run_in_period : int, optional
            How many periods to cut off at the beginning of the simulation. 
        save : str, optional
            filename to save the .txt output. 'mark0_covid' is prepended
        output_folder : str, default: 'output'
            where the output .txt is saved

        Returns
        -------
        output : pd.DataFrame
            time-series of the model simulation
        """
        # Adjust the hyperparameters if necessary
        if t_end:
            self.hyper_parameters['T'] = t_end
        if run_in_period:
            self.hyper_parameters['Teq'] = run_in_period
        if t_print:
            self.hyper_parameters['tprint'] = t_print
        if num_firms:
            self.hyper_parameters['N'] = num_firms
        if seed:
            self.hyper_parameters['seed'] = seed
        if seed_init:
            self.hyper_parameters['seed_init'] = seed_init

        # Decide whether the output .txt is temporary
        if save:
            save = self.model_name + '_' + save
        else:
            save = self.model_name + '_temp'
        
        # Name of save file and model call
        filename = os.sep.join([os.getcwd(), output_folder, save])
        function_call = self._model_call(filename)

        # Generate an output folder to store txt files
        if not os.path.isdir(output_folder):
            print('Created output folder')
            os.makedirs(output_folder)

        # Run the Mark-0 model
        os.system(function_call)

        # Retrieve output as a pandas DataFrame
        try:
            self.output = self._read_output(filename + '.txt')
        except pd.errors.EmptyDataError:
            self.output = pd.DataFrame(columns=self.variables)

        if save==self.model_name + '_temp':
            os.remove(filename + '.txt')
        
        return self.output.loc[:, self.variables]

    def phase_classifier(self) -> str:
        """ Classify the phases into four different sections sequentially, that is
        pick the first one that applies. This is admittedly ad hoc, and it
        should be noted that phases may change over time 
        (e.g. in the long-run EC->FU)
        
        Phases include:
        1. Full Employment (average employment in last cutoff periods < 10%)
        2. Full Unemployment (average unemployment in last cutoff periods > 90%)
        3. Endogenous Crises (std. dev. of employment in last cutoff periods > 0.1)
        4. Residual Unemployment (standard deviation of employment <)

        Returns
        -------
        phase : str (one of "FE", FU", "EC", "RU")
        """
        if self.output.loc[:, 'u'].mean() < 0.1:
            return "FE"
        elif self.output.loc[:, 'u'].mean() > 0.9:
            return "FU"
        elif self.output.loc[:, 'u'].std() > 0.1:
            return "EC"
        else:
            return "RU"

    def _read_output(self, path: str) -> pd.DataFrame:
        """ Read tab separated values from txt file generated by the C++ rules.

        Parameters
        ----------
        path : str
            Location of the tab-separated .txt file of output

        Returns
        -------
        df: pd.DataFrame
            DataFrame of the simulation time-series
        """
        df = pd.read_csv(path, sep="\t", header=None, index_col=False)
        
        assert df.shape[1] == len(self.output_series)
        df.columns = self.output_series
        df.set_index('t', inplace=True)

        return df

    def _model_call(self, output_file: str = '_temp') -> str:
        """ Generate the function call for Mark-0 Covid based on the given 
        parameters 

        Parameters
        ----------
        output_file : str, default: '_temp'
            Name of the output file that will be saved as .txt by program

        Returns
        -------
        call: str
            commandline call for Mark-0
        """
        order_hyper = ('seed', 'seed_init', 'N', 'T', 'Teq', 'tprint')
        items = [self.call_path,
                 *[f"{self.parameters[p]}" for p in self.parameter_names],
                 *[f"{self.hyper_parameters[v]}" for v in order_hyper],
                 output_file]
        return ' '.join(items)

    def _read_parameters(self, hyper: bool = False):
        """
        Read default parameters and hyper parameters from yaml files

        Parameters
        ----------
        hyper : bool, optional
            If true read hyper parameter, else read parameter. The default is False.

        Returns
        -------
        parameter: dict.
            Parameter dict from yaml file
        """
        path = os.path.dirname(os.path.abspath(__file__))

        if hyper:
            with open(path + '/default_hyper_parameters.yaml') as f:
                parameter = yaml.load(f)
        else:
            with open(path + '/default_parameters.yaml') as f:
                parameter = yaml.load(f)

        return parameter


if __name__ == "__main__":
    pass
