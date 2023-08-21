import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt
import seaborn as sns

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import dash
from dash import dcc, html, ctx, callback, Patch
from dash.dependencies import Input, Output, State

# Model definitions


def gene_expression_eq0(t, x, phi, s, delta, A):
    '''
    t - time span
    x - differentiated variable
    params - set of specific parameters for this equation!
    returns y = dx/dt
    '''

    #try:
    #    phi, s, delta, A = params
    #except ValueError as ve:
    #    print(ve, 'Check that your parameters satisfy Model 1 parameters!')

    return phi * s * A - x*delta

def gene_expression_eq1(t, x, phi, s, delta, n, m, A, E, S, E_trh, S_trh):
    '''
    t - time span
    x - differentiated variable
    params - set of specific parameters for this equation!
    returns y = dx/dt
    '''

    #try:
    #    phi, s, delta, n, m, A, E, S, E_trh, S_trh = params
    #except ValueError as ve:
    #    print(ve, 'Check that your parameters satisfy Model 1 parameters!')

    enhancer_term = np.power(E, n) / (np.power(E, n) + np.power(E_trh, n))
    silencer_term = np.power(S_trh, m) / (np.power(S, m) + np.power(E_trh, m))

    return phi * s * A * enhancer_term * silencer_term - x*delta


def gene_expression_eq2(t, x, phi, s, delta, A, E, S, Ke1, Ke2, Ks1, Ks2):
    '''
    t - time span
    x - differentiated variable
    params - set of specific parameters for this equation!
    returns y = dx/dt
    '''

    #try:
    #    phi, s, delta, A, E, S, Ek1, Ek2, Sk1, Sk2 = params
    #except ValueError as ve:
    #    print(ve, 'Check that your parameters satisfy Model 1 parameters!')

    return phi * s * A * ((E/Ke1 + 1)*(A + Ke2))/((A + Ks2)*(1+S/Ks1)) - x*delta


def solve_eq(eq, x0, params, t_span = (0, 10), nsteps=10):
    # solving the equation numerically
    solution  = solve_ivp(
        eq,
        t_span = t_span,
        y0 = x0,
        args = params,
        t_eval = np.linspace(t_span[0], t_span[1], nsteps)
        )
    return solution


# Generative functions

def sample_parameters(equation_num, params, n_genes,
                      promotor_on_percent=100,
                      only_enh_percent=100, only_sil_percent=0,
                      both_enh_sil_percent=0):
    '''
    equation_num: id of the DE equation, {0, 1, 2}
    params: params to sample, not mentioned params = default
    n_genes: numbe of genes to generate
    promotor_on_percent: percent of genes with A=1
    only_enh_percent: percent of genes with S=0, E=1 AMONG GENES WITH A=1
    only_sil_percent: percent of genes with S=1, E=0 AMONG GENES WITH A=1

    returns: np.array of size (number of params for eq, n_genes)
    '''

    if only_enh_percent+only_sil_percent+both_enh_sil_percent > 100:
        raise ValueError('Sum of percents should be < 100!')
    elif (n_genes < 2) & (
      {only_enh_percent, only_sil_percent, both_enh_sil_percent, promotor_on_percent} != {100, 0}
      ):
        if ({only_enh_percent, only_sil_percent, both_enh_sil_percent} == {0}):
            pass
        else:
            raise ValueError('You have 1 gene, correct percentage (all 0 or one 100)!')

    # common parameters
    delta = np.random.lognormal(mean=0.25, sigma=1.5, size=n_genes) if 'delta' in params else [1.0]*n_genes
    s = np.random.zipf(a=1+0.8, size=n_genes) if 's' in params else [1.0]*n_genes
    phi = np.random.choice([0.0, 0.5, 1.0, 10.0], n_genes) if 'phi' in params else [1.0]*n_genes

    # promotor active state
    # we should derive it from aux equation but here we'll use const
    on_genes = int(n_genes*(promotor_on_percent/100))
    A = np.array([1.0]*on_genes + [0.0]*(n_genes-on_genes))

    # enhancer and silencer
    both_genes = int(on_genes*(both_enh_sil_percent/100))
    enh_genes = both_genes + int(on_genes*(only_enh_percent/100)) # both + only enh
    sil_genes = both_genes + int(on_genes*(only_sil_percent/100)) # both + only sil
    E = np.array([1.0]*enh_genes + [0.0]*(n_genes-enh_genes))
    S = np.array([1.0]*sil_genes + [0.0]*(n_genes-sil_genes))

    match str(equation_num):
        # constants for equation 0
        case '0':
          # gathering all parameters in one array
            params_array = np.column_stack((phi, s, delta, A))

        # constants for equation 1
        case '1':
            # random pick from [0, 2, 4]
            n = np.random.choice([0.0, 2.0, 4.0], n_genes) if 'n' in params else [1.0]*n_genes
            m = np.random.choice([0.0, 2.0, 4.0], n_genes) if 'm' in params else [1.0]*n_genes
            # random pick from [0.01, 0.05, 0.1, 0.5, 1.0, 1.5]
            E_trh = np.random.choice([0.01, 0.05, 0.1, 0.5, 1.0, 1.5], n_genes) if 'E_trh' in params else [0.1]*n_genes
            S_trh = np.random.choice([0.01, 0.05, 0.1, 0.5, 1.0, 1.5], n_genes) if 'S_trh' in params else [0.5]*n_genes

            # gathering all parameters in one array
            params_array = np.column_stack((phi, s, delta, n, m, A, E, S, E_trh, S_trh))

        # constants for equation 2
        case '2':
            # random pick from [0, 0.1, 0.5, 1, 1.5, 5, 10]
            Ke1 = np.random.choice([0.01, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0], n_genes) if 'Ke1' in params else [0.1]*n_genes
            Ke2 = np.random.choice([0.0, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0], n_genes) if 'Ke2' in params else [0.1]*n_genes
            Ks1 = np.random.choice([0.01, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0], n_genes) if 'Ks1' in params else [0.1]*n_genes
            Ks2 = np.random.choice([0.0, 0.1, 0.5, 1.0, 1.5, 5.0, 10.0], n_genes) if 'Ks2' in params else [0.1]*n_genes

            # gathering all parameters in one array
            params_array = np.column_stack((phi, s, delta, A, E, S, Ke1, Ke2, Ks1, Ks2))

        case _:
            raise ValueError('Equation number should be 0, 1 or 2!')

    return params_array




# for mouse add=0.06, mult=0.1
def add_noise(xval, lvl_add_noise=0.08, lvl_mult_noise=0.1):
    '''
    xval: solutions of the equation, x
    lvl_add_noise: level for the additive noise
    lvl_mult_noise: level for the multiplicative noise

    returns: noised solutions in log scale, log x
    '''
    additive_noise = np.random.lognormal(mean=0, sigma=lvl_add_noise*abs(xval))
    multiplicative_noise = np.random.lognormal(mean=0, sigma=lvl_mult_noise)
    log_x = multiplicative_noise*(abs(xval) + additive_noise)

    if np.isinf(log_x):
        return xval
    else:
        return log_x
    
    




def generate_sample(equation, variable_params=[], common_parameters=True,
                    percents_dict= {'promotor_on_percent': 100,
                                    'only_enh_percent': 50,
                                    'only_sil_percent': 50,
                                    'both_enh_sil_percent': 0
                                    },
                    variable_x0=False,
                    size=1000, timepoints=50,
                    white_noise=False):
    '''
    equation: function for the DE equation
    variable_params: which parameters to sample (others will be set to default values)
    common_parameters: if genes within one sample should share parameters' values
    size: number of genes in the sample
    timepoints: number of time points to sample

    returns: np.array of size (size, timepoints)
    '''

    # checking var params
    if not set(variable_params).issubset(set(equation.__code__.co_varnames)):
        raise ValueError('Check that variable parameters match the equation!')

    # generating x0 initial values
    if variable_x0:
        x0 = np.random.lognormal(mean=0.73/(2**(-1.58)), sigma=2.53/(2**(0.73)), size=size)
    else:
        x0 = np.array([10]*size)

    # generating parameters and solving equation
    # sampling parameters for each gene within sample
    if not common_parameters:
        # parameters generation
        params = sample_parameters(equation_num = equation.__name__[-1],
                                   params = variable_params, n_genes=size,
                                   promotor_on_percent = percents_dict['promotor_on_percent'],
                                   only_enh_percent = percents_dict['only_enh_percent'],
                                   only_sil_percent = percents_dict['only_sil_percent'],
                                   both_enh_sil_percent = percents_dict['both_enh_sil_percent']
                                   )
        # solving and saving solutions to array
        def __helper_solv(params_row):
            return solve_eq(equation, x0, params_row, nsteps=timepoints)

        solutions_list =  np.apply_along_axis(__helper_solv, 1, params)
        time_vector = solutions_list[0].t
        solution_array = np.array([sol.y[0] for sol in solutions_list])


    # same parameters for all genes whithin sample
    # just sampling params for one gene --> solving --> replicating * n_genes
    else:
        # parameters generation
        promotor_on_percent = 100 if percents_dict['promotor_on_percent'] > 50 else 0
        if percents_dict['both_enh_sil_percent'] >= 50:
            both_enh_sil_percent, only_enh_percent, only_enh_percent = 100, 0, 0
        else:
            both_enh_sil_percent = 0
            only_enh_percent, only_sil_percent = (100, 0) if percents_dict['only_enh_percent'] > percents_dict['only_sil_percent'] else (0, 100)

        params = sample_parameters(equation_num = equation.__name__[-1],
                                   params = variable_params, n_genes=1,
                                   promotor_on_percent = promotor_on_percent,
                                   only_enh_percent = only_enh_percent,
                                   only_sil_percent = only_sil_percent,
                                   both_enh_sil_percent = both_enh_sil_percent
                                   )

        # solving and saving solutions to array
        solution_one_gene = solve_eq(equation, [x0[0]], params[0], nsteps=timepoints)
        solution_array = np.tile(solution_one_gene.y[0], (size, 1))
        time_vector = solution_one_gene.t

    # adding the noise
    if white_noise:
        gene_expression_levels_noised = np.vectorize(add_noise)(solution_array)
        return params, gene_expression_levels_noised, time_vector

    else:
        return params, solution_array, time_vector


# Sample generations

def generate_steady_states(common_parameters=True,
                           variable_params=[],
                           percents_dict= {'promotor_on_percent': 100,
                                           'only_enh_percent': 100,
                                           'only_sil_percent': 0,
                                           'both_enh_sil_percent': 0},
                           variable_x0=False,
                           size=50, timepoints=50,
                           white_noise=True,
                           store_params=True):
    '''
    This function takes same parameters as generate_sample function with
    addition of store_params flag that regulates if function returns dataframe
    of parameters for each gene of not.

    returns: df of solutions at the last time point for 2 samples
    for equations 0, 1, and 2
    '''

    df = pd.DataFrame()
    # can be super big if we have many genes --> store_params=False recommended
    if store_params:
        columns = {'sample': [], 'equation': [],
                   'phi': [], 's': [], 'delta': [],
                   'A': [], 'E': [], 'S': [],
                   'n': [], 'm': [], 'E_trh': [], 'S_trh': [],
                   'Ke1': [], 'Ke2': [], 'Ks1': [], 'Ks2': []}
        params_df = pd.DataFrame(columns)

    # generating 2 samples for each of 3 functions
    for eq_id, eq_func in enumerate([gene_expression_eq0,
                                  gene_expression_eq1,
                                  gene_expression_eq2]):
        for sample in range(2):
            variable_params = set(variable_params).intersection(set(eq_func.__code__.co_varnames))
            params, solutions, time = generate_sample(eq_func,
                                                    variable_params=variable_params,
                                                    common_parameters=common_parameters,
                                                    percents_dict=percents_dict,
                                                    variable_x0=variable_x0,
                                                    size=size, timepoints=timepoints,
                                                    white_noise=white_noise)

            df[f'eq{eq_id}_s{sample}'] = solutions[:, 1]

            # handling parameters
            if store_params:
                if common_parameters:
                    params = np.tile(params[0], (size, 1))
                # equation-specific parameters
                match eq_id:
                    case 0:
                        params_local = pd.DataFrame(params, columns = [
                          'phi', 's', 'delta', 'A'
                          ])
                    case 1:
                        params_local = pd.DataFrame(params, columns = [
                          'phi', 's', 'delta', 'n', 'm', 'A', 'E', 'S', 'E_trh', 'S_trh'
                          ])
                    case 2:
                        params_local = pd.DataFrame(params, columns = [
                          'phi', 's', 'delta', 'A', 'E', 'S', 'Ke1', 'Ke2', 'Ks1', 'Ks2'
                          ])
                # equation and sample info        
                params_local['sample'] = [sample]*size
                params_local['equation'] = [eq_id]*size
                params_df = pd.concat([params_df, params_local])

    return df if not store_params else (df, params_df)

# APP CREATION


# app creation
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(children=[
    # header with parameters for generation
    html.H4('Select parameters for samples generation:'),
    
    # parameters block 1
    html.Div(children=[
        html.Div([
            'Common parameters for genes within sample:',
            dcc.Dropdown(
                ['True', 'False'],
                'True',
                id='common_params'
            )],
            style={'width': '40%', 'float': 'left', 'display': 'inline-block'}
        ),
        html.Div([
            'Variable x0:',
            dcc.Dropdown(
                ['True', 'False'],
                'True',
                id='variable_x0'
            )],
            style={'width': '20%', 'display': 'inline-block'}
        ),
        html.Div([
            'Number of genes:',
            dcc.Dropdown(
                [20, 100, 1000, 10000],
                100,
                id='size'
            )],
            style={'width': '20%', 'float': 'right', 'display': 'inline-block'}
        ),
        html.Div([
            'Number of timepoints:',
            dcc.Dropdown(
                [10, 20, 50, 100],
                10,
                id='time'
            )],
            style={'width': '20%', 'float': 'right', 'display': 'inline-block'}
        ),
    ], style={'padding-top': '20px', 'padding-bottom': '20px'}),
    
    # parameters block 2
    html.Div(children=[
        html.Div([
            'Variable parameters:',
            dcc.Dropdown(
                ['phi', 's', 'n', 'm', 'E_trh', 'S_trh', 'Ke1', 'Ke2', 'Ks1', 'Ks2',],
                'None',
                multi=True,
                id='variable_params'
            )],
            style={'width': '40%', 'display': 'inline-block'}
        ),
        html.Div([
            'Add white noise:',
            dcc.Dropdown(
                ['True', 'False'],
                'True',
                id='noise'
            )],
            style={'width': '20%', 'display': 'inline-block'}
        ),
    ], style={'padding-bottom': '20px'}),
    
    # parameters block 3
    html.Div(children=[
        html.Div([
            '% genes with P=1 (promotor ON):',
            dcc.Dropdown(
                [0, 25, 50, 100],
                100,
                id='pA'
            )],
            style={'width': '45%', 'float': 'left', 'display': 'inline-block'}
        ),
        html.Div([
            '% genes with A=1, I=1 (both Enh and Sil ON):',
            dcc.Dropdown(
                [0, 25, 50, 100],
                0,
                id='bothES'
            )],
            style={'width': '45%', 'float': 'right', 'display': 'inline-block'}
        ),
    ], style={'padding-bottom': '20px'}),
    
    # parameters block 4
    html.Div(children=[
        html.Div([
            '% genes with A=1, I=0 (only enhancer ON):',
            dcc.Dropdown(
                [0, 25, 50, 100],
                100,
                id='pE'
            )],
            style={'width': '45%', 'float': 'left', 'display': 'inline-block'}
        ),
        html.Div([
            '% genes with A=0, I=1 (only silencer ON):',
            dcc.Dropdown(
                [0, 25, 50, 100],
                0,
                id='pS'
            )],
            style={'width': '45%', 'float': 'right', 'display': 'inline-block'}
        ),
    ], style={'padding-bottom': '20px'}),
     
    # add generation button
    html.Div(children=[
        html.Button(id='button', n_clicks=0, children='Generate samples!'),
    ], style= {'padding-top': '100px', 'padding-bottom': '50px'}),
    
    # the plot itself
    html.Div(children=[
        dcc.Graph(id='ode-graph'),
    ], style={'width': '100%'}),
    
    # add dropdowns for plotting
    html.Div(children=[
        html.Div(children=[
            'Color the points according to:',
            dcc.Dropdown(
                ['phi, sample 1', 'phi, sample 2',
                's, sample 1', 's, sample 2',
                'delta, sample 1', 'delta, sample 2', 
                'P, sample 1', 'P, sample 2',
                'A, sample 1', 'A, sample 2',
                'I, sample 1', 'I, sample 2',
                'E_trh, sample 1', 'E_trh, sample 2',
                'S_trh, sample 1', 'S_trh, sample 2',
                'n, sample 1', 'n, sample 2',
                'm, sample 1', 'm, sample 2',
                'Ke1, sample 1', 'Ke1, sample 2',
                'Ke2, sample 1', 'Ke2, sample 2',
                'Ks1, sample 1', 'Ks1, sample 2',
                'Ks2, sample 1', 'Ks2, sample 2'
                ], 
                id='color'
            )
        ], style={'width': '33.3%', 'float': 'left', 'display': 'inline-block'}),
        
        html.Div(children=[
            'Style the points according to:',
            dcc.Dropdown(
                ['phi, sample 1', 'phi, sample 2',
                's, sample 1', 's, sample 2',
                'delta, sample 1', 'delta, sample 2', 
                'P, sample 1', 'P, sample 2',
                'A, sample 1', 'A, sample 2',
                'I, sample 1', 'I, sample 2',
                'E_trh, sample 1', 'E_trh, sample 2',
                'S_trh, sample 1', 'S_trh, sample 2',
                'n, sample 1', 'n, sample 2',
                'm, sample 1', 'm, sample 2',
                'Ke1, sample 1', 'Ke1, sample 2',
                'Ke2, sample 1', 'Ke2, sample 2',
                'Ks1, sample 1', 'Ks1, sample 2',
                'Ks2, sample 1', 'Ks2, sample 2'
                ], 
                id='style'
            )
        ], style={'width': '33.3%', 'float': 'right', 'display': 'inline-block'}),
        
        html.Div(children=[
            'Size of the points according to:',
            dcc.Dropdown(
                ['phi, sample 1', 'phi, sample 2',
                's, sample 1', 's, sample 2',
                'delta, sample 1', 'delta, sample 2', 
                'P, sample 1', 'P, sample 2',
                'A, sample 1', 'A, sample 2',
                'I, sample 1', 'I, sample 2',
                'E_trh, sample 1', 'E_trh, sample 2',
                'S_trh, sample 1', 'S_trh, sample 2',
                'n, sample 1', 'n, sample 2',
                'm, sample 1', 'm, sample 2',
                'Ke1, sample 1', 'Ke1, sample 2',
                'Ke2, sample 1', 'Ke2, sample 2',
                'Ks1, sample 1', 'Ks1, sample 2',
                'Ks2, sample 1', 'Ks2, sample 2'
                ],
                id='point_size'
            )
        ], style={'width': '33.3%', 'float': 'right', 'display': 'inline-block'}),
    ]),
    
    html.Div(children=[
        html.Button(id='button-style', n_clicks=0, children='Change the style!'),
    ], style= {'padding-top': '100px', 'padding-bottom': '50px'}),
])

# interactive part

# data generation callback
'''
    common_parameters=False,
    variable_params=['phi', 's', 'Ke2', 'Ks2', 'n', 'm'],
    percents_dict= {'promotor_on_percent': 100,
       'only_enh_percent': 50,
       'only_sil_percent': 50,
       'both_enh_sil_percent': 0},
    variable_x0=True,
    size=100,
    timepoints=50,
    white_noise=False,
'''
@callback(
    Output('ode-graph', 'figure', allow_duplicate=True),
    Input('button', 'n_clicks'),
        State('common_params', 'value'),
        State('variable_x0', 'value'),
        State('size', 'value'),
        State('time', 'value'),
        State('variable_params', 'value'),
        State('noise', 'value'),
        State('pA', 'value'),
        State('bothES', 'value'),
        State('pE', 'value'),
        State('pS', 'value'),
    prevent_initial_call='initial_duplicate'
)
def update_graph(button, common_params, variable_x0, size, time, variable_params, noise, 
                 pA, bothES, pE, pS):
    # generating samples - basic options
    common_params, variable_x0, size, time, variable_params, noise, pA, bothES, pE, pS = (
        common_params == 'True', 
        variable_x0 == 'True',
        int(size), int(time), list(variable_params), 
        noise == 'True', int(pA), int(bothES), int(pE), int(pS)
    )
    
    fig = make_subplots(rows=1, cols=3)
    
    # regenerating if button is clicked
    if "button" == ctx.triggered_id:
        percents = {'promotor_on_percent': pA,
                    'only_enh_percent': pE,
                    'only_sil_percent': pS,
                    'both_enh_sil_percent': bothES}

        df, params = generate_steady_states(common_parameters=common_params,
                                            variable_params=variable_params,
                                            percents_dict=percents,
                                            variable_x0=variable_x0,
                                            size=size, timepoints=time,
                                            white_noise=noise,
                                            store_params=True)

    else:
        df, _ = generate_steady_states()
        
    for i in range(3):
        fig.add_scatter(x=df[f'eq{i}_s0'],
                        y = df[f'eq{i}_s1'], mode='markers', row=1, col=i+1,
                        name=f'Equation {i}')
        #upper_lim = max(max(df[f'eq{i}_s0']),
        #                max(df[f'eq{i}_s1']))
        #lower_lim = min(min(df[f'eq{i}_s0']),
        #                min(df[f'eq{i}_s1']))
        #fig.update_yaxes(range=[lower_lim*0.9, upper_lim*1.1], row=1, col=i+1)
        #fig.update_xaxes(range=[lower_lim*0.9, upper_lim*1.1], row=1, col=i+1)

    fig.update_layout(xaxis_title='sample 1', yaxis_title='sample 2',
                      title='Scatterplots of steady states in sample 1 vs sample 2:')
    return fig

    
@callback(
    Output('ode-graph', 'figure', allow_duplicate=True),
    Input('button-style', 'n_clicks'),
    State('color', 'value'),
    State('style', 'value'),
    State('point_size', 'value'), 
    prevent_initial_call=True
)
def update_styling(button, color, style, point_size):
    if "button-style" == ctx.triggered_id:
        # updating styling if its not none
        patched_figure = Patch()

        # if parameter is common = update everywhere, if specific = update only one subplot
        def _check_param(par):
            if par in ['Ke1', 'Ke2', 'Ks1', 'Ks2']:
                return 2
            elif par in ['n', 'm', 'E_trh', 'S_trh']:
                return 1
            else:
                return 0

        if color:
            # check the params df - how to make a np array from (2, n) and col=sample
            eq = _check_param(color[:color.find(',')])
            color_by = params[(params['sample'].values == color[-1]) & (params['equation'].values == eq)][color[:color.find(',')]]
            patched_figure['data'][0]['marker']['color'] = color_by
        if style:
            eq = _check_param(style[:style.find(',')])
            style_by = params[(params['sample'].values == style[-1]) & (params['equation'].values == eq)][style[:style.find(',')]]
            patched_figure['data'][0]['marker']['style'] = style_by
            patched_figure['data'][0]
        if point_size:
            eq = _check_param(point_size[:point_size.find(',')])
            size_by = params[(params['sample'].values == point_size[-1]) & (params['equation'].values == eq)][point_size[:point_size.find(',')]]
            patched_figure['data'][0]['marker']['size'] = size_by

        return patched_figure


# plotting callbacks

if __name__ == '__main__':
    app.run(debug=True, port='7801')  # by default host='127.0.0.1', port='8050'
