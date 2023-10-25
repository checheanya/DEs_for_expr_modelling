# -*- coding: utf-8 -*-
"""## Imports"""

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

import plotly.express as px
import plotly.graph_objects as go
import dash
from dash import dcc, html
from dash.dependencies import Input, Output

"""## Model definitions"""

def gene_expression_eq1(t, x, params):
    '''
    t - time span
    x - differentiated variable
    params - set of specific parameters for this equation!
    returns: x(t)
    '''

    try:
        phi, s, delta, n, m, A, E, S, E_trh, S_trh = params
    except ValueError as ve:
        print(ve, 'Check that your parameters satisfy Model 1 parameters!')

    enhancer_term = np.power(E, n) / (np.power(E, n) + np.power(E_trh, n))
    silencer_term = np.power(S_trh, m) / (np.power(S, m) + np.power(E_trh, m))
    x_t = phi * s * A * enhancer_term * silencer_term - x*delta

    return x_t


def gene_expression_eq2(t, x, params):
    '''
    t - time span
    x - differentiated variable
    params - set of specific parameters for this equation!
    returns: x(t)
    '''

    try:
        phi, s, delta, A, E, S, Ek1, Ek2, Sk1, Sk2 = params
    except ValueError as ve:
        print(ve, 'Check that your parameters satisfy Model 1 parameters!')


    x_t = phi * s * A * ((E/Ek1 + 1)*(A + Ek2))/((A + Sk2)*(1+S/Sk1)) - x*delta

    return x_t

def solve_eq(eq, x0, params, t_span = (0, 10), nsteps=10):
    '''
    eq - function for DE
    x0 - initial values (np.array)
    params - set of specific parameters for this equation (list)
    returns: y = dx/dt
    '''
    # solving the equation numerically
    solution  = solve_ivp(
        eq,
        t_span = t_span,
        y0 = x0,
        args = [params],
        t_eval = np.linspace(t_span[0], t_span[1], nsteps)
        )
    return solution

"""## Dash app - plotting two equations"""

app = dash.Dash(__name__)

app.layout = html.Div(children=[
       html.Div([
          dcc.Graph(id='ode-graph'),
          ]) ,

       # base parameters
       html.Div([
          html.H4('Common parameters:'),
          html.H6('Phi'),
        dcc.Slider(id='slider-phi', min=0.0, max=5.0, step=1.0, value=1.0),
          html.H6('s'),
        dcc.Slider(id='slider-s', min=0.0, max=5.0, step=1.0, value=1.0),
          html.H6('delta'),
        dcc.Slider(id='slider-delta', min=0.0, max=5.0, step=1.0, value=1.0),
          html.H6('A'),
        dcc.Slider(id='slider-A', min=0.0, max=2.0, step=0.5, value=1.0),
          html.H6('E'),
        dcc.Slider(id='slider-E', min=0.0, max=2.0, step=0.5, value=1.0),
          html.H6('S'),
        dcc.Slider(id='slider-S', min=0.0, max=2.0, step=0.5, value=1.0),

          html.H4('Parameters for equation 2:'),
          html.H6('Ke1'),
        dcc.Slider(id='slider-Ke1', min=0.1, max=3.6, step=0.5, value=1.0),
          html.H6('Ke2'),
        dcc.Slider(id='slider-Ke2', min=0.1, max=3.6, step=0.5, value=1.0),
          html.H6('Ks1'),
        dcc.Slider(id='slider-Ks1', min=0.1, max=3.6, step=0.5, value=1.0),
          html.H6('Ks2'),
        dcc.Slider(id='slider-Ks2', min=0.1, max=3.6, step=0.5, value=1.0),

          html.H4('Parameters for equation 1:'),
          html.H6('n'),
        dcc.Slider(id='slider-n', min=0, max=8, step=2, value=2),
          html.H6('m'),
        dcc.Slider(id='slider-m', min=0, max=8, step=2, value=2),
          html.H6('E_trh'),
        dcc.Slider(id='slider-E_trh', min=0.1, max=3.6, step=0.5, value=1.0),
          html.H6('S_trh'),
        dcc.Slider(id='slider-S_trh', min=0.1, max=3.6, step=0.5, value=1.0),
       ])
])


@app.callback(
    Output('ode-graph', 'figure'),
    [
        Input('slider-phi', 'value'),
        Input('slider-s', 'value'),
        Input('slider-delta', 'value'),
        Input('slider-A', 'value'),
        Input('slider-E', 'value'),
        Input('slider-S', 'value'),

        Input('slider-Ke1', 'value'),
        Input('slider-Ke2', 'value'),
        Input('slider-Ks1', 'value'),
        Input('slider-Ks2', 'value'),

        Input('slider-n', 'value'),
        Input('slider-m', 'value'),
        Input('slider-E_trh', 'value'),
        Input('slider-S_trh', 'value')
    ]
)
def update_graph(phi, s, delta, A, E, S,
                 Ek1, Ek2, Sk1, Sk2,
                 n, m, E_trh, S_trh):
    x0 = np.array([10])

    params1 = [phi, s, delta, n, m, A, E, S, E_trh, S_trh]
    params2 = [phi, s, delta, A, E, S, Ek1, Ek2, Sk1, Sk2]

    sol1 = solve_eq(gene_expression_eq1, x0, params1)
    sol2 = solve_eq(gene_expression_eq2, x0, params2)

    sol_df = pd.DataFrame({'solution, dx/td': sol1.y.tolist()[0] + sol2.y.tolist()[0],
                           'time, units': sol1.t.tolist() + sol2.t.tolist(),
                           'equation': [1]*len(sol1.t) + [2]*len(sol2.t)})

    #fig = go.Figure()
    #fig.add_trace(go.Scatter(x=selected_df['sol'], y=selected_df['time'], mode='lines', ))
    #fig.update_layout(title='ODE Solution', xaxis_title='Time', yaxis_title='Solution')
    fig = px.line(
        sol_df, x="time, units", y="solution, dx/td", color='equation')

    return fig

if __name__ == '__main__':
    app.run(debug=True)  # by default host='127.0.0.1', port='8050'
