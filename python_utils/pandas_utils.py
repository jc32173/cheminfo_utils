import pandas as pd
import random

import sys
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO


def read_str_into_df(string, **kwargs):
    """
    Read csv string into DataFrame

    >>> read_str_into_df("col1;col2;col3\\n1;4.4;99\\n2;4.5;200", sep=';')
       col1  col2  col3
    0     1   4.4    99
    1     2   4.5   200
    """

    data = StringIO(string)
    df = pd.read_csv(data, **kwargs)
    return df


def conv_df_to_str(df, **kwargs):
    """
    Convert DataFrame to csv string

    >>> df = pd.DataFrame(data=[[1, 4.4, 99], [2, 4.5, 200]], 
    ...                   columns=['col1', 'col2', 'col3'], 
    ...                   index=[0, 1])
    >>> conv_df_to_str(df)
    ',col1,col2,col3\\n0,1,4.4,99\\n1,2,4.5,200\\n'
    """

    string = StringIO()
    df.to_csv(string, **kwargs)
    return string.getvalue()


def choose_random_element(series, rand_seed=None):
    """
    Choose random element from lists in a Series or DataFrame column.
    """

    # Reproducible random choice:
    repro_random = random.Random(rand_seed)
    rand_choice = []
    for ls in series.apply(eval):
        rand_choice.append(repro_random.choice(ls)[1])
    return rand_choice


def expand_lists_to_columns(series, 
                            col_prefix='col_', 
                            pad_with_final_value=True):
    """
    Convert column of lists to separate columns and optionally pad final 
    columns from short lists with final values in list.
    
    >>> series = pd.Series([[1, 2, 3], [2, 3]])
    >>> series
    0    [1, 2, 3]
    1       [2, 3]
    dtype: object
    >>> expand_lists_to_columns(series, col_prefix='val_')
       val_0  val_1  val_2
    0    1.0    2.0    3.0
    1    2.0    3.0    3.0
    >>> expand_lists_to_columns(series, col_prefix='val_', 
    ...                         pad_with_final_value=False)
       val_0  val_1  val_2
    0    1.0    2.0    3.0
    1    2.0    3.0    NaN
    """

    series = series.squeeze()

    # Maintain order of lists in resulting columns:
    series = series.apply(lambda x: [[i, x_i] for i, x_i in enumerate(x)])

    df_expl = series.explode(ignore_index=False)\
                    .to_frame()

    ls_col = df_expl.columns[0]

    df_expl['order'] = df_expl[ls_col].apply(lambda x: x[0])
    df_expl['val'] = df_expl[ls_col].apply(lambda x: x[1])

    df_expl.set_index('order', append=True, inplace=True)

    df_expl.drop(columns=[ls_col], inplace=True)

    df_expl = df_expl.unstack(level=-1)\
                     .droplevel(0, axis=1)
    df_expl.columns = [col_prefix+str(i) for i in df_expl.columns]

    if pad_with_final_value:
        df_expl.fillna(method='ffill', axis=1, inplace=True)

    return df_expl
