import mordred
from mordred import Calculator, descriptors
import numpy as np


all_descriptors = Calculator(descriptors).descriptors

descriptors2D = [d for d in Calculator(descriptors).descriptors if not d.require_3D]

descriptors3D = [d for d in Calculator(descriptors).descriptors if d.require_3D]


def calc_mordred_desc(cmpd_ls, id_ls=[], ignore_3D=False, descriptor_ls=descriptors, rm_const_desc=True):
    """
    Function to calculate mordred descriptors.
    """

    if isinstance(cmpd_ls, str):
        cmpd_ls = [cmpd_ls]

    # Set up index:
    if len(id_ls) == 0:
        if isinstance(cmpd_ls[0], str):
            id_ls = cmpd_ls
        else:
            id_ls = range(len(cmpd_ls))

    calc = Calculator(descriptor_ls, ignore_3D=ignore_3D)
    df_desc = calc.pandas(cmpd_ls)
    df_desc.index = id_ls

    # Remove any column which contains missing values, based on entries of type: mordred.error.Missing
    # or error, based on mordred.error.Error:
    miss_descs = []
    err_descs = []
    obj_descs = []
    const_desc = []
    for col in df_desc.columns:
        # Check for mordred.error.Missing:
        if np.any(df_desc[col].apply(lambda row: type(row)) == mordred.error.Missing):
            miss_descs.append(col)
        # Check for mordred.error.Error:
        elif np.any(df_desc[col].apply(lambda row: type(row)) == mordred.error.Error):
            err_descs.append(col)
        # Check for other objects:
        elif df_desc[col].dtype == object:
            try:
                pd.to_numeric(df_desc[col])
            except TypeError:
                obj_descs.append(col)
        # Remove any descriptors which are constant for all molecules:
        elif rm_const_desc and np.all(df_desc[col].iloc[0] == df_desc[col]):
            const_desc.append(col)
    print('Descriptors with missing values: {}'.format(len(miss_descs)))
    print('Descriptors with errors: {}'.format(len(err_descs)))
    print('Descriptors which are not numeric: {}'.format(len(obj_descs)))
    print('Descriptors which are constant: {}'.format(len(const_desc)))
    df_desc.drop(list(set(miss_descs + err_descs + obj_descs + const_desc)), axis=1, inplace=True)
    return df_desc
