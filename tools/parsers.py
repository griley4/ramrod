"""
This module contains parsers for madx output files
"""
import pandas as pd




def read_twiss_file(filename):
    """
    Read the twiss datafile and fills in a dictionary

    Parameters
    ----------
    filename : string
        Name of the twiss file

    Returns
    -------
    header : dict
        Dictionnary of the twiss file header variables

    twiss : pandas.DataFrame
        dataframe of all the data present in the twiss file

    """


    with open(filename, 'r') as f:

        header = dict()
        count = 0

        #Parsing of the global variable
        while True:
            line = f.readline()
            count += 1
            if '@' not in line.strip()[0]:
                if '*' in line.strip()[0]:
                    break
                else:
                    raise SyntaxError('Twiss file missin *, label line :' + \
                                      line)

            words = line.split()
            header[words[1]] = line[line.find('"'):].strip(' "\n') \
                              if ('s' in words[2]) else float(words[3])


        columns = line.strip('* \n').split()
        # formats = f.readline().strip('$ \n').split()


    twiss = pd.read_csv(filename, skiprows=count+1, names=columns, sep=r'\s+')

    return header, twiss



def parse_trackone(filename, nmax=0):
    """
    Reads trackone file

    Returns :
    Dataframe witht the data

    Warnings
    --------
    This function is no stable
    """
    f = open(filename, 'r')
    df = pd.DataFrame(columns=['keyword', 'id', 'X', 'PX', 'Y', 'PY', 'PT'])

    while True:
        line = f.readline()
        if line.find('$') != -1:
            break


    while True:
        line = f.readline().split()
        if len(line) == 0:
            break
        keyword = line[5]
        nparts = int(line[3])
        parts = []
        for _ in range(nparts):
            line = f.readline().split()
            parts += [[keyword, int(line[0]),
                       float(line[2]), float(line[3]),
                       float(line[4]), float(line[5]),
                       float(line[7])]]

        dftemp = pd.DataFrame(parts, columns=df.columns)
        if nmax == 0:
            df = pd.concat([df, dftemp.sort_index()])
        else:
            df = pd.concat([df, dftemp.sample(min(len(dftemp), nmax)).sort_index()])
        del dftemp, parts

    f.close()
    return df

def parse_losses(filename):
    """
    Reads losses file

    Returns:
    DataFrame with the data

    Warnings
    --------
    This function is no stable
    """
    f = open(filename, 'r')

    while True:
        line = f.readline()
        if line.find('$') != -1:
            break

    parts = []
    while True:
        line = f.readline().split()
        if len(line) < 5:
            break
        parts += [[int(line[0]), line[10].strip('"'),
                   float(line[2]), float(line[3]),
                   float(line[4]), float(line[5]),
                   float(line[7]), float(line[8])]]

    df = pd.DataFrame(parts, columns=['id', 'keyword', 'X', 'PX', 'Y', 'PY', 'PT', 'S'])

    f.close()
    return df
