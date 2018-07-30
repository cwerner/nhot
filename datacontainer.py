import numpy as np
import pandas as pd

class Data(object):
    """Data class holding 1d array information. 
    """
    ONES = None
    ZEROS = None

    def __init__(self):
        self.df = pd.DataFrame()
    def add(self, data, colname='Col'):

        # ravel data to 1d format
        if len(data.shape) > 1:
            data = data.ravel()

        if colname not in self.df.columns:
            if len(self.df.columns) > 0:
                if len(self.df) == len(data):
                    self.df[colname] = data
            else:
                self.df[colname] = data
            
            # set default properties
            if Data.ONES is None:
                Data.ONES = np.ones(data.shape)
                Data.ZEROS = np.zeros(data.shape)

