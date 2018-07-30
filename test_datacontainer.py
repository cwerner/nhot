import pytest
import numpy as np

from datacontainer import Data


# constant


# fixures
@pytest.fixture
def data():
    return Data()

# helpers
def verify_answer(expected, answer, last_answer):
    assert expected == answer
    assert expected == last_answer

# test cases

def test_add(data):
    data.add(np.ones(10), colname='Test')
    assert len(data.df) == 10  
    assert 'Test' in data.df.columns
    
def test_add_twice(data):
    data.add(np.ones(10), colname='Test')
    data.add(np.ones(10), colname='Test')
    assert len(data.df.columns) == 1

def test_add_wrong_size(data):
    data.add(np.ones(10), colname='Test1')
    data.add(np.ones(11), colname='Test2')
    assert data.df.columns == ['Test1']
    assert len(data.df) == 10


     