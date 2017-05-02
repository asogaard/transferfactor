# In order to make classes available as e.g.
#   transferfactor.calculator
# instead of
#   transferfactor.calculator.calculator

__all__ = ['calculator', 'config', 'utils']

from calculator import calculator
from config     import *
from utils      import *