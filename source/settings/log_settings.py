import logging
import sys

stream=sys.stdout
level=logging.INFO
format='[%(asctime)s] {%(filename)s:%(lineno)d:%(funcName)s} %(levelname)s - %(message)s'