import logging
import sys




stream=sys.stdout
level=logging.INFO
format='[%(asctime)s] {%(filename)s:%(lineno)d:%(funcName)s} %(levelname)s - %(message)s'





logging.basicConfig(
        stream=stream,
        level=level,
        format=format
        )


logger = logging.getLogger(__name__)