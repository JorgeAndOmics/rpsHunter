import coloredlogs
import logging
import os

import defaults

# Logging
LEVEL_STYLES: dict = {
    'debug': {'color': 'white'},     # Standard debug level
    'info': {'color': 'cyan', 'bold': 'yes'},     # Standard info level
    'warning': {'color': 'yellow'},  # Standard warning level
    'error': {'color': 'red', 'bold': 'yes'},      # Standard error level
    'critical': {'color': 'black', 'bold': 'yes', 'background': 'red'},  # Standard critical level
}

FIELD_STYLES: dict = {
    'asctime': {'color': 'green'},
    'hostname': {'color': 'green'},
    'levelname': {'color': 'green'},
    'name': {'color': 'green'},
    'programname': {'color': 'green'},
    'username': {'color': 'green'},
    'process': {'color': 'green'},
    'thread': {'color': 'green'}
}


def colored_logging(log_file_name: str) -> None:
    """
    Sets up logging and configures coloredlogs with the custom fields and level styles

        Parameters
        ----------
            :param log_file_name: The name of the file to save the log in.
    """
    # Configure coloredlogs with the custom field and level styles
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(message)s', handlers=[
        logging.FileHandler(os.path.join(defaults.PATH_DICT['LOG_DIR'], log_file_name), mode='w'),
        logging.StreamHandler()
    ]
                        )

    coloredlogs.install(
        level='DEBUG',
        fmt='%(asctime)s - %(message)s',
        level_styles=LEVEL_STYLES,
        field_styles=FIELD_STYLES
    )
