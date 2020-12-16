LOGGING_CONFIG = {
    'version': 1,
    'disable_existing_loggers': True,
    'loggers': {
        '': {  # root logger
            'level': 'ERROR',
            'handlers': ['info_rotating_file_handler'],
            'propagate': True
        }
    },
    'handlers': {
        'info_rotating_file_handler': {
            'level': 'ERROR',
            'formatter': 'info',
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': 'logs/info.log',
            'mode': 'a',
            'maxBytes': 100000,
            'backupCount': 20
        }
    },
    'formatters': {
        'info': {
            'format': '[%(asctime)s] %(levelname)s [%(name)s::%(module)s.%(funcName)s:%(lineno)d] %(message)s',
            'datefmt': '%m-%d-%Y@%H:%M:%S'
        }
    }
}
