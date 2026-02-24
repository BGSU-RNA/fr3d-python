# New fr3d web server September 2024

import sys
import logging
import traceback

sys.path.insert(0, '/var/www/fr3d')
sys.path.insert(1, '/var/www/fr3d/app')
sys.path.insert(2, '/usr/local/pipeline/fr3d-python')
sys.path.insert(3, '/usr/local/pipeline/fr3d-python/fr3d')
sys.path.insert(4, '/usr/local/pipeline/fr3d-python/classifiers')

# logging for regular users
logging.basicConfig(filename='/var/www/fr3d/flask.log', level=logging.DEBUG, format='%(asctime)s %(message)s')

logging.info("==============================")
logging.info("Starting app import in wsgi.py")

from app import app as application  # Import your Flask app

logging.info("Starting wsgi run")

# these are apparently optional
application.config['ENV'] = 'production'
application.config['DEBUG'] = False

logging.info('End of wsgi.py')
